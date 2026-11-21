import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gtotree_worker as worker
import gtotree_report as report
from job_validation import read_accessions, read_metadata, validate_manifest, bound_resource_flags


def metadata(n=1):
    return {"version": 1, "source": "assembly-grasp", "created_at": "2026-11-20", "genomes": [
        {"assembly_accession": "GCF_%09d.1" % i, "organism_name": "Example bacterium", "biosample": "SAMN1"} for i in range(n)]}


class FakeS3:
    def __init__(self):
        self.objects = {}; self.events = []
    def get_paginator(self, name):
        parent = self
        class Paginator:
            def paginate(self, Bucket, Prefix=""):
                yield {"Contents": [{"Key": k, "Size": len(v)} for (b,k),v in parent.objects.items() if b==Bucket and k.startswith(Prefix)]}
        return Paginator()
    def head_object(self, Bucket, Key):
        if (Bucket, Key) not in self.objects:
            raise worker.ClientError({"Error": {"Code": "404"}}, "HeadObject")
        return {}
    def get_object(self, Bucket, Key):
        return {"Body": io.BytesIO(self.objects[Bucket,Key])}
    def put_object(self, Bucket, Key, Body, **kwargs):
        self.objects[Bucket,Key] = Body; self.events.append(Key)
    def download_file(self, bucket, key, target):
        Path(target).write_bytes(self.objects[bucket,key])
    def upload_file(self, filename, bucket, key, ExtraArgs=None):
        self.objects[bucket,key] = Path(filename).read_bytes(); self.events.append(key)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.root=Path(self.tmp.name)
        self.args=SimpleNamespace(results_bucket="results",input_bucket="inputs",work_root=str(self.root),force=False,
            retry_failed=False,max_upload_bytes=100*1024*1024,max_references=100,max_uploaded_genomes=100,
            max_total_genomes=200,max_threads=16,job_timeout=30,command_prefix="",gtotree_bin="GToTree",
            report_script="",contact_email="ark@midauthorbio.com",continue_on_error=False,clean=False)
    def tearDown(self): self.tmp.cleanup()
    def job(self):
        s3=FakeS3(); code="Example123"
        form="Name: Test run\nSCG set: Bacteria\nAccessions File: accessions.txt\nGenome File: mine.fna\nGenome Metadata File: selected-genomes.json\nGtotreeFlags: -t -n 999 -M 99 -j 99\n"
        for name,content in {"form-data.txt":form,"accessions.txt":"gcf_000000000.1\nGCF_000000000.1\n",
                             "mine.fna":">contig_1\nATGC\n","selected-genomes.json":json.dumps(metadata())}.items():
            s3.objects["inputs",code+"/"+name]=content.encode()
        return s3,code
    def fake_command(self,cmd,log_handle,timeout=30):
        if cmd[0]=="GToTree":
            out=Path(cmd[cmd.index("-o")+1]);out.mkdir(parents=True)
            (out/"GToTree_output.tre").write_text("(GCF_000000000.1:0.1,mine:0.2);")
            (out/"Aligned_SCGs.faa").write_text(">GCF_000000000.1\nMAAA\n>mine\nMAAT\n")
            (out/"Genomes_summary_info.tsv").write_text("genome\thits\nGCF_000000000.1\t10\nmine\t10\n")
            self.command=cmd
        else:
            subprocess.run([sys.executable,*cmd[1:]],stdout=log_handle,stderr=subprocess.STDOUT,check=True,timeout=timeout)
    def test_combined_job_report_publication_and_final_marker(self):
        s3,code=self.job()
        with patch.object(worker,"run_command",side_effect=self.fake_command): worker.process_job(self.args,s3,code)
        final=json.loads(s3.objects["results",code+"/status.json"])
        self.assertEqual(final["state"],"complete")
        manifest=json.loads(s3.objects["results",code+"/result.json"])
        self.assertEqual(manifest["input_summary"],{"reference_genomes":1,"uploaded_genomes":1})
        self.assertIn("-a",self.command);self.assertIn("-f",self.command)
        self.assertEqual(self.command[self.command.index("-n")+1],"16")
        html=s3.objects["results",code+"/gtotree-report.html"].decode()
        self.assertIn("Request a dashboard for my lab",html);self.assertIn(code,html)
        self.assertIn('"gene_set":"Bacteria"',html);self.assertIn("Example bacterium",html)
        self.assertIn("output/selected-genomes.json",manifest["files"])
        self.assertEqual(s3.events[-1],code+"/status.json")
        self.assertLess(s3.events.index(code+"/gtotree-report.html"),len(s3.events)-1)
        self.assertEqual(s3.objects["inputs",code+"/accessions.txt"],b"gcf_000000000.1\nGCF_000000000.1\n")
    def test_failure_is_preserved_and_only_selected_retry_runs(self):
        s3,code=self.job(); self.args.clean=True
        with patch.object(worker,"run_command",side_effect=RuntimeError("simulated failure")):
            with self.assertRaises(RuntimeError):worker.process_job(self.args,s3,code)
        self.assertTrue((self.root/code/"run.log").exists())
        with patch.object(worker,"run_command") as run:
            worker.process_job(self.args,s3,code);run.assert_not_called()
        self.args.retry_failed=True
        with patch.object(worker,"run_command",side_effect=self.fake_command):worker.process_job(self.args,s3,code)
        self.assertEqual(json.loads(s3.objects["results",code+"/status.json"])["state"],"complete")
        self.assertFalse((self.root/code).exists())
    def test_reference_limit_and_local_normalization(self):
        p=self.root/"accessions.txt";p.write_text("gcf_000000000.1\nGCF_000000000.1\n")
        self.assertEqual(read_accessions(p),["GCF_000000000.1"])
        p.write_text("\n".join(row["assembly_accession"] for row in metadata(101)["genomes"]))
        with self.assertRaisesRegex(ValueError,"100"):read_accessions(p)
    def test_metadata_cannot_claim_unsubmitted_genomes_and_report_escapes_text(self):
        p=self.root/"meta.json";d=metadata();d["genomes"][0]["organism_name"]="</script><img src=x onerror=alert(1)>";p.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"absent"):read_metadata(p,set())
        parsed=read_metadata(p,{"GCF_000000000.1"})
        html=report.build_html("test","Bacteria","(A,B);",None,None,None,parsed,"Example123")
        self.assertNotIn('</script><img',html);self.assertIn('<\\/script>',html)
        with self.assertRaises(ValueError):report.build_html("t","B","(A,B);",None,None,None,contact_email="x\nBcc:evil@example.com")
    def test_input_paths_and_s3_objects_cannot_escape_work_directory(self):
        s3=FakeS3();s3.objects["inputs","job/../../escape"]=b"secret"
        with self.assertRaises(ValueError):worker.download_prefix(s3,"inputs","job/",self.root/"input")
        s3=FakeS3();s3.objects["inputs","job/big.fna"]=b"12345"
        with self.assertRaisesRegex(ValueError,"upload-size"):worker.download_prefix(s3,"inputs","job/",self.root/"input",max_bytes=4)
        m=worker.JobManifest(slug="job",gene_set="Bacteria",genomes=["../secret.fna"])
        with self.assertRaises(ValueError):validate_manifest(m,self.root,self.args)
    def test_resource_product_fits_budget(self):
        flags=bound_resource_flags(["-t","-n","999","-M","16","-j","100"],16)
        n=int(flags[flags.index("-n")+1]);m=int(flags[flags.index("-M")+1]);j=int(flags[flags.index("-j")+1])
        self.assertLessEqual(max(n,m)*j,16)
        with self.assertRaises(ValueError):bound_resource_flags(["-j","0"],16)
    def test_lock_survives_stale_pid_file_and_excludes_second_worker(self):
        p=self.root/"worker.lock";p.write_text("stale-pid")
        first=worker.acquire_lock(p)
        try:
            with self.assertRaises(FileExistsError):worker.acquire_lock(p)
        finally:worker.release_lock(first)
        second=worker.acquire_lock(p);worker.release_lock(second)
    def test_pattern_preference_uses_mod_names_first(self):
        (self.root/"a.tre").write_text("(A,B);");(self.root/"z_mod_names.tre").write_text("(namedA,namedB);")
        self.assertEqual(worker.collect_outputs(self.root)["tree"].name,"z_mod_names.tre")
    def test_untargeted_force_is_rejected(self):
        with patch.object(sys,"argv",["worker","--force"]),self.assertRaises(SystemExit):worker.parse_args()


if __name__=="__main__":unittest.main()
