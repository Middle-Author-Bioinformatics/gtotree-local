import json
import fcntl
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_workflow as fixtures
import gtotree_worker as worker
from public_results import is_public_file, public_text, PUBLIC_ERROR
from cleanup_public_results import plan_job, apply_job
import cleanup_public_results as cleanup


class PublicResultTests(unittest.TestCase):
    def test_apply_refuses_the_running_workers_outer_lock(self):
        with tempfile.TemporaryDirectory() as td:
            lock_path=Path(td)/'cron.lock'
            with lock_path.open('a+') as held:
                fcntl.flock(held,fcntl.LOCK_EX|fcntl.LOCK_NB)
                with patch.object(sys,'argv',['cleanup','--apply','--lock-file',str(lock_path)]),patch.object(cleanup.boto3,'client') as client:
                    with self.assertRaises(SystemExit):cleanup.main()
                    client.assert_not_called()

    def test_absolute_paths_are_removed_without_changing_branches_sequences_or_citations(self):
        self.assertEqual(public_text("('/home/ark/private/mine.fna':0.123,B:0.456);"),"('mine.fna':0.123,B:0.456);")
        self.assertEqual(public_text(r'C:\Users\ark\private\mine.fna'),'mine.fna')
        self.assertEqual(public_text('file:///home/ark/private/data.tsv'),'data.tsv')
        for value in ['https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001.1/','AA/BB','0/1','</script>']:
            self.assertEqual(public_text(value),value)
        for name in ['run.log','output/gtotree-runlog.txt','raw-results.tar.gz','output/command.txt','result.json']:
            self.assertFalse(is_public_file(name),name)

    def test_worker_publishes_only_sanitized_scientific_outputs(self):
        f=fixtures.WorkflowTests();f.setUp()
        try:
            s3,code=f.job()
            def command(cmd,handle,timeout):
                f.fake_command(cmd,handle,timeout)
                if cmd[0]=='GToTree':
                    out=Path(cmd[cmd.index('-o')+1])
                    (out/'GToTree_output.tre').write_text('(/home/ark/private/mine.fna:0.1,B:0.2);')
                    (out/'Aligned_SCGs.faa').write_text('>/home/ark/private/mine.fna\nMAAA\n>B\nMAAT\n')
                    (out/'Genomes_summary_info.tsv').write_text('genome\tfile\n/home/ark/private/mine.fna\tC:\\Users\\ark\\private\\mine.fna\n')
                    (out/'gtotree-runlog.txt').write_text('/home/ark/private — PRIVATE_CREDENTIAL_TEST')
                    handle.write('/home/ark/private — PRIVATE_CREDENTIAL_TEST\n')
            with patch.object(worker,'run_command',side_effect=command):worker.process_job(f.args,s3,code)
            for (bucket,key),data in s3.objects.items():
                if bucket!='results':continue
                self.assertNotIn(b'/home/ark',data,key)
                self.assertNotIn(b'C:\\Users',data,key)
                self.assertNotIn(b'PRIVATE_CREDENTIAL_TEST',data,key)
                self.assertTrue(key.endswith('/status.json') or key.endswith('/result.json') or is_public_file(key.split('/',1)[1]),key)
            self.assertIn('/home/ark', (f.root/code/'GToTree_output'/'Genomes_summary_info.tsv').read_text())
            self.assertIn('PRIVATE_CREDENTIAL_TEST',(f.root/code/'run.log').read_text())
            self.assertEqual((f.root/code).stat().st_mode&0o777,0o700)
        finally:f.tearDown()

    def test_failure_detail_stays_private(self):
        f=fixtures.WorkflowTests();f.setUp()
        try:
            s3,code=f.job()
            with patch.object(worker,'run_command',side_effect=RuntimeError('/home/ark/private/command failed')):
                with self.assertRaises(RuntimeError):worker.process_job(f.args,s3,code)
            result=json.loads(s3.objects['results',code+'/status.json'])
            self.assertEqual(result['error'],PUBLIC_ERROR)
            self.assertNotIn(('results',code+'/run.log'),s3.objects)
        finally:f.tearDown()

    def legacy(self):
        s3=fixtures.FakeS3();slug='OldExample'
        data={'title':'Old example','gene_set':'Bacteria','newick':'(/home/ark/private/mine.fna:0.1,B:0.2);','alignment':None,'summary':{'columns':['genome','file'],'rows':[['mine',r'C:\Users\ark\private\mine.fna']]},'hits':None}
        objects={'run.log':b'PRIVATE_CREDENTIAL_TEST /home/ark/private',
                 'raw-results.tar.gz':b'private-archive',
                 'output/gtotree-runlog.txt':b'/home/ark/private',
                 'output/GToTree_output.tre':data['newick'].encode(),
                 'gtotree-report.html':('<script>window.GT_DATA = '+json.dumps(data)+';</script>').encode(),
                 'status.json':json.dumps({'state':'complete','has_report':True,'error':'/home/ark/private','files':['run.log','raw-results.tar.gz']}).encode(),
                 'result.json':json.dumps({'run_name':'/home/ark/private/analysis','has_report':True,'command':'/home/ark/command','files':['run.log']}).encode()}
        for name,body in objects.items():s3.objects['results',slug+'/'+name]=body
        return s3,slug,{name:{'Key':slug+'/'+name,'Size':len(body)} for name,body in objects.items()}

    def test_legacy_cleanup_backs_up_originals_and_removes_direct_diagnostic_objects(self):
        s3,slug,objects=self.legacy();originals=dict(s3.objects)
        updates,removals=plan_job(s3,'results',slug,objects)
        self.assertEqual(s3.events,[],'audit must not mutate')
        with tempfile.TemporaryDirectory() as td:
            apply_job(s3,'results',slug,objects,updates,removals,Path(td))
            index=json.loads((Path(td)/slug/'index.json').read_text())
            for key,filename in index.items():self.assertEqual((Path(td)/slug/filename).read_bytes(),originals['results',key])
            self.assertEqual((Path(td)/slug).stat().st_mode&0o777,0o700)
            self.assertEqual((Path(td)/slug/'index.json').stat().st_mode&0o777,0o600)
        for (bucket,key),body in s3.objects.items():
            self.assertNotIn(b'/home/ark',body,key)
            self.assertNotIn(b'PRIVATE_CREDENTIAL_TEST',body,key)
            self.assertNotIn(b'C:\\\\Users',body,key)
            self.assertFalse(key.endswith('run.log') or key.endswith('gtotree-runlog.txt') or key.endswith('tar.gz'))
        report=s3.objects['results',slug+'/gtotree-report.html'].decode()
        self.assertIn('mine.fna',report)
        manifest=json.loads(s3.objects['results',slug+'/result.json'])
        self.assertEqual(s3.events[-1],slug+'/status.json')
        self.assertNotIn('command',manifest)
        self.assertNotIn('run.log',manifest['files'])

    def test_backup_failure_prevents_all_remote_mutations(self):
        s3,slug,objects=self.legacy();originals=dict(s3.objects)
        updates,removals=plan_job(s3,'results',slug,objects)
        with tempfile.TemporaryDirectory() as td,patch.object(s3,'download_file',side_effect=OSError('backup unavailable')):
            with self.assertRaises(OSError):apply_job(s3,'results',slug,objects,updates,removals,Path(td))
        self.assertEqual(s3.objects,originals)
        self.assertEqual(s3.events,[])

    def test_reports_only_refresh_preserves_every_other_object(self):
        s3,slug,objects=self.legacy();originals=dict(s3.objects)
        updates,removals=plan_job(s3,'results',slug,objects,reports_only=True)
        self.assertEqual(set(updates),{'gtotree-report.html'})
        self.assertEqual(removals,[])
        report=updates['gtotree-report.html'].decode()
        self.assertIn('https://midauthorbio.com/#contactdash',report)
        self.assertIn('Final pricing depends on scope',report)
        self.assertIn('aria-label="Report views"',report)
        self.assertNotIn('class="brand"',report)
        self.assertNotIn('brand-logo',report)
        self.assertNotIn('mailto:',report)
        self.assertIn('mine.fna:0.1,B:0.2',report)
        with tempfile.TemporaryDirectory() as td:
            apply_job(s3,'results',slug,objects,updates,removals,Path(td))
            index=json.loads((Path(td)/slug/'index.json').read_text())
            self.assertEqual(set(index),{slug+'/gtotree-report.html'})
            self.assertEqual((Path(td)/slug/index[slug+'/gtotree-report.html']).read_bytes(),originals['results',slug+'/gtotree-report.html'])
        self.assertEqual(set(s3.objects),set(originals))
        for key,body in originals.items():
            if key!=('results',slug+'/gtotree-report.html'):
                self.assertEqual(s3.objects[key],body,key)

if __name__=='__main__':unittest.main()
