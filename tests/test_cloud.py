from pathlib import Path
import subprocess,sys,unittest
from unittest.mock import Mock,patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from cloud_runner import AuditError,run_audit,CALCULATION_GATE
from streamlit.testing.v1 import AppTest
ROOT=Path(__file__).resolve().parents[1]

class Cloud(unittest.TestCase):
    def test_real_manual_and_generated_names(self):
        for candidate in ['', 'ethanol', 'methanol']:
            result=run_audit('CCO',candidate)
            self.assertEqual(result['pin_status'],'unknown')
            if candidate:self.assertEqual(result['name'],candidate)
            if candidate=='methanol':self.assertEqual(result['overall'],'fail')
    def test_invalid_input(self):
        for smiles,candidate in [('', ''),('C'*3001,''),('CCO','x'*2001),('not a smiles','')]:
            with self.assertRaises(AuditError):run_audit(smiles,candidate)
    def test_deadline_cleans_entire_process_group(self):
        proc=Mock(pid=123);proc.communicate.side_effect=[subprocess.TimeoutExpired('worker',120),('',None)]
        with patch('cloud_runner.subprocess.Popen',return_value=proc),patch('cloud_runner.os.killpg') as kill:
            with self.assertRaisesRegex(AuditError,'120秒'):run_audit('CCO')
            kill.assert_called_once();self.assertEqual(kill.call_args.args[0],123)
    def test_bad_worker_reply(self):
        proc=Mock(returncode=0);proc.communicate.return_value=('not-json',None)
        with patch('cloud_runner.subprocess.Popen',return_value=proc):
            with self.assertRaises(AuditError):run_audit('CCO')
    def test_streamlit_real_result_and_invalidation(self):
        at=AppTest.from_file(str(ROOT/'streamlit_app.py'),default_timeout=40).run()
        self.assertFalse(at.exception)
        at.button[0].click().run();self.assertFalse(at.exception)
        self.assertEqual(at.session_state['report']['pin_status'],'unknown')
        self.assertTrue(at.get('download_button'))
        at.text_area[0].set_value('CC').run();self.assertFalse(at.exception)
        with self.assertRaises(KeyError):at.session_state['report']
        at.selectbox[0].select('α-ピネン型（立体未指定）').run()
        self.assertEqual(at.text_area[0].value,'CC1=CCC2CC1C2(C)C')
    def test_streamlit_busy_does_not_start_second_job(self):
        at=AppTest.from_file(str(ROOT/'streamlit_app.py'),default_timeout=40).run()
        gate=CALCULATION_GATE;gate.acquire()
        try:
            at.button[0].click().run();self.assertFalse(at.exception)
            self.assertTrue(any('別の検証' in x.value for x in at.warning))
        finally:gate.release()

if __name__=='__main__':unittest.main()
