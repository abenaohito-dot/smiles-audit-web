"""Bound each cloud audit and clean up its entire subprocess group on timeout."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading

ROOT = Path(__file__).resolve().parent
TIMEOUT_SECONDS = 120
CALCULATION_GATE = threading.BoundedSemaphore(1)

class AuditError(Exception):
    pass

def run_audit(smiles, candidate=''):
    if not isinstance(smiles, str) or not smiles.strip():
        raise AuditError('SMILESを入力してください。')
    if len(smiles) > 3000 or not isinstance(candidate, str) or len(candidate) > 2000:
        raise AuditError('SMILESは3,000文字、検査する名前は2,000文字以内にしてください。')
    proc = subprocess.Popen(
        [sys.executable, str(ROOT / 'cloud_worker.py')], cwd=ROOT,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        text=True, start_new_session=True,
    )
    try:
        stdout, _ = proc.communicate(json.dumps({'smiles': smiles, 'candidate': candidate}), timeout=TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise AuditError('処理が120秒を超えたため中止しました。小さい構造か、既存の名前の検査をお試しください。') from None
    if proc.returncode:
        raise AuditError('計算を完了できませんでした。しばらく待って再実行してください。')
    try:
        reply = json.loads(stdout)
    except (ValueError, TypeError):
        raise AuditError('計算結果を読み取れませんでした。再実行してください。') from None
    if 'error' in reply:
        raise AuditError(reply['error'])
    return reply['result']
