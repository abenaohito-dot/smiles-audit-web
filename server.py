"""Loopback-only server. No telemetry, remote calls, database, or input logging."""
import argparse
import json
import hashlib
import importlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import threading
from urllib.parse import urlsplit
import audit as audit_engine
import extended_rules
import carbon_rules
import natural_rules
import linkage_rules
import aromatic_rules
import functional_tree_rules
import oxygen_ring_rules
import ring_stereo_rules
import glycoside_analysis
import carbohydrate_rules
import enone_rules
import flavonoid_rules
import monocycle_rules
import fused_rules
import bridged_rules

ROOT = Path(__file__).resolve().parent
LOCK = threading.BoundedSemaphore(1)
RULE_FILES = ('extended_rules.py', 'carbon_rules.py', 'natural_rules.py', 'linkage_rules.py', 'aromatic_rules.py', 'functional_tree_rules.py', 'oxygen_ring_rules.py', 'glycoside_analysis.py', 'ring_stereo_rules.py', 'carbohydrate_rules.py', 'enone_rules.py', 'flavonoid_rules.py', 'monocycle_rules.py', 'fused_rules.py', 'bridged_rules.py', 'audit.py', 'OpsinBridge.java')


def rule_stamp():
    return tuple(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in RULE_FILES)


ENGINE_STAMP = rule_stamp()


def refresh_rules():
    global ENGINE_STAMP
    stamp = rule_stamp()
    if stamp != ENGINE_STAMP:
        # Calls are serialized by LOCK. Import dependency order matters.
        for module in (extended_rules, carbon_rules, natural_rules, linkage_rules, aromatic_rules, functional_tree_rules, oxygen_ring_rules, ring_stereo_rules, glycoside_analysis, carbohydrate_rules, enone_rules, flavonoid_rules, monocycle_rules, fused_rules, bridged_rules, audit_engine):
            # Avoid timestamp/size-based .pyc reuse after quick edits.
            if module.__cached__:
                Path(module.__cached__).unlink(missing_ok=True)
            importlib.invalidate_caches()
            importlib.reload(module)
        ENGINE_STAMP = stamp
    return ENGINE_STAMP

class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def send(self, status, body, content_type='application/json; charset=utf-8'):
        if isinstance(body, dict):
            body = json.dumps(body, ensure_ascii=False).encode()
        elif isinstance(body, str):
            body = body.encode()
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Referrer-Policy', 'no-referrer')
        self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; frame-ancestors 'none'; base-uri 'none'")
        self.end_headers()
        self.wfile.write(body)

    def allowed_host(self):
        return self.headers.get('Host') in {f'127.0.0.1:{self.server.server_port}', f'localhost:{self.server.server_port}'}

    def do_GET(self):
        if not self.allowed_host():
            return self.send(403, {'error': 'Host is not allowed.'})
        path = urlsplit(self.path).path
        assets = {'/': ('index.html', 'text/html; charset=utf-8'), '/app.js': ('app.js', 'text/javascript; charset=utf-8'), '/style.css': ('style.css', 'text/css; charset=utf-8'), '/coverage.json': ('coverage.json', 'application/json; charset=utf-8')}
        if path == '/health':
            if LOCK.acquire(blocking=False):
                try:
                    refresh_rules()
                except Exception:
                    return self.send(503, {'status':'rules_unavailable'})
                finally:
                    LOCK.release()
            return self.send(200, {'status': 'ok', 'version': audit_engine.APP_VERSION, 'network': 'offline'})
        if path not in assets:
            return self.send(404, {'error': 'Not found'})
        name, mime = assets[path]
        self.send(200, (ROOT / 'static' / name).read_bytes(), mime)

    def do_POST(self):
        if not self.allowed_host():
            return self.send(403, {'error': 'Host is not allowed.'})
        origin = self.headers.get('Origin')
        if origin and origin not in {f'http://127.0.0.1:{self.server.server_port}', f'http://localhost:{self.server.server_port}'}:
            return self.send(403, {'error': 'Cross-origin request rejected.'})
        if self.path != '/api/audit':
            return self.send(404, {'error': 'Not found'})
        if self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
            return self.send(415, {'error': 'JSONが必要です。'})
        try:
            length = int(self.headers.get('Content-Length', '0'))
            if not 0 < length <= 20000:
                return self.send(413, {'error': 'リクエストが大きすぎます。'})
            data = json.loads(self.rfile.read(length))
            if not isinstance(data, dict) or not isinstance(data.get('smiles'), str) or not isinstance(data.get('candidate', ''), str):
                raise ValueError('SMILESと名前は文字列で指定してください。')
        except (ValueError, UnicodeDecodeError):
            return self.send(400, {'error': '入力形式を確認してください。'})
        if not LOCK.acquire(blocking=False):
            return self.send(429, {'error': '別の検証が進行中です。完了してから再実行してください。'})
        try:
            stamp = refresh_rules()
            result = audit_engine.audit(data['smiles'], data.get('candidate', ''))
            if rule_stamp() != stamp:
                return self.send(409, {'error': '検証中に規則が更新されました。新しい規則で再実行してください。'})
            self.send(200, result)
        except ValueError as exc:
            self.send(400, {'error': str(exc)})
        except Exception:
            self.send(500, {'error': '内部エラーで検証が完了しませんでした。入力は保存されていません。'})
        finally:
            LOCK.release()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', type=int, default=8897)
    args = parser.parse_args()
    server = ThreadingHTTPServer(('127.0.0.1', args.port), Handler)
    print(f'SMILES Audit → http://127.0.0.1:{args.port}  (Ctrl+Cで終了)', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()
