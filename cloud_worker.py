"""One request per process. Inputs and errors are never written to application logs."""
import json
import sys

def main():
    from audit import audit
    try:
        data = json.load(sys.stdin)
        result = audit(data['smiles'], data.get('candidate', ''))
        reply = {'result': result}
    except ValueError as exc:
        reply = {'error': str(exc)}
    except Exception:
        reply = {'error': '計算を完了できませんでした。入力や実行環境を確認してください。'}
    print(json.dumps(reply, ensure_ascii=False))

if __name__ == '__main__':
    main()
