"""Isolated naming process: permits a hard timeout without stopping the server."""
import json
import sys
from nispo import smiles_to_iupac

if __name__ == '__main__':
    try:
        name = smiles_to_iupac(json.load(sys.stdin)['smiles'])
        print(json.dumps({'name': name or ''}))
    except Exception as exc:
        print(json.dumps({'name': '', 'error': type(exc).__name__}))
