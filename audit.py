"""Offline molecular naming audit. Rule coverage is deliberately narrow and explicit."""
import base64
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime, timezone
from rdkit import Chem, rdBase
from extended_rules import extended_rule_checks
from glycoside_analysis import analyze_glycoside
from carbohydrate_rules import analyze_carbohydrates
from flavonoid_rules import analyze_flavonoids, flavonoid_checks
from monocycle_rules import analyze_monocycles, monocycle_checks
from fused_rules import analyze_fused, fused_checks
from bridged_rules import analyze_bridged, bridged_checks
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

ROOT = Path(__file__).resolve().parent
JAR = ROOT / 'vendor' / 'opsin-cli-2.9.0.jar'
APP_VERSION = '0.30.0'
RULESET = 'bounded-organic-rules-0.30'
P1 = 'https://iupac.qmul.ac.uk/BlueBook/P1.html'
P2 = 'https://iupac.qmul.ac.uk/BlueBook/P2.html'
P6 = 'https://iupac.qmul.ac.uk/BlueBook/P6.html'
ALKANES = ('', 'methane', 'ethane', 'propane', 'butane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane', 'decane')
RULES = [
    ('parent', '親鎖の選択・炭素数', 'P-21.1.1 / P-21.2.1', P2),
    ('suffix', '主特性基と接尾辞', 'P-63.1.2 / P-21.2.1', P6),
    ('numbering', 'OH基の位置番号', 'P-14.4 / P-15.1.7.2.1', P1),
    ('spelling', '名前の組み立て・位置番号の表記', 'P-14.3.4.2 / P-15.1.7.2.1', P1),
    ('stereo', '立体・同位体・その他の命名規則', 'P-12.1（PIN全体の判定は対象外）', P1),
]


def check(key, title, status, detail, reference='', url=''):
    return dict(id=key, title=title, status=status, detail=detail, reference=reference, url=url)


def parse_smiles(text):
    if not isinstance(text, str) or not text.strip():
        raise ValueError('SMILESを入力してください。')
    if len(text) > 3000:
        raise ValueError('初版はSMILES 3,000文字以内です。')
    params = Chem.SmilesParserParams()
    params.parseName = False
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(text.strip(), params)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError('SMILESを解釈できません。構文・原子価を確認してください。')
    if mol.GetNumAtoms() > 200:
        raise ValueError('初版は200原子までです。')
    return mol


def normalized(mol, stereo=True):
    # Ignore only atom maps and removable ordinary explicit H. Never neutralize,
    # merge tautomers, drop fragments, isotopes, radicals, or stereo.
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    copy = Chem.RemoveHs(copy)
    if not stereo:
        Chem.RemoveStereochemistry(copy)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def picture(mol, highlights=None):
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 260)
    drawer.drawOptions().addAtomIndices = True
    if highlights:
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=sorted(highlights), highlightAtomColors=highlights)
    else:
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def structure_checks(original, restored):
    uncertain = bool(original.GetStereoGroups() or restored.GetStereoGroups())
    uncertain |= any(b.GetStereo() == Chem.BondStereo.STEREOANY for m in (original, restored) for b in m.GetBonds())
    # CXSMILES groups and coordinate/dative bonding need additional comparison rules.
    uncertain |= any(b.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC) for m in (original, restored) for b in m.GetBonds())
    topology = normalized(original, False) == normalized(restored, False)
    exact = normalized(original) == normalized(restored)
    return [
        check('constitution', '結合・電荷・同位体・成分', 'pass' if topology else 'fail',
              '原子マップを除き、通常の明示水素・芳香族表現を正規化して比較。互変異性体や塩は同一化しません。'),
        check('roundtrip', '立体化学を含む往復照合', 'unknown' if uncertain else ('pass' if exact else 'fail'),
              '拡張立体情報・未確定結合立体・特殊結合は比較保証の対象外です。' if uncertain else
              'RDKitのcanonical isomeric SMILESを照合。未指定の立体を勝手に補完しません。'),
    ]


def linear_scope(mol):
    mol = Chem.RemoveHs(Chem.Mol(mol))
    atoms = list(mol.GetAtoms())
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetRingInfo().NumRings():
        return None
    if mol.GetStereoGroups() or any(a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in atoms):
        return None
    if any(a.GetAtomicNum() not in (6, 8) for a in atoms):
        return None
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in mol.GetBonds()):
        return None
    carbons = [a for a in atoms if a.GetAtomicNum() == 6]
    oxygen = [a for a in atoms if a.GetAtomicNum() == 8]
    n = len(carbons)
    if not 1 <= n <= 10 or len(oxygen) > 1:
        return None
    adjacent = {a.GetIdx(): [b.GetIdx() for b in a.GetNeighbors() if b.GetAtomicNum() == 6] for a in carbons}
    if any(len(v) > 2 for v in adjacent.values()):
        return None
    if n > 1 and sum(map(len, adjacent.values())) != 2 * (n - 1):
        return None
    ends = [i for i, ns in adjacent.items() if len(ns) <= 1]
    if not ends:
        return None
    path, previous, current = [], None, ends[0]
    while current is not None:
        path.append(current)
        next_nodes = [i for i in adjacent[current] if i != previous]
        previous, current = current, next_nodes[0] if next_nodes else None
    if len(path) != n:
        return None
    if not oxygen:
        return dict(n=n, kind='alkane', locant=None, expected=ALKANES[n], path=path)
    oh = oxygen[0]
    if oh.GetDegree() != 1 or oh.GetTotalNumHs() != 1 or oh.GetNeighbors()[0].GetAtomicNum() != 6:
        return None
    p = path.index(oh.GetNeighbors()[0].GetIdx()) + 1
    locant = min(p, n + 1 - p)
    stem = ALKANES[n][:-1]
    expected = stem + ('ol' if n <= 2 else f'-{locant}-ol')
    return dict(n=n, kind='alcohol', locant=locant, expected=expected, path=path,
                alternatives=[p, n + 1 - p])


def parsed_simple_name(name):
    name = name.strip().lower()
    if name in ALKANES[1:]:
        return dict(n=ALKANES.index(name), kind='alkane', locant=None)
    for n, alkane in enumerate(ALKANES[1:], 1):
        stem = alkane[:-1]
        if name == stem + 'ol':
            return dict(n=n, kind='alcohol', locant=1 if n <= 2 else None)
        match = re.fullmatch(re.escape(stem) + r'-(\d+)-ol', name)
        if match and 1 <= int(match[1]) <= n:
            return dict(n=n, kind='alcohol', locant=int(match[1]))
    return None


def rule_checks(mol, name):
    extended = extended_rule_checks(mol, name)
    if extended is not None:
        return extended
    scope = linear_scope(mol)
    parsed = parsed_simple_name(name)
    reason = ('今回の規則検証の対象外です。指定立体・同位体・塩・ラジカル、複雑な置換基、'
              '実装上限を超える不飽和、未対応の官能基の組合せ、複素環・多環などは未実装です。対応範囲表を確認してください。')
    if scope is None or parsed is None:
        if scope and parsed is None:
            reason = '構造は対応範囲内ですが、この名前の構文・別名は規則検証器が未対応です。'
        return dict(scope=scope, checks=[check(k, t, 'unknown', reason, r, u) for k, t, r, u in RULES])
    rows = []
    n_ok = scope['n'] == parsed['n']
    kind_ok = scope['kind'] == parsed['kind']
    k, t, r, u = RULES[0]
    rows.append(check(k, t, 'pass' if n_ok else 'fail',
                      f"炭素鎖は無分岐で{scope['n']}原子。名称の親鎖は{parsed['n']}炭素。鎖の原子index: {scope['path']}。", r, u))
    k, t, r, u = RULES[1]
    rows.append(check(k, t, 'pass' if kind_ok else 'fail',
                      'OHが1個なら-ol、OHがなければ-ane。競合する官能基はこの対応範囲にはありません。', r, u))
    k, t, r, u = RULES[2]
    if scope['kind'] == parsed['kind'] == 'alcohol' and n_ok:
        rows.append(check(k, t, 'unknown' if parsed['locant'] is None else ('pass' if scope['locant'] == parsed['locant'] else 'fail'),
                          f"両端からのOH位置は{scope['alternatives']}。採用する最小位置番号は{scope['locant']}。" +
                          ('名称に位置番号がないため、番号の選択は未判定。省略の可否は表記の項目で検査します。' if parsed['locant'] is None else f"名称の位置番号は{parsed['locant']}。"), r, u))
    else:
        rows.append(check(k, t, 'na' if scope['kind'] == parsed['kind'] == 'alkane' else 'unknown',
                          '無置換アルカンにOH位置番号は不要です。' if scope['kind'] == parsed['kind'] == 'alkane' else '親鎖または接尾辞の不一致により、この番号付けは判定できません。', r, u))
    k, t, r, u = RULES[3]
    exact = name.strip().lower() == scope['expected']
    rows.append(check(k, t, 'pass' if exact else 'fail',
                      f"限定した規則セットによる期待表記：{scope['expected']}。methanol/ethanolは位置番号を省略。", r, u))
    k, t, r, u = RULES[4]
    rows.append(check(k, t, 'unknown', '指定立体・同位体は対象外。PIN全体の適合性は未判定です。', r, u))
    return dict(scope=scope, checks=rows)


def generate_name(smiles):
    proc = subprocess.run([sys.executable, str(ROOT / 'naming_worker.py')],
                          input=json.dumps({'smiles': smiles}), text=True, capture_output=True, timeout=25, cwd=ROOT)
    if proc.returncode:
        raise RuntimeError('NISPO命名プロセスが失敗しました。')
    result=json.loads(proc.stdout)
    if result.get('error'):raise RuntimeError('NISPO命名エラー：'+str(result['error']))
    return result['name']


def restore_name(name):
    if not JAR.exists():
        raise RuntimeError('OPSINのJARがありません。READMEのセットアップ手順を確認してください。')
    proc = subprocess.run(['java', '-Xmx256m', '-cp', str(JAR), str(ROOT / 'OpsinBridge.java')],
                          input=name + '\n', text=True, capture_output=True, timeout=25, cwd=ROOT)
    if proc.returncode:
        raise RuntimeError('OPSINを起動できません。Java 11以降のJDKを確認してください。')
    lines = proc.stdout.splitlines()
    if len(lines) < 3 or lines[0] not in ('SUCCESS', 'WARNING', 'FAILURE'):
        raise RuntimeError('OPSIN応答を解釈できません。')
    return dict(status=lines[0], smiles=base64.b64decode(lines[1]).decode(), message=base64.b64decode(lines[2]).decode())


def expected_name_check(mol, expected, candidate, parsed_candidate):
    """A suggested reference spelling must also reconstruct the input structure."""
    if not expected:
        return None, None
    try:
        restored = parsed_candidate if expected == candidate else restore_name(expected)
        if not restored or restored['status'] != 'SUCCESS' or not restored['smiles']:
            return check('expected_roundtrip', '期待表記の構造照合', 'unknown',
                         '期待表記の逆変換に警告・失敗があるか、処理を完了できませんでした。'), restored
        result = structure_checks(mol, parse_smiles(restored['smiles']))[-1]
        return check('expected_roundtrip', '期待表記の構造照合', result['status'],
                     '規則から導いた期待表記もOPSINで構造へ戻し、入力と照合しました。' +
                     ('期待表記が構造に一致しません。検証器の結果を要確認。' if result['status']=='fail' else '')), restored
    except (RuntimeError, ValueError, subprocess.TimeoutExpired):
        return check('expected_roundtrip', '期待表記の構造照合', 'unknown',
                     '期待表記の構造照合を完了できませんでした。'), None


def audit(smiles, candidate=''):
    mol = parse_smiles(smiles)
    original = normalized(mol)
    candidate = candidate.strip()
    if len(candidate) > 2000 or any(ord(c) < 32 for c in candidate):
        raise ValueError('名前は制御文字を含まない2,000文字以内にしてください。')
    data = dict(schema_version=1, app_version=APP_VERSION, ruleset=RULESET,
                created_at=datetime.now(timezone.utc).isoformat(), input_smiles=smiles.strip(),
                canonical_input=original, formula=rdMolDescriptors.CalcMolFormula(mol),
                name_source='manual' if candidate else 'NISPO', input_svg=picture(mol),
                versions={'Python': sys.version.split()[0], 'verifier_SHA256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'extended_rules_SHA256': hashlib.sha256((ROOT / 'extended_rules.py').read_bytes()).hexdigest(), 'carbon_rules_SHA256': hashlib.sha256((ROOT / 'carbon_rules.py').read_bytes()).hexdigest(), 'natural_rules_SHA256': hashlib.sha256((ROOT / 'natural_rules.py').read_bytes()).hexdigest(), 'linkage_rules_SHA256': hashlib.sha256((ROOT / 'linkage_rules.py').read_bytes()).hexdigest(), 'aromatic_rules_SHA256': hashlib.sha256((ROOT / 'aromatic_rules.py').read_bytes()).hexdigest(), 'functional_tree_rules_SHA256': hashlib.sha256((ROOT / 'functional_tree_rules.py').read_bytes()).hexdigest(), 'oxygen_ring_rules_SHA256': hashlib.sha256((ROOT / 'oxygen_ring_rules.py').read_bytes()).hexdigest(), 'glycoside_analysis_SHA256': hashlib.sha256((ROOT / 'glycoside_analysis.py').read_bytes()).hexdigest(), 'ring_stereo_rules_SHA256': hashlib.sha256((ROOT / 'ring_stereo_rules.py').read_bytes()).hexdigest(), 'enone_rules_SHA256': hashlib.sha256((ROOT / 'enone_rules.py').read_bytes()).hexdigest(), 'fused_rules_SHA256': hashlib.sha256((ROOT / 'fused_rules.py').read_bytes()).hexdigest(), 'bridged_rules_SHA256': hashlib.sha256((ROOT / 'bridged_rules.py').read_bytes()).hexdigest(), 'monocycle_rules_SHA256': hashlib.sha256((ROOT / 'monocycle_rules.py').read_bytes()).hexdigest(), 'flavonoid_rules_SHA256': hashlib.sha256((ROOT / 'flavonoid_rules.py').read_bytes()).hexdigest(), 'carbohydrate_rules_SHA256': hashlib.sha256((ROOT / 'carbohydrate_rules.py').read_bytes()).hexdigest(), 'OPSIN_bridge_SHA256': hashlib.sha256((ROOT / 'OpsinBridge.java').read_bytes()).hexdigest(), 'NISPO': importlib.metadata.version('nispo'), 'RDKit': rdBase.rdkitVersion,
                          'OPSIN': '2.9.0', 'OPSIN_SHA256': hashlib.sha256(JAR.read_bytes()).hexdigest() if JAR.exists() else None},
                comparison_policy='Atom maps ignored; ordinary H/aromatic representations normalized; no tautomer, protonation, fragment or stereo collapsing.',
                disclaimer='研究用の開発版。構造照合と限定規則の結果です。PIN全体を保証しません。')
    topology=analyze_glycoside(mol)
    data['glycoside_analysis']=topology
    data['carbohydrate_analysis']=analyze_carbohydrates(mol,topology)
    data['flavonoid_analysis']=analyze_flavonoids(mol)
    data['monocycle_analysis']=analyze_monocycles(mol)
    data['fused_analysis']=analyze_fused(mol)
    data['bridged_analysis']=analyze_bridged(mol)
    if topology['rings']:
        highlights={idx:(0.65,0.83,0.70) for ring in topology['rings'] for idx in ring['ring_atoms']}
        for link in topology['o_linkages']:
            for idx in link['bond_atoms']:highlights[idx]=(1.0,0.76,0.40)
        for link in topology['c_linkages']:
            for idx in link['bond_atoms']:highlights[idx]=(0.66,0.70,1.0)
        data['glycoside_svg']=picture(mol,highlights)
    naming_issue=None;prechecked_opsin=None
    data['generation']=dict(engine='NISPO',status='not_requested' if candidate else 'success')
    try:
        name = candidate or generate_name(original)
        if not name:
            data['generation'].update(status='empty')
            naming_issue=check('naming','NISPOによる命名','unknown','NISPOはこの構造の名前を生成できませんでした。')
    except (RuntimeError, subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
        name=''
        data['generation'].update(status='timeout' if isinstance(exc,subprocess.TimeoutExpired) else 'failed',error_type=type(exc).__name__,detail=str(exc))
        naming_issue=check('naming','NISPOによる命名','unknown','NISPO命名が完了しませんでした：'+type(exc).__name__+'。限定規則の期待名は独立して構造照合します。')
    if not name:
        reference_kind=next((kind for kind in ('carbohydrate','flavonoid','monocycle','fused','bridged') if data[kind+'_analysis']['expected_name']),'flavonoid')
        reference_label={'carbohydrate':'期待糖類名','flavonoid':'期待フラボノイド名','monocycle':'期待単環骨格名','fused':'期待縮合二環名','bridged':'期待架橋二環名'}[reference_kind]
        expected=data[reference_kind+'_analysis']['expected_name']
        reference,parsed=expected_name_check(mol,expected,'',None)
        if reference:
            reference.update(id=reference_kind+'_expected_roundtrip',title=reference_label+'の構造照合')
            data[reference_kind+'_expected_opsin']=parsed
            data[reference_kind+'_expected_roundtrip_status']=reference['status']
        if reference and reference['status']=='pass':
            # Promote only a reference whose complete input stereostructure was
            # independently restored. Do not hide the original naming failure.
            name=expected;prechecked_opsin=parsed
            data['name_source']='bounded_rules'
            data['fallback_reason']='NISPO naming unavailable; bounded reference name independently restored and matched.'
        else:
            rows=[naming_issue] if naming_issue else []
            if reference:rows.append(reference)
            data.update(name='',checks=rows,roundtrip_status='unknown',overall='fail' if any(c['status']=='fail' for c in rows) else 'unknown',pin_status='unknown')
            return data
    data['name']=name
    rules = rule_checks(mol, name)
    data['rule_scope'] = rules['scope']
    rows = [naming_issue] if naming_issue else []
    try:
        opsin = prechecked_opsin or restore_name(name)
        data['opsin'] = opsin
        restored = parse_smiles(opsin['smiles']) if opsin['smiles'] else None
        rows.append(check('opsin', 'OPSINによる名前の解釈', 'pass' if opsin['status'] == 'SUCCESS' else 'unknown',
                          opsin['message'] or opsin['status']))
        if restored is not None:
            data.update(canonical_restored=normalized(restored), restored_svg=picture(restored))
            matches = structure_checks(mol, restored)
            # A warning must never lead to an apparently clean round-trip result.
            if opsin['status'] != 'SUCCESS':
                for row in matches:
                    if row['status'] == 'pass':
                        row.update(status='unknown', detail='OPSINの警告・失敗があるため照合の合格を保留。' + row['detail'])
            rows.extend(matches)
        else:
            rows.append(check('roundtrip', '構造の往復照合', 'unknown', '名前を構造に戻せないため、検証不能です。誤った名前と確定したわけではありません。'))
    except (RuntimeError, ValueError, subprocess.TimeoutExpired) as exc:
        rows.append(check('roundtrip', '構造の往復照合', 'unknown', '照合を完了できません：' + str(exc)))
    expected = (rules['scope'] or {}).get('expected')
    reference_check, reference_opsin = expected_name_check(mol, expected, name, data.get('opsin'))
    carbohydrate_expected=data['carbohydrate_analysis']['expected_name']
    if carbohydrate_expected:
        carbohydrate_check,carbohydrate_opsin=expected_name_check(mol,carbohydrate_expected,name,data.get('opsin')) if carbohydrate_expected!=expected else (dict(reference_check) if reference_check else None,reference_opsin)
        data['carbohydrate_expected_opsin']=carbohydrate_opsin
        data['carbohydrate_expected_roundtrip_status']=carbohydrate_check['status'] if carbohydrate_check else 'unknown'
        if carbohydrate_check:
            carbohydrate_check.update(id='carbohydrate_expected_roundtrip',title='期待糖類名の構造照合')
            rows.append(carbohydrate_check)
    flavonoid_expected=data['flavonoid_analysis']['expected_name']
    if flavonoid_expected:
        flav_check,flav_opsin=expected_name_check(mol,flavonoid_expected,name,data.get('opsin'))
        data['flavonoid_expected_opsin']=flav_opsin
        data['flavonoid_expected_roundtrip_status']=flav_check['status']
        flav_check.update(id='flavonoid_expected_roundtrip',title='期待フラボノイド名の構造照合')
        rows.append(flav_check)
    rows.extend(flavonoid_checks(data['flavonoid_analysis'],name))
    monocycle_expected=data['monocycle_analysis']['expected_name']
    if monocycle_expected:
        mono_check,mono_opsin=expected_name_check(mol,monocycle_expected,name,data.get('opsin'))
        data['monocycle_expected_opsin']=mono_opsin
        data['monocycle_expected_roundtrip_status']=mono_check['status']
        mono_check.update(id='monocycle_expected_roundtrip',title='期待単環骨格名の構造照合')
        rows.append(mono_check)
    rows.extend(monocycle_checks(data['monocycle_analysis'],name))
    fused_expected=data['fused_analysis']['expected_name']
    if fused_expected:
        fused_check,fused_opsin=expected_name_check(mol,fused_expected,name,data.get('opsin'))
        data['fused_expected_opsin']=fused_opsin
        data['fused_expected_roundtrip_status']=fused_check['status']
        fused_check.update(id='fused_expected_roundtrip',title='期待縮合二環名の構造照合')
        rows.append(fused_check)
    rows.extend(fused_checks(data['fused_analysis'],name))
    bridged_expected=data['bridged_analysis']['expected_name']
    if bridged_expected:
        bridged_check,bridged_opsin=expected_name_check(mol,bridged_expected,name,data.get('opsin'))
        data['bridged_expected_opsin']=bridged_opsin
        data['bridged_expected_roundtrip_status']=bridged_check['status']
        bridged_check.update(id='bridged_expected_roundtrip',title='期待架橋二環名の構造照合')
        rows.append(bridged_check)
    rows.extend(bridged_checks(data['bridged_analysis'],name))
    if reference_check:
        rows.append(reference_check)
        data['expected_roundtrip_status'] = reference_check['status']
        data['expected_opsin'] = reference_opsin
    rows.extend(rules['checks'])
    data['checks'] = rows
    data['overall'] = 'fail' if any(r['status'] == 'fail' for r in rows) else 'unknown'
    data['pin_status'] = 'unknown'
    data['roundtrip_status'] = next((r['status'] for r in rows if r['id'] == 'roundtrip'), 'unknown')
    return data
