"""Bounded absolute R/S correspondence for fully specified oxygen-ring polyols.
CIP assignment is delegated to RDKit rdCIPLabeler, not an independent CIP proof.
Skeleton numbering uses the existing finite evaluator. Symmetry ties needing
stereochemical numbering, pseudoasymmetry, partial/enhanced stereo stay unknown.
"""
import re
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from extended_rules import eligible,row,unknown_result
from oxygen_ring_rules import oxygen_ring_checks,parse_name

NOTE='環炭素の絶対R/SをRDKit CIP計算と原子対応で限定照合します。CIP規則自体の独立証明、保持名の優先性、PIN全体は未判定です。糖類番号・糖種・D/L・α/βは別の糖類解析欄で限定照合します。'


def cip_inventory(mol):
    copy=Chem.Mol(mol)
    Chem.AssignStereochemistry(copy,cleanIt=True,force=True)
    for atom in copy.GetAtoms():
        if atom.HasProp('_CIPCode'):atom.ClearProp('_CIPCode')
    try:rdCIPLabeler.AssignCIPLabels(copy)
    except (RuntimeError,ValueError):return dict(status='unknown',centers=[],reason='CIP計算を完了できませんでした。')
    centers=[];other=False
    for info in Chem.FindPotentialStereo(copy):
        if info.type!=Chem.StereoType.Atom_Tetrahedral:other=True;continue
        atom=copy.GetAtomWithIdx(info.centeredOn)
        centers.append(dict(atom_index=info.centeredOn,cip=atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None,assigned=info.specified==Chem.StereoSpecified.Specified))
    return dict(status='calculated',centers=sorted(centers,key=lambda c:c['atom_index']),other_stereo=other)


def split_name(name):
    text=name.strip()
    if not text.startswith('('):return [],text,True
    match=re.fullmatch(r'\((\d+[RS](?:,\d+[RS])*)\)-(.+)',text)
    if not match:return [],text,False
    return [(int(token[:-1]),token[-1]) for token in match[1].split(',')],match[2],True


def ring_stereo_checks(mol,name):
    if not any(a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):return None
    base=Chem.Mol(mol);Chem.RemoveStereochemistry(base)
    if not eligible(base):return None
    skeleton=oxygen_ring_checks(base,'unrecognized')
    if skeleton is None:return None
    scope=dict(skeleton['scope'],kind='oxygen_ring_stereo',scope_note=NOTE)
    scope.pop('expected',None)
    info=cip_inventory(mol)
    if not skeleton['scope'].get('expected'):return unknown_result(scope,'骨格の置換・親骨格の選択が限定範囲外のため、R/Sと親骨格番号の対応を保留します。')
    if mol.GetStereoGroups() or info['status']!='calculated' or info.get('other_stereo'):
        return unknown_result(scope,'拡張立体・非四面体・二重結合立体・CIP計算不能は未対応です。')
    centers=info['centers'];path=skeleton['scope']['path']
    if not centers or any(not c['assigned'] or c['cip'] not in ('R','S') or c['atom_index'] not in path for c in centers):
        return unknown_result(scope,'環上の全ての四面体立体中心がR/Sで指定された場合のみ対応します。部分指定、環外立体、擬不斉r/sは保留します。')
    codes={c['atom_index']:c['cip'] for c in centers}
    def mapped(path):return [(i,codes[idx]) for i,idx in enumerate(path,1) if idx in codes]
    alternatives={tuple(mapped(c['path'])) for c in skeleton['scope']['preferred_candidates']}
    if len(alternatives)!=1:return unknown_result(scope,'骨格規則だけでは環の向きが同順位です。立体配置を使う番号選択は未実装のため期待名を保留します。')
    descriptors=mapped(path);prefix='('+','.join(str(i)+code for i,code in descriptors)+')-'
    expected=prefix+skeleton['scope']['expected']
    scope.update(expected=expected,stereo_numbering=[dict(locant=i,atom_index=path[i-1],cip=code) for i,code in descriptors],cip_engine='RDKit rdCIPLabeler')
    parsed,body,syntax=split_name(name)
    result=oxygen_ring_checks(base,body)
    checks=result['checks']
    parsed_body=parse_name(body) if syntax else None
    mappings=[]
    if parsed_body and parsed_body['parent']==scope['parent_name']:
        mappings=[mapped(c['path']) for c in scope['candidates'] if c['locants']==parsed_body['locants'] and c['groups']==parsed_body['groups']]
    complete=bool(mappings) and syntax
    same=any(sorted(parsed)==mapping for mapping in mappings)
    checks=[r for r in checks if r['id']!='stereo']
    for item in checks:
        if item['id']=='spelling':
            item.update(status='pass' if name.strip()==expected else ('different' if item['status']=='different' and same and parsed==sorted(parsed) and len({i for i,_ in parsed})==len(parsed) else 'fail'),detail='R/Sを含む限定期待表記：'+expected+'。糖類保持名・PINの優先性は未検証です。')
        if item['id']=='sugar_nomenclature':item['detail']=NOTE
    checks.extend([
      row('ring_rs','環炭素のR/Sと原子対応','unknown' if not complete else ('pass' if same else 'fail'),f"名称のR/Sは{parsed}。候補番号ごとの入力R/Sは{mappings}。指定の欠落・反転・余分な中心を区別し、入力にない立体は補完しません。",'P-93.5.1.1 / P-92','P9'),
      row('ring_rs_order','R/Sの列挙順・位置番号の重複','unknown' if not complete else ('pass' if parsed==sorted(parsed) and len({i for i,_ in parsed})==len(parsed) else 'fail'),'R/Sは位置番号昇順で列挙し、同じ位置への二重指定を許しません。','P-93.4.1.2 / P-93.5','P9'),
      row('stereo','その他の立体・同位体・PIN全体','unknown',NOTE,'P-93.5 / P-102','P10')])
    # Unsupported aliases/grammars stay unknown rather than being declared invalid.
    if parsed_body is None:
        for item in checks:
            if item['id'] not in ('stereo','sugar_nomenclature'):item['status']='unknown'
    return dict(scope=scope,checks=checks)
