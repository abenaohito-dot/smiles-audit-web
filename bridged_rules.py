"""Finite norbornane and pinane graph rules; no general polycycle or PIN claim."""
import re
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from extended_rules import group_rank

SOURCE='https://iupac.qmul.ac.uk/BlueBook/P2.html'
# Each registered graph carries the fixed von Baeyer numbering. Bond orders
# are checked separately; the pinane graph has three SymmSSSR rings but cycle rank 2.
SCAFFOLDS={
    '2.2.1': Chem.MolFromSmarts('[C:1]12~[C:2]~[C:3]~[C:4](~[C:7]~1)~[C:5]~[C:6]~2'),
    '3.1.1': Chem.MolFromSmarts('[C:1]12~[C:2]~[C:3]~[C:4]~[C:5](~[C:7]~1)~[C:6]~2'),
}
QUERY=SCAFFOLDS['2.2.1']
NOTE=('登録したbicyclo[2.2.1]heptane／bicyclo[3.1.1]heptaneの7炭素骨格、メチル0〜3個、環上OHまたはケトン0/1個だけを検査。'
      'ピナン型のみ、最長橋の内側2–3または3–4結合にC=Cを1本まで許可します。'
      '固定番号の候補を主特性基→二重結合→メチル位置番号の順に比較します。全指定の大文字R/Sまたは全未指定に限定。'
      '環外C=C（β-ピネン型）、橋頭C=C、他の側鎖・部分立体・r/s・他の架橋骨格と保持名の優先性・PIN全体は未判定。CIPはRDKitに依存します。')

def analyze_bridged(mol):
    out=dict(detected=False,scope='unsupported',expected_name=None,numbering={},bridgeheads={},substituents=[],ring_stereo=[],reason='',note=NOTE,source=SOURCE)
    mol=Chem.Mol(mol)
    for a in mol.GetAtoms():a.SetIntProp('_bridged_original',a.GetIdx())
    mol=Chem.RemoveHs(mol)
    rings=mol.GetRingInfo().AtomRings()
    if not rings:return out
    core=set().union(*map(set,rings))
    if len(core)!=7 or any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 or mol.GetAtomWithIdx(i).GetIsAromatic() for i in core):return out
    if sum(b.GetBeginAtomIdx() in core and b.GetEndAtomIdx() in core for b in mol.GetBonds())!=8:return out
    registrations=[(key,q) for key,q in SCAFFOLDS.items() if any(set(m)==core for m in mol.GetSubstructMatches(q,uniquify=False))]
    if len(registrations)!=1:return out
    scaffold,query=registrations[0];head=4 if scaffold=='2.2.1' else 5
    core_doubles=[b for b in mol.GetBonds() if b.GetBeginAtomIdx() in core and b.GetEndAtomIdx() in core and b.GetBondType()==Chem.BondType.DOUBLE]
    out.update(detected=True,scaffold=scaffold,reason='bicyclo['+scaffold+']heptane型候補。置換基・主特性基・立体を限定照合します。')
    if len(core_doubles)>(1 if scaffold=='3.1.1' else 0):return out
    if len(Chem.GetMolFrags(mol))!=1 or mol.GetStereoGroups():return out
    if any(a.GetAtomicNum() not in (6,8) or a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() or a.GetChiralTag() not in (Chem.ChiralType.CHI_UNSPECIFIED,Chem.ChiralType.CHI_TETRAHEDRAL_CW,Chem.ChiralType.CHI_TETRAHEDRAL_CCW) for a in mol.GetAtoms()):return out
    if any(b.GetStereo()!=Chem.BondStereo.STEREONONE or b.GetBondType() not in (Chem.BondType.SINGLE,Chem.BondType.DOUBLE) for b in mol.GetBonds()):return out
    oxy=[a for a in mol.GetAtoms() if a.GetAtomicNum()==8]
    kind=None;functional=None
    if len(oxy)>1:return out
    if oxy:
        o=oxy[0]
        if o.GetDegree()!=1 or o.GetNeighbors()[0].GetIdx() not in core:return out
        functional=o.GetNeighbors()[0].GetIdx();bond=mol.GetBondBetweenAtoms(o.GetIdx(),functional)
        if bond.GetBondType()==Chem.BondType.SINGLE and o.GetTotalNumHs()==1:kind='alcohol'
        elif bond.GetBondType()==Chem.BondType.DOUBLE and o.GetTotalNumHs()==0:kind='ketone'
        else:return out
    if any(b.GetBondType()!=Chem.BondType.SINGLE and b.GetIdx() not in {d.GetIdx() for d in core_doubles} and not (kind=='ketone' and {b.GetBeginAtomIdx(),b.GetEndAtomIdx()}=={oxy[0].GetIdx(),functional}) for b in mol.GetBonds()):return out
    used=core|{a.GetIdx() for a in oxy};sides=[]
    for idx in core:
        for a in mol.GetAtomWithIdx(idx).GetNeighbors():
            first=a.GetIdx()
            if first in used:continue
            if a.GetAtomicNum()!=6 or mol.GetBondBetweenAtoms(idx,first).GetBondType()!=Chem.BondType.SINGLE:return out
            if a.GetDegree()!=1:return out
            sides.append((idx,first));used.add(first)
    if len(sides)>3 or used!=set(range(mol.GetNumAtoms())):return out
    options=[]
    for match in mol.GetSubstructMatches(query,uniquify=False,maxMatches=64):
        if set(match)!=core:continue
        mapping={a.GetAtomMapNum():match[a.GetIdx()] for a in query.GetAtoms()}
        pos={idx:n for n,idx in mapping.items()}
        if set(pos)!=core or len(set(pos.values()))!=7:continue
        methyl=sorted(pos[ring] for ring,atom in sides)
        loc=pos[functional] if kind else None
        double=tuple(sorted(pos[i] for i in (core_doubles[0].GetBeginAtomIdx(),core_doubles[0].GetEndAtomIdx()))) if core_doubles else ()
        # Deliberately bounded to the interior of the long bridge: no Bredt,
        # compound-locant, exocyclic or small-ring alkene claims.
        if double and double not in ((2,3),(3,4)):continue
        rank=((loc,) if kind else (),double,tuple(methyl))
        options.append((rank,mapping,methyl,loc,double))
    # Enumerate all graph automorphisms consistent with the fixed numbering.
    uniq={tuple(sorted((n,i) for n,i in item[1].items())):item for item in options}
    options=list(uniq.values())
    if len(options)<2:return out
    best_rank=min(o[0] for o in options);best=[o for o in options if o[0]==best_rank]
    centers=dict(Chem.FindMolChiralCenters(mol,includeUnassigned=True,useLegacyImplementation=False))
    if any(i not in core for i in centers):return out
    for a in mol.GetAtoms():
        if a.HasProp('_CIPCode'):a.ClearProp('_CIPCode')
    cip_error=False
    try:rdCIPLabeler.AssignCIPLabels(mol)
    except (RuntimeError,ValueError):cip_error=True
    index={a.GetIdx():a.GetIntProp('_bridged_original') for a in mol.GetAtoms()}
    arrangements=[]
    for rank,mapping,methyl,loc,double in best:
        stereo=[];bad=cip_error
        for n,idx in sorted(mapping.items()):
            if idx not in centers:continue
            a=mol.GetAtomWithIdx(idx);cip=a.GetProp('_CIPCode') if a.HasProp('_CIPCode') and not cip_error else None
            if a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED and cip not in ('R','S'):bad=True
            stereo.append(dict(position=str(n),atom_index=index[idx],cip=cip if cip in ('R','S','r','s') else None))
        assigned=[s for s in stereo if s['cip']]
        partial=0<len(assigned)<len(stereo)
        descriptor='('+','.join(s['position']+s['cip'] for s in assigned)+')-' if assigned else ''
        name=((','.join(map(str,methyl))+'-methyl') if len(methyl)==1 else (','.join(map(str,methyl))+'-'+{2:'di',3:'tri'}[len(methyl)]+'methyl') if methyl else '')
        parent='bicyclo['+scaffold+']hept'
        if double:parent+='-'+str(double[0])+('-en-'+str(loc)+'-'+('ol' if kind=='alcohol' else 'one') if kind else '-ene')
        else:parent+='an'+(('-'+str(loc)+'-'+('ol' if kind=='alcohol' else 'one')) if kind else 'e')
        arrangements.append(dict(mapping=mapping,methyl=methyl,loc=loc,double=double,stereo=stereo,partial=partial,bad=bad,name=descriptor+name+parent))
    ambiguous=len({x['name'] for x in arrangements})>1
    selected=arrangements[0]
    scope='unsupported_stereo' if selected['bad'] or ambiguous else ('partial_stereo' if selected['partial'] else 'bounded')
    reason='登録骨格の固定番号候補で主特性基→二重結合→メチル位置番号を比較。'
    if selected['bad'] or ambiguous:reason+='r/s、CIP計算失敗、または立体番号の曖昧さがあるため期待名を保留。'
    elif selected['partial']:reason+='立体中心の一部が未指定のため期待名を保留。'
    else:reason+='全指定R/Sまたは全未指定を保持。未指定はラセミ体と断定しません。'
    out.update(scope=scope,kind=kind,functional_atom=index[functional] if functional is not None else None,functional_oxygen=index[oxy[0].GetIdx()] if oxy else None,functional_locant=selected['loc'],methyl_locants=selected['methyl'],ring_stereo=selected['stereo'],numbering={str(n):index[i] for n,i in sorted(selected['mapping'].items())},bridgeheads={str(n):index[selected['mapping'][n]] for n in (1,head)},double_bond_locants=list(selected['double']),double_bond_atoms=[index[i] for i in (core_doubles[0].GetBeginAtomIdx(),core_doubles[0].GetEndAtomIdx())] if core_doubles else [],expected_name=selected['name'] if scope=='bounded' else None,base_expected_name=selected['name'],reason=reason,orientations=[dict(functional=o[3],methyl=o[2],double=list(o[4]),numbering={str(n):index[i] for n,i in sorted(o[1].items())}) for o in options],substituents=[dict(position=str(next(n for n,i in o['mapping'].items() if i==ring)),kind='methyl',ring_atom=index[ring],atoms=[index[atom]]) for o in (selected,) for ring,atom in sides])
    return out

def parse_bridged_name(name):
    pattern=(r'(?:(\((?P<stereo>[0-9RS,]+)\)-))?'
             r'(?:(?P<methyl>\d+(?:,\d+)*)-(?P<count>di|tri)?methyl-?)?'
             r'(?:(?P<oxo>\d+)-oxo-)?bicyclo\[(?P<a>\d+)\.(?P<b>\d+)\.(?P<c>\d+)\]'
             r'(?P<root>hept|oct|non|dec)(?:ane|an-(?P<loc>\d+)-(?P<suffix>ol|one)|'
             r'-(?P<first>\d+)(?:\((?P<second>\d+)\))?-(?:ene|en-(?P<enloc>\d+)-(?P<ensuffix>ol|one)))')
    m=re.fullmatch(pattern,name.strip())
    if not m:return None
    g=m.groupdict();stereo={}
    if g['stereo']:
        for token in g['stereo'].split(','):
            x=re.fullmatch(r'(\d+)([RS])',token)
            if not x or x[1] in stereo:return None
            stereo[x[1]]=x[2]
    methyl=list(map(int,g['methyl'].split(','))) if g['methyl'] else []
    amount={'di':2,'tri':3,None:1}[g['count']] if methyl else 0
    suffix=g['suffix'] or g['ensuffix']
    kind=('alcohol' if suffix=='ol' else 'ketone') if suffix else ('oxo_prefix' if g['oxo'] else None)
    loc=g['oxo'] or g['loc'] or g['enloc']
    double=[int(g['first']),int(g['second']) if g['second'] else int(g['first'])+1] if g['first'] else []
    return dict(scaffold='.'.join(g[k] for k in ('a','b','c')),root=g['root'],methyl=methyl,multiplicity_ok=len(methyl)==amount,kind=kind,functional=int(loc) if loc else None,double=double,stereo=stereo)

def bridged_checks(analysis,name):
    if analysis['scope'] not in ('bounded','partial_stereo'):return []
    p=parse_bridged_name(name)
    def row(k,t,s,d):return dict(id='bridged_'+k,title=t,status=s,detail=d,reference='P-23.2.2 / P-23.2.3 / P-14.4 / P-31.1.4 / P-63',url=SOURCE)
    if p is None:return [row('grammar','架橋二環の限定系統名','unknown','候補名は限定構文外です。保持名や別表記の判断は保留します。')]
    parent=p['scaffold']==analysis['scaffold'] and p['root']=='hept'
    double_ok=p['double']==analysis['double_bond_locants']
    positions_same=p['methyl']==analysis['methyl_locants'] and p['functional']==analysis['functional_locant'] and double_ok
    same=parent and positions_same and p['kind']==analysis['kind']
    # Oxo-prefix names can represent the same positions but violate suffix policy.
    represented=parent and any(o['functional']==p['functional'] and o['methyl']==p['methyl'] and o['double']==p['double'] for o in analysis['orientations'])
    actual={s['position']:s['cip'] for s in analysis['ring_stereo'] if s['cip']}
    stereo_ok=actual==p['stereo']
    suffix_ok=p['kind']==analysis['kind']
    return [row('parent','架橋親骨格と橋長','pass' if parent else 'fail','登録骨格bicyclo['+analysis['scaffold']+']heptane、環骨格7炭素を対応。'),
            row('suffix','主特性基と接尾辞','pass' if suffix_ok else 'fail','OHはol、ケトンはone接尾辞を限定検査。入力の主特性基：'+str(analysis['kind'])+'。'),
            row('numbering','橋頭を基準にした位置番号','pass' if same else 'fail','採用位置：主特性基 '+str(analysis['functional_locant'])+'、二重結合 '+str(analysis['double_bond_locants'])+'、メチル '+str(analysis['methyl_locants'])+'。'),
            row('double','環内二重結合の位置','pass' if parent and double_ok else 'fail','採用結合番号 '+str(analysis['double_bond_locants'])+'。'),
            row('descriptor','置換位置と構造の対応','pass' if represented else 'fail','許容する固定番号候補と置換位置を照合。'),
            row('stereo','環立体のR/S','unknown' if analysis['scope']=='partial_stereo' else ('pass' if stereo_ok else 'fail'),'大文字R/Sのみを採用番号で比較し、未指定やr/sは補完しません。'),
            row('spelling','限定表記・置換基数','fail' if not p['multiplicity_ok'] else ('pass' if name.strip()==analysis['expected_name'] else 'different'),'限定期待表記：'+str(analysis['expected_name'])+'。PIN全体は未判定。')]
