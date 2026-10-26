"""Unspecified-stereo oxane/oxolane polyols: a carbohydrate-skeleton foothold.
This module does not identify sugars or assign carbohydrate stereodescriptors.
"""
import re
from rdkit import Chem
from extended_rules import MULT,SUBS,row,final_rows,unknown_result,group_rank
from natural_rules import TEMPLATES,SUGAR_NOMENCLATURE_NOTE
from linkage_rules import ALKOXY,flat_parse
from carbon_rules import HALO

PARENTS=('oxane','oxolane')
ARYLOXY={label:Chem.MolFromSmiles(smi) for label,smi in [('phenoxy','Oc1ccccc1'),('benzyloxy','OCc1ccccc1')]}
OXY_LABELS=set(ALKOXY.values())|set(ARYLOXY)
ARYL={label:Chem.MolFromSmiles(smi) for label,smi in [('phenyl','c1ccccc1'),('benzyl','Cc1ccccc1')]}
AROMATIC_LABELS=set(ARYLOXY)|set(ARYL)
ALKOXYMETHYL={n:label+'methyl' for n,label in ALKOXY.items()}
COMPOUND_LABELS={'hydroxymethyl'}|set(ALKOXYMETHYL.values())
LABELS={SUBS[n] for n in range(1,5)}|set(HALO.values())|set(ALKOXY.values())|{'amino','hydroxymethyl'}|set(ARYLOXY)|set(ARYL)|set(ALKOXYMETHYL.values())


def component(mol,start,excluded):
    seen=set();stack=[start]
    while stack:
        idx=stack.pop()
        if idx in seen or idx in excluded:continue
        seen.add(idx);stack.extend(a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors())
    return seen


def plain_alkyl(mol,indices,attachment):
    if not 1<=len(indices)<=4 or any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 for i in indices):return None
    degrees={i:sum(a.GetIdx() in indices for a in mol.GetAtomWithIdx(i).GetNeighbors()) for i in indices}
    if degrees[attachment]>1 or any(d>2 for d in degrees.values()):return None
    if any(b.GetBondType()!=Chem.BondType.SINGLE for b in mol.GetBonds() if b.GetBeginAtomIdx() in indices and b.GetEndAtomIdx() in indices):return None
    return SUBS[len(indices)]


def inventory(mol,path):
    ring=set(path);oh=[];groups=[]
    for loc,idx in enumerate(path,1):
        for atom in mol.GetAtomWithIdx(idx).GetNeighbors():
            first=atom.GetIdx()
            if first in ring:continue
            if mol.GetAtomWithIdx(idx).GetAtomicNum()!=6:return None
            if mol.GetBondBetweenAtoms(idx,first).GetBondType()!=Chem.BondType.SINGLE:return None
            side=component(mol,first,ring)
            if len(side)==1:
                if atom.GetAtomicNum()==8 and atom.GetTotalNumHs()==1:oh.append(loc);continue
                if atom.GetAtomicNum()==7 and atom.GetTotalNumHs()==2:groups.append((loc,'amino'));continue
                if atom.GetAtomicNum() in HALO:groups.append((loc,HALO[atom.GetAtomicNum()]));continue
            if atom.GetAtomicNum()==6:
                aryl=next((label for label,template in ARYL.items() if any(set(match)==side and match[0]==first for match in mol.GetSubstructMatches(template,uniquify=False))),None)
                if aryl:groups.append((loc,aryl));continue
                label=plain_alkyl(mol,side,first)
                if label:groups.append((loc,label));continue
                if len(side)==2 and atom.GetTotalNumHs()==2:
                    oxygen=next((mol.GetAtomWithIdx(i) for i in side if i!=first),None)
                    if (oxygen and oxygen.GetAtomicNum()==8 and oxygen.GetDegree()==1 and oxygen.GetTotalNumHs()==1
                            and mol.GetBondBetweenAtoms(first,oxygen.GetIdx()).GetBondType()==Chem.BondType.SINGLE):
                        groups.append((loc,'hydroxymethyl'));continue
            if atom.GetAtomicNum()==6 and atom.GetTotalNumHs()==2 and atom.GetDegree()==2:
                oxy=next((a for a in atom.GetNeighbors() if a.GetIdx() not in ring),None)
                if oxy and oxy.GetAtomicNum()==8 and oxy.GetDegree()==2 and mol.GetBondBetweenAtoms(first,oxy.GetIdx()).GetBondType()==Chem.BondType.SINGLE:
                    rest=side-{first,oxy.GetIdx()};anchor=next((a.GetIdx() for a in oxy.GetNeighbors() if a.GetIdx()!=first),None)
                    if anchor is not None and plain_alkyl(mol,rest,anchor) and mol.GetBondBetweenAtoms(oxy.GetIdx(),anchor).GetBondType()==Chem.BondType.SINGLE:
                        groups.append((loc,ALKOXYMETHYL[len(rest)]));continue
            if atom.GetAtomicNum()==8 and atom.GetDegree()==2:
                aromatic=next((label for label,template in ARYLOXY.items() if any(set(match)==side and match[0]==first for match in mol.GetSubstructMatches(template,uniquify=False))),None)
                if aromatic:groups.append((loc,aromatic));continue
                rest=side-{first};carbon=next((a.GetIdx() for a in atom.GetNeighbors() if a.GetIdx() not in ring),None)
                if carbon is not None and plain_alkyl(mol,rest,carbon):
                    if mol.GetBondBetweenAtoms(first,carbon).GetBondType()==Chem.BondType.SINGLE:
                        groups.append((loc,ALKOXY[len(rest)]));continue
            return None
    return dict(locants=sorted(oh),groups=sorted(groups))


def prefix_text(groups):
    grouped={}
    for loc,label in groups:grouped.setdefault(label,[]).append(loc)
    chunks=[]
    for label in sorted(grouped):
        locs=sorted(grouped[label]);lead=','.join(map(str,locs))+'-'
        if label in COMPOUND_LABELS:chunks.append(lead+('bis' if len(locs)==2 else '')+'('+label+')')
        else:chunks.append(lead+MULT[len(locs)]+label)
    return '-'.join(chunks)


def parse_name(name):
    text=name.strip().lower()
    for label in sorted(COMPOUND_LABELS,key=len,reverse=True):
        text=text.replace('bis('+label+')','di'+label).replace('('+label+')',label)
    for parent in PARENTS:
        for count in range(1,5):
            base=parent[:-1] if count==1 else parent
            ending=(MULT[count].removesuffix('a') if count>1 else '')+'ol'
            match=re.fullmatch(r'(.*?)'+base+r'(?:-(\d+(?:,\d+)*)-)?'+ending,text)
            if not match:continue
            parsed=flat_parse(match[1],LABELS)
            if parsed is None or any(isinstance(i,str) for i,_ in parsed['groups']):continue
            locants=list(map(int,match[2].split(','))) if match[2] else None
            parsed['multiplicity_ok'] &= locants is None or len(locants)==count
            parsed['groups']=sorted(parsed['groups'])
            return dict(parent=parent,locants=locants,**parsed)
    # Narrow legacy shape observed from NISPO. Keep its OH-as-prefix form
    # visible rather than rewriting the supplied name to the expected suffix.
    legacy=re.fullmatch(r'(.+?)(oxan|oxolan)',text)
    if legacy:
        parsed=flat_parse(legacy[1],LABELS|{'hydroxy'})
        if parsed is not None and all(isinstance(i,int) for i,_ in parsed['groups']):
            locants=sorted(i for i,label in parsed['groups'] if label=='hydroxy')
            if locants:
                parsed['groups']=sorted((i,label) for i,label in parsed['groups'] if label!='hydroxy')
                return dict(parent=legacy[2]+'e',locants=locants,hydroxy_prefix_only=True,**parsed)
    return None


def oxygen_ring_checks(mol,name):
    rings=mol.GetRingInfo().AtomRings()
    oxygen_rings=[set(r) for r in rings if len(r) in (5,6) and sum(mol.GetAtomWithIdx(i).GetAtomicNum()==8 for i in r)==1]
    if len(oxygen_rings)!=1 or len(rings)>2:return None
    ring=oxygen_rings[0]
    # A second ring is allowed only inside one complete phenoxy/benzyloxy side.
    for other in rings:
        if set(other)==ring:continue
        if set(other)&ring or len(other)!=6 or any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 or not mol.GetAtomWithIdx(i).GetIsAromatic() for i in other):return None
    if (sum(mol.GetAtomWithIdx(i).GetAtomicNum()==8 for i in ring)!=1
            or any(mol.GetAtomWithIdx(i).GetAtomicNum() not in (6,8) or mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)):return None
    # Keep existing bare-ring/alkyl-only behavior. At least one OH on ring C
    # is required; an external OH alone needs a different parent-selection rule.
    has_oh=any(a.GetAtomicNum()==8 and a.GetDegree()==1 and a.GetTotalNumHs()==1 and a.GetNeighbors()[0].GetIdx() in ring for a in mol.GetAtoms())
    if not has_oh:return None
    parent='oxolane' if len(ring)==5 else 'oxane'
    scope=dict(kind='oxygen_ring_polyol',parent_name=parent,n=len(ring),scope_note=SUGAR_NOMENCLATURE_NOTE,
               sugar_nomenclature_status='unknown',registry_numbering=Chem.MolToSmiles(TEMPLATES[parent]))
    candidates=[]
    for match in mol.GetSubstructMatches(TEMPLATES[parent],uniquify=False,maxMatches=1000):
        if set(match)!=ring:continue
        path=[None]*len(match)
        for idx,target in enumerate(match):path[TEMPLATES[parent].GetAtomWithIdx(idx).GetAtomMapNum()-1]=target
        if any(c['path']==path for c in candidates):continue
        inv=inventory(mol,path)
        if inv is None:return unknown_result(scope,'この環の置換には未対応です。羰基・環の不飽和・分岐や官能基を持つ側鎖・複雑なアセタール等は保留します。')
        candidates.append(dict(path=path,**inv))
    if not candidates:return unknown_result(scope,'飽和oxane/oxolaneの固定した環骨格に一致しないため保留します。')
    first=candidates[0];oh=first['locants'];groups=first['groups']
    if not 1<=len(oh)<=4 or len(set(oh))!=len(oh):return unknown_result(scope,'環上のOHは別々の炭素に1〜4個の限定対応です。gem-ジオールや5個以上のOHは保留します。')
    counts={label:sum(x==label for _,x in groups) for _,label in groups}
    if len(groups)>3 or counts.get('hydroxymethyl',0)>2 or sum(label in ALKOXYMETHYL.values() for _,label in groups)>1 or counts.get('amino',0)>1 or sum(label in OXY_LABELS for _,label in groups)>2 or sum(label in AROMATIC_LABELS for _,label in groups)>1:
        return unknown_result(scope,'環上の接頭辞は合計3個まで、hydroxymethyl最大2個、alkoxymethyl最大1個、amino最大1個、O結合側鎖は合計2個まで、phenoxy/benzyloxy/phenyl/benzylは合計1個までです。')
    candidates.sort(key=lambda c:(c['locants'],group_rank(c['groups'])))
    best=candidates[0];locants=best['locants'];groups=best['groups'];count=len(locants)
    suffix=(parent[:-1] if count==1 else parent)+'-'+','.join(map(str,locants))+'-'+('' if count==1 else MULT[count].removesuffix('a'))+'ol'
    expected=prefix_text(groups)+suffix
    scope.update(expected=expected,path=best['path'],ring_oh_count=count,hydroxymethyl_count=counts.get('hydroxymethyl',0),
                 candidates=candidates,preferred_candidates=[c for c in candidates if (c['locants'],group_rank(c['groups']))==(locants,group_rank(groups))],
                 selection_order=['環のOを1とする固定番号','環上OHの位置番号列','接頭辞の位置番号列','名称中の列挙順'])
    parsed=parse_name(name)
    if parsed is None:return unknown_result(scope,'骨格は限定範囲ですが、この名称構文・糖類保持名は未対応です。glucose等の糖名を骨格表記と同じものとして扱いません。')
    same_parent=parsed['parent']==parent;complete=parsed['locants'] is not None
    represented=same_parent and complete and any(c['locants']==parsed['locants'] and c['groups']==parsed['groups'] for c in candidates)
    compatible=same_parent and complete
    def conditional(ok,value):return 'unknown' if not ok else ('pass' if value else 'fail')
    checks=[
      row('parent','環状ポリオールの親骨格','pass' if same_parent else 'fail',f"飽和した{parent}環、O位置1。名称の親骨格は{parsed['parent']}。外側のOHがあっても、この範囲では環上OHを主特性基として扱います。",'P-22.2.2 / P-44.1','P4'),
      row('suffix','環上OHの接尾辞と期待表記の比較','different' if parsed.get('hydroxy_prefix_only') else 'pass','期待表記では環上OHをol/diol/triol/tetrol接尾辞で表します。hydroxy接頭辞による表記との差を示し、構造・立体の誤りやIUPAC名全体の無効を意味しません。','P-43 / P-63.1','P6'),
      row('descriptor','OH・ヒドロキシメチル等と環の構造対応',conditional(compatible,represented),f"名称のOH位置は{parsed['locants']}、接頭辞は{parsed['groups']}。固定番号を保つ両方向と照合。",'P-15.1.7'),
      row('numbering','環上OHを優先した番号付け',conditional(compatible,parsed['locants']==locants),f"採用OH位置は{locants}。環のO位置1を保ち、その後OH、接頭辞の順に比較。糖類命名の番号との一致は保証しません。",'P-14.4 / P-44.4.1.8','P4'),
      row('functional_prefix','環の接頭辞と列挙順',conditional(compatible,parsed['groups']==groups),f"採用する配置は{groups}。hydroxymethylはOH1個を持つCH2側鎖として構造確認。",'P-45.2 / P-63.1','P6'),
      row('multiplicity','OH・置換基数と倍数表記','pass' if parsed['multiplicity_ok'] else 'fail','diol/triol/tetrolと位置番号、接頭辞の倍数を照合。','P-14.2'),
      row('spelling','置換命名法による限定した骨格表記','pass' if name.strip().lower()==expected else ('different' if parsed.get('hydroxy_prefix_only') and parsed['locants']==locants and parsed['groups']==groups and same_parent else 'fail'),f"期待する骨格表記：{expected}。糖類保持名やPINとしての優先性は検査していません。",'P-22.2 / P-63.1 / P-102','P10'),
      row('sugar_nomenclature','糖の同定・D/L・α/β・糖類命名','unknown',SUGAR_NOMENCLATURE_NOTE,'P-102','P10'),
    ]
    return dict(scope=scope,checks=final_rows(checks))
