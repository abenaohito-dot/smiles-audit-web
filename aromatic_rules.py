"""A single benzene ring with explicitly bounded functional substituents."""
import re
from rdkit import Chem
from extended_rules import SUBS,MULT,row,final_rows,unknown_result,group_rank
from linkage_rules import ALKOXY,flat_parse
from carbon_rules import HALO

PRIORITY={'acid':0,'ester':1,'amide':2,'nitrile':3,'aldehyde':4,'alcohol':5,'amine':6}
LOWER={'nitrile':'cyano','aldehyde':'formyl','alcohol':'hydroxy','amine':'amino'}
PARENTS={'ester':'benzoate','amide':'benzamide','acid':'benzoic acid','nitrile':'benzonitrile','aldehyde':'benzaldehyde','alcohol':'phenol','amine':'aniline',None:'benzene'}
LABELS=set(SUBS[n] for n in range(1,5))|set(HALO.values())|set(ALKOXY.values())|set(LOWER.values())|{'isopropyl'}|{label+'carbonyl' for label in ALKOXY.values()}


def side_component(mol,start,ring):
    seen=set();stack=[start]
    while stack:
        idx=stack.pop()
        if idx in seen or idx in ring:continue
        seen.add(idx);stack.extend(a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors())
    return seen


def linear_alkyl(mol,indices,start):
    if not 1<=len(indices)<=4:return None
    if any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 for i in indices):return None
    if any(b.GetBondType()!=Chem.BondType.SINGLE for b in mol.GetBonds() if b.GetBeginAtomIdx() in indices and b.GetEndAtomIdx() in indices):return None
    degrees={i:sum(nb.GetIdx() in indices for nb in mol.GetAtomWithIdx(i).GetNeighbors()) for i in indices}
    if len(indices)==3 and degrees[start]==2:return 'propan-2-yl'
    if degrees[start]>1 or any(d>2 for d in degrees.values()):return None
    return SUBS[len(indices)]


def classify_side(mol,idx,first,ring):
    atom=mol.GetAtomWithIdx(first);indices=side_component(mol,first,ring)
    if mol.GetBondBetweenAtoms(idx,first).GetBondType()!=Chem.BondType.SINGLE:return None
    number=atom.GetAtomicNum()
    if len(indices)==1:
        if number in HALO:return ('prefix',HALO[number])
        if number==8 and atom.GetTotalNumHs()==1:return ('alcohol',None)
        if number==7 and atom.GetTotalNumHs()==2:return ('amine',None)
    if number==6:
        label=linear_alkyl(mol,indices,first)
        if label:return ('prefix',label)
        outside=[a for a in atom.GetNeighbors() if a.GetIdx() not in ring]
        if len(indices)==2 and len(outside)==1:
            nb=outside[0];bond=mol.GetBondBetweenAtoms(first,nb.GetIdx()).GetBondType()
            if nb.GetAtomicNum()==7 and bond==Chem.BondType.TRIPLE:return ('nitrile',None)
            if nb.GetAtomicNum()==8 and bond==Chem.BondType.DOUBLE and atom.GetTotalNumHs()==1:return ('aldehyde',None)
        derivative=classify_acyl(mol,first,indices,ring)
        if derivative is not None:return derivative
        if len(indices)==3 and len(outside)==2 and all(a.GetAtomicNum()==8 and a.GetDegree()==1 for a in outside):
            bonds=[mol.GetBondBetweenAtoms(first,a.GetIdx()).GetBondType() for a in outside]
            if sorted(map(str,bonds))==['DOUBLE','SINGLE'] and any(a.GetTotalNumHs()==1 for a in outside):return ('acid',None)
    if number==8 and atom.GetDegree()==2:
        others=[a for a in atom.GetNeighbors() if a.GetIdx() not in ring]
        if len(others)==1:
            alkyl=indices-{first};label=linear_alkyl(mol,alkyl,others[0].GetIdx())
            if label in SUBS.values() and mol.GetBondBetweenAtoms(first,others[0].GetIdx()).GetBondType()==Chem.BondType.SINGLE:return ('prefix',ALKOXY[len(alkyl)])
    return None


def classify_acyl(mol,carbonyl,indices,ring):
    outside=[a for a in mol.GetAtomWithIdx(carbonyl).GetNeighbors() if a.GetIdx() not in ring]
    if len(outside)!=2:return None
    oxygens=[a for a in outside if a.GetAtomicNum()==8 and a.GetDegree()==1 and mol.GetBondBetweenAtoms(carbonyl,a.GetIdx()).GetBondType()==Chem.BondType.DOUBLE]
    if len(oxygens)!=1:return None
    oxygen=oxygens[0].GetIdx();linker=next(a for a in outside if a.GetIdx()!=oxygen)
    if mol.GetBondBetweenAtoms(carbonyl,linker.GetIdx()).GetBondType()!=Chem.BondType.SINGLE:return None
    if linker.GetAtomicNum()==8 and linker.GetDegree()==2:kind='ester'
    elif linker.GetAtomicNum()==7 and 1<=linker.GetDegree()<=3:kind='amide'
    else:return None
    rest=[a for a in linker.GetNeighbors() if a.GetIdx()!=carbonyl]
    consumed={carbonyl,oxygen,linker.GetIdx()};alkyls=[];paths=[]
    for atom in rest:
        if mol.GetBondBetweenAtoms(linker.GetIdx(),atom.GetIdx()).GetBondType()!=Chem.BondType.SINGLE:return None
        # Do not traverse back through N/O into a second branch.
        seen=set();stack=[atom.GetIdx()]
        while stack:
            idx=stack.pop()
            if idx in seen or idx==linker.GetIdx():continue
            seen.add(idx)
            stack.extend(a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors() if a.GetIdx()!=linker.GetIdx())
        label=linear_alkyl(mol,seen,atom.GetIdx())
        if label not in {SUBS[n] for n in range(1,5)}:return None
        if seen & consumed:return None
        consumed|=seen;alkyls.append(label);paths.append(sorted(seen))
    if consumed!=indices:return None
    return kind,dict(organyl=alkyls[0] if kind=='ester' else None,
                     n_groups=[('N',label) for label in alkyls] if kind=='amide' else [],
                     linker_index=linker.GetIdx(),carbonyl_index=carbonyl,side_paths=paths)


def sort_groups(groups):
    return sorted(groups,key=lambda x:(str(x[0]),x[1]))


def aromatic_prefix(groups,omit_single=False):
    grouped={}
    for loc,label in groups:grouped.setdefault(label,[]).append(loc)
    chunks=[]
    for label in sorted(grouped):
        locs=sorted(grouped[label],key=lambda x:(0 if isinstance(x,str) else 1,x));lead='' if omit_single and len(groups)==1 else ','.join(map(str,locs))+'-'
        if label=='propan-2-yl':chunks.append(lead+'(propan-2-yl)')
        else:chunks.append(lead+MULT[len(locs)]+label)
    return '-'.join(chunks)


def ring_paths(mol,ring):
    allowed=set(ring);adj={i:[a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() in allowed] for i in ring}
    for start in ring:
        for nextidx in adj[start]:
            path=[start,nextidx]
            while len(path)<6:path.append(next(i for i in adj[path[-1]] if i!=path[-2]))
            yield path


def parse_aromatic(name):
    text=name.strip().lower().replace('(propan-2-yl)','isopropyl')
    # Inspect observed lower-priority substitutive forms without endorsing
    # them as the preferred functional-class ester spelling.
    for label in ALKOXY.values():text=text.replace('('+label+'carbonyl)',label+'carbonyl')
    organyl=None
    if text.endswith('benzoate'):
        if ' ' not in text:return None
        organyl,text=text.split(' ',1)
        if organyl not in {SUBS[n] for n in range(1,5)}:return None
    if text=='toluene':return dict(principal=None,locants=[],groups=[(1,'methyl')],multiplicity_ok=True)
    if text=='anisole':return dict(principal=None,locants=[],groups=[(1,'methoxy')],multiplicity_ok=True)
    match=re.fullmatch(r'1,([234])-xylene',text)
    if match:return dict(principal=None,locants=[],groups=[(1,'methyl'),(int(match[1]),'methyl')],multiplicity_ok=True)
    for kind,parent in PARENTS.items():
        if not text.endswith(parent) or (kind=='ester')!=(organyl is not None):continue
        prefix=text[:-len(parent)]
        parsed=flat_parse(prefix,LABELS)
        if parsed is None and kind is None and prefix in LABELS:parsed=dict(groups=[(1,prefix)],multiplicity_ok=True)
        if parsed is not None:
            n_groups=[(i,label) for i,label in parsed['groups'] if isinstance(i,str)]
            if n_groups and (kind!='amide' or any(label not in {SUBS[n] for n in range(1,5)} for _,label in n_groups)):continue
            parsed['groups']=sorted((i,'propan-2-yl' if label=='isopropyl' else label) for i,label in parsed['groups'] if not isinstance(i,str))
            return dict(principal=kind,locants=[1] if kind else [],organyl=organyl,n_groups=sort_groups(n_groups),**parsed)
    for kind,suffix in [('alcohol','ol'),('amine','amine')]:
        match=re.fullmatch(r'(.*?)benzene-([\d,]+)-(di|tri)'+suffix,text)
        if not match:continue
        parsed=flat_parse(match[1],LABELS)
        if parsed is None or any(isinstance(i,str) for i,_ in parsed['groups']):continue
        locants=list(map(int,match[2].split(',')))
        parsed['multiplicity_ok'] &= len(locants)==(2 if match[3]=='di' else 3)
        parsed['groups']=sorted((i,'propan-2-yl' if label=='isopropyl' else label) for i,label in parsed['groups'])
        return dict(principal=kind,locants=locants,**parsed)
    return None


def aromatic_rule_checks(mol,name):
    rings=mol.GetRingInfo().AtomRings()
    if len(rings)!=1 or len(rings[0])!=6:return None
    ring=set(rings[0])
    if not all(mol.GetAtomWithIdx(i).GetAtomicNum()==6 and mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):return None
    attachments=[]
    for idx in ring:
        for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
            if nb.GetIdx() in ring:continue
            side=classify_side(mol,idx,nb.GetIdx(),ring)
            if side is None:return unknown_result(dict(kind='functional_benzene'),'ベンゼン環ですが、この側鎖・官能基は未対応です。複雑な分岐、N置換アニリン、環を含む側鎖、芳香族アシル基の複雑なO側／N側置換は保留します。')
            attachments.append((idx,*side))
    if not attachments:return None  # Keep the established unsubstituted benzene checker.
    scope=dict(kind='functional_benzene',n=6)
    if len(attachments)>4 or sum(label=='propan-2-yl' for _,_,label in attachments)>1:return unknown_result(scope,'ベンゼン上は合計4置換まで、propan-2-ylは1個までの限定対応です。')
    inventory={kind:[idx for idx,k,_ in attachments if k==kind] for kind in PRIORITY}
    if any(len(inventory[kind])>1 for kind in ('acid','ester','amide','nitrile','aldehyde')) or any(len(inventory[kind])>3 for kind in ('alcohol','amine')):return unknown_result(scope,'酸・エステル・アミド・CHO・CNは各1個、OH／NH2は各3個までの限定対応です。')
    if sum(bool(inventory[k]) for k in ('acid','ester','amide'))>1:return unknown_result(scope,'酸・エステル・アミドの複合と多アシル誘導体は、下位アシル基の接頭辞規則を保留します。')
    principal=next((k for k in PRIORITY if inventory[k]),None)
    acyl=next((label for _,kind,label in attachments if kind==principal and kind in ('ester','amide')),None)
    n_groups=sort_groups(acyl['n_groups']) if acyl else []
    organyl=acyl['organyl'] if acyl else None
    candidates=[]
    for path in ring_paths(mol,ring):
        locants=sorted(path.index(i)+1 for i in inventory[principal]) if principal else []
        if principal and len(locants)==1 and locants!=[1]:continue
        groups=[]
        for idx,kind,label in attachments:
            if kind==principal:continue
            groups.append((path.index(idx)+1,label if kind=='prefix' else LOWER[kind]))
        candidates.append(dict(path=path,locants=locants,groups=sorted(groups)))
    candidates.sort(key=lambda c:(c['locants'],group_rank(c['groups'])))
    best=candidates[0];locants=best['locants'];groups=best['groups']
    parent=PARENTS[principal]
    if principal in ('alcohol','amine') and len(locants)>1:parent='benzene-'+','.join(map(str,locants))+'-'+MULT[len(locants)]+('ol' if principal=='alcohol' else 'amine')
    expected=(organyl+' ' if organyl else '')+aromatic_prefix(groups+n_groups,omit_single=principal is None and len(groups)==1)+parent
    if principal is None:
        if [label for _,label in groups]==['methyl']:expected='toluene'
        elif [label for _,label in groups]==['methyl','methyl']:expected='1,'+str(groups[1][0])+'-xylene'
        elif [label for _,label in groups]==['methoxy']:expected='anisole'
    if acyl:scope.update(acyl=acyl,side_paths=acyl['side_paths'],n_groups=n_groups,organyl=organyl)
    scope.update(expected=expected,path=best['path'],principal=principal,parent_name=parent,inventory={k:len(v) for k,v in inventory.items()},
                 candidates=candidates,preferred_candidates=[c for c in candidates if (c['locants'],group_rank(c['groups']))==(best['locants'],group_rank(groups))],
                 selection_order=['主特性基の選択','主特性基の位置番号列','接頭辞の位置番号列','名称中の列挙順'])
    parsed=parse_aromatic(name)
    if parsed is None:return unknown_result(scope,'構造は限定範囲ですが、この保持名・別名や芳香族の名称構文は未対応です。期待表記の構造照合は別に実行します。')
    compatible=parsed['principal']==principal
    attachment_ok=parsed.get('organyl')==organyl and parsed.get('n_groups',[])==n_groups
    represented=compatible and attachment_ok and any(c['locants']==parsed['locants'] and c['groups']==parsed['groups'] for c in candidates)
    checks=[
      row('parent','単一ベンゼン環と主特性基','pass' if compatible else 'fail',f"ベンゼン単環上の主特性基は{principal or 'なし'}、名称は{parsed['principal'] or 'なし'}。採用骨格名は{parent}。",'P-41 / P-44.1','P4'),
      row('descriptor','芳香環の置換基と構造の対応','unknown' if not compatible else ('pass' if represented else 'fail'),f"名称の主特性基位置は{parsed['locants']}、接頭辞は{parsed['groups']}。環の候補と照合。",'P-15.1.7'),
      row('numbering','主特性基と環の番号付け','unknown' if not compatible else ('pass' if parsed['locants']==locants and parsed['groups']==groups else 'fail'),f"採用する主特性基位置は{locants}、接頭辞配置は{groups}。主特性基を先に小さくします。",'P-14.4 / P-45.2'),
      row('multiplicity','官能基・置換基数と倍数表記','pass' if parsed['multiplicity_ok'] else 'fail','di/tri接尾辞・接頭辞と位置番号の個数を照合。','P-14.2'),
      row('spelling','芳香族の限定した優先表記','pass' if name.strip().lower()==expected.lower() else 'fail',f"期待表記：{expected}。phenol/aniline/benzoic acid等の限定保持名を含みます。全保持名・PIN全体は未判定です。",'P-22.1.3 / P-62.2 / P-63.1 / P-65.1 / P-66.6','P6'),
    ]
    if acyl:
        checks.insert(1,row('attachment','O側アルキル／アミドN置換の構造対応','unknown' if not compatible else ('pass' if attachment_ok else 'fail'),f"O側は{organyl}、N置換は{n_groups}。名称は{parsed.get('organyl')}、{parsed.get('n_groups',[])}。環の番号とN位置を区別して照合。",'P-65.6.3.3.1' if principal=='ester' else 'P-66.1.1.3.1','P6' if principal=='ester' else 'P6a'))
    return dict(scope=scope,checks=final_rows(checks))
