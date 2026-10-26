"""Bounded acyclic ether, monoester and monoamide nomenclature checks.
No naming-engine calls: structures and finite name grammars are inspected here.
"""
from collections import defaultdict
import re
from rdkit import Chem
from extended_rules import ALKANES,SUBS,MULT,row,final_rows,unknown_result,group_rank,prefix_name
from natural_rules import PRIORITY,PREFIX,function_suffix
from carbon_rules import HALO

ALKOXY={1:'methoxy',2:'ethoxy',3:'propoxy',4:'butoxy'}
LABELS=set(HALO.values())|{'amino','hydroxy','oxo'}|set(ALKOXY.values())


def carbon_components(mol):
    """Every carbon component must be a saturated path of 1--10 carbons."""
    remaining={a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==6}
    paths=[]
    while remaining:
        seen=set();stack=[min(remaining)]
        while stack:
            idx=stack.pop()
            if idx in seen:continue
            seen.add(idx)
            for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
                if nb.GetAtomicNum()==6:
                    if mol.GetBondBetweenAtoms(idx,nb.GetIdx()).GetBondType()!=Chem.BondType.SINGLE:return None
                    stack.append(nb.GetIdx())
        remaining-=seen
        adj={i:[nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() in seen] for i in seen}
        if len(seen)>10 or any(len(v)>2 for v in adj.values()):return None
        start=min(i for i,v in adj.items() if len(v)<=1)
        path=[];previous=None;current=start
        while current is not None:
            path.append(current)
            nexts=[i for i in adj[current] if i!=previous]
            previous,current=current,nexts[0] if nexts else None
        if len(path)!=len(seen):return None
        paths.append(path)
    return paths


def terminal_inventory(mol,path,excluded):
    inv={k:[] for k in ('aldehyde','ketone','alcohol','amine')};halo=[]
    for idx in path:
        for a in mol.GetAtomWithIdx(idx).GetNeighbors():
            if a.GetAtomicNum()==6 or a.GetIdx() in excluded:continue
            if a.GetDegree()!=1:return None
            bond=mol.GetBondBetweenAtoms(idx,a.GetIdx()).GetBondType()
            if a.GetAtomicNum() in HALO and bond==Chem.BondType.SINGLE:halo.append((idx,HALO[a.GetAtomicNum()]))
            elif a.GetAtomicNum()==8 and bond==Chem.BondType.SINGLE and a.GetTotalNumHs()==1:inv['alcohol'].append(idx)
            elif a.GetAtomicNum()==7 and bond==Chem.BondType.SINGLE and a.GetTotalNumHs()==2:inv['amine'].append(idx)
            elif a.GetAtomicNum()==8 and bond==Chem.BondType.DOUBLE:
                degree=sum(nb.GetAtomicNum()==6 for nb in mol.GetAtomWithIdx(idx).GetNeighbors())
                if degree==1 and mol.GetAtomWithIdx(idx).GetTotalNumHs()==1:inv['aldehyde'].append(idx)
                elif degree==2:inv['ketone'].append(idx)
                else:return None
            else:return None
    # Acids and geminal diols need additional rules; never treat a COOH as OH+oxo.
    carbonyl=set(inv['aldehyde']+inv['ketone'])
    if carbonyl & set(inv['alcohol']+inv['amine']):return None
    if len(inv['alcohol'])>3 or len(inv['amine'])>2 or len(carbonyl)>2 or len(halo)>3:return None
    if len(set(inv['alcohol']))!=len(inv['alcohol']):return None
    return inv,halo


def bare_alkyl(mol,path,attachment,linker):
    if not 1<=len(path)<=4 or attachment not in (path[0],path[-1]):return False
    return all(nb.GetIdx() in path or nb.GetIdx()==linker for idx in path for nb in mol.GetAtomWithIdx(idx).GetNeighbors())


def flat_parse(text,labels=LABELS):
    if not text:return dict(groups=[],multiplicity_ok=True)
    pattern=re.compile(r'((?:n|\d+)(?:,(?:n|\d+))*)-('+ '|'.join(sorted((x for x in MULT.values() if x),key=len,reverse=True))+r')?('+ '|'.join(sorted(labels,key=len,reverse=True))+r')')
    groups=[];offset=0;ok=True
    while offset<len(text):
        match=pattern.match(text,offset)
        if not match:return None
        locs=['N' if value=='n' else int(value) for value in match[1].split(',')]
        if any(isinstance(i,int) and (i<1 or i>30) for i in locs):return None
        amount=next(i for i,value in MULT.items() if value==(match[2] or ''))
        ok &= len(locs)==amount
        groups.extend((i,match[3]) for i in locs)
        offset=match.end()
        if offset<len(text):
            if text[offset]!='-':return None
            offset+=1
            if offset==len(text):return None
    return dict(groups=sorted(groups,key=lambda x:(str(x[0]),x[1])),multiplicity_ok=ok)


def combined_prefix(groups):
    """N-locants and numeric locants share alphabetical substituent citation."""
    grouped=defaultdict(list)
    for loc,label in groups:grouped[label].append(loc)
    chunks=[]
    for label in sorted(grouped):
        locs=sorted(grouped[label],key=lambda x:(0 if isinstance(x,str) else 1,x))
        chunks.append(','.join(map(str,locs))+'-'+MULT[len(locs)]+label)
    return '-'.join(chunks)


def acyl_parent(n,kind):
    if kind=='ester':return {1:'formate',2:'acetate'}.get(n,ALKANES[n][:-1]+'oate')
    return {1:'formamide',2:'acetamide'}.get(n,ALKANES[n][:-1]+'amide')


def parse_acyl(name):
    name=name.strip().lower()
    if name=='dimethylformamide':return dict(kind='amide',n=1,organyl=None,groups=[('N','methyl'),('N','methyl')],multiplicity_ok=True)
    organyl=None
    if ' ' in name:
        organyl,name=name.split(' ',1)
        if organyl not in SUBS.values():return None
    for kind in ('ester','amide'):
        if (kind=='ester')!=(organyl is not None):continue
        for n in range(1,11):
            # Systematic short-chain aliases can describe the structure but
            # retained preferred spellings remain separate checks.
            endings={acyl_parent(n,kind),ALKANES[n][:-1]+('oate' if kind=='ester' else 'amide')}
            for ending in endings:
                if not name.endswith(ending):continue
                parsed=flat_parse(name[:-len(ending)],set(HALO.values())|{'amino','hydroxy','oxo'}|set(SUBS.values()))
                if parsed is not None:return dict(kind=kind,n=n,organyl=organyl,**parsed)
    return None


def acyl_linkage(mol,name,paths):
    cores=[]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()!=6:continue
        oxo=[a for a in atom.GetNeighbors() if a.GetAtomicNum()==8 and a.GetDegree()==1 and mol.GetBondBetweenAtoms(atom.GetIdx(),a.GetIdx()).GetBondType()==Chem.BondType.DOUBLE]
        if len(oxo)!=1:continue
        for linker in atom.GetNeighbors():
            if mol.GetBondBetweenAtoms(atom.GetIdx(),linker.GetIdx()).GetBondType()!=Chem.BondType.SINGLE:continue
            if linker.GetAtomicNum()==8 and linker.GetDegree()==2:kind='ester'
            elif linker.GetAtomicNum()==7 and 1<=linker.GetDegree()<=3:kind='amide'
            else:continue
            others=[a for a in linker.GetNeighbors() if a.GetIdx()!=atom.GetIdx()]
            if any(a.GetAtomicNum()!=6 or mol.GetBondBetweenAtoms(a.GetIdx(),linker.GetIdx()).GetBondType()!=Chem.BondType.SINGLE for a in others):continue
            cores.append((kind,atom.GetIdx(),oxo[0].GetIdx(),linker.GetIdx(),[a.GetIdx() for a in others]))
    if len(cores)!=1:return None
    kind,carbonyl,oxygen,linker,attachments=cores[0]
    scope=dict(kind='acyclic_'+kind,linker_index=linker,carbonyl_index=carbonyl)
    parent=next(p for p in paths if carbonyl in p)
    if carbonyl not in (parent[0],parent[-1]):return None
    if parent[0]!=carbonyl:parent=list(reversed(parent))
    if len(parent)==1 and mol.GetAtomWithIdx(carbonyl).GetTotalNumHs()!=1:return None
    sides=[]
    for idx in attachments:
        side=next(p for p in paths if idx in p)
        if side==parent or not bare_alkyl(mol,side,idx,linker):return unknown_result(scope,'アルコキシ側／N置換基は末端結合の無置換直鎖C1–C4のみ対応します。')
        sides.append(side)
    covered=set(parent)|{oxygen,linker}|{i for p in sides for i in p}
    if len(paths)!=1+len(sides):return None
    inventory=terminal_inventory(mol,parent,{oxygen,linker})
    if inventory is None:return unknown_result(scope,'酸や他のアシル基との複合、側鎖の環・分岐・不飽和、未対応の官能基があるため保留します。')
    inv,halo=inventory
    for idx in parent:
        covered.update(a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors())
    if covered!=set(range(mol.GetNumAtoms())):return None
    groups=[(parent.index(i)+1,PREFIX[k]) for k,indices in inv.items() for i in indices]
    groups.extend((parent.index(i)+1,label) for i,label in halo)
    if len(parent)<=2 and groups:return unknown_result(scope,'置換されたC1/C2の保持アシル名は省略位置番号などの規則を保留します。N上の単純アルキル置換は対応します。')
    organyl=SUBS[len(sides[0])] if kind=='ester' else None
    if kind=='amide':groups.extend(('N',SUBS[len(p)]) for p in sides)
    expected=(organyl+' ' if organyl else '')+combined_prefix(groups)+acyl_parent(len(parent),kind)
    # Explicit preferred example in Blue Book P-66.1.1.3.1.1.
    if kind=='amide' and len(parent)==1 and sorted(groups)==[('N','methyl'),('N','methyl')]:expected='dimethylformamide'
    normalized_groups=sorted(groups,key=lambda x:(str(x[0]),x[1]))
    scope.update(n=len(parent),path=parent,groups=normalized_groups,side_paths=sides,organyl=organyl,expected=expected,
                 candidates=[dict(path=parent,groups=normalized_groups,locants=[1])],
                 selection_order=['エステル／アミド結合の認識','アシル炭素を含む親鎖','アシル炭素を1として番号付け','O側／N上の置換と列挙順'])
    parsed=parse_acyl(name)
    if parsed is None:return unknown_result(scope,'構造は限定範囲ですが、このエステル／アミドの名称構文・別名は未対応です。')
    chapter='P6' if kind=='ester' else 'P6a';reference='P-65.6.3.3.1' if kind=='ester' else 'P-66.1.1 / P-66.1.1.3.1'
    compatible=parsed['kind']==kind and parsed['n']==len(parent)
    checks=[
      row('suffix','エステル／アミド結合と名称形式','pass' if parsed['kind']==kind else 'fail',f"構造は{kind}、名称形式は{parsed['kind']}。エステルはO側アルキル名と酸由来名、アミドはamide接尾辞とN位置番号を検査。",reference,chapter),
      row('parent','アシル炭素を含む親鎖','pass' if parsed['n']==len(parent) else 'fail',f"アシル炭素を含む{len(parent)}炭素。名称は{parsed['n']}炭素。親鎖index {parent}。",'P-44.1 / P-44.3.2','P4'),
      row('attachment','O側アルキル／N置換と結合位置','unknown' if not compatible else ('pass' if parsed['organyl']==organyl and parsed['groups']==normalized_groups else 'fail'),f"O側アルキル: {organyl}、接頭辞とN置換: {normalized_groups}。名称の配置: {parsed['groups']}。",reference,chapter),
      row('numbering','アシル基を基準にした位置番号','unknown' if not compatible else ('pass' if parsed['groups']==normalized_groups else 'fail'),'カルボニル炭素は1。数値位置とN位置を区別して対応を照合します。',reference,chapter),
      row('multiplicity','置換基数と倍数接頭辞','pass' if parsed['multiplicity_ok'] else 'fail','N,N-diや数値位置の個数を倍数接頭辞と比較。','P-14.2'),
      row('spelling','アシル化合物の限定した系統表記','pass' if name.strip().lower()==expected.lower() else 'fail',f"期待表記：{expected}。全保持名・PIN全体の正しさは未判定です。",reference,chapter),
    ]
    return dict(scope=scope,checks=final_rows(checks))


def parse_ether(name):
    name=name.strip().lower()
    for n in range(1,11):
        endings=[(ALKANES[n],None,None)]
        for kind in ('aldehyde','ketone','alcohol','amine'):
            for count in (1,2,3):
                if kind=='aldehyde':
                    if count>2:continue
                    endings.append(((ALKANES[n][:-1] if count==1 else ALKANES[n])+('' if count==1 else MULT[count])+'al',kind,[1] if count==1 else [1,n]))
                else:
                    ending={'ketone':'one','alcohol':'ol','amine':'amine'}[kind]
                    endings.append(((ALKANES[n][:-1] if count==1 else ALKANES[n])+r'-([\d,]+)-'+MULT[count]+ending,kind,count))
        for suffix,kind,locs in endings:
            match=re.fullmatch(r'(.*?)'+suffix,name)
            if not match:continue
            if isinstance(locs,int):locants=list(map(int,match[2].split(',')));ok=len(locants)==locs
            else:locants=locs;ok=True
            text=match[1]
            parsed=flat_parse(text)
            if parsed is None and n<=2 and kind is None and text in ALKOXY.values():parsed=dict(groups=[(1,text)],multiplicity_ok=True)
            if parsed is not None:
                if any(isinstance(i,str) for i,_ in parsed['groups']):continue
                parsed['multiplicity_ok'] &= ok
                return dict(n=n,principal=kind,locants=locants,**parsed)
    return None


def ether_linkage(mol,name,paths):
    bridges=[a for a in mol.GetAtoms() if a.GetAtomicNum()==8 and a.GetDegree()==2 and all(nb.GetAtomicNum()==6 and mol.GetBondBetweenAtoms(a.GetIdx(),nb.GetIdx()).GetBondType()==Chem.BondType.SINGLE for nb in a.GetNeighbors())]
    if len(bridges)!=1 or len(paths)!=2:return None
    bridge=bridges[0];neighbors=[a.GetIdx() for a in bridge.GetNeighbors()]
    if all(len(p)>4 for p in paths):return None
    candidates=[]
    for parent in paths:
        side=next(p for p in paths if p is not parent)
        attachment=next((i for i in neighbors if i in side),None)
        if attachment is None or not bare_alkyl(mol,side,attachment,bridge.GetIdx()):continue
        inventory=terminal_inventory(mol,parent,{bridge.GetIdx()})
        if inventory is None:continue
        inv,halo=inventory
        present=[k for k,v in inv.items() if v]
        principal=min(present,key=PRIORITY.get) if present else None
        if principal=='aldehyde' and len(parent)<3:continue
        for path in (parent,list(reversed(parent))):
            groups=[(path.index(i)+1,PREFIX[k]) for k,indices in inv.items() if k!=principal for i in indices]
            groups.extend((path.index(i)+1,label) for i,label in halo)
            attached=next(i for i in neighbors if i in path)
            groups.append((path.index(attached)+1,ALKOXY[len(side)]))
            locants=sorted(path.index(i)+1 for i in inv[principal]) if principal else []
            candidates.append(dict(path=path,groups=sorted(groups),locants=locants,principal=principal,side_path=side))
    if not candidates:return None
    # A side that contains a characteristic group cannot be named as a bare
    # alkoxy group; those candidates were excluded by bare_alkyl above.
    candidates.sort(key=lambda c:(PRIORITY.get(c['principal'],9),-len(c['path']),c['locants'],group_rank(c['groups'])))
    best=candidates[0];n=len(best['path']);principal=best['principal']
    omit=n<=2 and not principal and len(best['groups'])==1
    suffix=function_suffix(n,principal,best['locants'],True) if principal else ALKANES[n]
    expected=prefix_name(best['groups'],omit_single=omit)+suffix
    scope=dict(kind='acyclic_ether',n=n,principal=principal,path=best['path'],expected=expected,linker_index=bridge.GetIdx(),
               candidates=candidates,preferred_candidates=[best],selection_order=['主特性基を含む炭素鎖','親鎖の炭素数','主特性基位置番号','アルコキシ等の接頭辞位置番号・列挙順'])
    parsed=parse_ether(name)
    if parsed is None:return unknown_result(scope,'このアルコキシ名称構文・別名は未対応です。直鎖の単一エーテルと主鎖上の限定官能基のみ検査します。')
    represented=any(parsed['n']==len(c['path']) and parsed['principal']==c['principal'] and (parsed['locants'] or [])==c['locants'] and parsed['groups']==c['groups'] for c in candidates)
    compatible=parsed['n']==n and parsed['principal']==principal
    checks=[
      row('parent','エーテルの親鎖選択','pass' if parsed['n']==n else 'fail',f"採用する炭素鎖は{n}炭素、名称は{parsed['n']}炭素。主特性基を含む鎖を優先し、次に鎖長を比較。",'P-44.1 / P-44.3.2','P4'),
      row('suffix','アルコキシと主特性基の選択','pass' if parsed['principal']==principal else 'fail',f"エーテル酸素はalkoxy接頭辞。主特性基は{principal or 'なし'}、名称は{parsed['principal'] or 'なし'}。",'P-63.2.4.1 / P-41','P6'),
      row('descriptor','エーテル結合位置と官能基の構造対応','pass' if represented else 'fail',f"名称の接頭辞配置は{parsed['groups']}。比較可能な親鎖・両方向の候補と照合。",'P-63.2.2 / P-63.2.4.1','P6'),
      row('numbering','主特性基・アルコキシ基の番号付け','unknown' if not compatible else ('pass' if (parsed['locants'] or [])==best['locants'] and parsed['groups']==best['groups'] else 'fail'),f"期待する主特性基位置は{best['locants']}、接頭辞位置は{best['groups']}。",'P-14.4'),
      row('multiplicity','位置番号と倍数表記','pass' if parsed['multiplicity_ok'] else 'fail','接頭辞・接尾辞の個数を検査。','P-14.2'),
      row('spelling','アルコキシ化合物の限定した系統表記','pass' if name.strip().lower()==expected else 'fail',f"期待表記：{expected}。全エーテル命名・PIN全体の保証ではありません。",'P-63.2.2 / P-63.2.4.1','P6'),
    ]
    return dict(scope=scope,checks=final_rows(checks))


def linkage_rule_checks(mol,name):
    if mol.GetRingInfo().NumRings():return None
    paths=carbon_components(mol)
    if not paths:return None
    for evaluator in (acyl_linkage,ether_linkage):
        result=evaluator(mol,name,paths)
        if result is not None:return result
    return None
