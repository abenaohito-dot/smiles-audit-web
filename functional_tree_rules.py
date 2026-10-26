"""Bounded branched/unsaturated acyclic functional-carbon trees.
Parent candidates are ranked before unsupported substituents are discarded.
"""
from functools import lru_cache
import re
from rdkit import Chem
from extended_rules import ALKANES,ROOTS,MULT,SUBS,all_paths,group_rank,prefix_name,row,final_rows,unknown_result
from carbon_rules import HALO,MAX_MULTIPLE,parent_text,parse_name as parse_carbon_name
from natural_rules import PRIORITY,PREFIX,function_suffix
from linkage_rules import flat_parse

LABELS=set(SUBS.values())|set(HALO.values())|set(PREFIX.values())


def inventory(mol):
    carbons={a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==6}
    if not carbons:return None
    adj={i:[nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetAtomicNum()==6] for i in carbons}
    if sum(map(len,adj.values()))!=2*(len(carbons)-1):return None
    oh=[];oxo=[];amino=[];halo=[]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==6:continue
        if atom.GetDegree()!=1 or atom.GetNeighbors()[0].GetAtomicNum()!=6:return None
        idx=atom.GetNeighbors()[0].GetIdx();bond=mol.GetBondBetweenAtoms(idx,atom.GetIdx()).GetBondType()
        if atom.GetAtomicNum()==8 and bond==Chem.BondType.SINGLE and atom.GetTotalNumHs()==1:oh.append(idx)
        elif atom.GetAtomicNum()==8 and bond==Chem.BondType.DOUBLE:oxo.append(idx)
        elif atom.GetAtomicNum()==7 and bond==Chem.BondType.SINGLE and atom.GetTotalNumHs()==2:amino.append(idx)
        elif atom.GetAtomicNum() in HALO and bond==Chem.BondType.SINGLE:halo.append((idx,HALO[atom.GetAtomicNum()]))
        else:return None
    acid=[];aldehyde=[];ketone=[]
    for idx in oxo:
        if len(adj[idx])==1 and idx in oh:acid.append(idx);oh.remove(idx)
        elif len(adj[idx])==1 and mol.GetAtomWithIdx(idx).GetTotalNumHs()==1:aldehyde.append(idx)
        elif len(adj[idx])==2:ketone.append(idx)
        else:return None
    inv=dict(acid=acid,aldehyde=aldehyde,ketone=ketone,alcohol=oh,amine=amino)
    present=[k for k in PRIORITY if inv[k]]
    if not present:return None
    if len(oh)>4 or len(amino)>2 or len(oxo)>2 or len(halo)>3:return None
    if len(set(oh))!=len(oh):return None
    multiple=[b for b in mol.GetBonds() if b.GetBeginAtomIdx() in carbons and b.GetEndAtomIdx() in carbons and b.GetBondType()!=Chem.BondType.SINGLE]
    if not multiple and not any(len(v)>2 for v in adj.values()):return None  # Preserve v0.6 straight saturated rules.
    return adj,inv,halo,present[0],multiple


def describe(mol,path,adj,inv,halo,principal):
    parent=set(path);locants=sorted(path.index(i)+1 for i in inv[principal] if i in parent)
    groups=[(path.index(i)+1,PREFIX[k]) for k,indices in inv.items() if k!=principal for i in indices if i in parent]
    groups.extend((path.index(i)+1,label) for i,label in halo if i in parent)
    locations=[loc for loc,_ in groups];supported=True
    decorated={i for indices in inv.values() for i in indices}|{i for i,_ in halo}
    for loc,idx in enumerate(path,1):
        for first in adj[idx]:
            if first in parent:continue
            locations.append(loc);seen=set();stack=[first]
            while stack:
                current=stack.pop()
                if current in parent or current in seen:continue
                seen.add(current);stack.extend(adj[current])
            degrees={i:sum(nb in seen for nb in adj[i]) for i in seen}
            simple=(1<=len(seen)<=10 and degrees[first]<=1 and all(d<=2 for d in degrees.values()) and not seen & decorated
                    and mol.GetBondBetweenAtoms(idx,first).GetBondType()==Chem.BondType.SINGLE
                    and all(b.GetBondType()==Chem.BondType.SINGLE for b in mol.GetBonds() if b.GetBeginAtomIdx() in seen and b.GetEndAtomIdx() in seen))
            if not simple:supported=False
            else:groups.append((loc,SUBS[len(seen)]))
    double=[];triple=[]
    for loc,(left,right) in enumerate(zip(path,path[1:]),1):
        bond=mol.GetBondBetweenAtoms(left,right).GetBondType()
        if bond==Chem.BondType.DOUBLE:double.append(loc)
        elif bond==Chem.BondType.TRIPLE:triple.append(loc)
    return dict(path=list(path),locants=locants,groups=sorted(groups) if supported else None,
                doubles=double,triples=triple,branch_count=len(locations),prefix_locants=sorted(locations))


def rank(candidate):
    d,t=candidate['doubles'],candidate['triples']
    return (-len(candidate['locants']),-len(candidate['path']),-len(d+t),-len(d),tuple(candidate['locants']),
            tuple(sorted(d+t)),tuple(d),-candidate['branch_count'],tuple(candidate['prefix_locants']))


def expected_name(best,principal):
    n=len(best['path']);locants=best['locants'];d=best['doubles'];t=best['triples']
    if not d and not t:suffix=function_suffix(n,principal,locants,bool(best['groups']))
    else:
        hydride=parent_text(n,d,t,substituted=True)
        count=len(locants);stem=hydride[:-1] if count==1 else hydride
        if principal=='acid':suffix=stem+MULT[count]+'oic acid'
        elif principal=='aldehyde':suffix=stem+MULT[count]+'al'
        else:suffix=stem+'-'+','.join(map(str,locants))+'-'+('tetr' if principal=='alcohol' and count==4 else MULT[count])+{'alcohol':'ol','ketone':'one','amine':'amine'}[principal]
    return prefix_name(best['groups'])+suffix


@lru_cache(maxsize=5)
def grammar(kind):
    patterns=[]
    for n in range(1,11):
        for dcount in range(MAX_MULTIPLE+1):
            for tcount in range(MAX_MULTIPLE+1-dcount):
                for count in range(1,5):
                    if kind in ('acid','aldehyde') and count>2:continue
                    if dcount or tcount:
                        first=dcount or tcount
                        core=ROOTS[n]+('a' if first>1 else '')
                        if dcount:core+=r'(?:-(?P<double>\d+(?:,\d+)*)-)?'+MULT[dcount]+'en'+('e' if not tcount and count>1 else '')
                        if tcount:core+=r'(?:-(?P<triple>\d+(?:,\d+)*)-)?'+MULT[tcount]+'yn'+('e' if count>1 else '')
                    else:core=ALKANES[n] if count>1 else ALKANES[n][:-1]
                    if kind in ('acid','aldehyde'):
                        end=MULT[count]+('oic acid' if kind=='acid' else 'al')
                    else:end=r'(?:-(?P<func>\d+(?:,\d+)*)-)?'+('tetr' if kind=='alcohol' and count==4 else MULT[count])+{'alcohol':'ol','ketone':'one','amine':'amine'}[kind]
                    patterns.append((re.compile(r'(?P<prefix>.*?)'+core+end),n,dcount,tcount,count))
    return patterns


def parse_name(name):
    text=name.strip().lower()
    for kind in PRIORITY:
        if not text.endswith({'acid':'oic acid','aldehyde':'al','ketone':'one','alcohol':'ol','amine':'amine'}[kind]):continue
        for pattern,n,dc,tc,count in grammar(kind):
            match=pattern.fullmatch(text)
            if not match:continue
            parsed=flat_parse(match['prefix'],LABELS)
            if parsed is None or any(isinstance(i,str) for i,_ in parsed['groups']):continue
            parts=match.groupdict()
            doubles=list(map(int,parts['double'].split(','))) if dc and parts.get('double') else ([] if not dc else None)
            triples=list(map(int,parts['triple'].split(','))) if tc and parts.get('triple') else ([] if not tc else None)
            locants=([1] if count==1 else [1,n]) if kind in ('acid','aldehyde') else (list(map(int,parts['func'].split(','))) if parts.get('func') else None)
            ok=parsed['multiplicity_ok'] and (doubles is None or len(doubles)==dc) and (triples is None or len(triples)==tc) and (locants is None or len(locants)==count)
            return dict(n=n,principal=kind,locants=locants,doubles=doubles,triples=triples,groups=sorted(parsed['groups']),multiplicity_ok=ok)
    # Narrow observed older saturated alcohol spelling: 2-methyl-2-propanol.
    legacy=re.fullmatch(r'(?P<prefix>.*?)(?P<locs>\d+(?:,\d+)*)-(?P<root>'+ '|'.join(sorted(ROOTS[1:],key=len,reverse=True))+r')an(?P<end>ol|one|amine)',text)
    if legacy:
        prefix=legacy['prefix'].removesuffix('-')
        modern=prefix+legacy['root']+'an-'+legacy['locs']+'-'+legacy['end']
        if modern!=text:
            result=parse_name(modern)
            if result is not None:result['legacy']=True;return result
    # NISPO can express characteristic groups as prefixes on a hydrocarbon.
    def functional_prefix(text,n,doubles,triples,cyclic):
        result=flat_parse(text,LABELS)
        if result is None or any(isinstance(i,str) for i,_ in result['groups']):return None
        result['groups']=sorted(result['groups']);return result
    parent=parse_carbon_name(text,prefix_parser=functional_prefix)
    if parent is not None:return dict(parent,principal=None,locants=[])
    return None


def functional_tree_checks(mol,name):
    if mol.GetRingInfo().NumRings():return None
    found=inventory(mol)
    if found is None:return None
    adj,inv,halo,principal,multiple=found
    scope=dict(kind='functional_carbon_tree',principal=principal,inventory={k:len(v) for k,v in inv.items()})
    if len(multiple)>MAX_MULTIPLE:return unknown_result(scope,'C=C/C#Cは合計4個までの限定対応です。')
    paths=all_paths(adj)
    # Retain every orientation at the maximum suffix coverage before considering
    # nameability. A shorter nameable chain must never hide a senior parent.
    count=max(sum(i in path for i in inv[principal]) for path in paths)
    scope['selection_summary']=dict(principal_count=len(inv[principal]),max_in_one_path=count,unrestricted_longest=max(map(len,paths)),paths_examined=len(paths))
    if count!=len(inv[principal]):return unknown_result(scope,'すべての主特性基を1本の炭素鎖に含められません。鎖外の主特性基を複雑な置換基として命名する規則は保留します。')
    paths=[p for p in paths if sum(i in p for i in inv[principal])==count]
    longest=max(map(len,paths));scope['n']=longest
    scope['selection_summary']['selected_length']=longest
    if longest>10:return unknown_result(scope,'主特性基を含む親鎖はC10までの限定対応です。')
    paths=[p for p in paths if len(p)==longest]
    candidates=[describe(mol,p,adj,inv,halo,principal) for p in paths]
    first=min(map(rank,candidates));competitive=[c for c in candidates if rank(c)==first]
    scope.update(candidates=candidates,selection_order=['主特性基を含む個数','親鎖の炭素数','C=C/C#Cの個数','C=Cの個数','主特性基の位置番号列','不飽和位置番号列・C=C位置','接頭辞の個数','接頭辞の位置番号列・列挙順'])
    if any(c['groups'] is None for c in competitive):return unknown_result(scope,'親鎖選択に残る候補の側鎖が、不飽和・官能基付き・複雑な分岐を含みます。別の名前付け可能な親鎖に置き換えず保留します。')
    secondary=min(group_rank(c['groups'])[1:] for c in competitive)
    preferred=[c for c in competitive if group_rank(c['groups'])[1:]==secondary]
    best=preferred[0]
    if len(best['path'])<3 and (best['doubles'] or best['triples']):return unknown_result(scope,'C2以下の不飽和官能基名は、保持名と位置番号省略の規則を保留します。')
    expected=expected_name(best,principal)
    scope.update(expected=expected,path=best['path'],preferred_candidates=preferred)
    parsed=parse_name(name)
    if parsed is None:return unknown_result(scope,'構造は限定範囲ですが、この系統名構文・保持名・別名は未対応です。')
    scope['name_syntax']='legacy' if parsed.get('legacy') else 'systematic'
    complete=all(parsed[key] is not None for key in ('locants','doubles','triples'))
    same_parent=parsed['n']==longest;same_kind=parsed['principal']==principal
    represented=False
    if complete and same_kind:
        for path in all_paths(adj):
            if len(path)!=parsed['n']:continue
            candidate=describe(mol,path,adj,inv,halo,principal)
            if candidate['groups'] is not None and all(candidate[k]==parsed[k] for k in ('locants','doubles','triples','groups')):represented=True;break
    compatible=complete and same_parent and same_kind
    def conditional(ok,value):return 'unknown' if not ok else ('pass' if value else 'fail')
    checks=[
      row('parent','主特性基を含む親鎖の選択','pass' if same_parent else 'fail',f"主特性基{count}個を含む最長候補は{longest}炭素。名称は{parsed['n']}炭素。単なる最長炭素鎖より、主特性基を含む個数を先に比較。",'P-44.1.1 / P-44.3.2','P4'),
      row('suffix','主特性基と接尾辞','pass' if same_kind else 'fail',f"主特性基は{principal}、名称は{parsed['principal']}。検出した官能基数は{scope['inventory']}。",'P-41 / P-43','P4'),
      row('descriptor','官能基・側鎖・不飽和結合の構造対応',conditional(complete and same_kind,represented),'名称の分解が、入力構造のいずれかの炭素鎖・向きと一致するかを照合。','P-15.1.7'),
      row('numbering','主特性基を優先した位置番号',conditional(compatible,parsed['locants']==best['locants']),f"期待位置は{best['locants']}、名称は{parsed['locants']}。主特性基を不飽和・接頭辞の番号より先に比較。",'P-44.4.1.8','P4'),
      row('unsaturation','C=C/C#Cの位置と同点解消',conditional(compatible,parsed['doubles']==best['doubles'] and parsed['triples']==best['triples']),f"期待C=Cは{best['doubles']}、C#Cは{best['triples']}。主特性基の番号を先に確定。",'P-44.4.1.10','P4'),
      row('functional_prefix','下位官能基・側鎖・ハロゲン接頭辞',conditional(compatible,parsed['groups']==best['groups']),f"期待配置は{best['groups']}、名称は{parsed['groups']}。側鎖は末端結合の飽和直鎖アルキルのみ。",'P-45.2 / P-60.2','P4'),
      row('multiplicity','接頭辞・不飽和・接尾辞の倍数表記','pass' if parsed['multiplicity_ok'] else 'fail','各位置番号の個数とdi/tri等を照合。','P-14.2'),
      row('spelling','分岐・不飽和官能基の限定した系統表記','pass' if name.strip().lower()==expected else 'fail',f"期待表記：{expected}。保持名の優先性・立体命名・PIN全体は未判定です。",'P-61 / P-62 / P-63 / P-64 / P-65','P6'),
    ]
    if not multiple:checks[4].update(status='na',detail='C=C/C#Cのない飽和炭素骨格です。')
    return dict(scope=scope,checks=final_rows(checks))
