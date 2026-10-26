"""Bounded carbon-chain unsaturation and halo-substitution rules, independent of NISPO."""
from collections import Counter
import re
from rdkit import Chem
from extended_rules import (ALKANES, ROOTS, MULT, SUBS, all_paths, row, final_rows,
                            group_rank, unknown_result)

HALO = {9:'fluoro',17:'chloro',35:'bromo',53:'iodo'}
PREFIXES = tuple(SUBS.values()) + tuple(HALO.values())
MAX_MULTIPLE = 4


def prefix_parse(text, n, doubles, triples, cyclic=False):
    if not text:
        return dict(groups=[], multiplicity_ok=True, omitted=False)
    # Unnumbered halo prefixes are supported only when the positions are derivable
    # without consulting the input structure: methane, monosubstituted ethane/
    # ethene/ethyne/cycloalkanes, or complete identical halo substitution.
    multipliers = '|'.join(sorted((p for p in MULT.values() if p),key=len,reverse=True))
    halo = re.fullmatch('('+multipliers+')?('+'|'.join(HALO.values())+')',text)
    if halo:
        count = next(k for k,v in MULT.items() if v==(halo[1] or ''))
        if n==1:
            groups=[(1,halo[2])]*count
        elif count==1 and (cyclic or n==2):
            groups=[(1,halo[2])]
        elif doubles is not None and triples is not None:
            capacities=[2 if cyclic else 4 for _ in range(n)]
            if not cyclic:
                for i in range(n-1):
                    capacities[i]-=1;capacities[i+1]-=1
            for loc in doubles:
                if not 1<=loc<n:return None
                capacities[loc-1]-=1;capacities[loc]-=1
            for loc in triples:
                if not 1<=loc<n:return None
                capacities[loc-1]-=2;capacities[loc]-=2
            if count!=sum(capacities) or min(capacities)<0:
                return None
            groups=[(i+1,halo[2]) for i,c in enumerate(capacities) for _ in range(c)]
        else:
            return None
        return dict(groups=groups,multiplicity_ok=True,omitted=True)
    if n==1:
        pattern=re.compile('('+multipliers+')?('+'|'.join(HALO.values())+')')
        offset=0; groups=[]
        while offset<len(text):
            match=pattern.match(text,offset)
            if not match:return None
            count=next(k for k,v in MULT.items() if v==(match[1] or ''))
            groups.extend([(1,match[2])]*count)
            offset=match.end()
        return dict(groups=sorted(groups),multiplicity_ok=True,omitted=True)
    # A single plain alkyl substituent on an otherwise unsubstituted ring.
    if cyclic and text in SUBS.values():
        return dict(groups=[(1,text)],multiplicity_ok=True,omitted=True)
    pattern=re.compile(r'(\d+(?:,\d+)*)-('+multipliers+r')?('+'|'.join(PREFIXES)+')')
    offset=0; groups=[]; count_ok=True
    while offset<len(text):
        match=pattern.match(text,offset)
        if not match:return None
        locants=list(map(int,match[1].split(',')))
        if any(i<1 or i>30 for i in locants):return None
        count=next(k for k,v in MULT.items() if v==(match[2] or ''))
        count_ok &= count==len(locants)
        groups.extend((loc,match[3]) for loc in locants)
        offset=match.end()
        if offset<len(text):
            if text[offset]!='-':return None
            offset+=1
            if offset==len(text):return None
    return dict(groups=sorted(groups),multiplicity_ok=count_ok,omitted=False)


def parent_text(n,doubles,triples,substituted=False,cyclic=False):
    if not doubles and not triples:
        return ('cyclo' if cyclic else '')+ALKANES[n]
    if cyclic:
        return None  # This release covers only saturated halo-substituted rings.
    if n==2 and not substituted and not doubles and triples==[1]:
        return 'acetylene'
    if (n==2 or n==3 and not substituted) and len(doubles)+len(triples)==1:
        return ROOTS[n]+('ene' if doubles else 'yne')
    first_count=len(doubles) if doubles else len(triples)
    result=ROOTS[n]+('a' if first_count>1 else '')
    if doubles:
        result+='-'+','.join(map(str,doubles))+'-'+MULT[len(doubles)]+'en'+('' if triples else 'e')
    if triples:
        result+='-'+','.join(map(str,triples))+'-'+MULT[len(triples)]+'yne'
    return result


def parse_name(name,cyclic=False,prefix_parser=prefix_parse):
    text=name.strip().lower()
    if not cyclic and text=='acetylene':
        return dict(n=2, doubles=[],triples=[1],groups=[],multiplicity_ok=True,omitted=False)
    for n in range(1,11):
        for dcount in range(MAX_MULTIPLE+1):
            for tcount in range(MAX_MULTIPLE+1-dcount):
                if cyclic and (dcount or tcount):continue
                if not dcount and not tcount:
                    parent=('cyclo' if cyclic else '')+ALKANES[n]
                    if not text.endswith(parent):continue
                    prefix=text[:-len(parent)]
                    double=[];triple=[]
                    syntax_count=True
                else:
                    root=ROOTS[n]+('a' if (dcount or tcount)>1 else '')
                    regex=re.escape(root)
                    if dcount:
                        regex+=r'(?:-(?P<double>\d+(?:,\d+)*)-)?'+MULT[dcount]+'en'+('' if tcount else 'e')
                    if tcount:
                        regex+=r'(?:-(?P<triple>\d+(?:,\d+)*)-)?'+MULT[tcount]+'yne'
                    match=re.fullmatch(r'(?P<prefix>.*?)'+regex,text)
                    if not match:continue
                    prefix=match['prefix']
                    double=list(map(int,match['double'].split(','))) if dcount and match['double'] else ([] if not dcount else None)
                    triple=list(map(int,match['triple'].split(','))) if tcount and match['triple'] else ([] if not tcount else None)
                    # The omitted single-bond locant is unambiguous for any
                    # dinuclear chain, or unsubstituted trinuclear ene/yne.
                    if n==2 or n==3 and not prefix:
                        if double is None and dcount==1:double=[1]
                        if triple is None and tcount==1:triple=[1]
                    syntax_count=(double is None or len(double)==dcount) and (triple is None or len(triple)==tcount)
                parsed=prefix_parser(prefix,n,double,triple,cyclic)
                if parsed is not None:
                    parsed.update(n=n,doubles=double,triples=triple,multiplicity_ok=parsed['multiplicity_ok'] and syntax_count)
                    return parsed
    if not cyclic:
        for parent in ALKANES[1:]:
            if text.endswith('-'+parent):
                prefix=text[:-len(parent)-1]
                if prefix and not prefix.endswith('-'):
                    parsed=parse_name(prefix+parent,prefix_parser=prefix_parser)
                    if parsed is not None:
                        parsed['legacy']=True
                        return parsed
        # Observed NISPO spellings: '2-methyl-1,3-butadien', '1-butene'.
        # Interpret the narrow older locant placement independently, but retain
        # the original text for preferred-spelling comparison.
        old=re.fullmatch(r'(?P<prefix>.*?)(?P<locs>\d+(?:,\d+)*)-(?P<root>'+ '|'.join(sorted(ROOTS[1:],key=len,reverse=True)) + r')(?P<a>a?)(?P<mult>di|tri|tetra)?(?P<end>ene?|yne?)',text)
        if old:
            count=next(k for k,v in MULT.items() if v==(old['mult'] or ''))
            if bool(old['a']) == (count>1):
                prefix=old['prefix'][:-1] if old['prefix'].endswith('-') else old['prefix']
                modern=prefix+old['root']+old['a']+'-'+old['locs']+'-'+(old['mult'] or '')+('ene' if old['end'].startswith('en') else 'yne')
                if modern!=text:
                    parsed=parse_name(modern,prefix_parser=prefix_parser)
                    if parsed is not None:
                        parsed['legacy']=True
                        return parsed
    return None


def describe(mol,path):
    parent=set(path)
    groups=[];locations=[];supported=True
    for loc,idx in enumerate(path,1):
        for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
            first=neighbor.GetIdx()
            if first in parent:continue
            locations.append(loc)
            if neighbor.GetAtomicNum() in HALO:
                groups.append((loc,HALO[neighbor.GetAtomicNum()]));continue
            if mol.GetBondBetweenAtoms(idx,first).GetBondType()!=Chem.BondType.SINGLE:
                supported=False
            seen=set();stack=[first]
            while stack:
                current=stack.pop()
                if current in parent or current in seen:continue
                seen.add(current)
                stack.extend(a.GetIdx() for a in mol.GetAtomWithIdx(current).GetNeighbors())
            if not 1<=len(seen)<=10:
                supported=False;continue
            degrees={i:sum(nb.GetIdx() in seen for nb in mol.GetAtomWithIdx(i).GetNeighbors()) for i in seen}
            if degrees[first]>1 or any(degree>2 for degree in degrees.values()):
                supported=False
            if any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 for i in seen):
                supported=False
            if any(b.GetBeginAtomIdx() in seen and b.GetEndAtomIdx() in seen and b.GetBondType()!=Chem.BondType.SINGLE for b in mol.GetBonds()):
                supported=False
            groups.append((loc,SUBS[len(seen)]))
    doubles=[];triples=[]
    for loc,(a,b) in enumerate(zip(path,path[1:]),1):
        bond=mol.GetBondBetweenAtoms(a,b)
        if bond.GetBondType()==Chem.BondType.DOUBLE:doubles.append(loc)
        if bond.GetBondType()==Chem.BondType.TRIPLE:triples.append(loc)
    return dict(path=list(path),groups=sorted(groups) if supported else None,
                doubles=doubles,triples=triples,branch_count=len(locations),prefix_locants=sorted(locations))


def primary_rank(candidate):
    d,t=candidate['doubles'],candidate['triples']
    return (-len(candidate['path']),-len(d+t),-len(d),tuple(sorted(d+t)),tuple(d),
            -candidate['branch_count'],tuple(candidate['prefix_locants']))


def expected_text(mol,best,cyclic=False):
    grouped=Counter(name for _,name in best['groups'])
    if any(count not in MULT for count in grouped.values()):return None
    n=len(best['path']);groups=best['groups']; omit=n==1
    if groups and len(grouped)==1:
        only=next(iter(grouped))
        # P-14.3.4: mono-substitution of symmetrical parents, mononuclear
        # parents, and complete identical substitution of all hydrogens.
        complete=(only in HALO.values() and len(groups)==sum(4-sum(int(b.GetBondTypeAsDouble()) for b in mol.GetAtomWithIdx(i).GetBonds() if b.GetOtherAtomIdx(i) in set(best['path'])) for i in best['path']))
        omit=(n==1 or len(groups)==1 and (cyclic or n==2) or complete)
    parts=[]
    for label in sorted(grouped):
        locants=[loc for loc,name in groups if name==label]
        lead='' if omit else ','.join(map(str,locants))+'-'
        parts.append(lead+MULT[len(locants)]+label)
    parent=parent_text(n,best['doubles'],best['triples'],bool(groups),cyclic)
    return (('' if omit else '-').join(parts))+parent if parent else None


def carbon_decorations(mol,name):
    """Neutral, unspecified-stereo C/H/F/Cl/Br/I trees and saturated monocyclic parents."""
    if any(a.GetAtomicNum() not in (6,*HALO) for a in mol.GetAtoms()):return None
    halos=[a for a in mol.GetAtoms() if a.GetAtomicNum() in HALO]
    multiple=[b for b in mol.GetBonds() if b.GetBondType() in (Chem.BondType.DOUBLE,Chem.BondType.TRIPLE)]
    if not halos and not multiple:return None  # Preserve the v0.2 saturated-carbon paths.
    if any(a.GetDegree()!=1 or a.GetNeighbors()[0].GetAtomicNum()!=6 or mol.GetBondBetweenAtoms(a.GetIdx(),a.GetNeighbors()[0].GetIdx()).GetBondType()!=Chem.BondType.SINGLE for a in halos):return None
    carbons=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==6]
    if not carbons:return None
    cyclic=bool(mol.GetRingInfo().NumRings())
    if cyclic:
        rings=mol.GetRingInfo().AtomRings()
        if len(rings)!=1 or multiple or any(b.GetIsAromatic() for b in mol.GetBonds()) or not 3<=len(rings[0])<=10:return None
        cycle=list(rings[0]);paths=[]
        for direction in (cycle,list(reversed(cycle))):
            paths.extend(direction[i:]+direction[:i] for i in range(len(cycle)))
        scope=dict(kind='halo_cycloalkane',n=len(cycle))
    else:
        adj={i:[a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetAtomicNum()==6] for i in carbons}
        if sum(map(len,adj.values()))!=2*(len(carbons)-1):return None
        paths=all_paths(adj)
        longest=max(map(len,paths))
        scope=dict(kind='decorated_carbon_chain',n=longest)
        if longest>10:return unknown_result(scope,'親鎖C10まで対応しています。長い鎖の名称は未実装です。')
        paths=[p for p in paths if len(p)==longest]
    if len(multiple)>MAX_MULTIPLE:
        return unknown_result(scope,'今回の不飽和結合の検証は4個までです。')
    candidates=[describe(mol,path) for path in paths]
    best_primary=min(map(primary_rank,candidates))
    competitive=[c for c in candidates if primary_rank(c)==best_primary]
    if any(c['groups'] is None for c in competitive):
        return unknown_result(scope,'親鎖の選択に残る候補に、不飽和・ハロゲン置換・分岐した側鎖があります。その側鎖名の比較は未実装です。')
    best_secondary=min(group_rank(c['groups'])[1:] for c in competitive)
    preferred=[c for c in competitive if group_rank(c['groups'])[1:]==best_secondary]
    best=preferred[0]
    if not cyclic and len(best['path'])>=3 and multiple and len(set(label for _,label in best['groups']))==1 and best['groups'] and best['groups'][0][1] in HALO.values() and all(mol.GetAtomWithIdx(i).GetTotalNumHs()==0 for i in carbons):
        return unknown_result(scope,'完全ハロゲン置換された不飽和鎖の位置番号省略は、この版では未判定です。')
    expected=expected_text(mol,best,cyclic)
    if expected is None:return unknown_result(scope,'接頭辞・不飽和語尾の個数が今回の実装上限を超えています。')
    scope.update(expected=expected,path=best['path'],candidates=candidates,preferred_candidates=preferred,
                 selection_order=['親骨格の長さ','不飽和結合数','二重結合数','不飽和位置番号列','二重結合位置番号','置換基数','置換基位置番号列','名称中の列挙順'])
    parsed=parse_name(name,cyclic)
    if parsed is None:return unknown_result(scope,'構造は対応範囲内ですが、この保持名・複雑な名称構文は未対応です。')
    scope['name_syntax']='legacy' if parsed.get('legacy') else 'systematic'
    same_n=parsed['n']==len(best['path'])
    missing=parsed['doubles'] is None or parsed['triples'] is None
    same_counts=not missing and len(parsed['doubles'])==len(best['doubles']) and len(parsed['triples'])==len(best['triples'])
    same_unsat=not missing and parsed['doubles']==best['doubles'] and parsed['triples']==best['triples']
    same_groups=parsed['groups']==best['groups']
    # Descriptor comparison includes nonpreferred orientations. For chains, use
    # the supplied parent length to recognize an equivalent, shorter parent.
    source_paths=paths if cyclic else [p for p in all_paths(adj) if len(p)==parsed['n']]
    descriptors=[describe(mol,p) for p in source_paths]
    match=not missing and any(c['groups'] is not None and c['groups']==parsed['groups'] and c['doubles']==parsed['doubles'] and c['triples']==parsed['triples'] for c in descriptors)
    parsed_all=None if missing else sorted(parsed['doubles']+parsed['triples'])
    expected_all=sorted(best['doubles']+best['triples'])
    def conditional(condition, ok):return 'unknown' if not condition else ('pass' if ok else 'fail')
    checks=[
        row('parent','親骨格の選択・炭素数','pass' if same_n else 'fail',f"採用する親骨格は{len(best['path'])}炭素、名称は{parsed['n']}炭素。鎖では長さを不飽和結合数より先に比較します。",'P-44.3.2','P4'),
        row('suffix','二重・三重結合の種類と数',conditional(not missing,same_counts),f"期待する二重結合位置は{best['doubles']}、三重結合位置は{best['triples']}。名称はそれぞれ{parsed['doubles']}、{parsed['triples']}。",'P-31.1.1','P3'),
        row('descriptor','名称と構造の分解の対応',conditional(not missing,match),'名称中の親骨格・結合位置・置換基が入力構造の候補と対応するかを検査します。','P-15.1.7'),
        row('numbering','不飽和結合を先に番号付け',conditional(same_n and same_counts and match,parsed_all==expected_all),f"不飽和位置番号列を最初の相違点で比較。採用列は{expected_all}、名称は{parsed_all}。置換基の番号より先に判定します。",'P-14.4(e) / P-31.1.1.1'),
        row('double_tiebreak','同点なら二重結合の番号を優先',conditional(same_n and same_counts and match and parsed_all==expected_all,same_unsat),f"二重結合の採用位置は{best['doubles']}、名称は{parsed['doubles']}。",'P-31.1.1.1','P3'),
        row('prefix_numbering','置換基の位置番号と列挙順',conditional(same_n and same_unsat and match,same_groups),f"採用する置換基配置は{best['groups']}、名称は{parsed['groups']}。不飽和の番号付けを確定してから比較します。",'P-45.2 / P-14.5','P4'),
        row('multiplicity','接頭辞・語尾の個数','pass' if parsed['multiplicity_ok'] else 'fail','di/tri等と位置番号の個数を照合します。','P-14.2'),
        row('spelling','不飽和・ハロゲン置換の限定した優先表記','pass' if name.strip().lower()==expected else 'fail',f"期待表記：{expected}。{'旧式の位置番号・語尾表記として解析しました。' if parsed.get('legacy') else ''}eneをyneより先に記し、連結時の母音省略と接頭辞の列挙順を検査します。",'P-31.1.1 / P-61.3','P3'),
    ]
    if not best['doubles'] and not best['triples'] and same_counts:
        checks[3].update(status='na',detail='飽和骨格のため、不飽和結合の位置番号はありません。')
    if not best['doubles'] or not best['triples']:
        checks[4].update(status='na',detail='二重・三重結合が混在しないため、この同点解消は不要です。')
    if halos:
        checks.insert(2,row('halogen','ハロゲンの接頭辞と個数','pass' if Counter(label for _,label in parsed['groups'] if label in HALO.values())==Counter(HALO[a.GetAtomicNum()] for a in halos) else 'fail','F/Cl/Br/Iをfluoro/chloro/bromo/iodoとして、名前の中の個数と照合します。','P-61.3.1','P6'))
    return dict(scope=scope,checks=final_rows(checks))
