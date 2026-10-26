"""Independent bounded rule evaluator. Never imports NISPO or uses its decisions."""
from collections import defaultdict
import re
from rdkit import Chem

BASE = 'https://iupac.qmul.ac.uk/BlueBook/'
ALKANES = ('', 'methane', 'ethane', 'propane', 'butane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane', 'decane')
ROOTS = ('', 'meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non', 'dec')
MULT = {1:'', 2:'di', 3:'tri', 4:'tetra', 5:'penta', 6:'hexa', 7:'hepta', 8:'octa', 9:'nona', 10:'deca'}
SUBS = {n: ALKANES[n][:-3] + 'yl' for n in range(1, 11)}


def row(key, title, status, detail, reference, chapter='P1'):
    return dict(id=key, title=title, status=status, detail=detail, reference=reference, url=BASE+chapter+'.html')


def final_rows(rows):
    return rows + [row('stereo', '立体・同位体・PIN全体', 'unknown',
                       'この規則セットは構造・表記の限定検証です。立体命名・同位体命名・PIN全体の適合性は未判定です。', 'P-12.1')]


def eligible(mol):
    return (len(Chem.GetMolFrags(mol)) == 1 and mol.GetNumAtoms() <= 30 and not mol.GetStereoGroups()
            and not any(a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons()
                        or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms())
            and not any(b.GetStereo() != Chem.BondStereo.STEREONONE or b.GetBondType() not in (Chem.BondType.SINGLE,Chem.BondType.DOUBLE,Chem.BondType.TRIPLE,Chem.BondType.AROMATIC) for b in mol.GetBonds()))


def all_paths(adj):
    # Carbon trees only; at most 30 nodes. Include shorter paths for wrong-parent explanations.
    result = []
    def walk(node, path):
        result.append(tuple(path))
        for nb in adj[node]:
            if nb not in path:
                walk(nb, path + [nb])
    for start in adj:
        walk(start, [start])
    return result


def side_groups(mol, path):
    """Return simple terminal unbranched alkyl groups, or None for unsupported sides."""
    parent = set(path)
    groups = []
    for pos, atomid in enumerate(path, 1):
        for neighbor in mol.GetAtomWithIdx(atomid).GetNeighbors():
            first = neighbor.GetIdx()
            if first in parent:
                continue
            seen, stack = set(), [first]
            while stack:
                current = stack.pop()
                if current in seen or current in parent:
                    continue
                seen.add(current)
                stack.extend(a.GetIdx() for a in mol.GetAtomWithIdx(current).GetNeighbors())
            if not 1 <= len(seen) <= 10:
                return None
            # Attachment must be at the end of a linear substituent.
            degrees = {i: sum(nb.GetIdx() in seen for nb in mol.GetAtomWithIdx(i).GetNeighbors()) for i in seen}
            if degrees[first] > 1 or any(d > 2 for d in degrees.values()):
                return None
            groups.append((pos, SUBS[len(seen)]))
    return sorted(groups)


def group_rank(groups):
    named = sorted(groups, key=lambda x:(x[1], x[0]))
    return (tuple(sorted(x[0] for x in groups)), tuple(x[0] for x in named), tuple(x[1] for x in named))


def prefix_name(groups, omit_single=False):
    grouped = defaultdict(list)
    for locant, name in groups:
        grouped[name].append(locant)
    chunks = []
    for name in sorted(grouped):
        locants = sorted(grouped[name])
        if len(locants) not in MULT:
            return None
        lead = '' if omit_single and len(groups) == 1 else ','.join(map(str, locants)) + '-'
        chunks.append(lead + MULT[len(locants)] + name)
    return '-'.join(chunks)


def parse_prefix(text):
    """Only flat simple alkyl prefixes. Retained/compound substituents stay unknown."""
    if not text:
        return [], True
    tokens = re.compile(r'(\d+(?:,\d+)*)-(' + '|'.join(sorted((x for x in MULT.values() if x), key=len, reverse=True)) + r')?(' + '|'.join(SUBS.values()) + ')')
    groups, offset, multiplicity_ok = [], 0, True
    while offset < len(text):
        match = tokens.match(text, offset)
        if not match:
            return None
        locants = list(map(int, match[1].split(',')))
        if any(x < 1 or x > 30 for x in locants):
            return None
        expected_count = next((n for n, p in MULT.items() if p == (match[2] or '')), 1)
        multiplicity_ok &= expected_count == len(locants)
        groups.extend((loc, match[3]) for loc in locants)
        offset = match.end()
        if offset < len(text):
            if text[offset] != '-':
                return None
            offset += 1
            if offset == len(text):
                return None
    return sorted(groups), multiplicity_ok


def parse_hydrocarbon(name, cyclic=False):
    name = name.strip().lower()
    for n in range(1,11):
        suffix = ('cyclo' if cyclic else '') + ALKANES[n]
        if name.endswith(suffix):
            prefix = name[:-len(suffix)]
            # Single ring substituent omits its locant.
            if cyclic and prefix in SUBS.values():
                return dict(n=n, groups=[(1,prefix)], multiplicity_ok=True)
            parsed = parse_prefix(prefix)
            if parsed is not None:
                return dict(n=n, groups=parsed[0], multiplicity_ok=parsed[1])
    return None


def unknown_result(scope, detail):
    return dict(scope=scope, checks=final_rows([row('scope', '今回の規則の適用範囲', 'unknown', detail, 'P-12.1')]))


def hydrocarbon_rows(scope, name, parsed, candidates):
    if parsed is None:
        return unknown_result(scope, '構造は対応範囲内ですが、この名称構文・別名・複雑な置換基名は未対応です。')
    best = candidates[0]
    n = scope['n']
    groups = parsed['groups']
    locset = sorted(x[0] for x in groups)
    best_locset = sorted(x[0] for x in best['groups'])
    same_n = parsed['n'] == n
    signature_match = any(groups == c['groups'] for c in candidates)
    # Matching a non-preferred orientation/path may still describe the same structure.
    valid_match = any(groups == c['groups'] for c in scope.pop('_all_descriptors', candidates) if len(c['path']) == parsed['n'])
    same_count = len(groups) == len(best['groups'])
    compatible = same_n and same_count and valid_match
    rows = [
        row('parent', '親骨格の選択・炭素数', 'pass' if same_n else 'fail',
            f"採用する親骨格は{n}炭素、名称の親骨格は{parsed['n']}炭素。選択した原子index: {best['path']}。", 'P-44.3.2 / P-15.1.7.1.3', 'P4'),
        row('substituent_count', '同長候補の置換基数', ('pass' if same_count else 'fail') if same_n else 'unknown',
            f"対象候補の最大置換基数は{len(best['groups'])}、名称の置換基数は{len(groups)}。", 'P-45.2.1', 'P4'),
        row('descriptor', '親骨格・置換基の構造対応', 'pass' if valid_match else 'fail',
            '名称中の位置と置換基が入力構造の分解の一つに対応します。' if valid_match else '対応する親骨格・置換基の分解が入力構造に見つかりません。', 'P-15.1.7.1'),
        row('numbering', '置換基の位置番号・最初の相違点', ('pass' if locset == best_locset else 'fail') if compatible else 'unknown',
            f"候補の最小位置番号列は{best_locset}、名称は{locset}。合計ではなく、昇順に並べて最初の相違点で比較します。", 'P-45.2.2 / P-14.3.5', 'P4'),
        row('citation_order', '列挙順による同点の解消', ('pass' if signature_match else 'fail') if compatible and locset == best_locset else 'unknown',
            f"採用する配置は{best['groups']}。位置番号列が同じなら、名称中の列挙順で番号を比較します。", 'P-45.2.3 / P-45.5', 'P4'),
        row('multiplicity', '置換基の個数と倍数接頭辞', 'pass' if parsed['multiplicity_ok'] else 'fail',
            '位置番号の個数とdi/tri等の数を照合します。', 'P-14.2 / P-15.1.7.1.4'),
        row('spelling', '置換基の列挙・限定した優先表記', 'pass' if name.strip().lower() == scope['expected'] else 'fail',
            f"この対応範囲の期待表記：{scope['expected']}。単純置換基は倍数接頭辞を除いてアルファベット順に列挙します。", 'P-14.5 / P-15.1.7.1.5'),
    ]
    return dict(scope=scope, checks=final_rows(rows))


def branched_alkane(mol, name):
    if mol.GetRingInfo().NumRings() or any(a.GetAtomicNum()!=6 for a in mol.GetAtoms()) or any(b.GetBondType()!=Chem.BondType.SINGLE for b in mol.GetBonds()):
        return None
    if not any(a.GetDegree()>2 for a in mol.GetAtoms()):
        return None  # Preserve the original straight-chain evaluator.
    adj = {a.GetIdx():[b.GetIdx() for b in a.GetNeighbors()] for a in mol.GetAtoms()}
    paths = all_paths(adj)
    longest = max(map(len, paths))
    scope = dict(kind='branched_alkane', n=longest, total_carbons=mol.GetNumAtoms())
    if longest > 10:
        return unknown_result(scope, '親鎖C10までに対応しています。より長い親鎖名は未実装です。')
    longest_paths = [p for p in paths if len(p)==longest]
    # Do not discard an unsupported competing parent: that could falsely select a name.
    decomposed = [dict(path=list(p), groups=side_groups(mol,p), branch_count=sum(nb.GetIdx() not in p for i in p for nb in mol.GetAtomWithIdx(i).GetNeighbors())) for p in longest_paths]
    count = max(c['branch_count'] for c in decomposed)
    competitive = [c for c in decomposed if c['branch_count']==count]
    if any(c['groups'] is None for c in competitive):
        return unknown_result(scope, '最長鎖候補に、内部結合・分岐した置換基があります。複雑な置換基同士の選択規則は未実装です。')
    candidates = sorted(competitive, key=lambda c:group_rank(c['groups']))
    best_rank = group_rank(candidates[0]['groups'])
    preferred = [c for c in candidates if group_rank(c['groups']) == best_rank]
    prefix = prefix_name(preferred[0]['groups'])
    if prefix is None:
        return unknown_result(scope, '同一置換基の個数が接頭辞の実装上限を超えています。')
    expected = prefix + ALKANES[longest]
    # Collect descriptors for the supplied parent length, including shorter, nonpreferred chains.
    parsed = parse_hydrocarbon(name)
    requested = parsed['n'] if parsed else longest
    all_desc = [dict(path=list(p),groups=g) for p in paths if len(p)==requested and (g:=side_groups(mol,p)) is not None]
    scope.update(expected=expected, path=preferred[0]['path'], candidates=decomposed, preferred_candidates=preferred, _all_descriptors=all_desc)
    return hydrocarbon_rows(scope, name, parsed, preferred)


def simple_ring(mol, name):
    rings = mol.GetRingInfo().AtomRings()
    if len(rings)!=1 or any(a.GetAtomicNum()!=6 for a in mol.GetAtoms()):
        return None
    # Only unsubstituted benzene is covered among aromatic rings.
    if all(a.GetIsAromatic() for a in mol.GetAtoms()) and mol.GetNumAtoms()==6:
        scope = dict(kind='benzene', n=6, expected='benzene')
        if name.strip().lower()!='benzene':
            return unknown_result(scope, 'ベンゼンの保持名のみを検査します。別の環命名構文は未対応です。')
        return dict(scope=scope, checks=final_rows([row('retained_name','環の保持名','pass','芳香族炭素6員環の保持名はbenzene。','P-22.1.3','P2')]))
    if any(b.GetBondType()!=Chem.BondType.SINGLE for b in mol.GetBonds()):
        return None
    cycle = list(rings[0]); n = len(cycle)
    if not 3<=n<=10:
        return None
    candidates = []
    for direction in (cycle,list(reversed(cycle))):
        for offset in range(n):
            path = direction[offset:] + direction[:offset]
            groups = side_groups(mol,path)
            if groups is None:
                return unknown_result(dict(kind='cycloalkane',n=n), '環の側鎖は末端で結合した直鎖アルキルのみ対応しています。')
            candidates.append(dict(path=path,groups=groups))
    candidates.sort(key=lambda c:group_rank(c['groups']))
    best_rank = group_rank(candidates[0]['groups'])
    preferred = [c for c in candidates if group_rank(c['groups'])==best_rank]
    prefix = prefix_name(preferred[0]['groups'],omit_single=True)
    if prefix is None:
        return None
    scope = dict(kind='cycloalkane', n=n, expected=prefix+'cyclo'+ALKANES[n], path=preferred[0]['path'], candidates=candidates, preferred_candidates=preferred, _all_descriptors=candidates)
    result = hydrocarbon_rows(scope,name,parse_hydrocarbon(name,True),preferred)
    for r in result['checks']:
        if r['id']=='parent':
            r.update(title='単環の親骨格・環員数', detail=f"炭素{n}員の単環を親骨格に選択。原子index: {scope['path']}。", reference='P-22.1.1 / P-44.1.2.2', url=BASE+'P2.html')
        if r['id']=='substituent_count':
            r.update(status='na', detail='環の向きによって置換基数は変わりません。')
        if r['id']=='numbering' and not preferred[0]['groups']:
            r.update(status='na',detail='無置換の単環には置換基の位置番号はありません。')
        if r['id']=='spelling':
            r['detail'] += ' 単一置換基だけなら環の位置番号1を省略します。'
    return result


def straight_carbon_path(mol):
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==6]
    if not 1<=len(carbons)<=10 or mol.GetRingInfo().NumRings():
        return None
    adj = {i:[nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetAtomicNum()==6] for i in carbons}
    if any(len(v)>2 for v in adj.values()) or sum(map(len,adj.values()))!=2*(len(carbons)-1):
        return None
    start = next((i for i,v in adj.items() if len(v)<=1),None)
    if start is None:
        return None
    path, previous, current = [], None, start
    while current is not None:
        path.append(current)
        nexts = [v for v in adj[current] if v!=previous]
        previous,current = current, nexts[0] if nexts else None
    return path if len(path)==len(carbons) else None


def unbranched_unsaturated(mol,name):
    if any(a.GetAtomicNum()!=6 for a in mol.GetAtoms()):
        return None
    path = straight_carbon_path(mol)
    if path is None:
        return None
    multiple = [b for b in mol.GetBonds() if b.GetBondType() in (Chem.BondType.DOUBLE,Chem.BondType.TRIPLE)]
    if len(multiple)!=1:
        return None
    bond = multiple[0]; n = len(path)
    pos = min(path.index(bond.GetBeginAtomIdx()),path.index(bond.GetEndAtomIdx()))+1
    preferred = min(pos,n-pos)
    ending = 'ene' if bond.GetBondType()==Chem.BondType.DOUBLE else 'yne'
    expected = ROOTS[n] + (ending if n<=3 else f'-{preferred}-{ending}')
    scope = dict(kind='mono_unsaturated_chain',n=n,expected=expected,path=path,bond_locants=[pos,n-pos])
    parsed = None
    for count in range(2,11):
        for suffix in ('ene','yne'):
            match = re.fullmatch(re.escape(ROOTS[count])+r'(?:-(\d+)-)?'+suffix,name.strip().lower())
            if match:
                loc = int(match[1]) if match[1] else (1 if count<=3 else None)
                parsed = dict(n=count,ending=suffix,locant=loc)
    if parsed is None:
        return unknown_result(scope,'単一の二重結合または三重結合を持つ無分岐鎖の系統名のみ対応しています。')
    loc = parsed['locant']
    rows = [
        row('parent','親鎖の炭素数','pass' if n==parsed['n'] else 'fail',f"構造は{n}炭素。名称は{parsed['n']}炭素。",'P-44.3.2','P4'),
        row('suffix','不飽和結合と語尾','pass' if ending==parsed['ending'] else 'fail',f"結合に対応する語尾は-{ending}。",'P-31.1.1','P3'),
        row('numbering','不飽和結合の位置番号','unknown' if loc is None or n!=parsed['n'] else ('pass' if loc==preferred else 'fail'),f"両端からの位置番号は{[pos,n-pos]}。採用値は{preferred}、名称の値は{loc}。",'P-14.4 / P-15.1.7.2.2'),
        row('spelling','不飽和鎖の限定した優先表記','pass' if name.strip().lower()==expected else 'fail',f"期待表記：{expected}。ethene/propene/ethyne/propyneは位置番号を省略。",'P-14.3 / P-31.1.1','P3'),
    ]
    return dict(scope=scope,checks=final_rows(rows))


def parse_functional(name):
    text = name.strip().lower()
    if text in ('formic acid','acetic acid'):
        return dict(n=1 if text=='formic acid' else 2,kind='acid',locants=[1],hydroxy=[])
    for n in range(1,11):
        alkane, stem = ALKANES[n], ALKANES[n][:-1]
        prefix = r'(?:(\d+)-hydroxy)?'
        for ending,kind in ((r'oic acid','acid'),(r'al','aldehyde')):
            match = re.fullmatch(prefix+re.escape(stem)+ending,text)
            if match:
                return dict(n=n,kind=kind,locants=[1],hydroxy=[int(match[1])] if match[1] else [])
        for ending,kind in (('one','ketone'),('amine','amine'),('ol','alcohol')):
            match = re.fullmatch(prefix+re.escape(stem)+r'(?:-(\d+)-)?'+ending,text)
            if match:
                loc = int(match[2]) if match[2] else (1 if kind in ('amine','alcohol') and n<=2 else None)
                return dict(n=n,kind=kind,locants=[loc] if loc is not None else None,hydroxy=[int(match[1])] if match[1] else [])
        for count in (2,3):
            match = re.fullmatch(re.escape(alkane)+r'(?:-([\d,]+)-)?'+MULT[count]+'ol',text)
            if match:
                locants = list(map(int,match[1].split(','))) if match[1] else None
                return dict(n=n,kind='polyol',locants=locants,hydroxy=[],multiplicity_ok=locants is None or len(locants)==count)
    return None


def functional_chain(mol,name):
    path = straight_carbon_path(mol)
    if path is None or any(a.GetAtomicNum() not in (6,7,8) for a in mol.GetAtoms()):
        return None
    if any(mol.GetBondBetweenAtoms(a,b).GetBondType()!=Chem.BondType.SINGLE for a,b in zip(path,path[1:])):
        return None
    oh, carbonyl, amine = [], [], []
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z==6:
            continue
        if atom.GetDegree()!=1 or atom.GetNeighbors()[0].GetAtomicNum()!=6:
            return None
        carbon = atom.GetNeighbors()[0].GetIdx()
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(),carbon)
        if z==8 and bond.GetBondType()==Chem.BondType.SINGLE and atom.GetTotalNumHs()==1:
            oh.append(carbon)
        elif z==8 and bond.GetBondType()==Chem.BondType.DOUBLE:
            carbonyl.append(carbon)
        elif z==7 and bond.GetBondType()==Chem.BondType.SINGLE and atom.GetTotalNumHs()==2:
            amine.append(carbon)
        else:
            return None
    n = len(path)
    if len(carbonyl)==1 and not amine:
        center = carbonyl[0]
        carbon_degree = sum(a.GetAtomicNum()==6 for a in mol.GetAtomWithIdx(center).GetNeighbors())
        if center in oh and carbon_degree<=1:
            kind='acid'; oh.remove(center)
            if n<=2 and oh:
                return None  # Retained-name substitution conventions not implemented.
        elif carbon_degree==1 and n>=3:
            kind='aldehyde'
        elif carbon_degree==2:
            kind='ketone'
        else:
            return None
        if len(oh)>1:
            return None
        principal=[center]
    elif not carbonyl and not amine and 2<=len(oh)<=3 and len(set(oh))==len(oh) and n>=2:
        kind='polyol'; principal=oh; oh=[]
    elif not carbonyl and len(amine)==1 and not oh:
        kind='amine'; principal=amine
    else:
        return None
    orientations=[]
    for direction in (path,list(reversed(path))):
        locants=sorted(direction.index(i)+1 for i in principal)
        hydroxyl=sorted(direction.index(i)+1 for i in oh)
        orientations.append(dict(path=direction,locants=locants,hydroxy=hydroxyl))
    orientations.sort(key=lambda c:(c['locants'],c['hydroxy']))
    best=orientations[0]; stem=ALKANES[n][:-1]
    if kind=='acid':
        suffix={1:'formic acid',2:'acetic acid'}.get(n,stem+'oic acid')
    elif kind=='aldehyde':
        suffix=stem+'al'
    elif kind=='ketone':
        suffix=stem+f"-{best['locants'][0]}-one"
    elif kind=='amine':
        suffix=stem+('amine' if n<=2 else f"-{best['locants'][0]}-amine")
    else:
        suffix=ALKANES[n]+'-'+','.join(map(str,best['locants']))+'-'+MULT[len(principal)]+'ol'
    expected=(f"{best['hydroxy'][0]}-hydroxy" if best['hydroxy'] else '')+suffix
    scope=dict(kind='functional_'+kind,n=n,expected=expected,path=best['path'],orientations=orientations,principal=kind)
    parsed=parse_functional(name)
    if parsed is None:
        return unknown_result(scope,'この名称構文・保持名・別名は官能基検証器が未対応です。')
    same_n=parsed['n']==n; same_kind=parsed['kind']==kind
    rows=[
        row('parent','官能基を含む親鎖の炭素数','pass' if same_n else 'fail',f"構造は{n}炭素、名称は{parsed['n']}炭素。鎖の原子index: {best['path']}。",'P-44.3.2','P4'),
        row('suffix','主特性基と接尾辞の選択','pass' if same_kind else 'fail',f"構造の主特性基は{kind}、名称の接尾辞は{parsed['kind']}。この範囲ではカルボン酸／アルデヒド／ケトンをOHより優先します。",'P-41 / P-43（主特性基の選択）','P4'),
        row('numbering','主特性基を優先した位置番号','unknown' if not same_n or not same_kind or parsed['locants'] is None else ('pass' if parsed['locants']==best['locants'] else 'fail'),f"両方向の候補：{orientations}。採用する主特性基位置は{best['locants']}、名称は{parsed['locants']}。",'P-14.4 / P-15.1.7.2'),
        row('functional_prefix','下位のOH基と接頭辞','pass' if parsed['hydroxy']==best['hydroxy'] else 'fail',f"主特性基を番号付けしてからOHをhydroxy接頭辞で表します。期待位置は{best['hydroxy']}、名称は{parsed['hydroxy']}。",'P-63.1.4','P6'),
        row('spelling','官能基名の限定した優先表記','pass' if name.strip().lower()==expected else 'fail',f"期待表記：{expected}。未対応の保持名・別名を一般的に無効とは判定しません。",'P-62.2 / P-63.1 / P-64 / P-65','P6'),
    ]
    if 'multiplicity_ok' in parsed:
        rows.append(row('multiplicity','OH基の数と倍数接尾辞','pass' if parsed['multiplicity_ok'] else 'fail','diol/triolの数と位置番号の個数を照合。','P-14.2'))
    return dict(scope=scope,checks=final_rows(rows))


def extended_rule_checks(mol,name):
    mol=Chem.RemoveHs(Chem.Mol(mol))
    from carbohydrate_rules import carbohydrate_name_checks
    carbohydrate=carbohydrate_name_checks(mol,name)
    if carbohydrate is not None:return carbohydrate
    from ring_stereo_rules import ring_stereo_checks
    stereo_ring=ring_stereo_checks(mol,name)
    if stereo_ring is not None:return stereo_ring
    if not eligible(mol):
        return None
    from oxygen_ring_rules import oxygen_ring_checks
    oxygen_ring=oxygen_ring_checks(mol,name)
    if oxygen_ring is not None:return oxygen_ring
    from aromatic_rules import aromatic_rule_checks
    aromatic=aromatic_rule_checks(mol,name)
    if aromatic is not None:return aromatic
    from linkage_rules import linkage_rule_checks
    linked=linkage_rule_checks(mol,name)
    if linked is not None:return linked
    from functional_tree_rules import functional_tree_checks
    functional=functional_tree_checks(mol,name)
    if functional is not None:return functional
    from natural_rules import natural_rule_checks
    result=natural_rule_checks(mol,name)
    if result is not None:return result
    from carbon_rules import carbon_decorations
    for evaluator in (carbon_decorations,branched_alkane,unbranched_unsaturated,functional_chain,simple_ring):
        result=evaluator(mol,name)
        if result is not None:
            # Private working keys are never part of a saved report.
            if result['scope']:
                result['scope'].pop('_all_descriptors',None)
            return result
    return None
