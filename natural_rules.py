"""Bounded mixed functional groups and an explicitly numbered heterocycle registry."""
import re
from rdkit import Chem
from extended_rules import (ALKANES,MULT,SUBS,BASE,row,final_rows,unknown_result,
                            straight_carbon_path,group_rank,prefix_name)
from carbon_rules import HALO,prefix_parse

SUGAR_NOMENCLATURE_NOTE = 'この項目は置換命名法による骨格表記の限定検証です。糖類番号・糖種・D/L・α/βの解析結果は糖類解析欄を参照してください。保持名の優先性・PIN全体は未判定です。'

# Atom-map numbers are IUPAC parent locants, not RDKit atom indices.
PARENTS = {
 'pyridine':'[n:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1',
 'pyridazine':'[n:1]1[n:2][cH:3][cH:4][cH:5][cH:6]1',
 'pyrimidine':'[n:1]1[cH:2][n:3][cH:4][cH:5][cH:6]1',
 'pyrazine':'[n:1]1[cH:2][cH:3][n:4][cH:5][cH:6]1',
 'furan':'[o:1]1[cH:2][cH:3][cH:4][cH:5]1',
 'thiophene':'[s:1]1[cH:2][cH:3][cH:4][cH:5]1',
 'oxolane':'[O:1]1[CH2:2][CH2:3][CH2:4][CH2:5]1',
 'oxane':'[O:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1',
 'pyrrolidine':'[NH:1]1[CH2:2][CH2:3][CH2:4][CH2:5]1',
 'piperidine':'[NH:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1',
 'morpholine':'[O:1]1[CH2:2][CH2:3][NH:4][CH2:5][CH2:6]1',
 'piperazine':'[NH:1]1[CH2:2][CH2:3][NH:4][CH2:5][CH2:6]1',
 '1,3-oxazole':'[o:1]1[cH:2][n:3][cH:4][cH:5]1',
 '1,2-oxazole':'[o:1]1[n:2][cH:3][cH:4][cH:5]1',
 '1,3-thiazole':'[s:1]1[cH:2][n:3][cH:4][cH:5]1',
 '1,2-thiazole':'[s:1]1[n:2][cH:3][cH:4][cH:5]1',
 '1,3-dioxolane':'[O:1]1[CH2:2][O:3][CH2:4][CH2:5]1',
 '1,4-dioxane':'[O:1]1[CH2:2][CH2:3][O:4][CH2:5][CH2:6]1',
 '1,3-thiazolidine':'[S:1]1[CH2:2][NH:3][CH2:4][CH2:5]1',
 'thiomorpholine':'[S:1]1[CH2:2][CH2:3][NH:4][CH2:5][CH2:6]1',
}
TEMPLATES={name:Chem.MolFromSmiles(smi) for name,smi in PARENTS.items()}


def hetero_groups(mol,path):
    ring=set(path);groups=[]
    for loc,idx in enumerate(path,1):
        atom=mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            first=neighbor.GetIdx()
            if first in ring:continue
            if atom.GetAtomicNum()!=6:
                if atom.GetAtomicNum()!=7 or atom.GetIsAromatic() or neighbor.GetAtomicNum()!=6:return None
                # Neutral saturated ring N may carry a simple alkyl group;
                # its position is the registered numeric heterocycle locant.
            if mol.GetBondBetweenAtoms(idx,first).GetBondType()!=Chem.BondType.SINGLE:return None
            if neighbor.GetAtomicNum() in HALO:
                groups.append((loc,HALO[neighbor.GetAtomicNum()]));continue
            seen=set();stack=[first]
            while stack:
                current=stack.pop()
                if current in seen or current in ring:continue
                seen.add(current);stack.extend(a.GetIdx() for a in mol.GetAtomWithIdx(current).GetNeighbors())
            if not 1<=len(seen)<=4 or any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 for i in seen):return None
            degrees={i:sum(nb.GetIdx() in seen for nb in mol.GetAtomWithIdx(i).GetNeighbors()) for i in seen}
            if degrees[first]>1 or any(d>2 for d in degrees.values()):return None
            if any(b.GetBeginAtomIdx() in seen and b.GetEndAtomIdx() in seen and b.GetBondType()!=Chem.BondType.SINGLE for b in mol.GetBonds()):return None
            groups.append((loc,SUBS[len(seen)]))
    return sorted(groups)


def parse_hetero(name):
    name=name.strip().lower()
    # NISPO 0.1.7 was observed to emit '2-methylpyridin'.
    # Preserve the submitted spelling in the spelling check; parse only this
    # narrow legacy form so its structural locants can still be inspected.
    if name.endswith('pyridin'):name+='e'
    for parent in sorted(PARENTS,key=len,reverse=True):
        if name.endswith(parent):
            prefix=name[:-len(parent)]
            if prefix and not prefix[0].isdigit():continue
            if prefix and parent[0].isdigit():
                if not prefix.endswith('-'):continue
                prefix=prefix[:-1]
            parsed=prefix_parse(prefix,len(TEMPLATES[parent].GetAtoms()),[],[])
            if parsed is not None:return dict(parent=parent,**parsed)
    return None


def heterocycle(mol,name):
    rings=mol.GetRingInfo().AtomRings()
    if len(rings)!=1:return None
    ring=set(rings[0])
    if not any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 for i in ring):return None
    for parent,template in TEMPLATES.items():
        if len(ring)!=template.GetNumAtoms():continue
        paths=[]
        for match in mol.GetSubstructMatches(template,uniquify=False,maxMatches=1000):
            if set(match)!=ring:continue
            numbered=[None]*len(match)
            for idx,target in enumerate(match):numbered[template.GetAtomWithIdx(idx).GetAtomMapNum()-1]=target
            if numbered not in paths:paths.append(numbered)
        if not paths:continue
        scope=dict(kind='registered_heterocycle',n=len(ring),parent_name=parent,registry_size=len(PARENTS),
                   registry_reference=BASE+'P2.html',registry_numbering=PARENTS[parent])
        candidates=[dict(path=p,groups=hetero_groups(mol,p)) for p in paths]
        if any(c['groups'] is None for c in candidates):
            return unknown_result(scope,'登録した環骨格ですが、O/Sへの置換、N上のアルキル以外の置換、OH等の官能基、不飽和・分岐した側鎖は未対応です。')
        if len(candidates[0]['groups'])>3:
            return unknown_result(scope,'複素環の置換基は3個までの限定対応です。')
        candidates.sort(key=lambda c:group_rank(c['groups']))
        best=candidates[0]; rank=group_rank(best['groups'])
        preferred=[c for c in candidates if group_rank(c['groups'])==rank]
        expected=prefix_name(best['groups'])+('-' if best['groups'] and parent[0].isdigit() else '')+parent
        scope.update(expected=expected,path=best['path'],candidates=candidates,preferred_candidates=preferred,
                     selection_order=['登録した親骨格の固定番号','等価な向きの置換基位置番号列','名称中の列挙順'])
        parsed=parse_hetero(name)
        if parsed is None:return unknown_result(scope,'この保持名の別名・名称構文は未対応です。登録名と数値の置換基位置番号のみ検査します。')
        same_parent=parsed['parent']==parent
        compatible=same_parent and any(c['groups']==parsed['groups'] for c in candidates)
        correct=parsed['groups']==best['groups']
        checks=[
          row('parent','登録した複素環の親骨格','pass' if same_parent else 'fail',f"入力は{parent}骨格。名称の親骨格は{parsed['parent']}。{len(PARENTS)}種類の有限の登録骨格のみを照合します。",'P-22.2.1 / P-22.2.2','P2'),
          row('hetero_numbering','ヘテロ原子に基づく固定番号','pass' if same_parent else 'unknown',f"位置番号順の原子index: {best['path']}。登録した番号付き骨格に一致する向きだけを比較します。",'P-22.2 / P-14.4(a)','P2'),
          row('descriptor','置換位置・側鎖と環骨格の対応','pass' if compatible else 'fail',f"名称の置換基は{parsed['groups']}。実際の環の各候補と照合しました。",'P-15.1.7'),
          row('numbering','等価な向きの位置番号・列挙順','unknown' if not compatible else ('pass' if correct else 'fail'),f"採用配置は{best['groups']}、名称は{parsed['groups']}。固定番号を保つ等価な向きの中で比較。",'P-14.4 / P-45.2'),
          row('multiplicity','置換基数と倍数接頭辞','pass' if parsed['multiplicity_ok'] else 'fail','位置番号の個数とdi/tri等を照合。','P-14.2'),
          row('spelling','登録骨格の限定した系統表記','pass' if name.strip().lower()==expected else 'fail',f"期待表記：{expected}。全複素環の命名・PINの検証ではありません。",'P-22.2 / P-14.5','P2'),
        ]
        return dict(scope=scope,checks=final_rows(checks))
    return None


PRIORITY={'acid':0,'aldehyde':1,'ketone':2,'alcohol':3,'amine':4}
PREFIX={'aldehyde':'oxo','ketone':'oxo','alcohol':'hydroxy','amine':'amino'}


def mixed_inventory(mol):
    path=straight_carbon_path(mol)
    if path is None or len(path)<2:return None
    if any(mol.GetBondBetweenAtoms(a,b).GetBondType()!=Chem.BondType.SINGLE for a,b in zip(path,path[1:])):return None
    oh=[];oxo=[];amino=[];halo=[]
    for a in mol.GetAtoms():
        if a.GetAtomicNum()==6:continue
        if a.GetDegree()!=1 or a.GetNeighbors()[0].GetAtomicNum()!=6:return None
        idx=a.GetNeighbors()[0].GetIdx();bond=mol.GetBondBetweenAtoms(a.GetIdx(),idx).GetBondType()
        if a.GetAtomicNum()==8 and bond==Chem.BondType.SINGLE and a.GetTotalNumHs()==1:oh.append(idx)
        elif a.GetAtomicNum()==8 and bond==Chem.BondType.DOUBLE:oxo.append(idx)
        elif a.GetAtomicNum()==7 and bond==Chem.BondType.SINGLE and a.GetTotalNumHs()==2:amino.append(idx)
        elif a.GetAtomicNum() in HALO and bond==Chem.BondType.SINGLE:halo.append((idx,HALO[a.GetAtomicNum()]))
        else:return None
    acid=[];aldehyde=[];ketone=[]
    for idx in oxo:
        degree=sum(a.GetAtomicNum()==6 for a in mol.GetAtomWithIdx(idx).GetNeighbors())
        if idx in oh and degree==1:acid.append(idx);oh.remove(idx)
        elif degree==1 and mol.GetAtomWithIdx(idx).GetTotalNumHs()==1:aldehyde.append(idx)
        elif degree==2:ketone.append(idx)
        else:return None
    inventory={'acid':acid,'aldehyde':aldehyde,'ketone':ketone,'alcohol':oh,'amine':amino}
    if len(oh)>6 or len(amino)>2 or len(oxo)>2 or len(halo)>3:return None
    if len(set(oh))!=len(oh) or len(set(amino))!=len(amino):return None
    present=[k for k,v in inventory.items() if v]
    if not present:return None
    # Simple cases stay with v0.3; this evaluator handles mixed groups and
    # multiple carbonyl/amino groups without changing legacy diagnostics.
    extended=(len(present)>1 or len(ketone)>1 or len(aldehyde)>1 or len(amino)>1 or len(acid)>1 or len(oh)>3 or bool(halo))
    if not extended:return None
    principal=min(present,key=PRIORITY.get)
    return path,inventory,halo,principal


def function_suffix(n,kind,locants,substituted):
    stem=ALKANES[n][:-1];count=len(locants); loc=','.join(map(str,locants))
    if kind=='acid':
        if count==1:return stem+'oic acid'
        return ALKANES[n]+MULT[count]+'oic acid'
    if kind=='aldehyde':return stem+'al' if count==1 else ALKANES[n]+MULT[count]+'al'
    if count==1:
        ending={'ketone':'one','alcohol':'ol','amine':'amine'}[kind]
        if n<=2 and not substituted and kind in ('alcohol','amine'):return stem+ending
        return stem+'-'+loc+'-'+ending
    return ALKANES[n]+'-'+loc+'-'+(MULT[count].removesuffix('a') if kind=='alcohol' else MULT[count])+{'ketone':'one','alcohol':'ol','amine':'amine'}[kind]


def parse_mixed(name):
    labels='amino|hydroxy|oxo|fluoro|chloro|bromo|iodo'
    multiple='|'.join(sorted((x for x in MULT.values() if x),key=len,reverse=True))
    for n in range(2,11):
        stem=ALKANES[n][:-1]
        for kind in PRIORITY:
            for count in range(1,7):
                if kind in ('acid','aldehyde'):
                    if count>2:continue
                    ending=(stem if count==1 else ALKANES[n])+(('' if count==1 else MULT[count])+('oic acid' if kind=='acid' else 'al'))
                    if not name.endswith(ending):continue
                    text=name[:-len(ending)];locants=[1] if count==1 else [1,n]
                else:
                    ending=(stem if count==1 else ALKANES[n])+r'(?:-([\d,]+)-)?'+(MULT[count].removesuffix('a') if kind=='alcohol' else MULT[count])+{'ketone':'one','alcohol':'ol','amine':'amine'}[kind]
                    match=re.fullmatch(r'(.*?)'+ending,name)
                    if not match:continue
                    text=match[1];locants=list(map(int,match[2].split(','))) if match[2] else None
                groups=[];offset=0;ok=locants is None or len(locants)==count
                pattern=re.compile(r'(\d+(?:,\d+)*)-('+multiple+r')?('+labels+')')
                while offset<len(text):
                    match=pattern.match(text,offset)
                    if not match:break
                    locs=list(map(int,match[1].split(',')))
                    amount=next(k for k,v in MULT.items() if v==(match[2] or ''))
                    ok &= len(locs)==amount
                    groups.extend((i,match[3]) for i in locs)
                    offset=match.end()
                    if offset<len(text):
                        if text[offset]!='-':break
                        offset+=1
                        if offset==len(text):offset=-1;break
                if offset!=len(text):continue
                return dict(n=n,principal=kind,locants=locants,groups=sorted(groups),multiplicity_ok=ok)
    return None


def mixed_functions(mol,name):
    found=mixed_inventory(mol)
    if found is None:return None
    path,inventory,halo,principal=found;n=len(path)
    scope=dict(kind='mixed_functional_chain',n=n,principal=principal)
    if n<=2 and principal in ('acid','aldehyde'):
        return unknown_result(scope,'短鎖の酸・アルデヒドの保持名と置換規則は今回の複合官能基検証では保留します。')
    orientations=[]
    for direction in (path,list(reversed(path))):
        suffix=sorted(direction.index(i)+1 for i in inventory[principal])
        groups=[(direction.index(i)+1,PREFIX[k]) for k,indices in inventory.items() if k!=principal for i in indices]
        groups.extend((direction.index(i)+1,label) for i,label in halo)
        orientations.append(dict(path=direction,locants=suffix,groups=sorted(groups)))
    orientations.sort(key=lambda c:(c['locants'],group_rank(c['groups'])))
    best=orientations[0]
    expected=prefix_name(best['groups'])+function_suffix(n,principal,best['locants'],bool(best['groups']))
    scope.update(expected=expected,path=best['path'],orientations=orientations,
                 inventory={k:len(v) for k,v in inventory.items()},
                 selection_order=['主特性基の種類','主特性基の位置番号列','下位官能基・ハロゲンの位置番号列','名称中の列挙順'])
    if len(inventory['alcohol'])>=4:
        scope.update(scope_note=SUGAR_NOMENCLATURE_NOTE,sugar_nomenclature_status='unknown')
    parsed=parse_mixed(name.strip().lower())
    if parsed is None:return unknown_result(scope,'構造は対応範囲内ですが、この名称構文・保持名は未対応です。有限の接尾辞とamino/hydroxy/oxo/ハロゲン接頭辞のみ検査します。')
    same_parent=parsed['n']==n; same_kind=parsed['principal']==principal
    compatible=same_parent and same_kind and parsed['locants'] is not None
    represented=compatible and any(c['locants']==parsed['locants'] and c['groups']==parsed['groups'] for c in orientations)
    checks=[
      row('parent','複合官能基を含む親鎖','pass' if same_parent else 'fail',f"無分岐の親鎖は{n}炭素、名称は{parsed['n']}炭素。",'P-44.3.2','P4'),
      row('suffix','主特性基の優先順位','pass' if same_kind else 'fail',f"検出した基の個数は{scope['inventory']}。採用する主特性基は{principal}、名称は{parsed['principal']}。酸→アルデヒド→ケトン→OH→アミンの順で選択。",'P-41 / P-43','P4'),
      row('numbering','主特性基を優先した番号付け','unknown' if not compatible else ('pass' if parsed['locants']==best['locants'] else 'fail'),f"採用する位置番号列は{best['locants']}、名称は{parsed['locants']}。",'P-14.4(c)'),
      row('functional_prefix','下位官能基と接頭辞','unknown' if not compatible else ('pass' if parsed['groups']==best['groups'] else 'fail'),f"採用する接頭辞配置は{best['groups']}、名称は{parsed['groups']}。",'P-60.2 / P-63.1.4','P6'),
      row('descriptor','主特性基・接頭辞の構造対応','unknown' if not compatible else ('pass' if represented else 'fail'),'位置番号と官能基が入力構造のいずれかの向きに対応するかを照合。','P-15.1.7'),
      row('multiplicity','官能基数と倍数表記','pass' if parsed['multiplicity_ok'] else 'fail','接頭辞・接尾辞の倍数と位置番号の個数を照合。','P-14.2'),
      row('spelling','複合官能基の限定した系統表記','pass' if name.strip().lower()==expected else 'fail',f"期待表記：{expected}。アミノ酸等の保持名の優先性やPIN全体は未判定です。",'P-60.2 / P-62 / P-63 / P-64 / P-65','P6'),
    ]
    if scope.get('scope_note'):
        checks.append(row('sugar_nomenclature','糖の同定・糖類命名','unknown',SUGAR_NOMENCLATURE_NOTE,'P-102','P10'))
    return dict(scope=scope,checks=final_rows(checks))


def natural_rule_checks(mol,name):
    for evaluator in (mixed_functions,heterocycle):
        result=evaluator(mol,name)
        if result is not None:return result
    return None
