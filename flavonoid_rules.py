"""Bounded flavonoid aglycones, with an extensible, separate scaffold registry.

Numbering: IUPAC flavonoid recommendations Flv-2.1.1.2, Flv-3.2.1.
The registry detects substructures; only complete neutral OH/OMe registered cyclic aglycones
are eligible for a whole-molecule expected name. No trivial-name lookup.
"""
import re
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from enone_rules import CHALCONE_CORE,AURONE_CORE,analyze_enone,enone_checks

SOURCE='https://iupac.qmul.ac.uk/flavonoid/'
CORE='[O:40]=[c:4]1[c:3][c:2](-[c:101]2[c:102][c:103][c:104][c:105][c:106]2)[o:1][c:9]2[c:8][c:7][c:6][c:5][c:10]12'
SATURATED_CORE=CORE.replace('[c:4]','[C:4]').replace('[c:3]','[C:3]').replace('[c:2]','[CH:2]').replace('[o:1]','[O:1]')
B_RING='[c:101]2[c:102][c:103][c:104][c:105][c:106]2'
ISO_CORE='[O:40]=[c:4]1[c:3](-'+B_RING+')[c:2][o:1][c:9]2[c:8][c:7][c:6][c:5][c:10]12'
NEO_CORE='[O:40]=[c:2]1[o:1][c:9]2[c:8][c:7][c:6][c:5][c:10]2[c:4](-'+B_RING+')[c:3]1'
PARENTS={'flavone':(2,4,4),'flavanone':(2,4,4),'isoflavone':(3,4,4),'neoflavone':(4,2,2)}
FAMILIES={
 'flavanone':dict(label='フラバノン系',smarts=SATURATED_CORE,naming=True),
 'flavone':dict(label='フラボン系',smarts=CORE,naming=True),
 'isoflavone':dict(label='イソフラボン系',smarts=ISO_CORE,naming=True),
 'neoflavone':dict(label='ネオフラボン系',smarts=NEO_CORE,naming=True),
 'chalcone':dict(label='カルコン系',smarts=CHALCONE_CORE,naming=True),
 'aurone':dict(label='オーロン系',smarts=AURONE_CORE,naming=True),
}
QUERIES={key:Chem.MolFromSmarts(value['smarts']) for key,value in FAMILIES.items()}
LOCANTS={1:'1',2:'2',3:'3',4:'4',5:'5',6:'6',7:'7',8:'8',9:'8a',10:'4a',**{100+i:str(i)+"'" for i in range(1,7)}}
ALLOWED={3,5,6,7,8,102,103,104,105,106}
MULT={1:'',2:'di',3:'tri',4:'tetra',5:'penta',6:'hexa',7:'hepta',8:'octa',9:'nona',10:'deca'}
NOTE='限定した非配糖体のフラボン系・フラバノン系・イソフラボン・ネオフラボン骨格とOH・メトキシ置換を検査します。絶対R/SはRDKitのCIP計算に依存し、規則自体の独立証明ではありません。慣用名の推定、他の置換基、配糖体、他骨格の全体命名、PIN全体は対象外です。'


def prefix(groups):
 parts=[]
 for label in ('hydroxy','methoxy'):
  locs=sorted(groups.get(label,[]))
  if locs:parts.append(','.join(map(str,locs))+'-'+MULT[len(locs)]+label)
 return '-'.join(parts)


def expected_name(a,b,attachment=2,carbonyl=4,indicated_h=4):
 terms=[]
 for label in ('hydroxy','methoxy'):
  locs=sorted(a.get(label,[]))
  if locs:terms.append((label,','.join(map(str,locs))+'-'+MULT[len(locs)]+label))
 bp=prefix(b);phenyl=str(attachment)+'-('+bp+'phenyl)' if bp else str(attachment)+'-phenyl'
 terms.append(('phenyl',phenyl))
 return '-'.join(term for _,term in sorted(terms))+f'-{indicated_h}H-chromen-{carbonyl}-one'


def substitution(mol,mapping,family='flavone'):
 attachment,carbonyl,indicated_h=PARENTS[family]
 allowed=({2,3,4,5,6,7,8}-{attachment,carbonyl})|{102,103,104,105,106}
 core=set(mapping.values());used=set(core);groups={'hydroxy':[],'methoxy':[]};details=[]
 for num,idx in mapping.items():
  extra=[a for a in mol.GetAtomWithIdx(idx).GetNeighbors() if a.GetIdx() not in core]
  if not extra:continue
  if num not in allowed or len(extra)!=1:return None
  o=extra[0];oi=o.GetIdx()
  if o.GetAtomicNum()!=8 or o.GetIsAromatic() or mol.GetBondBetweenAtoms(idx,oi).GetBondType()!=Chem.BondType.SINGLE:return None
  branch=[a for a in o.GetNeighbors() if a.GetIdx()!=idx]
  if not branch and o.GetTotalNumHs()==1:
   kind='hydroxy';atoms=[oi]
  elif len(branch)==1 and o.GetTotalNumHs()==0:
   c=branch[0]
   if c.GetAtomicNum()!=6 or c.GetDegree()!=1 or c.GetTotalNumHs()!=3 or mol.GetBondBetweenAtoms(oi,c.GetIdx()).GetBondType()!=Chem.BondType.SINGLE:return None
   kind='methoxy';atoms=[oi,c.GetIdx()]
  else:return None
  if used.intersection(atoms):return None
  used.update(atoms);groups[kind].append(num)
  details.append(dict(position=LOCANTS[num],kind=kind,ring_atom=idx,atoms=atoms))
 if used!=set(range(mol.GetNumAtoms())):return None
 a={k:sorted(n for n in ns if n<100) for k,ns in groups.items()}
 b={k:sorted(n-100 for n in ns if n>100) for k,ns in groups.items()}
 rank=(tuple(sorted(b['hydroxy']+b['methoxy'])),tuple(b['hydroxy']),tuple(b['methoxy']))
 return dict(a=a,b=b,substituents=sorted(details,key=lambda d:(int(d['position'].rstrip("'")),"'" in d['position'])),rank=rank,expected=expected_name(a,b,attachment,carbonyl,indicated_h))


def analyze_aglycones(mol):
 # Work on a copy; ordinary explicit H does not change scaffold eligibility.
 mol=Chem.Mol(mol)
 for atom in mol.GetAtoms():atom.SetIntProp('_flv_original_idx',atom.GetIdx())
 mol=Chem.RemoveHs(mol)
 hits=[]
 for family,q in QUERIES.items():
  matches=mol.GetSubstructMatches(q,uniquify=False,maxMatches=256)
  cores={frozenset(m) for m in matches}
  if cores:hits.append(dict(family=family,label=FAMILIES[family]['label'],core_count=len(cores)))
 result=dict(detected=bool(hits),candidates=hits,family=None,label=None,scope='unsupported',expected_name=None,numbering={},equivalent_numberings=[],substituents=[],reason='',note=NOTE,source=SOURCE)
 if len(hits)!=1 or hits[0]['core_count']!=1:
  result['reason']='複数の骨格候補があるため全体名を保留。' if hits else '登録した骨格候補なし。'
  return result
 family=hits[0]['family'];result.update(family=family,label=hits[0]['label'])
 if family in ('chalcone','aurone'):return analyze_enone(mol,family,result)
 if family not in PARENTS:
  result['reason']='骨格候補のみ。番号・置換・全体名の検証は次段階。';return result
 if len(Chem.GetMolFrags(mol))!=1 or mol.GetStereoGroups() or any(a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() or (family!='flavanone' and a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED) for a in mol.GetAtoms()) or any(b.GetStereo()!=Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
  result['reason']='塩・同位体・電荷・指定立体などがあるため全体名を保留。';return result
 q=QUERIES[family];options=[]
 for match in mol.GetSubstructMatches(q,uniquify=False,maxMatches=256):
  mapping={a.GetAtomMapNum():match[a.GetIdx()] for a in q.GetAtoms()}
  found=substitution(mol,mapping,family)
  if found:options.append((found,mapping))
 if not options:
  result['reason']='OH・メトキシ以外の置換、追加環、配糖体などは全体命名の範囲外。';return result
 rank=min(v['rank'] for v,_ in options);options=[(v,m) for v,m in options if v['rank']==rank]
 found,mapping=min(options,key=lambda pair:tuple(pair[1][n] for n in sorted(pair[1])))
 result.update(scope='bounded',family='flavonol' if 3 in found['a']['hydroxy'] else 'flavone',label='フラボノール' if 3 in found['a']['hydroxy'] else 'フラボン',expected_name=found['expected'],numbering={LOCANTS[n]:idx for n,idx in mapping.items() if n in LOCANTS},carbonyl_oxygen=mapping[40],substituents=found['substituents'],a_groups=found['a'],b_groups=found['b'],reason='中性・非配糖体、置換基はOH・メトキシのみ。')
 attachment,carbonyl,indicated_h=PARENTS[family]
 result.update(parent_attachment=attachment,parent_carbonyl=carbonyl,parent_indicated_h=indicated_h,parent_name=f'{attachment}-phenyl'+('chroman' if family=='flavanone' else 'chromen')+f'-{carbonyl}-one')
 if family in ('isoflavone','neoflavone'):result.update(family=family,label=FAMILIES[family]['label'].removesuffix('系'))
 if family=='flavanone':
  result.update(family='flavanonol' if 3 in found['a']['hydroxy'] else 'flavanone',label='フラバノノール' if 3 in found['a']['hydroxy'] else 'フラバノン',saturation='saturated')
  centers=dict(Chem.FindMolChiralCenters(mol,includeUnassigned=True,useLegacyImplementation=False))
  permitted={mapping[2],mapping[3]}
  if any(idx not in permitted for idx in centers):
   result.update(scope='unsupported',expected_name=None,reason='環外の立体中心は対象外。');return result
  for atom in mol.GetAtoms():
   if atom.HasProp('_CIPCode'):atom.ClearProp('_CIPCode')
  try:rdCIPLabeler.AssignCIPLabels(mol)
  except (RuntimeError,ValueError):
   result.update(scope='unsupported',expected_name=None,numbering={},substituents=[],reason='CIP計算を完了できないため全体名を保留。',stereo_status='unavailable');return result
  stereo=[]
  for n in (2,3):
   idx=mapping[n]
   if idx in centers:
    atom=mol.GetAtomWithIdx(idx);cip=atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
    if cip not in ('R','S'):cip=None
    if atom.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED and not cip:
     result.update(scope='unsupported',expected_name=None,numbering={},substituents=[],reason='R/Sを確定できない指定立体は対象外。',stereo_status='unavailable');return result
    stereo.append(dict(position=str(n),atom_index=idx,cip=cip))
  assigned=sum(c['cip'] is not None for c in stereo)
  result['ring_stereo']=stereo
  result['stereo_status']='complete' if assigned==len(stereo) else ('unspecified' if not assigned else 'partial')
  if 0<assigned<len(stereo):
   result.update(scope='partial_stereo',expected_name=None,reason='一部の環立体が未指定のため全体期待名を保留。未指定部分は補いません。')
  else:
   descriptor='('+','.join(c['position']+c['cip'] for c in stereo)+')-' if assigned else ''
   result['expected_name']=descriptor+found['expected'].replace('-4H-chromen-4-one','chroman-4-one')
   result['reason']+=' / '+('全環立体を絶対R/S指定。' if assigned else '環立体は未指定。ラセミ体とは断定しません。')
 else:result['saturation']='unsaturated'
 result['equivalent_numberings']=[{LOCANTS[n]:idx for n,idx in m.items() if n in LOCANTS} for _,m in options]
 index={a.GetIdx():a.GetIntProp('_flv_original_idx') for a in mol.GetAtoms()}
 result['numbering']={n:index[i] for n,i in result['numbering'].items()}
 result['equivalent_numberings']=[{n:index[i] for n,i in m.items()} for m in result['equivalent_numberings']]
 result['carbonyl_oxygen']=index[result['carbonyl_oxygen']]
 for center in result.get('ring_stereo',[]):center['atom_index']=index[center['atom_index']]
 for group in result['substituents']:
  group['ring_atom']=index[group['ring_atom']];group['atoms']=[index[i] for i in group['atoms']]
 return result


def normalize_name(name):
 return re.sub(r",\s+",",",name.strip().lower().replace('′',"'").replace('’',"'")).replace('2h-1-benzopyran','chromen').replace('2h-chromen','chromen').replace('4h-1-benzopyran','chromen').replace('1-benzopyran','chromen').replace('4h-chromen','chromen')


class UnsupportedName(ValueError):
 """A name outside the finite grammar, not evidence of incorrect chemistry."""


def parse_groups(text,primed=False,allowed=None):
 groups={'hydroxy':[],'methoxy':[]};errors=[];labels=[];seen=set()
 # Full consumption is required. Unknown words and punctuation never disappear.
 pattern=re.compile(r"(\d+'?(?:,\d+'?)*)-(deca|nona|octa|hepta|hexa|penta|tetra|tri|di)?(hydroxy|methoxy)(?:-|$)")
 while text:
  match=pattern.match(text)
  if not match:raise UnsupportedName()
  raw=match[1].split(',');label=match[3];labels.append(label)
  if len(raw)!=next(n for n,prefix in MULT.items() if prefix==(match[2] or '')):errors.append('倍数接頭辞と位置番号の個数が一致しません。')
  if groups[label]:errors.append('同じ置換基は位置番号をまとめて記載します。')
  locs=[]
  for value in raw:
   if value.endswith("'") and not primed:raise UnsupportedName()
   number=int(value.rstrip("'"));loc=(number,value.endswith("'"))
   if allowed is not None and loc not in allowed:errors.append('この親骨格・置換基では使用できない位置番号です：'+value)
   if loc in seen:errors.append('置換位置が重複しています：'+value)
   seen.add(loc);locs.append(loc)
  if locs!=sorted(locs):errors.append('位置番号が昇順に並んでいません。')
  groups[label].extend(locs);text=text[match.end():]
 return groups,errors,labels


def parse_flavonoid_name(name):
 """Finite systematic and flavone-parent grammars. No trivial-name inference."""
 if re.search(r'(?<!-)[24]h-',name.strip().lower()):return None
 raw_h=re.search(r'(?:^|-)([24])h-(?:1-benzopyran|chromen)-',name.strip().lower())
 indicated_h=int(raw_h[1]) if raw_h else None
 stereo={}
 if name.strip().startswith('('):
  match=re.match(r'^\(([23][RS](?:,[23][RS])?)\)-',name.strip())
  if not match:return None
  for token in match[1].split(','):
   if token[0] in stereo:return None
   stereo[token[0]]=token[1]
  name=name.strip()[match.end():]
 name=normalize_name(name)
 if '--' in name or name.startswith('-'):return None
 try:
  semi=next((parent for parent in ('flavan-4-one','isoflavone','neoflavone','flavone') if name.endswith(parent)),None)
  if semi:
   attachment,carbonyl,h=PARENTS[{'flavan-4-one':'flavanone'}.get(semi,semi)]
   allowed={(n,False) for n in ({2,3,4,5,6,7,8}-{attachment,carbonyl})}|{(n,True) for n in range(2,7)}
   groups,errors,labels=parse_groups(name[:-len(semi)],True,allowed)
   a={k:sorted(n for n,p in ns if not p) for k,ns in groups.items()}
   b={k:sorted(n for n,p in ns if p) for k,ns in groups.items()}
   return dict(style=semi,attachment=attachment,carbonyl=carbonyl,indicated_h=h,a=a,b=b,errors=errors,ordered=labels==sorted(labels),saturation='saturated' if semi=='flavan-4-one' else 'unsaturated',stereo=stereo)
  suffix=re.search(r'(chroman|(?:2,3-dihydro-)?chromen)-([234])-one$',name)
  if not suffix:return None
  head=name[:suffix.start()].removesuffix('-')
  phenyl=re.search(r'(?:^|-)([234])-(?:\(([^()]*)phenyl\)|phenyl)(?=-|$)',head)
  if not phenyl:return None
  rest='-'.join(part.strip('-') for part in (head[:phenyl.start()],head[phenyl.end():]) if part.strip('-'))
  a,ae,al=parse_groups(rest,False,{(n,False) for n in range(2,9)})
  b,be,bl=parse_groups(phenyl[2] or '',False,{(n,False) for n in range(2,7)})
  # Track order independently of structural correspondence. A reordered name
  # can still denote exactly the same molecule; it is not a structure failure.
  ordered=al==sorted(al) and bl==sorted(bl)
  if head[phenyl.end():].strip('-'):ordered=False
  return dict(style='systematic',attachment=int(phenyl[1]),carbonyl=int(suffix[2]),indicated_h=indicated_h,a={k:sorted(n for n,_ in ns) for k,ns in a.items()},b={k:sorted(n for n,_ in ns) for k,ns in b.items()},errors=ae+be,ordered=ordered,saturation='saturated' if suffix[1]!='chromen' else 'unsaturated',stereo=stereo)
 except UnsupportedName:
  return None


def flavonoid_checks(analysis,name):
 if analysis['scope']=='bounded_glycoside':return flavonoid_glycoside_checks(analysis,name)
 if analysis['scope']!='bounded':return []
 if analysis['family'] in ('chalcone','aurone'):return enone_checks(analysis,name)
 parsed=parse_flavonoid_name(name)
 def row(key,title,status,detail):return dict(id=key,title=title,status=status,detail=detail,reference='Flv-3 / Flv-3.1.1.3 / Flv-3.2.1–3.2.3',url=SOURCE)
 expected='限定期待名：'+analysis['expected_name']
 if parsed is None:
  reason='この候補名は有限構文の範囲外です。構造の往復照合は別欄を参照。'
  return [row('flavonoid_parent','フラボノイドの親骨格・番号','unknown',reason),row('flavonoid_substituents','OH・メトキシ基の位置','unknown',expected),row('flavonoid_suffix','環カルボニルとOHの表記','unknown',reason),row('flavonoid_spelling','位置番号・倍数・接頭辞の表記','unknown',reason)]
 if parsed['errors']:
  return [row('flavonoid_spelling','位置番号・倍数・接頭辞の表記','fail',' / '.join(dict.fromkeys(parsed['errors']))),row('flavonoid_parent','フラボノイドの親骨格・番号','unknown','位置番号等の表記に不整合があるため判定を保留。'),row('flavonoid_substituents','OH・メトキシ基の位置','unknown',expected),row('flavonoid_suffix','環カルボニルとOHの表記','unknown','不整合な表記のため判定を保留。')]
 parent=parsed['attachment']==analysis.get('parent_attachment',2) and parsed['carbonyl']==analysis.get('parent_carbonyl',4) and parsed['saturation']==analysis.get('saturation','unsaturated') and parsed.get('indicated_h') in (None,analysis.get('parent_indicated_h',4))
 same=parsed['a']==analysis['a_groups'] and parsed['b']==analysis['b_groups']
 rows=[row('flavonoid_parent','フラボノイドの親骨格・番号','pass' if parent else 'fail','入力の親骨格：'+analysis.get('parent_name','2-phenylchromen-4-one')+'。候補の飽和度：'+parsed['saturation']+'、phenyl結合位置：'+str(parsed['attachment'])+'、環カルボニル位置：'+str(parsed['carbonyl'])+'、指示H：'+str(parsed.get('indicated_h') or '省略')+'。'),row('flavonoid_substituents','OH・メトキシ基の位置',('pass' if same else 'fail') if parent else 'unknown',expected+' / 候補の非プライム位置：'+str(parsed['a'])+' / B環位置：'+str(parsed['b'])),row('flavonoid_suffix','環カルボニルとOHの表記','pass' if parent else 'fail','この骨格では'+str(analysis.get('parent_carbonyl',4))+'-oneの環カルボニルを保持し、OHをhydroxy接頭辞で示します。ネオフラボンの2-oneはラクトン骨格です。PIN全体の証明ではありません。'),row('flavonoid_spelling','位置番号・倍数・接頭辞の表記','pass' if parsed['ordered'] else 'different','位置番号の昇順・倍数・重複を確認。'+('接頭辞の有限な列挙順とも一致。' if parsed['ordered'] else '構造の対応と接頭辞の記載順の差を分けて表示します。'))]

 if analysis.get('saturation')=='saturated' or parsed['stereo']:
  input_stereo={c['position']:c['cip'] for c in analysis.get('ring_stereo',[]) if c['cip']}
  same_stereo=input_stereo==parsed['stereo']
  rows.append(row('flavonoid_stereo','環立体の絶対R/S','pass' if same_stereo and parent else ('fail' if parent else 'unknown'),'入力：'+str(input_stereo or '未指定')+' / 候補：'+str(parsed['stereo'] or '未指定')+'。未指定はラセミ体と断定せず、R/Sを補完しません。CIP割当てはRDKitに依存。'))
 return rows


GLYCOSIDE_NOTE='イソフラボン・フラバノン／フラバノノールのA環5–8位のO配糖体も対応。フラバノンの環立体R/Sは未切断の配糖体で再計算し、全指定または全未指定のみ。 フラボン／フラボノール（O: C環3位・A環5–8位、C: A環6・8位）に未修飾グルコピラノシル基が1個結合する限定範囲です。D/L・α/βは4つの局所立体参照と照合し、未指定部分を補いません。他糖種・複数糖・B環配糖体・他骨格・PIN全体は対象外です。'


def analyze_flavonoids(mol):
 result=analyze_aglycones(mol)
 if result['expected_name'] or not result['detected']:return result
 # Keep atom identity through H normalization, subgraph extraction and capping.
 original=Chem.Mol(mol)
 for atom in original.GetAtoms():atom.SetIntProp('_glyco_original_idx',atom.GetIdx())
 original=Chem.RemoveHs(original)
 if len(Chem.GetMolFrags(original))!=1 or original.GetStereoGroups() or any(a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in original.GetAtoms()):return result
 from carbohydrate_rules import REGISTRY,label
 options=[]
 for key,ref in REGISTRY.items():
  if ref['family']!='glucose':continue
  # Replace the anomeric OH query in place, preserving bond order and the
  # C1 tetrahedron. C-glycosyl identity never comes from a CIP-letter shortcut.
  query=Chem.RWMol(Chem.MolFromSmarts(key))
  roles={a.GetAtomMapNum():a.GetIdx() for a in query.GetAtoms()}
  query.ReplaceAtom(roles[11],Chem.AtomFromSmarts('[#6,#8]'))
  query=query.GetMol()
  for match in original.GetSubstructMatches(query,useChirality=True,uniquify=False,maxMatches=256):
   mapping={n:match[i] for n,i in roles.items()};anchor=mapping[11]
   atom=original.GetAtomWithIdx(anchor);kind='O' if atom.GetAtomicNum()==8 else 'C'
   sugar=set(match)-{anchor}
   if kind=='O':
    if atom.GetDegree()!=2 or atom.GetTotalNumHs()!=0:continue
    targets=[a.GetIdx() for a in atom.GetNeighbors() if a.GetIdx() not in sugar]
    if len(targets)!=1:continue
    target=targets[0]
   else:target=anchor
   if not original.GetAtomWithIdx(target).GetIsAromatic():continue
   # Each sugar atom must have no unaccounted substituents. This rejects
   # acylated/extended sugars and bridges back into the aglycon.
   if any(a.GetIdx() not in set(match) for i in sugar for a in original.GetAtomWithIdx(i).GetNeighbors()):continue
   copy=Chem.RWMol(original)
   for a in copy.GetAtoms():a.SetIntProp('_glyco_idx',a.GetIdx())
   if kind=='O':copy.GetAtomWithIdx(anchor).SetNumExplicitHs(1)
   for i in sorted(sugar,reverse=True):copy.RemoveAtom(i)
   try:Chem.SanitizeMol(copy)
   except (RuntimeError,ValueError):continue
   ag=analyze_aglycones(copy)
   if ag['scope']!='bounded' or ag['family'] not in ('flavone','flavonol','isoflavone','flavanone','flavanonol'):continue
   if kind=='C' and ag['family'] not in ('flavone','flavonol'):continue
   index={a.GetIdx():a.GetIntProp('_glyco_idx') for a in copy.GetAtoms()}
   num={n:index[i] for n,i in ag['numbering'].items()}
   position=next((n for n,i in num.items() if i==target),None)
   allowed=('3','5','6','7','8') if ag['family'] in ('flavone','flavonol') else ('5','6','7','8')
   if position not in (allowed if kind=='O' else ('6','8')):continue
   a={k:list(v) for k,v in ag['a_groups'].items()};b=ag['b_groups']
   if kind=='O':a['hydroxy'].remove(int(position))
   base=expected_name(a,b,ag['parent_attachment'],ag['parent_carbonyl'],ag['parent_indicated_h'])
   stereo=[];descriptor=''
   if ag.get('saturation')=='saturated':
    # CIP belongs to the intact glycoside, not the capped private aglycon.
    full=Chem.Mol(original)
    for atom in full.GetAtoms():
     if atom.HasProp('_CIPCode'):atom.ClearProp('_CIPCode')
    try:rdCIPLabeler.AssignCIPLabels(full)
    except (RuntimeError,ValueError):continue
    good=True
    for center in ag.get('ring_stereo',[]):
     idx=index[center['atom_index']];atom=full.GetAtomWithIdx(idx)
     cip=atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
     if atom.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED and cip not in ('R','S'):good=False;break
     stereo.append(dict(position=center['position'],atom_index=idx,cip=cip if cip in ('R','S') else None))
    if not good:continue
    assigned=[c for c in stereo if c['cip']]
    if assigned and len(assigned)!=len(stereo):continue
    descriptor='('+','.join(c['position']+c['cip'] for c in assigned)+')-' if assigned else ''
    base=base.replace('-4H-chromen-4-one','chroman-4-one')
   sugar_label=label(ref,'osyl')
   whole=descriptor+position+'-('+sugar_label+('oxy' if kind=='O' else '')+')-'+base
   out=dict(ag)
   out.update(a_groups=a,base_expected_name=descriptor+base,expected_name=whole,label=ag['label']+' '+kind+'配糖体',scope='bounded_glycoside',note=GLYCOSIDE_NOTE,reason='糖残基とアグリコンを全原子について照合。全体名の立体構造照合は別判定。',numbering=num,equivalent_numberings=[{n:index[i] for n,i in m.items()} for m in ag['equivalent_numberings']])
   out['carbonyl_oxygen']=index[ag['carbonyl_oxygen']]
   if ag.get('saturation')=='saturated':out['ring_stereo']=stereo
   out['substituents']=[dict(g,ring_atom=index[g['ring_atom']],atoms=[index[i] for i in g['atoms']]) for g in ag['substituents'] if not(kind=='O' and g['position']==position)]
   out['glycoside']=dict(kind=kind,position=position,sugar_label=sugar_label,**ref,carbon_atoms={str(n):mapping[n] for n in range(1,7)},ring_oxygen=mapping[15],anomeric_atom=mapping[1],aglycon_atom=target,bridge_oxygen=anchor if kind=='O' else None)
   # Export original input atom indices, including when ordinary H was explicit.
   restore=lambda i:original.GetAtomWithIdx(i).GetIntProp('_glyco_original_idx')
   out['numbering']={n:restore(i) for n,i in out['numbering'].items()}
   out['equivalent_numberings']=[{n:restore(i) for n,i in m.items()} for m in out['equivalent_numberings']]
   out['carbonyl_oxygen']=restore(out['carbonyl_oxygen'])
   for center in out.get('ring_stereo',[]):center['atom_index']=restore(center['atom_index'])
   for g in out['substituents']:g.update(ring_atom=restore(g['ring_atom']),atoms=[restore(i) for i in g['atoms']])
   g=out['glycoside'];g['carbon_atoms']={n:restore(i) for n,i in g['carbon_atoms'].items()}
   for n in ('ring_oxygen','anomeric_atom','aglycon_atom','bridge_oxygen'):
    if g[n] is not None:g[n]=restore(g[n])
   options.append(out)
 if len(options)==1:return options[0]
 if options:result['reason']='複数の配糖体の番号・立体候補があるため全体名を保留。'
 return result


def flavonoid_glycoside_checks(analysis,name):
 def row(status,detail):return dict(id='flavonoid_glycoside',title='糖残基・O/C結合・結合位置',status=status,detail=detail,reference='Flv-3.2.1.2 / Flv-3.2.1.3',url=SOURCE)
 stereo=re.match(r'^\(([23][RS](?:,[23][RS])?)\)-',name.strip())
 descriptor=stereo[0] if stereo else ''
 text=name.strip()[len(descriptor):]
 match=re.fullmatch(r'(\d+)-\((alpha|beta)-(D|L)-glucopyranosyl(oxy)?\)-(.+)',text)
 if not match:return [row('unknown','この候補名は限定したグルコピラノシル接頭辞の構文外です。全体の構造照合は別判定。')]
 g=analysis['glycoside']
 same=(match[1],match[2],match[3],bool(match[4]))==(g['position'],g['anomer'],g['series'],g['kind']=='O')
 base=dict(analysis,scope='bounded',expected_name=analysis.get('base_expected_name',analysis['expected_name'].split(')-',1)[1]))
 return [row('pass' if same else 'fail','入力：'+g['sugar_label']+' / '+g['kind']+'結合 / アグリコン'+g['position']+'位。糖の局所立体参照と原子対応を照合。全体PINの証明ではありません。')]+flavonoid_checks(base,descriptor+match[5])
