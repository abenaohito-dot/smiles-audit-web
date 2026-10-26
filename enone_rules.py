"""Bounded OH/OMe chalcone and aurone aglycones; E/Z is never inferred.
IUPAC Flv-3.3/3.4. CIP labelling is delegated to RDKit, not independently proved.
"""
import re
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
SOURCE='https://iupac.qmul.ac.uk/flavonoid/'
RING_B='[c:101]2[c:102][c:103][c:104][c:105][c:106]2'
RING_A='[c:201]1[c:202][c:203][c:204][c:205][c:206]1'
CHALCONE_CORE='[O:40]=[C;!R:1](-'+RING_A+')-[CH;!R:2]=[CH;!R:3]-'+RING_B
AURONE_CORE='[O:40]=[C:3]1[C:2](=[CH:20]-'+RING_B+')[O:1][c:9]2[c:7][c:6][c:5][c:4][c:8]21'
QUERIES={'chalcone':Chem.MolFromSmarts(CHALCONE_CORE),'aurone':Chem.MolFromSmarts(AURONE_CORE)}


def rank(groups):
 return (tuple(sorted(groups['hydroxy']+groups['methoxy'])),tuple(groups['hydroxy']),tuple(groups['methoxy']))


def expected(family,a,b,ez):
 from flavonoid_rules import prefix
 pre='(2'+ez+')-' if ez else ''
 if family=='chalcone':
  if not any(a.values()) and not any(b.values()):return pre+'1,3-diphenylprop-2-en-1-one'
  def phenyl(n,g):
   p=prefix(g);return str(n)+'-('+p+'phenyl)' if p else str(n)+'-phenyl'
  return pre+phenyl(1,a)+'-'+phenyl(3,b)+'prop-2-en-1-one'
 p=prefix(b);benzyl='2-('+p+'benzylidene)' if p else '2-benzylidene'
 ag=prefix(a)
 return pre+benzyl+('-'+ag if ag else '')+'-1-benzofuran-3(2H)-one'


def analyze_enone(mol,family,result):
 # Caller provides a private, H-normalized copy with original-index properties.
 result.update(note='OH・メトキシのみの中性非配糖体に限定。E/ZはRDKitのCIP計算へ委譲し、規則自体の独立証明とは区別します。未指定をE/Zに補完しません。慣用名、他官能基、配糖体、全体PINは対象外。')
 if len(Chem.GetMolFrags(mol))!=1 or mol.GetStereoGroups() or any(a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() or a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
  result['reason']='塩・同位体・電荷・四面体/拡張立体は対象外。';return result
 options=[];q=QUERIES[family]
 for match in mol.GetSubstructMatches(q,uniquify=False,maxMatches=256):
  mapping={a.GetAtomMapNum():match[a.GetIdx()] for a in q.GetAtoms()};core=set(mapping.values());used=set(core)
  allowed=set(range(202,207))|set(range(102,107)) if family=='chalcone' else {4,5,6,7}|set(range(102,107))
  a={'hydroxy':[],'methoxy':[]};b={'hydroxy':[],'methoxy':[]};details=[];valid=True
  for n,idx in mapping.items():
   extra=[x for x in mol.GetAtomWithIdx(idx).GetNeighbors() if x.GetIdx() not in core]
   if not extra:continue
   if n not in allowed or len(extra)!=1:valid=False;break
   o=extra[0];oi=o.GetIdx()
   if o.GetAtomicNum()!=8 or o.GetIsAromatic() or mol.GetBondBetweenAtoms(idx,oi).GetBondType()!=Chem.BondType.SINGLE:valid=False;break
   branch=[x for x in o.GetNeighbors() if x.GetIdx()!=idx]
   if not branch and o.GetTotalNumHs()==1:kind='hydroxy';atoms=[oi]
   elif len(branch)==1 and o.GetTotalNumHs()==0:
    c=branch[0]
    if c.GetAtomicNum()!=6 or c.GetDegree()!=1 or c.GetTotalNumHs()!=3 or mol.GetBondBetweenAtoms(oi,c.GetIdx()).GetBondType()!=Chem.BondType.SINGLE:valid=False;break
    kind='methoxy';atoms=[oi,c.GetIdx()]
   else:valid=False;break
   if used.intersection(atoms):valid=False;break
   used.update(atoms)
   if family=='chalcone':
    side='a' if n>200 else 'b';pos=n-(200 if side=='a' else 100);label=str(pos)+("'" if side=='a' else '')
   else:side='b' if n>100 else 'a';pos=n-100 if side=='b' else n;label=str(pos)+("'" if side=='b' else '')
   (a if side=='a' else b)[kind].append(pos);details.append(dict(position=label,kind=kind,ring_atom=idx,atoms=atoms))
  if valid and used==set(range(mol.GetNumAtoms())):
   for g in (a,b):
    for k in g:g[k].sort()
   options.append((rank(a)+rank(b),mapping,a,b,details))
 if not options:result['reason']='OH・メトキシ以外の置換、追加環、配糖体などは全体命名範囲外。';return result
 best=min(o[0] for o in options);options=[o for o in options if o[0]==best]
 _,mapping,a,b,details=min(options,key=lambda o:tuple(o[1][n] for n in sorted(o[1])))
 endpoints=(mapping[2],mapping[3] if family=='chalcone' else mapping[20]);bond=mol.GetBondBetweenAtoms(*endpoints)
 if bond.GetStereo()==Chem.BondStereo.STEREOANY or any(x.GetIdx()!=bond.GetIdx() and x.GetStereo()!=Chem.BondStereo.STEREONONE for x in mol.GetBonds()):
  result['reason']='未確定E/Zや他の結合立体は対象外。';return result
 for x in mol.GetBonds():
  if x.HasProp('_CIPCode'):x.ClearProp('_CIPCode')
 try:rdCIPLabeler.AssignCIPLabels(mol)
 except (RuntimeError,ValueError):result['reason']='E/ZのCIP計算を完了できません。';return result
 ez=bond.GetProp('_CIPCode') if bond.HasProp('_CIPCode') else None
 if ez not in ('E','Z'):ez=None
 if bond.GetStereo()!=Chem.BondStereo.STEREONONE and not ez:result['reason']='指定された結合立体をE/Zとして確定できません。';return result
 index={x.GetIdx():x.GetIntProp('_flv_original_idx') for x in mol.GetAtoms()}
 def numbers(m):
  if family=='chalcone':
   return {('鎖C'+str(n) if n<4 else (str(n-200)+"'" if n>200 else str(n-100))):index[i] for n,i in m.items() if n!=40}
  labels={1:'1',2:'2',3:'3',4:'4',5:'5',6:'6',7:'7',8:'3a',9:'7a',20:'環外CH',**{100+i:str(i)+"'" for i in range(1,7)}}
  return {labels[n]:index[i] for n,i in m.items() if n in labels}
 for d in details:d['ring_atom']=index[d['ring_atom']];d['atoms']=[index[i] for i in d['atoms']]
 result.update(scope='bounded',label='カルコン骨格' if family=='chalcone' else 'オーロン骨格',expected_name=expected(family,a,b,ez),numbering=numbers(mapping),equivalent_numberings=[numbers(o[1]) for o in options],carbonyl_oxygen=index[mapping[40]],a_groups=a,b_groups=b,substituents=details,alkene_stereo=ez,alkene_atoms=[index[i] for i in endpoints],stereo_status='specified' if ez else 'unspecified',parent_name='1,3-diphenylprop-2-en-1-one' if family=='chalcone' else '2-benzylidene-1-benzofuran-3(2H)-one',reason='中性・非配糖体、OH・メトキシのみ。 / '+('E/Z指定あり。' if ez else 'E/Z未指定。配置は補完しません。'))
 if family=='chalcone':result['note']+=' カルボニル側芳香環はプライム付き、鎖C3側芳香環はプライムなし。鎖C1〜C3は芳香環番号と別です。'
 return result


def parse_enone_name(name):
 from flavonoid_rules import parse_groups,UnsupportedName,normalize_name
 raw=name.strip();ez=None;loc=None
 if raw.startswith('('):
  m=re.match(r'^\((\d+)([EZ])\)-',raw)
  if not m:return None
  loc=int(m[1]);ez=m[2];raw=raw[m.end():]
 text=normalize_name(raw)
 if '--' in text:return None
 try:
  if text.endswith('prop-2-en-1-one'):
   head=text[:-len('prop-2-en-1-one')];parts={};errors=[];order=[]
   if head=='1,3-diphenyl':parts={1:'',3:''}
   else:
    pattern=re.compile(r'([13])-(?:\(([^()]*)phenyl\)|phenyl)(?:-|$)')
    while head:
     m=pattern.match(head)
     if not m:return None
     n=int(m[1]);order.append(n)
     if n in parts:errors.append('同じ鎖位置のphenylが重複しています。')
     parts[n]=m[2] or '';head=head[m.end():]
   if set(parts)!={1,3}:return None
   a,ae,al=parse_groups(parts[1],False,{(n,False) for n in range(2,7)})
   b,be,bl=parse_groups(parts[3],False,{(n,False) for n in range(2,7)})
   family='chalcone';ordered=al==sorted(al) and bl==sorted(bl) and order in ([],[1,3])
  else:
   suffix=re.search(r'(?:1-)?benzofuran-3\(2h\)-one$',text)
   if not suffix:return None
   head=text[:suffix.start()].removesuffix('-');m=re.match(r'2-(?:\(([^()]*)benzylidene\)|benzylidene)(?:-|$)',head)
   if not m:return None
   a,ae,al=parse_groups(head[m.end():],False,{(n,False) for n in range(4,8)})
   b,be,bl=parse_groups(m[1] or '',False,{(n,False) for n in range(2,7)})
   family='aurone';errors=[];ordered=al==sorted(al) and bl==sorted(bl)
  return dict(family=family,a={k:sorted(n for n,_ in v) for k,v in a.items()},b={k:sorted(n for n,_ in v) for k,v in b.items()},ez=ez,ez_locant=loc,errors=errors+ae+be,ordered=ordered)
 except UnsupportedName:return None


def enone_checks(analysis,name):
 parsed=parse_enone_name(name)
 def row(k,t,s,d):return dict(id=k,title=t,status=s,detail=d,reference='Flv-3.3 / Flv-3.4',url=SOURCE)
 titles=[('flavonoid_parent','カルコン・オーロンの親骨格'),('flavonoid_substituents','OH・メトキシと両芳香環の位置'),('flavonoid_alkene_stereo','二重結合のE/Z'),('flavonoid_spelling','位置番号・倍数の表記')]
 if not parsed:return [row(k,t,'unknown','候補名は今回の有限系統名構文の範囲外です。構造往復照合は別に表示。') for k,t in titles]
 if parsed['errors']:return [row(k,t,'fail' if k=='flavonoid_spelling' else 'unknown',' / '.join(dict.fromkeys(parsed['errors']))) for k,t in titles]
 parent=parsed['family']==analysis['family'];same=parsed['a']==analysis['a_groups'] and parsed['b']==analysis['b_groups']
 stereo=parsed['ez']==analysis['alkene_stereo'] and parsed['ez_locant'] in (None,2)
 return [row('flavonoid_parent',titles[0][1],'pass' if parent else 'fail','入力：'+analysis['parent_name']+' / 候補骨格：'+parsed['family']),row('flavonoid_substituents',titles[1][1],('pass' if same else 'fail') if parent else 'unknown','両側を区別して構造から比較。限定期待名：'+analysis['expected_name']),row('flavonoid_alkene_stereo',titles[2][1],('pass' if stereo else 'fail') if parent else 'unknown','入力：'+str(analysis['alkene_stereo'] or '未指定')+' / 候補：'+str(parsed['ez'] or '未指定')+' / E/Z位置：'+str(parsed['ez_locant'] or 'なし')+'。割当てはRDKit依存。未指定を補完しません。'),row('flavonoid_spelling',titles[3][1],'pass' if parsed['ordered'] else 'different','位置番号の昇順・重複・倍数を照合。記載順の差は構造対応と別に表示。')]
