"""Finite C5/C6 monocycle evaluator; not a terpenoid identity classifier.
Full atom accounting, explicit numbering choices, no naming-engine oracle.
"""
import re
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from extended_rules import group_rank

SOURCE='https://iupac.qmul.ac.uk/BlueBook/P2.html'
NOTE='5/6員炭素単環、環内二重結合0/1本、環上OHまたはケトン0/1個、methyl/ethyl/propyl/propan-2-yl/prop-1-en-2-yl側鎖3個まで（分岐C3側鎖は各1個まで）の限定検証です。天然物・テルペノイドの分類やPIN全体の証明ではありません。部分立体・多環・芳香環・環外OH/カルボニル・追加官能基は保留。R/SはRDKit CIPへ委譲。'
SIDE={'methyl':('C',0),'ethyl':('CC',0),'propyl':('CCC',0),'propan-2-yl':('CCC',1),'prop-1-en-2-yl':('C=C(C)',1)}


def side_label(mol,atoms,first):
 if any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 for i in atoms):return None
 for label,(smi,root) in SIDE.items():
  template=Chem.MolFromSmiles(smi)
  # Require the entire substituent, including its attachment vertex.
  for match in mol.GetSubstructMatches(template,uniquify=False):
   if set(match)==atoms and match[root]==first:return label
 return None


def prefixes(groups):
 grouped=defaultdict(list)
 for loc,label in groups:grouped[label].append(loc)
 terms=[]
 for label,locs in sorted(grouped.items(),key=lambda pair:re.sub(r'^\d+-','',pair[0])):
  if '-' in label:
   if len(locs)!=1:return None
   text='('+label+')'
  else:text={1:'',2:'di',3:'tri'}[len(locs)]+label
  terms.append(','.join(map(str,sorted(locs)))+'-'+text)
 return '-'.join(terms)


def base_name(n,kind,functional,doubles,groups):
 root={5:'cyclopent',6:'cyclohex'}[n]
 if kind:
  parent=root+('-'+str(doubles[0])+'-en' if doubles else 'an')+'-'+str(functional)+'-'+('ol' if kind=='alcohol' else 'one')
 else:parent=root+('-'+str(doubles[0])+'-ene' if doubles else 'ane')
 prefix=prefixes(groups)
 return (prefix or '')+parent if prefix is not None else None


def analyze_ring_monocycles(mol):
 result=dict(detected=False,scope='unsupported',expected_name=None,numbering={},substituents=[],ring_stereo=[],reason='',note=NOTE,source=SOURCE)
 mol=Chem.Mol(mol)
 for a in mol.GetAtoms():a.SetIntProp('_mc_original',a.GetIdx())
 mol=Chem.RemoveHs(mol)
 rings=mol.GetRingInfo().AtomRings()
 if len(rings)!=1 or len(rings[0]) not in (5,6) or any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 or mol.GetAtomWithIdx(i).GetIsAromatic() for i in rings[0]):return result
 result['detected']=True;result['reason']='単環候補。追加官能基や立体等が限定範囲外の場合は全体名を保留。'
 if len(Chem.GetMolFrags(mol))!=1 or mol.GetNumAtoms()>20 or mol.GetStereoGroups() or any(a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() or a.GetAtomicNum() not in (6,8) for a in mol.GetAtoms()) or any(b.GetBondType() not in (Chem.BondType.SINGLE,Chem.BondType.DOUBLE) or b.GetStereo()!=Chem.BondStereo.STEREONONE for b in mol.GetBonds()):return result
 ring=set(rings[0]);oxygen=[a for a in mol.GetAtoms() if a.GetAtomicNum()==8]
 if len(oxygen)>1:return result
 kind=None;functional=None
 if oxygen:
  o=oxygen[0]
  if o.GetDegree()!=1 or o.GetNeighbors()[0].GetIdx() not in ring:return result
  functional=o.GetNeighbors()[0].GetIdx();bond=mol.GetBondBetweenAtoms(o.GetIdx(),functional)
  if bond.GetBondType()==Chem.BondType.SINGLE and o.GetTotalNumHs()==1:kind='alcohol'
  elif bond.GetBondType()==Chem.BondType.DOUBLE and o.GetTotalNumHs()==0:kind='ketone'
  else:return result
 side=[];used=ring|{o.GetIdx() for o in oxygen}
 for idx in ring:
  for a in mol.GetAtomWithIdx(idx).GetNeighbors():
   first=a.GetIdx()
   if first in used:continue
   if mol.GetBondBetweenAtoms(idx,first).GetBondType()!=Chem.BondType.SINGLE:return result
   atoms=set();stack=[first]
   while stack:
    current=stack.pop()
    if current in ring or current in atoms:continue
    atoms.add(current);stack.extend(nb.GetIdx() for nb in mol.GetAtomWithIdx(current).GetNeighbors())
   if atoms.intersection(used) or not 1<=len(atoms)<=3:return result
   label=side_label(mol,atoms,first)
   if not label:return result
   side.append((idx,label,sorted(atoms)));used.update(atoms)
 if len(side)>3 or used!=set(range(mol.GetNumAtoms())):return result
 cycle=list(rings[0]);n=len(cycle);options=[]
 for direction in (cycle,list(reversed(cycle))):
  for offset in range(n):
   path=direction[offset:]+direction[:offset]
   pos={idx:i+1 for i,idx in enumerate(path)}
   doubles=sorted(min(pos[b.GetBeginAtomIdx()],pos[b.GetEndAtomIdx()]) if abs(pos[b.GetBeginAtomIdx()]-pos[b.GetEndAtomIdx()])==1 else n for b in mol.GetBonds() if b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring and b.GetBondType()==Chem.BondType.DOUBLE)
   if len(doubles)>1:return result
   groups=sorted((pos[idx],label) for idx,label,_ in side)
   loc=pos[functional] if functional is not None else None
   text=base_name(n,kind,loc,doubles,groups)
   if text is None:return result
   rank=((loc,) if kind else (),tuple(doubles),group_rank(groups))
   options.append(dict(path=path[:],functional=loc,doubles=doubles,groups=groups,base=text,rank=rank))
 # Every rotation, in both directions, must have been enumerated exactly once.
 best_rank=min(c['rank'] for c in options);best=[c for c in options if c['rank']==best_rank]
 centers=dict(Chem.FindMolChiralCenters(mol,includeUnassigned=True,useLegacyImplementation=False))
 if any(i not in ring for i in centers):return result
 for a in mol.GetAtoms():
  if a.HasProp('_CIPCode'):a.ClearProp('_CIPCode')
 try:rdCIPLabeler.AssignCIPLabels(mol)
 except (RuntimeError,ValueError):return result
 alternatives=[]
 for c in best:
  stereo=[]
  for loc,idx in enumerate(c['path'],1):
   if idx not in centers:continue
   a=mol.GetAtomWithIdx(idx);cip=a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None
   if a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED and cip not in ('R','S'):return result
   stereo.append(dict(position=str(loc),atom_index=idx,cip=cip if cip in ('R','S') else None))
  assigned=[x for x in stereo if x['cip']]
  descriptor='('+','.join(x['position']+x['cip'] for x in assigned)+')-' if assigned else ''
  alternatives.append((c,stereo,descriptor))
 if len({d for _,_,d in alternatives})!=1:
  result['reason']='同順位の番号付けで立体記述子が変わるため全体名を保留。';return result
 c,stereo,descriptor=alternatives[0];index={a.GetIdx():a.GetIntProp('_mc_original') for a in mol.GetAtoms()}
 partial=0<sum(x['cip'] is not None for x in stereo)<len(stereo)
 result.update(scope='partial_stereo' if partial else 'bounded',kind=kind,ring_size=n,expected_name=None if partial else descriptor+c['base'],base_expected_name=c['base'],numbering={str(i+1):index[idx] for i,idx in enumerate(c['path'])},double_locants=c['doubles'],functional_locant=c['functional'],groups=c['groups'],ring_stereo=[dict(x,atom_index=index[x['atom_index']]) for x in stereo],stereo_status='partial' if partial else ('complete' if descriptor else 'unspecified'),reason='環の主特性基 → 環内二重結合 → 側鎖の位置番号を順に比較。'+('一部の立体が未指定のため全体名を保留。' if partial else '未指定立体を補完しません。'),orientations=[dict(path=[index[i] for i in o['path']],functional=o['functional'],doubles=o['doubles'],groups=o['groups']) for o in options])
 result['substituents']=[dict(position=str(c['path'].index(idx)+1),kind=label,ring_atom=index[idx],atoms=[index[i] for i in atoms]) for idx,label,atoms in side]
 return result


def parse_name(name):
 text=name.strip();stereo={}
 if text.startswith('('):
  m=re.match(r'^\((\d+[RS](?:,\d+[RS])*)\)-',text)
  if not m:return None
  for loc,cip in re.findall(r'(\d+)([RS])',m[1]):
   if loc in stereo:return None
   stereo[loc]=cip
  text=text[m.end():]
 parent=re.search(r'cyclo(pent|hex)(an-(\d+)-(ol|one)|-(\d+)-en-(\d+)-(ol|one)|ane|-(\d+)-ene)$',text)
 if not parent:return None
 n=5 if parent[1]=='pent' else 6
 kind=('alcohol' if (parent[4] or parent[7])=='ol' else 'ketone') if parent[4] or parent[7] else None
 loc=int(parent[3] or parent[6]) if kind else None
 doubles=[int(parent[5] or parent[8])] if parent[5] or parent[8] else []
 prefix=text[:parent.start()];groups=[];errors=[];seen_labels=set()
 pattern=re.compile(r'(\d+(?:,\d+)*)-(?:(di|tri)?(hydroxymethyl|hydroxy|oxo|methyl|ethyl|propyl)|\((2-hydroxypropan-2-yl|propan-2-yl|prop-1-en-2-yl)\))(?:-|$)')
 while prefix:
  m=pattern.match(prefix)
  if not m:return None
  locs=list(map(int,m[1].split(',')));label=m[3] or m[4];amount={None:1,'di':2,'tri':3}[m[2]]
  if label in seen_labels:errors.append('同一置換基の位置をまとめて記載します')
  seen_labels.add(label)
  # Geminal substituents may repeat a locant.
  if len(locs)!=amount or locs!=sorted(locs):errors.append('倍数または昇順の不整合')
  groups.extend((i,label) for i in locs);prefix=prefix[m.end():]
 return dict(n=n,kind=kind,functional=loc,doubles=doubles,groups=sorted(groups),stereo=stereo,errors=errors)


def ring_monocycle_checks(analysis,name,parsed=None):
 if analysis['scope'] not in ('bounded','partial_stereo'):return []
 def row(key,title,status,detail):return dict(id='monocycle_'+key,title=title,status=status,detail=detail,reference='P-14.4 / P-22.1 / P-63.1 / P-64.1',url=SOURCE)
 p=parsed or parse_name(name)
 if p is None:return [row('grammar','単環骨格の限定系統名','unknown','有限構文外です。慣用名や別表記を推定しません。構造照合は別判定。')]
 parent=p['n']==analysis['ring_size'] and p['kind']==analysis['kind']
 same=p['functional']==analysis['functional_locant'] and p['doubles']==analysis['double_locants'] and p['groups']==analysis['groups']
 represented=parent and any(p['functional']==o['functional'] and p['doubles']==o['doubles'] and p['groups']==o['groups'] for o in analysis['orientations'])
 actual={x['position']:x['cip'] for x in analysis['ring_stereo'] if x['cip']}
 return [row('parent','単環の親骨格・主特性基','pass' if parent else 'fail','環炭素数 '+str(analysis['ring_size'])+' / 主特性基 '+str(analysis['kind'])+'。'),row('numbering','主特性基・不飽和・側鎖の番号','pass' if parent and same else ('fail' if parent else 'unknown'),'採用番号：主特性基 '+str(analysis['functional_locant'])+' / 二重結合 '+str(analysis['double_locants'])+' / 側鎖 '+str(analysis['groups'])),row('descriptor','単環と置換位置の構造対応','pass' if represented else 'fail','全方向の環番号と候補名の位置を照合。'),row('stereo','単環の絶対R/S','unknown' if analysis['scope']=='partial_stereo' else ('pass' if actual==p['stereo'] and parent and same else ('fail' if parent and same else 'unknown')),'入力 '+str(actual or '未指定')+' / 候補 '+str(p['stereo'] or '未指定')+'。部分立体は補完しません。'),row('spelling','単環名の有限表記','fail' if p['errors'] else ('pass' if name.strip()==analysis['expected_name'] else 'different'),'期待名：'+str(analysis['expected_name'])+'。表記差と構造不一致を区別。')]


EXTERNAL_NOTE='環に直結するCH2OHまたはC(OH)(CH3)2を1個追加した限定範囲。環上OH/ケトンは0/1個です。主特性基が環外OHのみならmethanol/propan-2-olを親とし、環上OHとの同順位または環上ケトンがある場合は環を親にします。全指定/全未指定のみ、環外カルボニル・長い含酸素側鎖・複数の環外OH・多環は保留。PIN全体は未判定。'


def analyze_monocycles(mol):
 result=analyze_ring_monocycles(mol)
 if not result['detected'] or result['scope']!='unsupported':return result
 original=Chem.Mol(mol)
 for a in original.GetAtoms():a.SetIntProp('_side_original',a.GetIdx())
 original=Chem.RemoveHs(original)
 if len(Chem.GetMolFrags(original))!=1 or original.GetStereoGroups() or any(a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() or a.GetAtomicNum() not in (6,8) or a.GetChiralTag() not in (Chem.ChiralType.CHI_UNSPECIFIED,Chem.ChiralType.CHI_TETRAHEDRAL_CW,Chem.ChiralType.CHI_TETRAHEDRAL_CCW) for a in original.GetAtoms()) or any(b.GetStereo()!=Chem.BondStereo.STEREONONE for b in original.GetBonds()):return result
 ring=set(original.GetRingInfo().AtomRings()[0]);external=[]
 for o in original.GetAtoms():
  if o.GetAtomicNum()!=8 or o.GetDegree()!=1 or o.GetTotalNumHs()!=1:continue
  c=o.GetNeighbors()[0]
  if c.GetIdx() in ring:continue
  if original.GetBondBetweenAtoms(o.GetIdx(),c.GetIdx()).GetBondType()!=Chem.BondType.SINGLE:continue
  anchors=[a.GetIdx() for a in c.GetNeighbors() if a.GetIdx() in ring]
  if len(anchors)!=1:continue
  other=[a for a in c.GetNeighbors() if a.GetIdx() not in ring and a.GetIdx()!=o.GetIdx()]
  if not other and c.GetTotalNumHs()==2:label='hydroxymethyl';chain='methanol';cap='methyl'
  elif len(other)==2 and c.GetDegree()==4 and all(a.GetAtomicNum()==6 and a.GetDegree()==1 and a.GetTotalNumHs()==3 and original.GetBondBetweenAtoms(c.GetIdx(),a.GetIdx()).GetBondType()==Chem.BondType.SINGLE for a in other):label='2-hydroxypropan-2-yl';chain='propan-2-ol';cap='propan-2-yl'
  else:continue
  external.append(dict(oxygen=o.GetIdx(),carbon=c.GetIdx(),anchor=anchors[0],label=label,chain=chain,cap=cap))
 if len(external)!=1:return result
 ext=external[0];copy=Chem.RWMol(original)
 for a in copy.GetAtoms():a.SetIntProp('_side_idx',a.GetIdx())
 c=copy.GetAtomWithIdx(ext['carbon']);c.SetNumExplicitHs(0);c.SetNoImplicit(False)
 copy.RemoveAtom(ext['oxygen'])
 try:Chem.SanitizeMol(copy)
 except (RuntimeError,ValueError):return result
 Chem.RemoveStereochemistry(copy)
 base=analyze_ring_monocycles(copy)
 if base['scope']!='bounded':return result
 index={a.GetIdx():a.GetIntProp('_side_idx') for a in copy.GetAtoms()}
 mode='ring' if base['kind'] else 'sidechain'
 centers=dict(Chem.FindMolChiralCenters(original,includeUnassigned=True,useLegacyImplementation=False))
 if any(i not in ring for i in centers) or any(a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED and a.GetIdx() not in centers for a in original.GetAtoms()):return result
 for a in original.GetAtoms():
  if a.HasProp('_CIPCode'):a.ClearProp('_CIPCode')
 try:rdCIPLabeler.AssignCIPLabels(original)
 except (RuntimeError,ValueError):return result
 options=[]
 for candidate in base['orientations']:
  path=[index[i] for i in candidate['path']];position=path.index(ext['anchor'])+1
  if mode=='sidechain' and position!=1:continue
  groups=list(candidate['groups']);groups.remove((position,ext['cap']))
  if mode=='ring':groups.append((position,ext['label']))
  groups.sort();doubles=candidate['doubles'];functional=candidate['functional'] if mode=='ring' else (1 if ext['chain']=='methanol' else 2)
  rank=((functional,) if mode=='ring' else (),tuple(doubles),group_rank([(i,re.sub(r'^\d+-','',label)) for i,label in groups]))
  stereo=[]
  for loc,idx in enumerate(path,1):
   if idx not in centers:continue
   a=original.GetAtomWithIdx(idx);cip=a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None
   if a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED and cip not in ('R','S'):return result
   stereo.append(dict(position=str(loc),atom_index=idx,cip=cip if cip in ('R','S') else None))
  assigned=[c for c in stereo if c['cip']];descriptor='('+','.join(c['position']+c['cip'] for c in assigned)+')-' if assigned else ''
  if mode=='ring':text=descriptor+base_name(base['ring_size'],base['kind'],functional,doubles,groups)
  else:
   root={5:'cyclopent',6:'cyclohex'}[base['ring_size']];prefix=prefixes(groups)
   if prefix is None:return result
   ring_name=descriptor+(prefix or '')+root+('-'+str(doubles[0])+'-en-1-yl' if doubles else 'an-1-yl')
   bracket='['+ring_name+']' if descriptor else '('+ring_name+')'
   text=('2-' if ext['chain']=='propan-2-ol' else '')+bracket+ext['chain']
  options.append(dict(path=path,groups=groups,doubles=doubles,functional=functional,rank=rank,stereo=stereo,text=text))
 if not options:return result
 best_rank=min(c['rank'] for c in options);best=[c for c in options if c['rank']==best_rank]
 if len({c['text'] for c in best})!=1:
  result['reason']='同順位の番号で立体記述子等が変わるため全体名を保留。';return result
 chosen=best[0];partial=0<sum(c['cip'] is not None for c in chosen['stereo'])<len(chosen['stereo'])
 restore=lambda i:original.GetAtomWithIdx(i).GetIntProp('_side_original')
 result.update(scope='partial_stereo' if partial else 'bounded',expected_name=None if partial else chosen['text'],parent_mode=mode,kind=base['kind'] if mode=='ring' else 'sidechain_alcohol',ring_size=base['ring_size'],numbering={str(i+1):restore(idx) for i,idx in enumerate(chosen['path'])},functional_locant=chosen['functional'],double_locants=chosen['doubles'],groups=chosen['groups'],ring_stereo=[dict(c,atom_index=restore(c['atom_index'])) for c in chosen['stereo']],stereo_status='partial' if partial else ('complete' if any(c['cip'] for c in chosen['stereo']) else 'unspecified'),orientations=[dict(functional=c['functional'],doubles=c['doubles'],groups=c['groups']) for c in options],external_hydroxyl={k:restore(ext[k]) for k in ('oxygen','carbon','anchor')},chain_parent=ext['chain'] if mode=='sidechain' else None,note=EXTERNAL_NOTE,reason=('環上の主特性基を接尾辞、環外OHを接頭辞にします。' if mode=='ring' else '主特性基を持つ環外鎖を親とし、環を置換基にします。')+('部分立体のため全体期待名を保留。' if partial else '環立体は元の全分子で計算。'),substituents=[])
 for g in base['substituents']:
  idx=index[g['ring_atom']];pos=str(chosen['path'].index(idx)+1)
  if ext['carbon'] in [index[i] for i in g['atoms']]:
   atoms=[restore(index[i]) for i in g['atoms']]+[restore(ext['oxygen'])];label=ext['label']
  else:atoms=[restore(index[i]) for i in g['atoms']];label=g['kind']
  result['substituents'].append(dict(position=pos,kind=label,ring_atom=restore(idx),atoms=atoms))
 return result


def parse_sidechain_name(name):
 match=re.fullmatch(r'(2-)?(?:\[(.+)\]|\((.+)\))(methanol|propan-2-ol)',name.strip())
 if not match or bool(match[1])!=(match[4]=='propan-2-ol'):return None
 inner=match[2] or match[3]
 if not inner.endswith('an-1-yl') and not re.search(r'-\d+-en-1-yl$',inner):return None
 text=re.sub(r'-(\d+)-en-1-yl$',r'-\1-ene',inner)
 text=re.sub(r'an-1-yl$','ane',text)
 parsed=parse_name(text)
 if parsed is None or parsed['kind']:return None
 parsed.update(kind='sidechain_alcohol',functional=2 if match[4]=='propan-2-ol' else 1,chain_parent=match[4],parent_mode='sidechain')
 return parsed


def monocycle_checks(analysis,name):
 if not analysis.get('external_hydroxyl'):return ring_monocycle_checks(analysis,name)
 parsed=parse_sidechain_name(name)
 if parsed is None:
  parsed=parse_name(name)
  if parsed:parsed.update(parent_mode='ring',chain_parent=None)
 if parsed is None:return ring_monocycle_checks(analysis,name)
 same=parsed['parent_mode']==analysis['parent_mode'] and parsed.get('chain_parent')==analysis.get('chain_parent')
 row=dict(id='monocycle_parent_selection',title='環外OHと親骨格の選択',status='pass' if same else 'fail',detail='入力の親：'+analysis['parent_mode']+' / 候補：'+parsed['parent_mode']+'。環上ケトンをOHより優先、OHが同順位なら環を選択。全体PINの証明ではありません。',reference='P-44.1.1 / P-44.1.2.2',url='https://iupac.qmul.ac.uk/BlueBook/P4.html')
 return [row]+ring_monocycle_checks(analysis,name,parsed)
