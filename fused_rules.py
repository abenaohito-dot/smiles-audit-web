"""Registered bicyclo[4.4.0]decane skeleton with bounded unsaturation; no general polycycle naming claim."""
import re
from collections import defaultdict
from extended_rules import group_rank
from monocycle_rules import side_label
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

SOURCE='https://iupac.qmul.ac.uk/BlueBook/P2.html'
CORE='[C:1]~1~2~[C:2]~[C:3]~[C:4]~[C:5]~[C:6]~1~[C:7]~[C:8]~[C:9]~[C:10]~2'
QUERY=Chem.MolFromSmarts(CORE)
NOTE='bicyclo[4.4.0]decane型骨格、メチル0〜4個とethyl/propan-2-ylのどちらか0/1個、環内二重結合0/1本、環上OHまたはケトン0/1個に限定。橋頭1・6と二つの4原子経路を固定し、4方向を比較します。全指定の大文字R/Sまたは全未指定、無置換デカリンの融合cis/transに対応。その他の小文字r/s、部分/相対/拡張立体、他側鎖、環外官能基、二重結合2本以上、結合E/Z指定、他の多環、保持名の優先性とPIN全体は保留。CIPはRDKitに依存。'
MULT={1:'',2:'di',3:'tri',4:'tetra'}


# Independent cis/trans references: the outer ten-membered perimeter is planar;
# explicit bridgehead H atoms lie on the same/opposite side of that plane.
# RDKit converts these signed tetrahedra to graph tags. No naming engine labels
# a reference, and no energy-minimized conformer is assumed.
def fusion_references():
 base=Chem.MolFromSmiles('C12CCCCC1CCCC2');mapping={i+1:i for i in range(10)}
 coords={1:(0,0,0),2:(-1,-.5,0),3:(-2,0,0),4:(-2,2,0),5:(-1,2.5,0),6:(0,2,0),7:(1,2.5,0),8:(2,2,0),9:(2,0,0),10:(1,-.5,0)}
 refs={}
 for label,z in [('cis',1),('trans',-1)]:
  m=Chem.AddHs(base,onlyOnAtoms=[mapping[1],mapping[6]]);conf=Chem.Conformer(m.GetNumAtoms())
  for n,idx in mapping.items():conf.SetAtomPosition(idx,coords[n])
  for n,hz in ((1,1),(6,z)):
   h=next(a.GetIdx() for a in m.GetAtomWithIdx(mapping[n]).GetNeighbors() if a.GetAtomicNum()==1)
   x,y,_=coords[n];conf.SetAtomPosition(h,(x,y,hz))
  m.AddConformer(conf);Chem.AssignAtomChiralTagsFromStructure(m)
  for a in m.GetAtoms():
   a.SetAtomMapNum(0)
   if a.GetIdx() not in (mapping[1],mapping[6]):a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
  refs[label]=Chem.RemoveHs(m)
 return refs

FUSION_REFERENCES=fusion_references()


def double_text(pair):
 if not pair:return None
 a,b=pair
 return str(a) if b-a==1 else str(a)+'('+str(b)+')'


def double_rank(pair):
 # P-31.1.4.2: minimize compound locants before comparing their leading
 # numbers; bracketed numbers are considered only at the final tie-break.
 if not pair:return (0,(),())
 a,b=pair
 return (int(b-a!=1),(a,),tuple(pair))


def assemble(groups,kind,locant,double_bond=None):
 if groups and isinstance(groups[0],int):groups=[(n,'methyl') for n in groups]
 grouped=defaultdict(list)
 for n,label in groups:grouped[label].append(n)
 terms=[]
 for label,locs in sorted(grouped.items()):
  if label=='propan-2-yl':
   if len(locs)!=1:return None
   term='(propan-2-yl)'
  else:term=MULT[len(locs)]+label
  terms.append(','.join(map(str,sorted(locs)))+'-'+term)
 prefix='-'.join(terms)
 if double_bond:
  parent='bicyclo[4.4.0]dec-'+double_text(double_bond)+('-en-'+str(locant)+'-'+('ol' if kind=='alcohol' else 'one') if kind else '-ene')
 else:parent='bicyclo[4.4.0]deca'+('n-'+str(locant)+'-'+('ol' if kind=='alcohol' else 'one') if kind else 'ne')
 return prefix+parent


def analyze_fused(mol):
 result=dict(detected=False,scope='unsupported',expected_name=None,numbering={},bridgeheads={},fusion_bond=[],substituents=[],ring_stereo=[],reason='',note=NOTE,source=SOURCE)
 mol=Chem.Mol(mol)
 for a in mol.GetAtoms():a.SetIntProp('_fused_original',a.GetIdx())
 mol=Chem.RemoveHs(mol);rings=mol.GetRingInfo().AtomRings()
 if len(rings)!=2 or any(len(r)!=6 for r in rings):return result
 shared=set(rings[0])&set(rings[1]);core=set(rings[0])|set(rings[1])
 if len(shared)!=2 or len(core)!=10 or mol.GetBondBetweenAtoms(*sorted(shared)) is None or any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 or mol.GetAtomWithIdx(i).GetIsAromatic() for i in core):return result
 result.update(detected=True,reason='デカリン型候補。側鎖・環内二重結合・環上OH/ケトンの限定範囲を確認します。')
 if len(Chem.GetMolFrags(mol))!=1 or mol.GetStereoGroups() or any(a.GetAtomicNum() not in (6,8) or a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() or a.GetChiralTag() not in (Chem.ChiralType.CHI_UNSPECIFIED,Chem.ChiralType.CHI_TETRAHEDRAL_CW,Chem.ChiralType.CHI_TETRAHEDRAL_CCW) for a in mol.GetAtoms()) or any(b.GetStereo()!=Chem.BondStereo.STEREONONE or b.GetBondType() not in (Chem.BondType.SINGLE,Chem.BondType.DOUBLE) for b in mol.GetBonds()):return result
 core_doubles=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds() if b.GetBeginAtomIdx() in core and b.GetEndAtomIdx() in core and b.GetBondType()==Chem.BondType.DOUBLE]
 if len(core_doubles)>1:return result
 oxy=[a for a in mol.GetAtoms() if a.GetAtomicNum()==8]
 if len(oxy)>1:return result
 kind=None;functional=None
 if oxy:
  o=oxy[0]
  if o.GetDegree()!=1 or o.GetNeighbors()[0].GetIdx() not in core:return result
  functional=o.GetNeighbors()[0].GetIdx();bond=mol.GetBondBetweenAtoms(o.GetIdx(),functional)
  if bond.GetBondType()==Chem.BondType.SINGLE and o.GetTotalNumHs()==1:kind='alcohol'
  elif bond.GetBondType()==Chem.BondType.DOUBLE and o.GetTotalNumHs()==0:kind='ketone'
  else:return result
 sides=[];used=core|{a.GetIdx() for a in oxy}
 for idx in core:
  for a in mol.GetAtomWithIdx(idx).GetNeighbors():
   first=a.GetIdx()
   if first in used:continue
   if mol.GetBondBetweenAtoms(idx,first).GetBondType()!=Chem.BondType.SINGLE:return result
   atoms=set();stack=[first]
   while stack:
    current=stack.pop()
    if current in core or current in atoms:continue
    atoms.add(current);stack.extend(nb.GetIdx() for nb in mol.GetAtomWithIdx(current).GetNeighbors())
   if atoms.intersection(used) or not 1<=len(atoms)<=3:return result
   label=side_label(mol,atoms,first)
   if label not in ('methyl','ethyl','propan-2-yl'):return result
   sides.append((idx,label,sorted(atoms)));used.update(atoms)
 if sum(label=='methyl' for _,label,_ in sides)>4 or sum(label!='methyl' for _,label,_ in sides)>1 or used!=set(range(mol.GetNumAtoms())):return result
 options=[]
 for match in mol.GetSubstructMatches(QUERY,uniquify=False,maxMatches=256):
  mapping={a.GetAtomMapNum():match[a.GetIdx()] for a in QUERY.GetAtoms()}
  if set(match)!=core or {mapping[1],mapping[6]}!=shared:continue
  pos={idx:n for n,idx in mapping.items()};groups=sorted((pos[i],label) for i,label,_ in sides);loc=pos[functional] if kind else None
  double=tuple(sorted(pos[i] for i in core_doubles[0])) if core_doubles else None
  options.append(dict(mapping=mapping,groups=groups,functional=loc,double=double,base=assemble(groups,kind,loc,double),rank=((loc,) if kind else (),double_rank(double),group_rank(groups))))
 if len(options)!=4:
  result['reason']='登録した4通りの番号付けを完了できないため全体名を保留。';return result
 rank=min(c['rank'] for c in options);best=[c for c in options if c['rank']==rank]
 centers=dict(Chem.FindMolChiralCenters(mol,includeUnassigned=True,useLegacyImplementation=False))
 if any(i not in core for i in centers):return result
 for a in mol.GetAtoms():
  if a.HasProp('_CIPCode'):a.ClearProp('_CIPCode')
 cip_error=False
 try:rdCIPLabeler.AssignCIPLabels(mol)
 except (RuntimeError,ValueError):cip_error=True
 index={a.GetIdx():a.GetIntProp('_fused_original') for a in mol.GetAtoms()}
 alternatives=[]
 for c in best:
  stereo=[];unsupported=cip_error
  for n,idx in sorted(c['mapping'].items()):
   a=mol.GetAtomWithIdx(idx)
   if idx not in centers:
    if a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:unsupported=True
    continue
   cip=a.GetProp('_CIPCode') if a.HasProp('_CIPCode') and not cip_error else None
   if a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED and cip not in ('R','S'):unsupported=True
   stereo.append(dict(position=str(n),atom_index=index[idx],cip=cip if cip in ('R','S','r','s') else None))
  assigned=[s for s in stereo if s['cip']]
  partial=0<len(assigned)<len(stereo)
  descriptor='('+','.join(s['position']+s['cip'] for s in assigned)+')-' if assigned else ''
  alternatives.append((c,stereo,descriptor,unsupported,partial))
 c,stereo,descriptor,unsupported,partial=alternatives[0]
 fusion=None
 if not sides and not kind and not core_doubles and not cip_error:
  fusion=next((label for label,ref in FUSION_REFERENCES.items() if mol.HasSubstructMatch(ref,useChirality=True)),None)
  if fusion:descriptor=fusion+'-';unsupported=False;partial=False
 ambiguous=len({d+c['base'] for c,_,d,_,_ in alternatives})!=1
 if fusion:ambiguous=False
 scope='unsupported_stereo' if unsupported or ambiguous else ('partial_stereo' if partial else 'bounded')
 reason='橋長4・4・0、橋頭1・6を固定した4方向で、主特性基→複合位置番号→二重結合→側鎖の番号を比較。'
 if unsupported:reason+='小文字r/s、未確定の指定立体またはCIP計算失敗のため全体名を保留。r/sをR/Sへ読み替えません。'
 elif ambiguous:reason+='同順位で立体番号が変わるため全体名を保留。'
 elif partial:reason+='一部の立体が未指定のため全体名を保留。'
 else:reason+=('橋頭Hの同側/反対側の独立参照と照合した'+fusion+'配置。' if fusion else '全指定R/Sまたは全未指定を保持。未指定はラセミ体と断定しません。')
 result.update(scope=scope,fusion_configuration=fusion,kind=kind,functional_atom=index[functional] if functional is not None else None,functional_oxygen=index[oxy[0].GetIdx()] if oxy else None,functional_locant=c['functional'],methyl_locants=[n for n,label in c['groups'] if label=='methyl'],side_groups=c['groups'],double_bond_locants=list(c['double']) if c['double'] else [],double_bond_atoms=[index[i] for i in core_doubles[0]] if core_doubles else [],double_bond_label=double_text(c['double']),expected_name=descriptor+c['base'] if scope=='bounded' else None,base_expected_name=c['base'],numbering={str(n):index[i] for n,i in sorted(c['mapping'].items())},bridgeheads={str(n):index[c['mapping'][n]] for n in (1,6)},fusion_bond=[index[c['mapping'][n]] for n in (1,6)],ring_stereo=stereo,stereo_status='fusion_relative' if fusion else ('unsupported' if unsupported or ambiguous else ('partial' if partial else ('complete' if descriptor else 'unspecified'))),reason=reason,orientations=[dict(numbering={str(n):index[i] for n,i in sorted(o['mapping'].items())},functional=o['functional'],methyl=[n for n,label in o['groups'] if label=='methyl'],side_groups=o['groups'],double_bond=list(o['double']) if o['double'] else []) for o in options])
 inverse={idx:n for n,idx in c['mapping'].items()}
 result['substituents']=[dict(position=str(inverse[i]),kind=label,ring_atom=index[i],atoms=[index[j] for j in atoms]) for i,label,atoms in sides]
 return result


def parse_fused_name(name):
 text=name.strip();stereo={};errors=[];fusion=None
 for label in ('cis','trans'):
  if text.startswith(label+'-'):fusion=label;text=text[len(label)+1:];break
 if text.startswith('('):
  m=re.match(r'^\((\d+[RS](?:,\d+[RS])*)\)-',text)
  if not m:return None
  tokens=re.findall(r'(\d+)([RS])',m[1])
  for loc,cip in tokens:
   if loc in stereo:errors.append('立体位置の重複')
   stereo[loc]=cip
  if [int(n) for n,_ in tokens]!=sorted(int(n) for n,_ in tokens):errors.append('立体位置が昇順ではありません')
  text=text[m.end():]
 parent=re.search(r'bicyclo\[(\d+)\.(\d+)\.(\d+)\](oct|non|dec|undec)(.+)$',text)
 if not parent:return None
 prefix=text[:parent.start()];groups=[];seen=set()
 pattern=re.compile(r'(\d+(?:,\d+)*)-(?:(di|tri|tetra)?(methyl|ethyl)|\((propan-2-yl)\))(?:-|$)')
 while prefix:
  m=pattern.match(prefix)
  if not m:return None
  locs=list(map(int,m[1].split(',')));label=m[3] or m[4]
  amount=next(n for n,v in MULT.items() if v==(m[2] or ''))
  if len(locs)!=amount or locs!=sorted(locs):errors.append('側鎖の倍数または位置の昇順が不一致')
  if label in seen:errors.append('同種の側鎖は位置をまとめて記載します')
  seen.add(label);groups.extend((n,label) for n in locs);prefix=prefix[m.end():]
 tail=parent[5];kind=None;functional=None;double=[]
 saturated=re.fullmatch(r'an-(\d+)-(ol|one)',tail)
 unsaturated=re.fullmatch(r'-(\d+)(?:\((\d+)\))?-(ene|en-(\d+)-(ol|one))',tail)
 if tail=='ane':pass
 elif saturated:kind='alcohol' if saturated[2]=='ol' else 'ketone';functional=int(saturated[1])
 elif unsaturated:
  first=int(unsaturated[1]);second=int(unsaturated[2]) if unsaturated[2] else first+1;double=[first,second]
  if unsaturated[2] and second-first<=1:errors.append('複合位置番号は非連続の端点を昇順に示します')
  if unsaturated[4]:kind='alcohol' if unsaturated[5]=='ol' else 'ketone';functional=int(unsaturated[4])
 else:return None
 return dict(fusion=fusion,bridges=tuple(map(int,parent.group(1,2,3))),root=parent[4],kind=kind,functional=functional,methyl=sorted(n for n,label in groups if label=='methyl'),side_groups=sorted(groups),double_bond=double,stereo=stereo,errors=errors)


def fused_checks(analysis,name):
 if analysis['scope'] not in ('bounded','partial_stereo'):return []
 def row(key,title,status,detail):return dict(id='fused_'+key,title=title,status=status,detail=detail,reference='P-23.2.2 / P-23.2.3 / P-31.1.4.2 / P-14.4 / P-63 / P-64',url=SOURCE)
 p=parse_fused_name(name)
 if p is None:return [row('grammar','縮合二環の限定系統名','unknown','この候補名は有限構文外です。decalin等の保持名や別表記は推定せず、構造照合を別判定します。')]
 parent=p['bridges']==(4,4,0) and p['root']=='dec';kind=p['kind']==analysis['kind'];same=p['functional']==analysis['functional_locant'] and p['side_groups']==analysis['side_groups'] and p['double_bond']==analysis['double_bond_locants']
 represented=parent and kind and any(p['functional']==o['functional'] and p['side_groups']==o['side_groups'] and p['double_bond']==o['double_bond'] for o in analysis['orientations'])
 actual={s['position']:s['cip'] for s in analysis['ring_stereo'] if s['cip']}
 stereo_same=(p['fusion']==analysis.get('fusion_configuration') and not p['stereo']) if analysis.get('fusion_configuration') else (actual==p['stereo'] and not p['fusion'])
 return [row('parent','縮合二環の橋長・炭素数','pass' if parent else 'fail','登録親骨格はbicyclo[4.4.0]decane型。二つの4原子経路と直接の融合結合、計10環炭素を対応。'),row('suffix','縮合二環の主特性基・接尾辞','pass' if parent and kind else ('fail' if parent else 'unknown'),'入力の主特性基：'+str(analysis['kind'])+'。OHはol、ケトンはone。'),row('numbering','橋頭を固定した位置番号','pass' if parent and kind and same else ('fail' if parent and kind else 'unknown'),'採用：主特性基 '+str(analysis['functional_locant'])+' / 二重結合 '+str(analysis['double_bond_label'] or 'なし')+' / 側鎖 '+str(analysis['side_groups'])+'。'),row('descriptor','縮合二環と置換位置の構造対応','pass' if represented else 'fail','4通りの許容番号と候補名を照合。構造が一致する番号と優先する番号を区別。'),row('stereo','融合部の相対配置' if analysis.get('fusion_configuration') else '融合部を含む全環R/S','unknown' if analysis['scope']=='partial_stereo' or not(parent and kind and same) else ('pass' if stereo_same else 'fail'),'入力 '+str(analysis.get('fusion_configuration') or actual or '未指定')+' / 候補 '+str(p['fusion'] or p['stereo'] or '未指定')+'。r/sや未指定をR/Sへ補完しません。'),row('spelling','縮合二環名の有限表記','fail' if p['errors'] else ('pass' if name.strip()==analysis['expected_name'] else 'different'),' / '.join(p['errors']) if p['errors'] else '限定期待名：'+str(analysis['expected_name'])+'。PIN全体は未判定。')]
