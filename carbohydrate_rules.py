"""Finite aldohexopyranose stereochemical references and O-linkage naming.
No NISPO or OPSIN decisions enter the reference registry. References derive
from Fischer OH patterns (IUPAC P-102 Table 10.2) and a mapped beta-D-glucose
pyranose, mirrored for L and inverted at C1 for alpha. Compare local tetrahedral
stereochemistry, never infer sugar identity/anomer from CIP labels alone.
"""
import re
from rdkit import Chem
from extended_rules import row,unknown_result

SOURCE='https://iupac.qmul.ac.uk/BlueBook/P10.html'
FISCHER={'allose':'RRRR','altrose':'LRRR','glucose':'RLRR','mannose':'LLRR','gulose':'RRLR','idose':'LRLR','galactose':'RLLR','talose':'LLLR'}
DEOXY_PARENTS={'quinovose':'glucose','rhamnose':'mannose','fucose':'galactose'}
AMINO_PARENTS={'glucosamine':'glucose','galactosamine':'galactose','mannosamine':'mannose'}
STEMS={name:name[:-3] for name in (*FISCHER,*DEOXY_PARENTS)}  # glucose -> gluc
BASE='[OH:16][CH2:6][C@H:5]1[O:15][C@@H:1]([OH:11])[C@H:2]([OH:12])[C@@H:3]([OH:13])[C@@H:4]1[OH:14]'
NOTE='限定したアルドヘキソピラノース8種・6-デオキシ糖3種・C2アミノ糖3種（遊離NH2・N-acetyl）の番号・局所立体参照との照合です。糖残基の識別と名称全体の適合性は別です。範囲外の糖、未指定立体、D/L・α/βの補完、PIN全体の保証は行いません。'


def core_key(mol,roles,ignore_c1=False):
    copy=Chem.RWMol(mol);keep=set(roles.values())
    for atom in copy.GetAtoms():atom.SetAtomMapNum(0)
    for label,idx in roles.items():
        atom=copy.GetAtomWithIdx(idx);atom.SetAtomMapNum(label)
        if label in (11,12,13,14,16):atom.SetNumExplicitHs(2 if atom.GetAtomicNum()==7 else 1);atom.SetNoImplicit(True)
        if label==1 and ignore_c1:atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for idx in reversed(range(copy.GetNumAtoms())):
        if idx not in keep:copy.RemoveAtom(idx)
    core=copy.GetMol();Chem.SanitizeMol(core)
    return Chem.MolToSmiles(core,canonical=True,isomericSmiles=True)


def references():
    result={};families={}
    for family,pattern in FISCHER.items():
        for series in ('D','L'):
            for anomer in ('alpha','beta'):
                mol=Chem.MolFromSmiles(BASE);roles={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()}
                for position in (2,3,4):
                    if pattern[position-2]!=FISCHER['glucose'][position-2]:mol.GetAtomWithIdx(roles[position]).InvertChirality()
                if series=='L':
                    for position in range(1,6):mol.GetAtomWithIdx(roles[position]).InvertChirality()
                if anomer=='alpha':mol.GetAtomWithIdx(roles[1]).InvertChirality()
                value=dict(family=family,series=series,anomer=anomer,stem=STEMS[family])
                result[core_key(mol,roles)]=value
                families[core_key(mol,roles,True)]=dict(family=family,series=series,stem=STEMS[family])
    return result,families

REGISTRY,FAMILY_REGISTRY=references()


def deoxy_references():
    # P-102.5.3.2: 6-deoxy glucose/mannose/galactose retain their own
    # names. C6 is not a stereocenter; removing its terminal O leaves all
    # C1-C5 tetrahedra unchanged. No naming engine supplies these references.
    result={};families={}
    for name,parent in DEOXY_PARENTS.items():
        for key,ref in REGISTRY.items():
            if ref['family']!=parent:continue
            mol=Chem.RWMol(Chem.MolFromSmiles(key))
            roles={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()}
            mol.GetAtomWithIdx(roles[6]).SetNumExplicitHs(3)
            mol.RemoveAtom(roles[16]);Chem.SanitizeMol(mol)
            roles={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()}
            value=dict(ref,family=name,stem=STEMS[name])
            result[core_key(mol,roles)]=value
            families[core_key(mol,roles,True)]=dict(family=name,series=ref['series'],stem=STEMS[name])
    return result,families

DEOXY_REGISTRY,DEOXY_FAMILY_REGISTRY=deoxy_references()


def amino_references():
    # P-102.5.4: replace the C2 OH by NH2 without changing the C2
    # tetrahedron. Only free, neutral NH2 and oxygenated C6 are referenced.
    result={};families={}
    for name,parent in AMINO_PARENTS.items():
        for key,ref in REGISTRY.items():
            if ref['family']!=parent:continue
            mol=Chem.MolFromSmiles(key)
            roles={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()}
            nitrogen=mol.GetAtomWithIdx(roles[12]);nitrogen.SetAtomicNum(7)
            nitrogen.SetNumExplicitHs(2);nitrogen.SetNoImplicit(True);Chem.SanitizeMol(mol)
            value=dict(ref,family=name,stem=STEMS[parent])
            result[core_key(mol,roles)]=value
            families[core_key(mol,roles,True)]=dict(family=name,series=ref['series'],stem=STEMS[parent])
    return result,families

AMINO_REGISTRY,AMINO_FAMILY_REGISTRY=amino_references()


def label(unit,ending='ose'):
    if not unit.get('family') or not unit.get('series'):return None
    prefix=(unit['anomer']+'-' if unit.get('anomer') else '')+unit['series']+'-'
    amino=unit['family'] in AMINO_PARENTS
    parent=AMINO_PARENTS.get(unit['family'],unit['family'])
    substituent='2-acetamido-2-deoxy-' if unit.get('amino_state')=='N_acetyl' else '2-amino-2-deoxy-'
    return (substituent if amino else '')+prefix+STEMS[parent]+'opyran'+ending


def external_oxygen(mol,idx,ring):
    candidates=[a for a in mol.GetAtomWithIdx(idx).GetNeighbors() if a.GetIdx() not in ring]
    if len(candidates)!=1:return None
    atom=candidates[0]
    if atom.GetAtomicNum()!=8 or atom.GetFormalCharge() or atom.GetNumRadicalElectrons():return None
    if any(b.GetBondType()!=Chem.BondType.SINGLE for b in atom.GetBonds()):return None
    if atom.GetDegree()==1 and atom.GetTotalNumHs()==1:return atom.GetIdx()
    if atom.GetDegree()==2 and all(a.GetAtomicNum()==6 for a in atom.GetNeighbors()):return atom.GetIdx()
    return None


def amino_substitution(mol,nitrogen,source):
    n=mol.GetAtomWithIdx(nitrogen)
    if n.GetAtomicNum()!=7 or n.GetFormalCharge() or n.GetNumRadicalElectrons() or n.IsInRing():return None
    if any(b.GetBondType()!=Chem.BondType.SINGLE for b in n.GetBonds()):return None
    others=[a for a in n.GetNeighbors() if a.GetIdx()!=source]
    if n.GetDegree()==1 and n.GetTotalNumHs()==2 and not others:return dict(kind='NH2',acetyl_atoms={})
    if n.GetDegree()!=2 or n.GetTotalNumHs()!=1 or len(others)!=1:return None
    carbon=others[0]
    if carbon.GetAtomicNum()!=6 or carbon.IsInRing() or carbon.GetDegree()!=3 or carbon.GetTotalNumHs()!=0 or carbon.GetFormalCharge() or carbon.GetNumRadicalElectrons():return None
    oxygen=[];methyl=[]
    for a in carbon.GetNeighbors():
        if a.GetIdx()==nitrogen:continue
        b=mol.GetBondBetweenAtoms(carbon.GetIdx(),a.GetIdx())
        if a.GetFormalCharge() or a.GetNumRadicalElectrons() or a.IsInRing():return None
        if a.GetAtomicNum()==8 and a.GetDegree()==1 and a.GetTotalNumHs()==0 and b.GetBondType()==Chem.BondType.DOUBLE:oxygen.append(a.GetIdx())
        elif a.GetAtomicNum()==6 and a.GetDegree()==1 and a.GetTotalNumHs()==3 and b.GetBondType()==Chem.BondType.SINGLE:methyl.append(a.GetIdx())
        else:return None
    if len(oxygen)!=1 or len(methyl)!=1:return None
    return dict(kind='N_acetyl',acetyl_atoms=dict(carbonyl_carbon=carbon.GetIdx(),carbonyl_oxygen=oxygen[0],methyl_carbon=methyl[0]))


def external_amino(mol,idx,ring):
    atoms=[a for a in mol.GetAtomWithIdx(idx).GetNeighbors() if a.GetIdx() not in ring]
    if len(atoms)!=1:return None
    a=atoms[0]
    return a.GetIdx() if amino_substitution(mol,a.GetIdx(),idx) else None


def number_ring(mol,ring):
    members=set(ring['ring_atoms'])
    if ring['parent']!='oxane':return None
    oxygen=ring['ring_oxygen']
    for first in mol.GetAtomWithIdx(oxygen).GetNeighbors():
        if first.GetIdx() not in members:continue
        path=[first.GetIdx()];prior=oxygen
        while len(path)<5:
            nxt=[a.GetIdx() for a in mol.GetAtomWithIdx(path[-1]).GetNeighbors() if a.GetIdx() in members and a.GetIdx()!=prior and a.GetIdx()!=oxygen]
            if len(nxt)!=1:break
            prior,current=path[-1],nxt[0];path.append(current)
        if len(path)!=5:continue
        if any(mol.GetAtomWithIdx(i).GetTotalNumHs()!=1 for i in path):continue
        roles={i:idx for i,idx in enumerate(path,1)};roles[15]=oxygen
        good=True
        for i in range(1,5):
            oxy=external_oxygen(mol,roles[i],members)
            if oxy is None and i==2:oxy=external_amino(mol,roles[i],members)
            if oxy is None:good=False;break
            roles[10+i]=oxy
        if not good:continue
        side=[a for a in mol.GetAtomWithIdx(path[4]).GetNeighbors() if a.GetIdx() not in members]
        if len(side)!=1 or side[0].GetAtomicNum()!=6 or side[0].IsInRing():continue
        roles[6]=side[0].GetIdx()
        if side[0].GetTotalNumHs()==2 and side[0].GetDegree()==2:
            oxy=external_oxygen(mol,roles[6],{roles[5]})
            if oxy is None:continue
            roles[16]=oxy
        elif not (side[0].GetTotalNumHs()==3 and side[0].GetDegree()==1):continue
        if len(set(roles.values()))!=len(roles):continue
        return roles
    return None


def aglycon_name(mol,oxygen,source,excluded):
    other=[a.GetIdx() for a in mol.GetAtomWithIdx(oxygen).GetNeighbors() if a.GetIdx()!=source]
    if len(other)!=1:return None
    stack=[other[0]];seen=set()
    while stack:
        idx=stack.pop()
        if idx in excluded or idx==oxygen or idx in seen:continue
        seen.add(idx);stack.extend(a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors())
    if any(mol.GetAtomWithIdx(i).GetAtomicNum()!=6 for i in seen):return None
    if len(seen)<=4 and all(mol.GetBondBetweenAtoms(i,a.GetIdx()).GetBondType()==Chem.BondType.SINGLE for i in seen for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() in seen):
        degrees={i:sum(a.GetIdx() in seen for a in mol.GetAtomWithIdx(i).GetNeighbors()) for i in seen}
        if degrees.get(other[0],99)<=1 and all(v<=2 for v in degrees.values()):return {1:'methyl',2:'ethyl',3:'propyl',4:'butyl'}[len(seen)]
    for name,smi in [('phenyl','c1ccccc1'),('benzyl','Cc1ccccc1')]:
        if any(set(match)==seen and match[0]==other[0] for match in mol.GetSubstructMatches(Chem.MolFromSmiles(smi),uniquify=False)):return name
    return None


def analyze_carbohydrates(mol,topology=None):
    if topology is None:
        from glycoside_analysis import analyze_glycoside
        topology=analyze_glycoside(mol)
    units=[]
    unsafe=bool(mol.GetStereoGroups()) or any(a.GetIsotope() or a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    unsafe |= any(b.GetBondType() not in (Chem.BondType.SINGLE,Chem.BondType.DOUBLE,Chem.BondType.TRIPLE,Chem.BondType.AROMATIC) for b in mol.GetBonds())
    for ring in topology['rings']:
        roles=number_ring(mol,ring)
        if roles is None:continue
        unit=dict(ring_id=ring['id'],carbon_atoms={str(i):roles[i] for i in range(1,7)},oxygen_atoms={str(i):roles[10+i] for i in (1,2,3,4,6) if 10+i in roles and mol.GetAtomWithIdx(roles[10+i]).GetAtomicNum()==8},amino_atoms={'2':roles[12]} if mol.GetAtomWithIdx(roles[12]).GetAtomicNum()==7 else {},deoxy_positions=([2] if mol.GetAtomWithIdx(roles[12]).GetAtomicNum()==7 else [])+([] if 16 in roles else [6]),ring_oxygen=roles[15],identity_status='unknown',anomer_status='unknown',family=None,series=None,anomer=None)
        substitution=amino_substitution(mol,roles[12],roles[2]) if unit['amino_atoms'] else None
        unit['amino_state']=substitution['kind'] if substitution else None
        unit['n_acetyl_atoms']=substitution['acetyl_atoms'] if substitution else {}
        if not unsafe:
            registry,families=(AMINO_REGISTRY,AMINO_FAMILY_REGISTRY) if unit['amino_atoms'] else ((REGISTRY,FAMILY_REGISTRY) if 16 in roles else (DEOXY_REGISTRY,DEOXY_FAMILY_REGISTRY))
            match=registry.get(core_key(mol,roles))
            family=families.get(core_key(mol,roles,True))
            if family:unit.update(family=family['family'],series=family['series'],identity_status='matched')
            if match:unit.update(anomer=match['anomer'],anomer_status='matched')
        unit['label']=label(unit);units.append(unit)
    indexed={u['ring_id']:u for u in units};links=[];aglycons=[]
    for bridge in topology['o_linkages']:
        directions=[]
        for d in bridge['directions']:
            source=indexed.get(d['source_ring_id'])
            if not source or source['carbon_atoms']['1']!=d['source_atom'] or source['oxygen_atoms']['1']!=bridge['bridge_oxygen']:continue
            if not d['targets']:
                aglycons.append(dict(source_ring_id=source['ring_id'],source_position=1,source_atom=d['source_atom'],bridge_oxygen=bridge['bridge_oxygen'],aglycon_atom=d['target_atom'],sugar_label=label(source,'osyl'),anomer=source['anomer']))
            for target in units:
                if target['ring_id']==source['ring_id']:continue
                loc=next((int(k) for k,v in target['carbon_atoms'].items() if v==d['target_atom']),None)
                if loc not in (1,2,3,4,6) or target['oxygen_atoms'].get(str(loc))!=bridge['bridge_oxygen']:continue
                directions.append(dict(source_ring_id=source['ring_id'],target_ring_id=target['ring_id'],source_position=1,target_position=loc,anomer=source['anomer'],series=source['series'],family=source['family']))
        if directions:links.append(dict(bridge_oxygen=bridge['bridge_oxygen'],bond_atoms=bridge['bond_atoms'],directions=directions))
    expected=None;name_kind=None
    if len(Chem.GetMolFrags(mol))==1 and len(units)==1:
        u=units[0];core=set(u['carbon_atoms'].values())|set(u['oxygen_atoms'].values())|set(u['amino_atoms'].values())|set(u['n_acetyl_atoms'].values())|{u['ring_oxygen']}
        others_free=all(mol.GetAtomWithIdx(idx).GetDegree()==1 for i,idx in u['oxygen_atoms'].items() if i!='1')
        if u['anomer'] and others_free:
            o1=u['oxygen_atoms']['1']
            if mol.GetAtomWithIdx(o1).GetDegree()==1 and len(core)==mol.GetNumAtoms():expected=label(u);name_kind='monosaccharide'
            else:
                aglycon=aglycon_name(mol,o1,u['carbon_atoms']['1'],core)
                if aglycon:expected=aglycon+' '+label(u,'oside');name_kind='simple_O_glycoside'
    if len(Chem.GetMolFrags(mol))==1 and len(units)==2 and len(links)==1:
        directions=links[0]['directions']
        all_core=set()
        for u in units:
            all_core.update(u['carbon_atoms'].values());all_core.update(u['oxygen_atoms'].values());all_core.update(u['amino_atoms'].values());all_core.update(u['n_acetyl_atoms'].values());all_core.add(u['ring_oxygen'])
        free=all(idx==links[0]['bridge_oxygen'] or mol.GetAtomWithIdx(idx).GetDegree()==1 for u in units for idx in u['oxygen_atoms'].values())
        complete=all(u['anomer'] for u in units) and free and len(all_core)==mol.GetNumAtoms()
        if complete and len(directions)==1:
            d=directions[0];source=indexed[d['source_ring_id']];target=indexed[d['target_ring_id']]
            if d['target_position'] in (2,3,4,6):
                if source['amino_atoms'] or target['amino_atoms']:
                    expected=label(source,'osyl')+'-(1->'+str(d['target_position'])+')-'+label(target)
                else:expected=str(d['target_position'])+'-O-'+label(source,'osyl')+'-'+label(target)
                name_kind='reducing_disaccharide'
        # Symmetric 1<->1 pairs only: choosing either residue as glycoside
        # gives the identical full name. No asymmetric parent priority guessed.
        if complete and len(directions)==2 and all(d['target_position']==1 for d in directions):
            first,second=units
            if all(first[k]==second[k] for k in ('family','series','anomer','amino_state')):
                expected=label(first,'osyl')+' '+label(second,'oside');name_kind='nonreducing_disaccharide'
    return dict(units=units,o_linkages=links,aglycon_o_linkages=aglycons,expected_name=expected,name_kind=name_kind,reference_count=len(REGISTRY)+len(DEOXY_REGISTRY)+len(AMINO_REGISTRY),reference_source=SOURCE,note=NOTE)


def normalize_carbohydrate_name(name):
    # This finite grammar normalizes the two P-102.7.1.2 linkage spellings,
    # not arbitrary synonyms, parent priorities, or general IUPAC syntax.
    canonical=name.strip().replace('α','alpha').replace('β','beta').replace('→','->')
    alias=re.fullmatch(r'N-acetyl-(alpha|beta)-(D|L)-(gluco|galacto|manno)samine',canonical)
    if alias:canonical='2-acetamido-2-deoxy-'+alias[1]+'-'+alias[2]+'-'+alias[3]+'pyranose'
    arrow=re.fullmatch(r'(.+?pyranosyl)-\(1->([1-6])\)-(.+?pyranose)',canonical)
    if arrow:canonical=arrow[2]+'-O-'+arrow[1]+'-'+arrow[3]
    grouped=re.fullmatch(r'([1-6])-O-\((.+?pyranosyl)\)-(.+?pyranose)',canonical)
    if grouped:canonical=grouped[1]+'-O-'+grouped[2]+'-'+grouped[3]
    return canonical


def carbohydrate_name_checks(mol,name):
    canonical=normalize_carbohydrate_name(name)
    sugar=r'(2-(?:amino|acetamido)-2-deoxy-)?(alpha|beta)-(D|L)-('+ '|'.join(stem+'o' for stem in STEMS.values())+r')pyran(ose|oside|osyl)'
    nonreducing=re.fullmatch(sugar.replace('(ose|oside|osyl)','osyl')+' '+sugar.replace('(ose|oside|osyl)','oside'),canonical)
    if not (nonreducing or re.fullmatch(sugar,canonical) or re.fullmatch(r'(methyl|ethyl|propyl|butyl|phenyl|benzyl) '+sugar,canonical) or re.fullmatch(r'[1-6]-O-'+sugar+'-'+sugar,canonical)):return None
    info=analyze_carbohydrates(mol)
    scope=dict(kind='carbohydrate_retained',scope_note=NOTE,carbohydrate_analysis=info)
    if not info['expected_name']:return unknown_result(scope,'この糖類名の限定検証には、対応する6員環の全立体指定・単純な単糖/O配糖体/還元性二糖、または同じ糖種・D/L・α/β同士の1↔1非還元性二糖が必要です。三糖以上、異なる糖種・D/L・α/β・N状態の1↔1、N-acetyl以外のN置換誘導体や不完全立体は保留します。')
    scope['expected']=info['expected_name'];same=canonical==normalize_carbohydrate_name(info['expected_name'])
    actual_tokens=re.findall(sugar,canonical);expected_tokens=re.findall(sugar,normalize_carbohydrate_name(info['expected_name']))
    def component_status(index):return 'pass' if [t[index] for t in actual_tokens]==[t[index] for t in expected_tokens] else 'fail'
    named_position=re.match(r'([1-6])-O-',canonical)
    expected_position=re.match(r'([1-6])-O-',normalize_carbohydrate_name(info['expected_name']))
    linkage_ok=bool(named_position and expected_position and named_position[1]==expected_position[1])
    disaccharide=info['name_kind'] in ('reducing_disaccharide','nonreducing_disaccharide')
    if info['name_kind']=='nonreducing_disaccharide':linkage_ok=bool(nonreducing)
    deoxy_stems={STEMS[f]+'o' for f in DEOXY_PARENTS}
    named_deoxy=[t[3] in deoxy_stems for t in actual_tokens];expected_deoxy=[t[3] in deoxy_stems for t in expected_tokens]
    deoxy_status=('pass' if named_deoxy==expected_deoxy else 'fail') if any(named_deoxy+expected_deoxy) else 'na'
    named_amino=[t[0] for t in actual_tokens];expected_amino=[t[0] for t in expected_tokens]
    amino_status=('pass' if named_amino==expected_amino else 'fail') if any(named_amino+expected_amino) else 'na'
    named_acetyl=['acetamido' in t[0] for t in actual_tokens];expected_acetyl=['acetamido' in t[0] for t in expected_tokens]
    acetyl_status=('pass' if named_acetyl==expected_acetyl else 'fail') if any(named_acetyl+expected_acetyl) else 'na'
    return dict(scope=scope,checks=[
      row('carbohydrate_name','糖種・D/L・α/β・結合位置を含む糖類表記','pass' if same else 'fail','限定期待糖類名：'+info['expected_name']+'。局所参照との立体照合と糖類C1〜C6の対応に基づきます。','P-102.3 / P-102.5 / P-102.6 / P-102.7','P10'),
      row('carbohydrate_family','糖の立体骨格（糖種）',component_status(3),str(actual_tokens)+' / 期待 '+str(expected_tokens),'P-102 Table 10.2 / P-102.5.3.2','P10'),
      row('carbohydrate_series','糖残基のD/L',component_status(2),'各糖残基のD/Lを局所立体参照と照合。旋光の符号やCIPの一つのR/Sからは推測しません。','P-102.3.2','P10'),
      row('carbohydrate_anomer','糖残基のα/β',component_status(1),'各糖残基のC1と配置参照原子を含む局所立体照合。未指定のアノマーは補完しません。','P-102.3.4','P10'),
      row('carbohydrate_deoxy','C6の6-デオキシ骨格',deoxy_status,'C6がCH3の残基：'+str([u['ring_id'] for u in info['units'] if 6 in u['deoxy_positions']])+'。ラムノース・フコース・キノボースの保持名を照合。その他の脱酸素位置や糖種は補完しません。','P-102.5.3.1–3','P10'),
      row('carbohydrate_amino','C2の窒素置換（NH2／N-acetyl）',amino_status,'C2–Nの原子対応：'+str([u['amino_atoms'] for u in info['units'] if u['amino_atoms']])+'。中性の遊離NH2とN-acetylを区別します。その他のN置換や塩は補完しません。','P-102.5.4','P10'),
      row('carbohydrate_n_acetyl','C2のN-アセチル基',acetyl_status,'N–C(=O)–CH3の原子対応：'+str([u['n_acetyl_atoms'] for u in info['units'] if u['n_acetyl_atoms']])+'。N上H1、末端CH3、単一カルボニルOを構造確認。O-acetyl、N-formyl/propionyl等は区別します。','P-102.5.6.4','P10'),
      row('carbohydrate_numbering','糖類番号と原子の対応','pass',str([u['carbon_atoms'] for u in info['units']]),'P-102.2','P10'),
      row('carbohydrate_links','糖残基間のO結合位置','pass' if disaccharide and linkage_ok else ('fail' if disaccharide else 'na'),str(info['o_linkages']),'P-102.7.1.1 / P-102.7.1.2','P10'),
      row('stereo','範囲外の糖類命名・PIN全体','unknown',NOTE,'P-102','P10')])
