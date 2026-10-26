"""Conservative input-graph annotations, independent of name/rule validation.
These are sugar-like motifs and O-linkage candidates, never sugar identities,
carbohydrate locants, anomer assignments, or a claim of complete detection.
"""
from rdkit import Chem
from ring_stereo_rules import cip_inventory

NOTE='入力構造の部分構造候補です。この候補表示だけでは糖類番号・糖種・D/L・α/β・名称との対応は決めません。糖類解析の結果は別欄に表示します。候補なしは配糖体ではないという意味ではありません。'


def analyze_glycoside(mol):
    all_rings=mol.GetRingInfo().AtomRings();rings=[]
    for raw in all_rings:
        atoms=set(raw)
        if len(raw) not in (5,6):continue
        if any(atoms&set(other) and atoms!=set(other) for other in all_rings):continue
        oxygen=[i for i in atoms if mol.GetAtomWithIdx(i).GetAtomicNum()==8]
        if len(oxygen)!=1:continue
        if any(mol.GetAtomWithIdx(i).GetAtomicNum() not in (6,8) or mol.GetAtomWithIdx(i).GetFormalCharge() or mol.GetAtomWithIdx(i).GetNumRadicalElectrons() or mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms):continue
        bonds=[b for b in mol.GetBonds() if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms]
        if len(bonds)!=len(atoms) or any(b.GetBondType()!=Chem.BondType.SINGLE for b in bonds):continue
        oxygenated=[];anomeric=[]
        for idx in sorted(atoms-{oxygen[0]}):
            atom=mol.GetAtomWithIdx(idx)
            exo=[a.GetIdx() for a in atom.GetNeighbors() if a.GetIdx() not in atoms and a.GetAtomicNum()==8 and not a.GetFormalCharge() and mol.GetBondBetweenAtoms(idx,a.GetIdx()).GetBondType()==Chem.BondType.SINGLE]
            if exo:oxygenated.append(idx)
            if len(exo)==1 and mol.GetBondBetweenAtoms(idx,oxygen[0]) is not None:anomeric.append(idx)
        if len(oxygenated)<2:continue
        rings.append(dict(ring_atoms=sorted(atoms),ring_oxygen=oxygen[0],parent='oxolane' if len(raw)==5 else 'oxane',oxygenated_carbon_atoms=oxygenated,anomeric_candidate_atoms=anomeric))
    stereo=cip_inventory(mol)
    for ring in rings:
        ring['stereo_centers']=[c for c in stereo['centers'] if c['atom_index'] in ring['ring_atoms']]
    rings.sort(key=lambda r:r['ring_atoms'])
    for i,ring in enumerate(rings,1):ring['id']=i
    links=[]
    for oxygen in mol.GetAtoms():
        if oxygen.GetAtomicNum()!=8 or oxygen.GetFormalCharge() or oxygen.GetDegree()!=2 or oxygen.IsInRing():continue
        neighbors=list(oxygen.GetNeighbors())
        if any(a.GetAtomicNum()!=6 or a.GetFormalCharge() for a in neighbors):continue
        if any(b.GetBondType()!=Chem.BondType.SINGLE for b in oxygen.GetBonds()):continue
        # Anomeric O-acyl / carbonate groups are outside this ether-linkage motif.
        if any(any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetAtomicNum() in (7,8,16) for b in a.GetBonds()) for a in neighbors):continue
        directions=[]
        for source,target in ((neighbors[0],neighbors[1]),(neighbors[1],neighbors[0])):
            for ring in rings:
                if source.GetIdx() not in ring['anomeric_candidate_atoms']:continue
                targets=[dict(ring_id=r['id'],attachment_atom=target.GetIdx(),type='ring_carbon') for r in rings if r['id']!=ring['id'] and target.GetIdx() in r['ring_atoms']]
                if not targets and target.GetTotalNumHs()==2 and target.GetDegree()==2:
                    for r in rings:
                        if r['id']==ring['id']:continue
                        for atom in target.GetNeighbors():
                            if atom.GetIdx() in r['ring_atoms'] and atom.GetAtomicNum()==6 and mol.GetBondBetweenAtoms(target.GetIdx(),atom.GetIdx()).GetBondType()==Chem.BondType.SINGLE:
                                targets.append(dict(ring_id=r['id'],attachment_atom=target.GetIdx(),ring_attachment_atom=atom.GetIdx(),type='exocyclic_methylene'))
                directions.append(dict(source_ring_id=ring['id'],source_atom=source.GetIdx(),target_atom=target.GetIdx(),targets=targets,target_type='sugar_like_ring' if targets else 'external_carbon'))
        if directions:links.append(dict(bridge_oxygen=oxygen.GetIdx(),bond_atoms=[neighbors[0].GetIdx(),oxygen.GetIdx(),neighbors[1].GetIdx()],directions=directions))
    # C-linkage candidates have no unique sugar C1 assignment here. A normal
    # terminal CH2OH is excluded; deoxy-sugar alkyl termini can still be ambiguous.
    carbon_links={};ambiguous_links={};side_chains=[]
    for ring in rings:
        members=set(ring['ring_atoms'])
        for source in mol.GetAtomWithIdx(ring['ring_oxygen']).GetNeighbors():
            if source.GetIdx() not in members or source.GetAtomicNum()!=6:continue
            for target in source.GetNeighbors():
                if target.GetIdx() in members or target.GetAtomicNum()!=6 or target.GetFormalCharge() or target.GetNumRadicalElectrons():continue
                bond=mol.GetBondBetweenAtoms(source.GetIdx(),target.GetIdx())
                if bond.GetBondType()!=Chem.BondType.SINGLE or bond.IsInRing():continue
                if any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(target).GetAtomicNum() in (7,8,16) for b in target.GetBonds()):continue
                side_oxygen=next((a for a in target.GetNeighbors() if a.GetAtomicNum()==8 and not a.GetFormalCharge() and mol.GetBondBetweenAtoms(target.GetIdx(),a.GetIdx()).GetBondType()==Chem.BondType.SINGLE),None)
                if target.GetTotalNumHs()==2 and target.GetDegree()==2 and not target.IsInRing() and side_oxygen is not None:
                    side_chains.append(dict(ring_id=ring['id'],ring_attachment_atom=source.GetIdx(),carbon_atom=target.GetIdx(),oxygen_atom=side_oxygen.GetIdx(),kind='hydroxymethyl' if side_oxygen.GetDegree()==1 and side_oxygen.GetTotalNumHs()==1 else 'o_substituted_methylene'))
                    continue
                is_terminal_methyl=target.GetTotalNumHs()==3 and target.GetDegree()==1
                if is_terminal_methyl:side_chains.append(dict(ring_id=ring['id'],ring_attachment_atom=source.GetIdx(),carbon_atom=target.GetIdx(),kind='methyl'))
                targets=[dict(ring_id=r['id'],attachment_atom=target.GetIdx(),type='ring_carbon') for r in rings if r['id']!=ring['id'] and target.GetIdx() in r['ring_atoms']]
                has_exo_oxygen=any(a.GetAtomicNum()==8 and a.GetIdx() not in members for a in source.GetNeighbors())
                collection=ambiguous_links if is_terminal_methyl else carbon_links
                entry=collection.setdefault(bond.GetIdx(),dict(bond_atoms=sorted([source.GetIdx(),target.GetIdx()]),directions=[]))
                entry['directions'].append(dict(source_ring_id=ring['id'],source_atom=source.GetIdx(),target_atom=target.GetIdx(),targets=targets,source_oxygen_retained=has_exo_oxygen))
    c_links=sorted(carbon_links.values(),key=lambda x:x['bond_atoms'])
    return dict(status='candidates' if rings else 'not_detected',rings=rings,o_linkages=links,c_linkages=c_links,ambiguous_c_linkages=sorted(ambiguous_links.values(),key=lambda x:x['bond_atoms']),ring_side_chains=side_chains,
                c_linkage_note='環Oに隣接する炭素からの環外C–C単結合の候補です。CH2OH・CH2–O側鎖は別表示とし、末端CH3のC–C結合は判別保留へ分けます。その他の側鎖と糖類のC1の同定は未判定です。C-糖の確定判定ではありません。',
                sugar_identity_status='unknown',anomer_status='unknown',nomenclature_status='unknown',note=NOTE)
