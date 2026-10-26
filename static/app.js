'use strict';
const $ = id => document.getElementById(id);
const labels = {pass:'適合・一致', fail:'違反・不一致', unknown:'未判定', different:'期待表記との差', na:'該当なし'};
let report = null;
function setImage(id, svg) {
  if(!svg) { $(id).hidden = true; $(id).removeAttribute('src'); return; }
  // SVG comes only from local RDKit; render as an image, never inject markup into DOM.
  const url = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
  $(id).src = url; $(id).hidden = false;
}
function render(data) {
  report = data;
  $('empty').hidden = true;
  $('empty').style.display = 'none';
  $('result').hidden = false;
  $('name').textContent = data.name || '命名できませんでした';
  $('source').textContent = data.name_source === 'manual' ? 'SUPPLIED NAME / 手入力名の検証' : (data.name_source === 'bounded_rules' ? 'BOUNDED RULES / 限定規則で命名・構造一致（NISPO命名は未完了）' : 'NISPO / GENERATED NAME');
  $('formula').textContent = `${data.formula} · ${data.canonical_input}`;
  const flav=data.flavonoid_analysis;
  $('flavonoid-result').hidden = !flav?.detected;
  $('flavonoid-summary').textContent = flav?.detected ? (flav.label || flav.candidates.map(c=>c.label).join(' / ')) + ' · ' + flav.reason : '';
  $('flavonoid-expected').textContent = '限定期待名：' + (flav?.expected_name || '未判定（骨格候補の識別のみ）');
  $('flavonoid-roundtrip').hidden = !flav?.expected_name;
  $('flavonoid-roundtrip').textContent = '期待フラボノイド名の往復照合：' + (labels[data.flavonoid_expected_roundtrip_status] || '未判定');
  $('flavonoid-alkene').hidden = !flav?.alkene_atoms?.length;
  $('flavonoid-alkene').textContent = '二重結合：原子 ' + (flav?.alkene_atoms || []).join(' ＝ 原子 ') + ' / E/Z：' + (flav?.alkene_stereo || '未指定（補完しません）');
  $('flavonoid-parent').hidden = !flav?.parent_name;
  $('flavonoid-parent').textContent = '親骨格：' + (flav?.parent_name || '') + (flav?.parent_carbonyl===2 ? '（ラクトン骨格）' : '');
  $('flavonoid-stereo').hidden = !flav?.ring_stereo?.length;
  $('flavonoid-stereo').textContent = '環立体：' + (flav?.ring_stereo || []).map(c=>c.position+'位 '+(c.cip || '未指定')+'（原子 '+c.atom_index+'）').join(' / ') + (flav?.stereo_status==='partial' ? ' / 部分指定のため全体期待名は保留' : '');
  $('flavonoid-map').textContent = '骨格番号 → 原子index：' + Object.entries(flav?.numbering || {}).map(([n,idx])=>n+' → '+idx).join(' / ');
  if((flav?.equivalent_numberings?.length || 0)>1) $('flavonoid-map').textContent += '（対称性により同等の原子対応あり）';
  $('flavonoid-groups').replaceChildren();
  for(const group of flav?.substituents || []) {
    const line=document.createElement('p'); line.textContent=group.position+'位 · '+group.kind+' / 骨格原子 '+group.ring_atom+' / 置換基の原子index ['+group.atoms.join(', ')+']';
    $('flavonoid-groups').append(line);
  }
  const flavGlyco=flav?.glycoside;
  $('flavonoid-glycoside').hidden=!flavGlyco;
  $('flavonoid-glycoside').textContent=flavGlyco ? flavGlyco.kind+'配糖体 · '+flavGlyco.sugar_label+' / 糖C1（原子 '+flavGlyco.anomeric_atom+'）'+(flavGlyco.kind==='O' ? ' → O（原子 '+flavGlyco.bridge_oxygen+'）' : '')+' → アグリコン'+flavGlyco.position+'位（原子 '+flavGlyco.aglycon_atom+'） / 糖の番号：'+Object.entries(flavGlyco.carbon_atoms).map(([n,i])=>'C'+n+' → '+i).join(' / ') : '';
  $('flavonoid-note').textContent=flav?.note || '';
  const mono=data.monocycle_analysis;
  $('monocycle-result').hidden=!mono?.detected;
  $('monocycle-summary').textContent=mono?.reason || '';
  $('monocycle-expected').textContent='限定期待名：'+(mono?.expected_name || '未判定（対応範囲・立体指定を確認してください）');
  $('monocycle-roundtrip').hidden=!mono?.expected_name;
  $('monocycle-roundtrip').textContent='期待単環骨格名の往復照合：'+(labels[data.monocycle_expected_roundtrip_status] || '未判定');
  $('monocycle-map').textContent='環の番号 → 原子index：'+Object.entries(mono?.numbering || {}).map(([n,i])=>n+' → '+i).join(' / ');
  $('monocycle-stereo').textContent='環立体：'+((mono?.ring_stereo || []).map(c=>c.position+'位 '+(c.cip || '未指定')+'（原子 '+c.atom_index+'）').join(' / ') || '指定なし・未判定')+(mono?.scope==='partial_stereo' ? ' / 部分指定のため全体期待名は保留' : '');
  $('monocycle-groups').replaceChildren();
  for(const group of mono?.substituents || []) {
    const line=document.createElement('p');line.textContent=group.position+'位 · '+group.kind+' / 環原子 '+group.ring_atom+' / 置換基の原子index ['+group.atoms.join(', ')+']';$('monocycle-groups').append(line);
  }
  $('monocycle-parent').hidden=!mono?.external_hydroxyl;
  $('monocycle-parent').textContent=mono?.external_hydroxyl ? '主骨格：'+(mono.parent_mode==='sidechain' ? mono.chain_parent+'（環は置換基）' : '環（環外OHは接頭辞）')+' / 環外OHのO原子 '+mono.external_hydroxyl.oxygen+' → C原子 '+mono.external_hydroxyl.carbon+' → 環原子 '+mono.external_hydroxyl.anchor : '';
  $('monocycle-note').textContent=mono?.note || '';
  const bicy=data.fused_analysis;
  $('fused-result').hidden=!bicy?.detected;
  $('fused-summary').textContent=bicy?.reason || '';
  $('fused-expected').textContent='限定期待名：'+(bicy?.expected_name || '未判定（対応範囲・立体指定を確認してください）');
  $('fused-roundtrip').hidden=!bicy?.expected_name;
  $('fused-roundtrip').textContent='期待縮合二環名の往復照合：'+(labels[data.fused_expected_roundtrip_status] || '未判定');
  $('fused-map').textContent='環の番号 → 原子index：'+Object.entries(bicy?.numbering || {}).map(([n,i])=>n+' → '+i).join(' / ');
  $('fused-fusion').textContent='橋頭原子（1・6位）：'+(bicy?.fusion_bond || []).join(' ↔ ')+' / 融合部：'+(bicy?.fusion_configuration || 'cis/transの判定対象外・未指定')+(bicy?.kind ? ' / 主特性基：'+(bicy.kind==='alcohol' ? 'OH' : 'ケトン')+'（'+bicy.functional_locant+'位、C原子 '+bicy.functional_atom+' → O原子 '+bicy.functional_oxygen+'）' : '');
  $('fused-stereo').textContent='CIP表示：'+((bicy?.ring_stereo || []).map(c=>c.position+'位 '+(c.cip || '未指定')+'（原子 '+c.atom_index+'）').join(' / ') || '未判定')+(bicy?.scope==='partial_stereo' ? ' / 部分指定のため全体期待名は保留' : '');
  $('fused-groups').replaceChildren();
  for(const group of bicy?.substituents || []) {
    const line=document.createElement('p');line.textContent=group.position+'位 · '+group.kind+' / 環原子 '+group.ring_atom+' / 置換基の原子index ['+group.atoms.join(', ')+']';$('fused-groups').append(line);
  }
  $('fused-double').textContent=bicy?.double_bond_atoms?.length ? '環内二重結合：'+bicy.double_bond_atoms.join('＝')+' / 骨格位置 '+(bicy.double_bond_label || '未判定') : '';
  $('fused-note').textContent=bicy?.note || '';
  const bridge=data.bridged_analysis;
  $('bridged-result').hidden=!bridge?.detected;
  $('bridged-summary').textContent=bridge?.reason || '';
  $('bridged-expected').textContent='限定期待名：'+(bridge?.expected_name || '未判定（対応範囲・立体指定を確認してください）');
  $('bridged-roundtrip').hidden=!bridge?.expected_name;
  $('bridged-roundtrip').textContent='期待架橋二環名の往復照合：'+(labels[data.bridged_expected_roundtrip_status] || '未判定');
  $('bridged-map').textContent='骨格番号 → 原子index：'+Object.entries(bridge?.numbering || {}).map(([n,i])=>n+' → '+i).join(' / ');
  $('bridged-fusion').textContent='橋頭原子（'+Object.keys(bridge?.bridgeheads || {}).join('・')+'位）：'+Object.values(bridge?.bridgeheads || {}).join(' ↔ ')+(bridge?.kind ? ' / 主特性基：'+(bridge.kind==='alcohol' ? 'OH' : 'ケトン')+'（'+bridge.functional_locant+'位、C原子 '+bridge.functional_atom+' → O原子 '+bridge.functional_oxygen+'）' : '');
  $('bridged-double').textContent='環内二重結合：'+(bridge?.double_bond_locants?.length ? bridge.double_bond_locants.join('–')+'位 / 原子 '+bridge.double_bond_atoms.join(' ↔ ') : 'なし');
  $('bridged-stereo').textContent='CIP表示：'+((bridge?.ring_stereo || []).map(c=>c.position+'位 '+(c.cip || '未指定')+'（原子 '+c.atom_index+'）').join(' / ') || '未判定')+(bridge?.scope==='partial_stereo' ? ' / 部分指定のため全体期待名は保留' : '');
  $('bridged-groups').replaceChildren();
  for(const group of bridge?.substituents || []) { const line=document.createElement('p');line.textContent=group.position+'位 · methyl / 骨格原子 '+group.ring_atom+' / メチル原子index ['+group.atoms.join(', ')+']';$('bridged-groups').append(line); }
  $('bridged-note').textContent=bridge?.note || '';
  const carb=data.carbohydrate_analysis;
  $('carbohydrate-result').hidden = !carb?.units?.length;
  $('carbohydrate-expected').textContent = flavGlyco ? '配糖体全体の期待名・構造照合はフラボノイド欄に表示しています。' : '限定期待糖類名：' + (carb?.expected_name || '未判定（置換・立体・構造が命名範囲外）');
  $('carbohydrate-roundtrip').hidden = !carb?.expected_name;
  $('carbohydrate-roundtrip').textContent = '期待糖類名の往復照合：' + (labels[data.carbohydrate_expected_roundtrip_status] || '未判定');
  $('carbohydrate-note').textContent = carb?.note || '';
  $('carbohydrate-units').replaceChildren(); $('carbohydrate-links').replaceChildren();
  for(const unit of carb?.units || []) {
    const title=document.createElement('p'); title.textContent='環 ' + unit.ring_id + ' · 糖立体骨格：' + (unit.label || '未判定') + ' / D/L ' + (unit.series || '未判定') + ' / α/β ' + (unit.anomer || '未指定・未判定');
    const mapping=document.createElement('p'); mapping.className='carbohydrate-map'; mapping.textContent='糖類番号 → 原子index：' + Object.entries(unit.carbon_atoms).map(([n,idx])=>'C' + n + ' → ' + idx).join(' / ');
    if(unit.amino_atoms?.['2']!==undefined) {
      mapping.textContent += ' / ' + (unit.amino_state==='N_acetyl' ? 'C2–NH–C(=O)CH3' : 'C2–NH2') + '：N原子 ' + unit.amino_atoms['2'];
      if(unit.amino_state==='N_acetyl') mapping.textContent += ' / アセチル基：C=OのC ' + unit.n_acetyl_atoms.carbonyl_carbon + '・O ' + unit.n_acetyl_atoms.carbonyl_oxygen + '・CH3のC ' + unit.n_acetyl_atoms.methyl_carbon;
      const names={glucosamine:'グルコサミン',galactosamine:'ガラクトサミン',mannosamine:'マンノサミン'};
      if(names[unit.family]) title.textContent += ' / ' + (unit.amino_state==='N_acetyl' ? 'N-アセチル' : '') + names[unit.family];
    }
    if(unit.deoxy_positions?.includes(6)) mapping.textContent += ' / C6 = CH3（6-デオキシ）';
    $('carbohydrate-units').append(title,mapping);
  }
  for(const bridge of carb?.o_linkages || []) {
    if(bridge.directions.length===2 && bridge.directions.every(d=>d.target_position===1)) {
      const d=bridge.directions[0],line=document.createElement('p');
      line.textContent='1↔1 O結合（1本） · 環 ' + d.source_ring_id + ' のC1 ↔ 環 ' + d.target_ring_id + ' のC1 / 橋Oの原子index ' + bridge.bridge_oxygen + '。両側のアノマー炭素が結合しています。';
      $('carbohydrate-links').append(line);
      continue;
    }
    for(const d of bridge.directions) {
      const source=carb.units.find(u=>u.ring_id===d.source_ring_id),target=carb.units.find(u=>u.ring_id===d.target_ring_id);
      const line=document.createElement('p'); line.textContent='糖間O結合 · 環 ' + d.source_ring_id + ' のC1 → 環 ' + d.target_ring_id + ' のC' + d.target_position + ' / ' + (source.label || '糖種未判定') + ' → ' + (target.label || '糖種未判定') + ' / 橋Oの原子index ' + bridge.bridge_oxygen;
      $('carbohydrate-links').append(line);
    }
  }
  for(const link of carb?.aglycon_o_linkages || []) {
    const line=document.createElement('p');line.textContent='アグリコンへのO結合 · 環 ' + link.source_ring_id + ' のC1 → O原子 ' + link.bridge_oxygen + ' → アグリコン原子 ' + link.aglycon_atom + ' / ' + (link.sugar_label || '糖種未判定');$('carbohydrate-links').append(line);
  }
  const glyco=data.glycoside_analysis;
  $('glycoside-result').hidden = !glyco?.rings?.length;
  $('glycoside-summary').textContent = glyco ? '糖様環候補 ' + glyco.rings.length + '個 / O結合候補 ' + glyco.o_linkages.length + '本 / C–C結合候補 ' + (glyco.c_linkages || []).length + '本 / 末端CH3の判別保留 ' + (glyco.ambiguous_c_linkages || []).length + '本（名称の検証とは独立）' : '';
  $('glycoside-note').textContent = glyco?.note || '';
  setImage('glycoside-structure',data.glycoside_svg);
  $('glycoside-links').replaceChildren();
  for(const ring of glyco?.rings || []) {
    const line=document.createElement('p'); line.textContent='環候補 ' + ring.id + ' · ' + ring.parent + ' / 原子index [' + ring.ring_atoms.join(', ') + '] / 環O ' + ring.ring_oxygen;
    const stereo=(ring.stereo_centers || []).map(c=>'原子 ' + c.atom_index + ': ' + (c.cip || '未指定')).join(' / ');
    if(stereo)line.textContent += ' / R/S ' + stereo;
    $('glycoside-links').append(line);
  }
  for(const link of glyco?.o_linkages || []) {
    const line=document.createElement('p');
    line.textContent='O結合候補 · 原子index [' + link.bond_atoms.join(' → ') + '] / ' + link.directions.map(d=>'環候補 ' + d.source_ring_id + ' から ' + (d.targets.length ? d.targets.map(t=>'環候補 ' + t.ring_id + (t.type==='exocyclic_methylene' ? ' の環外CH2' : ' の環炭素')).join('・') : '外部の炭素')).join(' ／ ');
    $('glycoside-links').append(line);
  }
  for(const link of glyco?.c_linkages || []) {
    const line=document.createElement('p');
    line.textContent='C–C結合候補 · 原子index [' + link.bond_atoms.join(' → ') + '] / ' + link.directions.map(d=>'環候補 ' + d.source_ring_id + ' の環O隣接炭素 ' + d.source_atom + ' から炭素 ' + d.target_atom + (d.source_oxygen_retained ? '（外側Oも保持）' : '')).join(' ／ ');
    $('glycoside-links').append(line);
  }
  for(const side of glyco?.ring_side_chains || []) {
    const line=document.createElement('p'); const kind={hydroxymethyl:'CH2OH',o_substituted_methylene:'CH2–O側鎖',methyl:'末端CH3'};
    line.textContent='環O隣接位置の側鎖 · 環候補 ' + side.ring_id + ' / ' + kind[side.kind] + ' / 原子index ' + side.ring_attachment_atom + ' → ' + side.carbon_atom + (side.oxygen_atom!==undefined ? ' → O ' + side.oxygen_atom : '');
    $('glycoside-links').append(line);
  }
  $('ambiguous-c-note').hidden = !glyco?.ambiguous_c_linkages?.length;
  $('ambiguous-c-note').textContent = '末端CH3のC–C結合は判別保留：' + (glyco?.ambiguous_c_linkages || []).map(link=>'原子index [' + link.bond_atoms.join(' → ') + ']').join('、') + '。デオキシ糖の末端かC-糖の置換基かは、この部分構造だけでは決めません。糖類解析欄で糖骨格との対応を確認できます。';
  $('c-linkage-note').hidden = !glyco?.c_linkages?.length;
  $('c-linkage-note').textContent = glyco?.c_linkage_note || '';
  const scope = data.rule_scope;
  $('scope-result').hidden = !(scope?.expected || scope?.candidates?.length || scope?.orientations?.length || scope?.selection_summary || scope?.scope_note);
  $('expected').textContent = scope?.expected || '未判定';
  $('expected-check').hidden = !scope?.expected;
  $('expected-check').textContent = '期待表記の往復照合：' + (labels[data.expected_roundtrip_status] || '未判定');
  $('scope-note').textContent = scope?.scope_note || '';
  $('scope-note').hidden = !scope?.scope_note;
  const family = {carbohydrate_retained:'限定した糖類名・O配糖体・二糖',oxygen_ring_stereo:'環状ポリオールの絶対R/S',oxygen_ring_polyol:'環状ポリオール・糖様骨格',functional_carbon_tree:'分岐・不飽和の官能基鎖',acyclic_ether:'エーテル・アルコキシ鎖',acyclic_ester:'非環式モノエステル',acyclic_amide:'非環式モノアミド',functional_benzene:'官能基を持つベンゼン環',mixed_functional_chain:'複合官能基の直鎖',registered_heterocycle:'登録した複素環',branched_alkane:'分岐アルカン',decorated_carbon_chain:'不飽和鎖・ハロゲン置換',halo_cycloalkane:'ハロゲン置換炭素単環',mono_unsaturated_chain:'単一不飽和鎖',functional_polyol:'ジオール・トリオール',functional_ketone:'ケトン',functional_acid:'カルボン酸',functional_aldehyde:'アルデヒド',functional_amine:'第一級アミン',cycloalkane:'炭素単環',benzene:'ベンゼン',alcohol:'単価アルコール',alkane:'直鎖アルカン'};
  $('family').textContent = 'RULE SCOPE / ' + (family[scope?.kind] || '限定検証');
  const candidates = scope?.candidates || scope?.orientations || [];
  $('trace-details').hidden = !(candidates.length || scope?.selection_summary);
  $('trace').replaceChildren();
  if(scope?.selection_order) { const order = document.createElement('p'); order.textContent = '選択順：' + scope.selection_order.join(' → '); $('trace').append(order); }
  if(scope?.selection_summary) { const s = scope.selection_summary; const line = document.createElement('p'); line.textContent = '全炭素の最長鎖：' + s.unrestricted_longest + '炭素 / 主特性基：' + s.principal_count + '個 / 1本の鎖に含められる最大数：' + s.max_in_one_path + '個 / 採用する親鎖：' + (s.selected_length ? s.selected_length + '炭素' : '未決定'); $('trace').append(line); }
  if(scope?.stereo_numbering) { const line=document.createElement('p'); line.textContent='R/Sの対応：' + scope.stereo_numbering.map(c=>'環位置 ' + c.locant + ' → 原子index ' + c.atom_index + ' → ' + c.cip).join(' / '); $('trace').append(line); }
  if(scope?.side_paths) { const sides = document.createElement('p'); sides.textContent = 'O側／N上のアルキル鎖の原子index：' + JSON.stringify(scope.side_paths); $('trace').append(sides); }
  for (const candidate of candidates) {
    const line = document.createElement('p');
    const preferred = scope?.preferred_candidates?.some(x => JSON.stringify(x.path) === JSON.stringify(candidate.path)) || scope?.orientations?.[0] === candidate;
    line.textContent = (preferred ? '採用候補 · ' : '比較候補 · ') + '骨格の原子index [' + candidate.path.join(', ') + '] / ' + (('groups' in candidate) ? '置換基 ' + JSON.stringify(candidate.groups) + (candidate.branch_count !== undefined ? ' / 個数 ' + candidate.branch_count : '') : '主特性基 ' + JSON.stringify(candidate.locants) + ' / OH ' + JSON.stringify(candidate.hydroxy));
    if (candidate.groups && candidate.locants) line.textContent += ' / 主特性基 ' + JSON.stringify(candidate.locants);
    if (candidate.doubles) line.textContent += ' / 二重結合 ' + JSON.stringify(candidate.doubles) + ' / 三重結合 ' + JSON.stringify(candidate.triples);
    if (candidate.groups === null) line.textContent += '（側鎖の命名は未対応。置換基数による候補比較のみ）';
    $('trace').append(line);
  }
  const hasDifference=data.checks.some(r=>r.status==='different');
  $('overall').className = `badge ${data.overall==='fail' ? 'fail' : hasDifference ? 'different' : 'unknown'}`;
  $('overall').textContent = data.overall === 'fail' ? '違反・不一致あり' : hasDifference ? '表記差あり・未判定' : '全体：未判定';
  setImage('input-structure', data.input_svg);
  setImage('restored-structure', data.restored_svg);
  $('no-structure').hidden = !!data.restored_svg;
  $('checks').replaceChildren();
  for(const row of data.checks) {
    const box = document.createElement('div'); box.className = 'check';
    const top = document.createElement('div'); top.className = 'check-top';
    const title = document.createElement('h4'); title.textContent = row.title;
    const badge = document.createElement('span'); badge.className = `badge ${row.status}`; badge.textContent = labels[row.status];
    top.append(title, badge); box.append(top);
    const detail = document.createElement('p'); detail.textContent = row.detail; box.append(detail);
    if(row.url) { const link = document.createElement('a'); link.href = row.url; link.textContent = row.reference; link.target = '_blank'; link.rel = 'noopener noreferrer'; box.append(link); }
    $('checks').append(box);
  }
  $('metadata').textContent = JSON.stringify({実行日時:data.created_at, 規則セット:data.ruleset, バージョン:data.versions, 入力の正規化:data.canonical_input, 再構築の正規化:data.canonical_restored || null, 比較方針:data.comparison_policy, 規則検証の根拠:data.rule_scope || null, 糖様環とO結合の候補:data.glycoside_analysis || null, 糖類番号と立体識別:data.carbohydrate_analysis || null, 期待糖類名の逆変換:data.carbohydrate_expected_opsin || null, フラボノイド骨格と番号:data.flavonoid_analysis || null, 単環炭素骨格:data.monocycle_analysis || null, 縮合二環骨格:data.fused_analysis || null, 架橋二環骨格:data.bridged_analysis || null, 期待縮合二環名の逆変換:data.fused_expected_opsin || null, 期待架橋二環名の逆変換:data.bridged_expected_opsin || null, 期待単環名の逆変換:data.monocycle_expected_opsin || null, 期待フラボノイド名の逆変換:data.flavonoid_expected_opsin || null, NISPO生成状態:data.generation || null, OPSIN:data.opsin || null, 期待表記の逆変換:data.expected_opsin || null}, null, 2);
  $('copy').disabled = !data.name;
}
$('form').addEventListener('submit', async event => {
  event.preventDefault(); $('error').textContent = ''; $('run').disabled = true;
  $('smiles').readOnly = $('candidate').readOnly = true;
  for(const button of document.querySelectorAll('[data-smiles]')) button.disabled = true;
  // Remove stale results while a different input is being checked.
  $('result').hidden = true; report = null;
  $('overall').textContent = '処理中'; $('overall').className = 'badge unknown';
  $('progress').textContent = 'ローカルで処理中… 命名と各名称の照合はそれぞれ最大25秒です。';
  try {
    const response = await fetch('/api/audit', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify({smiles:$('smiles').value, candidate:$('candidate').value})});
    const data = await response.json();
    if(!response.ok) throw new Error(data.error || '処理に失敗しました。');
    render(data); $('progress').textContent = '処理が完了しました。結果はブラウザ内にのみ保持しています。';
  } catch(error) { $('error').textContent = error.message; $('progress').textContent = ''; $('overall').textContent = '未完了'; }
  finally { $('run').disabled = false; $('smiles').readOnly = $('candidate').readOnly = false; for(const button of document.querySelectorAll('[data-smiles]')) button.disabled = false; }
});
function invalidateResult() {
  if (!report) return;
  report = null; $('result').hidden = true; $('overall').textContent = '入力変更・未実行'; $('overall').className = 'badge unknown'; $('progress').textContent = '入力が変わりました。再度検証してください。';
}
$('smiles').addEventListener('input', invalidateResult);
$('candidate').addEventListener('input', invalidateResult);
for(const button of document.querySelectorAll('[data-smiles]')) button.addEventListener('click', () => { invalidateResult(); $('smiles').value = button.dataset.smiles; $('candidate').value = button.dataset.candidate || ''; if(button.dataset.candidate) $('candidate').closest('details').open = true; $('smiles').focus(); });
$('copy').addEventListener('click', async () => {
  try { await navigator.clipboard.writeText(report.name); $('copy').textContent = 'コピー済み'; setTimeout(() => $('copy').textContent = 'コピー', 1600); }
  catch { $('error').textContent = 'コピーできませんでした。表示名を選択してコピーしてください。'; }
});
$('export').addEventListener('click', () => {
  if(!report) return;
  const {input_svg, restored_svg, ...data} = report;
  const url = URL.createObjectURL(new Blob([JSON.stringify(data, null, 2)], {type:'application/json'}));
  const anchor = document.createElement('a'); anchor.href = url; anchor.download = `smiles-audit-${Date.now()}.json`; anchor.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});

fetch('/coverage.json').then(r => {if(!r.ok) throw new Error('coverage'); return r.json();}).then(data => {
  $('coverage').replaceChildren();
  for(const item of data.families) {
    const box = document.createElement('div'); box.className = 'coverage-item';
    const title = document.createElement('h4'); title.textContent = item.name + ' · ' + item.status;
    const detail = document.createElement('p'); detail.textContent = item.detail;
    box.append(title, detail); $('coverage').append(box);
  }
  const limit = document.createElement('p'); limit.textContent = data.global_limit; $('coverage').append(limit);
}).catch(() => { $('coverage').textContent = '対応範囲を取得できません。READMEの対応表を確認してください。'; });
