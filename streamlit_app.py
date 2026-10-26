"""Small research beta, sharing the existing naming engine without changing its rules."""
import base64
import json
from pathlib import Path
import shutil

import streamlit as st
from cloud_runner import AuditError, run_audit, CALCULATION_GATE

ROOT = Path(__file__).resolve().parent
WEB_VERSION = '0.31.0'
LABELS = {'pass': '一致・適合', 'fail': '不一致・違反', 'different': '表記差', 'unknown': '未判定', 'na': '該当なし'}
EXAMPLES = {
    '自由入力': '',
    'エタノール': 'CCO',
    'フラボン': 'O=c1cc(-c2ccccc2)oc2ccccc12',
    'アピゲニン型': 'O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12',
    'メントール型': 'CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O',
    'カンファー型（立体未指定）': 'CC1(C2CCC1(C(=O)C2)C)C',
    'α-ピネン型（立体未指定）': 'CC1=CCC2CC1C2(C)C',
}

def use_example():
    choice = st.session_state.example
    if EXAMPLES[choice]:
        st.session_state.smiles = EXAMPLES[choice]
        st.session_state.candidate = ''
    st.session_state.pop('report', None)

def invalidate():
    st.session_state.pop('report', None)

st.set_page_config(page_title='SMILES → IUPAC', page_icon='🧪', layout='wide')
st.title('SMILES → IUPAC')
st.caption(f'研究用β版 · Web {WEB_VERSION} · 構造の往復照合と、限定した命名規則を別々に確認')
st.info('入力はこのアプリのサーバーへ送信して計算します。アプリは入力をファイル・履歴・共有キャッシュへ保存しません。機密の構造にはローカル版をお使いください。')
with st.sidebar:
    st.subheader('判定の読み方')
    st.write('構造一致：名前から元の構造に戻るか。')
    st.write('限定規則：対応している骨格・番号・官能基・立体についての個別検査。')
    st.warning('構造が一致しても、命名規則への全面適合や優先IUPAC名（PIN）は保証しません。')
    st.caption('同時に1件を計算します。混雑時は少し待って再実行してください。1件の上限は120秒です。')
    with st.expander('対応範囲'):
        coverage = json.loads((ROOT / 'static' / 'coverage.json').read_text())
        for family in coverage['families']:
            st.markdown(f"**{family['name']}** — {family['status']}")
            st.caption(family['detail'])

st.selectbox('例を選ぶ', list(EXAMPLES), key='example', on_change=use_example)
st.session_state.setdefault('smiles', 'CCO')
st.text_area('SMILES', key='smiles', max_chars=3000, height=100, on_change=invalidate, help='200原子まで。指定されていない立体は補完しません。')
with st.expander('既存の名前を検査する（任意）'):
    st.text_input('検査する英語名', key='candidate', max_chars=2000, on_change=invalidate, help='空欄ならNISPOで生成します。')

if st.button('名前を生成・検証', type='primary', width='stretch'):
    st.session_state.pop('report', None)
    if not shutil.which('java'):
        st.error('Javaが見つかりません。公開設定のpackages.txtを確認してください。')
    else:
        gate = CALCULATION_GATE
        if not gate.acquire(blocking=False):
            st.warning('別の検証を実行中です。少し待ってから再実行してください。')
        else:
            try:
                with st.spinner('命名・構造復元・限定規則を確認しています…'):
                    result = run_audit(st.session_state.smiles, st.session_state.candidate)
                    result['web_interface_version'] = WEB_VERSION
                    st.session_state.report = result
            except AuditError as exc:
                st.error(str(exc))
            finally:
                gate.release()

report = st.session_state.get('report')
if report:
    st.divider()
    source = {'manual': '入力された名前', 'bounded_rules': '限定規則による期待名（構造復元一致）'}.get(report.get('name_source'), 'NISPOの生成名')
    st.subheader(source)
    st.code(report.get('name') or '命名できませんでした', language=None)
    st.caption(f"{report['formula']} · {report['canonical_input']}")
    cols = st.columns(3)
    checks = {row['id']: row['status'] for row in report['checks']}
    cols[0].metric('構造の往復照合', LABELS.get(checks.get('roundtrip'), '未判定'))
    cols[1].metric('限定規則の違反', str(sum(row['status'] == 'fail' for row in report['checks'] if row['id'] not in ('constitution', 'roundtrip', 'opsin')))+'件')
    cols[2].metric('PIN全体', '未判定')
    for kind, label in [('carbohydrate', '糖類'), ('flavonoid', 'フラボノイド'), ('monocycle', '単環'), ('fused', '縮合二環'), ('bridged', '架橋二環')]:
        analysis = report.get(kind + '_analysis') or {}
        if analysis.get('expected_name'):
            st.markdown(f"**{label}の限定期待名**")
            st.code(analysis['expected_name'], language=None)
            st.caption('期待名の構造復元：'+LABELS.get(report.get(kind+'_expected_roundtrip_status'), '未判定'))
    cols = st.columns(2)
    for col, key, caption in [(cols[0], 'input_svg', '入力構造（原子index）'), (cols[1], 'restored_svg', '名前から復元した構造')]:
        with col:
            st.markdown('**'+caption+'**')
            if report.get(key):
                # Only a base64-encoded SVG produced by RDKit enters the img source.
                encoded = base64.b64encode(report[key].encode('utf-8')).decode('ascii')
                st.markdown(f'<img src="data:image/svg+xml;base64,{encoded}" alt="{caption}" style="width:100%">', unsafe_allow_html=True)
            else:
                st.caption('構造を復元できませんでした。')
    st.subheader('検証の内訳')
    for row in report['checks']:
        with st.expander(f"{LABELS.get(row['status'], row['status'])} · {row['title']}"):
            st.write(row['detail'])
            if row.get('url'):
                st.link_button(row.get('reference') or '根拠', row['url'])
    export = {k: v for k, v in report.items() if not k.endswith('_svg')}
    st.download_button('検証レポートを保存（JSON）', json.dumps(export, ensure_ascii=False, indent=2), file_name='smiles-audit.json', mime='application/json', width='stretch')
    with st.expander('再現性・バージョン情報'):
        st.json({'engine_version': report['app_version'], 'web_version': WEB_VERSION, 'ruleset': report['ruleset'], 'versions': report['versions']})

st.caption('NISPO × OPSIN × RDKit · 互変異性体・塩・未指定立体を無条件に同一化しません。')
