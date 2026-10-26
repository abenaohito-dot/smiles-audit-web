# SMILES → IUPAC — Streamlit β版

SMILESからNISPOで命名し、OPSINで構造へ戻して照合します。対応範囲の独立した命名規則も個別に検査します。構造一致はIUPAC規則全体への適合やPINを保証しません。

**最短で公開：** このフォルダの内容をGitHubリポジトリへ置き、[Streamlit Community Cloud](https://share.streamlit.io/)から`streamlit_app.py`を指定。Advanced settingsでPython **3.13** を選んでDeployします。APIキー・データベース・Dockerは不要です。

詳しい手順と動作条件は[DEPLOY.md](DEPLOY.md)。Python依存はrequirements.txt、Javaはpackages.txtで宣言しています。同梱OPSIN JARも含めてアップロードしてください。

```bash
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m streamlit run streamlit_app.py
```

- Web UI：v0.31.0。命名エンジン：検証済みv0.30.0から変更なし。
- 同時1件、最大120秒、200原子。入力はホスティング先へ送信されます。
- アプリによる入力のファイル保存・履歴・共有キャッシュなし。結果は利用者の一時セッションだけで保持。
- フラボノイド、糖、単環・デカリン・ノルボルナン・ピナン型などの有限範囲。未対応は未判定。
- ローカル専用HTML版は`python server.py`（localhostのみ）で利用可能。

NISPO・RDKitはBSD-3-Clause、OPSINはMIT。OPSINライセンスはvendor/OPSIN-LICENSE.txt、同梱依存のライセンスはJAR内、Python依存のライセンスは各インストールパッケージに含まれます。
