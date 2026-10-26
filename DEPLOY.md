# Streamlit Community Cloudへの公開

Web版 v0.31.0。命名・独立規則エンジンはv0.30.0のままです。新しい骨格の追加は止め、公開に必要なUIと実行上限を用意しました。

## 公開する

1. 公開用ZIPを展開し、中のファイルをGitHubリポジトリのルートに置きます。`.venv`や過去の検証データは不要です。`vendor/opsin-cli-2.9.0.jar`と`.streamlit/config.toml`も含めてください。
2. https://share.streamlit.io/ にGitHubでサインインし、Create appを選びます。
3. リポジトリ・ブランチを選び、Entrypointを **streamlit_app.py** に設定します。
4. Advanced settingsで **Python 3.13** を選び、Deployします。Secretsの設定は不要です。
5. 起動後、エタノール・フラボン・カンファー型の例を試し、JSONレポートを保存できればβ公開の確認は完了です。

Javaは`packages.txt`の`default-jdk-headless`から導入します。依存関係はルートのrequirements.txtから読み込みます。無料環境の速度・メモリは実際のデプロイ後に確認してください。ここまでの動作確認はローカルmacOS/Python 3.13です。

## ローカルでWeb版を動かす

Python 3.13とJava 11以降のJDKを用意して実行します。

```bash
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m streamlit run streamlit_app.py
```

元のローカル専用HTML版も残しています。

```bash
.venv/bin/python server.py
```

## β版の動作

- 全利用者を通して計算は同時1件。混雑時は待って再実行します。
- 1件120秒で計算プロセスと子プロセスを停止。SMILESは3,000文字・200原子まで、候補名は2,000文字まで。
- 入力と結果をアプリのファイル・DB・履歴・共有キャッシュへ保存しません。結果はStreamlitの利用者ごとの一時セッションに保持します。ブラウザから入力を変更すると古い結果を非表示にします。
- Web版の入力はホスティング先へ送信されます。ローカル処理とは異なります。ホスティング事業者のログ・運用はその事業者の規定に従います。
- NISPOの生成名や入力名を期待名で上書きしません。構造往復一致と限定規則、PIN全体の未判定を区別します。
- 機密構造には元のローカル版を使ってください。

## 不具合への対応

- Javaがない：packages.txtがリポジトリのルートにあるか確認し、Cloud上でReboot。
- 依存関係が解決しない：Advanced settingsのPythonを3.13へ合わせます。
- 処理が遅い：まず単純構造／既存の名前の検査で確認。継続的な混雑やメモリ不足が観測されたら、同じエンジンをDocker対応ホスティングへ移します。

公開を優先し、通常は変更箇所の試験＋代表例の実行を使います。命名規則やエンジンに変更があるときだけ、対応する回帰試験と公開前の全体試験を行います。

公式資料：
- https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/app-dependencies
- https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy
- https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app
