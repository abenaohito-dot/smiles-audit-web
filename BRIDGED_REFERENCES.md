# v0.29 架橋二環の根拠と限界

- IUPAC Blue Book P-23.2.2 / P-23.2.3: https://iupac.qmul.ac.uk/BlueBook/P2.html
- IUPAC Blue Book P-63 (suffixes): https://iupac.qmul.ac.uk/BlueBook/P6.html
- IUPAC Blue Book P-14.4 (locant selection): https://iupac.qmul.ac.uk/BlueBook/P1.html
- PubChem (-)-camphor record and displayed IUPAC name: https://pubchem.ncbi.nlm.nih.gov/compound/l-Camphor

対象グラフは bicyclo[2.2.1]heptane の7炭素骨格のみ。橋頭1と4の間の二つの2原子経路、短い1原子橋を登録SMARTSと全骨格原子の対応で確認する。分子のSMILES原子順から番号は導かない。対応する番号候補の主特性基 locant、続いて置換基 locant を比較し、候補名の構造対応と採用番号は別の判定にする。

全炭素骨格にメチル0〜3個と、環上OHまたはケトン0/1個だけを許可する。OHは ol、ケトンは one の接尾辞を使う。RDKitのCIPラベルが全て大文字R/Sとして確定し、指定が全てそろう場合、または立体指定が全て無い場合だけ全体期待名を作る。全未指定はラセミ体とみなさない。部分指定、r/s、相対配置、C=C、複数官能基、酸化窒素、塩、他の側鎖・架橋・縮合・スピロは未対応。

各期待名はOPSINで独立に全分子へ戻し、立体を含めて元入力と照合する。NISPOの非空の生成名は置き換えない。実測したカンファー出力 `(1S,4S)-4,7,7-trimethyl-3-oxo-bicyclo[2.2.1]heptane` は、その逆変換が一致することと、限定した主特性基の接尾辞・位置番号規則に一致しないことを別表示する。これは有限範囲の規則比較であり、慣用名の有効性全般やPIN全体を判定するものではない。

## v0.30 追加範囲と一次資料

- IUPAC Blue Book P-23.2.2 / P-23.2.3：親骨格と橋の番号。https://iupac.qmul.ac.uk/BlueBook/P2.html#P-23.2.2
- P-14.4 / P-15.1.7.2.3：固定番号の下で主特性基、次いで不飽和、次いで接頭辞へ低い番号。https://iupac.qmul.ac.uk/BlueBook/P1.html#P-14.4
- P-31.1.4.1：von Baeyer二環の二重結合番号。https://iupac.qmul.ac.uk/BlueBook/P3.html#P-31.1.4

資料確認日：2026-10-26。ピナン型は橋頭1と5を、1–2–3–4–5、5–6–1、5–7–1で接続した7炭素・8結合の登録グラフ。独立試験の構築子はこの辺リストを直接使い、解析SMARTSもOPSINの出力も参照構造の正解にしない。4つのグラフ対応に限って位置番号を比較する。独立性は有限グラフ規則についてのものであり、CIPと逆変換はそれぞれRDKit、OPSINへ依存する。

最長橋の内部結合2–3／3–4だけにC=Cを1本許可する。橋頭・短い橋・環外二重結合を受理せず、compound locantやBredt則の一般判定は実装しない。OH/ケトン0/1個、メチル0〜3個、全指定R/Sまたは全未指定に限定する。名称中の橋長と炭素数も照合し、異なる親骨格を合格にしない。既存ノルボルナンの範囲は飽和のまま。
