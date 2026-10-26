# v0.25の根拠と適用範囲

- IUPAC/IUBMB Flavonoid Nomenclature: https://iupac.qmul.ac.uk/flavonoid/
- Flv-3.2.2.2 Example 28のgenistin全体系統名をリテラル参照として試験。Flv-3.1.1.3のflavan-4-one/chroman-4-oneと環立体を使用。糖置換は既存の4つのglucose D/L×α/β局所参照へキラル照合。
- C配糖体とO配糖体の区別はFlv-3.2.1.2–3、フラボノール3-OはFlv-3.2.1.4。今回の新骨格はA環のOのみ。糖カタログは増やさない。
- 単環炭素親骨格はIUPAC Blue Book P-22.1: https://iupac.qmul.ac.uk/BlueBook/P2.html
- 番号の順位はP-14.4: https://iupac.qmul.ac.uk/BlueBook/P1.html
- 環上OH・ケトンのol/one接尾辞はP-63.1 / P-64.1: https://iupac.qmul.ac.uk/BlueBook/P6.html
- 絶対配置はRDKit rdCIPLabelerへ委譲。糖種・D/L・α/βをR/S文字だけから推定しない。糖を結合した元分子と切り出し分子のCIPは同じと仮定しない。

新しい単環検証は5/6員炭素単環・環内二重結合0/1本・環上OHまたはケトン0/1個・有限側鎖3個まで（分岐C3側鎖は各1個）で全原子を説明する範囲。番号・主特性基・置換位置・R/Sと表記差を分ける。慣用名の優先性、保持名、同順位の複雑な立体番号選択、相対立体、PIN全体は証明しない。全体期待名はOPSINで元の全分子へ戻して照合する。NISPOやOPSINを独立規則の正解ラベル生成には使わない。


## v0.26

P-44.1.1（主特性基の接尾辞に対応する基を持つ親）とP-44.1.2.2（同じクラスの環・鎖の選択）：https://iupac.qmul.ac.uk/BlueBook/P4.html 。今回の限定方針は環外OHのみなら鎖親、OHが同順位なら環親、環上ケトンがあればケトンを接尾辞にする環親。環外OHは直結CH2OHまたはC(OH)(CH3)2の1個に絞る。親選択方針の適合と構造一致は別であり、同じ構造を表す別名があり得る。全命名法やPIN全体は判定しない。

環上OH/ケトン0/1個との共存を含め、元の未切断分子でCIPを計算。環が置換基なら局所的な環番号のR/Sをその環名内に記載。一般の環外OH鎖の優先順位・長鎖・環外カルボニル、相対立体、複数環外OHや多環へは外挿しない。
