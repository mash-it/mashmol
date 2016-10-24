#mashmol document

## Residue class
list の派生クラスで Atom object のリスト。

`getCa(self)`
その残基におけるα炭素の atom object を取得する。

## GoProtein class

`class GoProtein(filename)`
PDB file を開いて Go モデルの Protein として解釈する。

### コンストラクタで呼び出されるメソッド

`readAtoms()`
PDB file に含まれる `ATOM` 行をすべて読んで `self.atoms` プロパティに入れる。

`classifyResidues()`
各原子を残基ごとに分類して `self.residues` に格納する。`self.residues` は `int resSeq -> Residue` の dict である。

`readAtomLine(line)`
string を Atom object に変換して返す。Atom object はそういうクラスではなく単なる dict であり、`serial`, `atomname`, `altLoc` といった PDB の基本的な要素のほか、座標を `numpy.array` 化した `atom['pos']` が含まれる。

### その他のメソッド

`getBoxSize()`
原子が存在する範囲を返す。

`getNativeBondLengt(resSeq a, resSeq b)`
2つの `resSeq` を指定すると、対応する2残基間の距離を返す。

`getNativeAngle(resSeq a, resSeq b, resSeq c)`
3つの `resSeq` を指定すると、対応する3残基のなす角を返す。値域は 0 < x < pi.


`getNativeDihedral(resSeq a, resSeq b, resSeq c, resSeq d)`
4つの `resSeq` を指定すると、対応する4残基の二面角を返す。値域は -pi < x < pi で、右ねじ方向が正。

`shift(distance, direction)`
構造全体を平行移動させる。

`getPdbText()`
構造全体を PDB ファイルとして出力する。
