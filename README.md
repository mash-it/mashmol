#mashmol

Coarse-Grained Go-model protein MD simluation based on OpenMM written in C++.

## Environment and Dependency

* C++11
* Python 3.x and numpy
* [OpenMM](http://openmm.org/)
* [json.hpp](https://github.com/nlohmann/json)

## Installation


```
$ make
```

and pray

## Usage 


```bash
$ python pdb2force.py foo.pdb > input.json
$ ./mashmol input.json
```

and see `output.pdb` and `output.ts`

入力ファイル (JSON) を読ませる形式。入力ファイルは PDB から pdb2force.py で自動生成できる。

結果は output.pdb および output.ts として生成される。

あとは `pdb2force.py` の中を見て、適当にパラメータとかいじる。

## issue

* 計算が遅すぎる
    * 並列化するとかえって遅い
    * 何か致命的な部分があるのか、それとも単に OpenMM が数百粒子というスケールに合わないのか



