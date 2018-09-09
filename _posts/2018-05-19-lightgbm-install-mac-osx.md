---
layout: article
date: 2018-05-19 18:28
title: "MacにLightGBMを再インストール"
tags: [python, mac, LightGBM]
comments: true
published: true
---

LightGBMを使おうと思いpyhtonのnotebook上でimportしようとしたら、エラーが出た。

```python
from LightGBM import LGBMClassifier
```

Output (Excerpt):

```
OSError: dlopen(/usr/local/lib/python3.6/site-packages/LightGBM/lib_LightGBM.so, 6): Library not loaded: /usr/local/opt/gcc/lib/gcc/7/libgomp.1.dylib
  Referenced from: /usr/local/lib/python3.6/site-packages/LightGBM/lib_LightGBM.so
  Reason: image not found
```

この前までは正常に使えていたと思ったのだが。エラーを見るとgccの7系のライブラリを見に行っているようだが、自分の環境には7系はなかった。じゃあなんで今まで動いていたんだろう？8にアップデートしたのが最近だったんだろうか？

と思ってbrewのlogを見てみたら、5/7にgccのバージョンを7.3.0_1 -> 8.1.0に上げていた。なるほど。

エラーの原因もわかったので、動くようにする。一旦LightGBMをアンインストールして再度インストールすることにした。

[LightGBMのインストールガイド](https://github.com/Microsoft/LightGBM/blob/master/docs/Installation-Guide.rst)を見ると、LightGBMはコンパイルにOpenMPを使っており、これはApple Clangじゃサポートしてないからgcc(の7系)を使ってコンパイルしてくれと言っている。gccは8系なので、7系をインストールする。

```sh
$ brew install gcc@7
```

環境:

```sh
$ sw_vers
ProductName:  Mac OS X
ProductVersion: 10.13.4
BuildVersion: 17E202
```

インストールガイドに沿って、cloneしてビルド、インストールする。

```sh
$ cd /tmp
$ git clone --recursive https://github.com/Microsoft/LightGBM
$ cd LightGBM
$ export CXX=g++-7 CC=gcc-7
$ mkdir build
$ cd build
$ cmake ..
$ make -j4
```

LightGBMのインストールはできたが、元々のエラーは`/usr/local/opt/gcc/lib/gcc/7/libgomp.1.dylib`がないということだったので、パスを確認してみる。

```sh
$ ls /usr/local/opt
lrwxr-xr-x 19 fhiyo  7 May 12:00 gcc -> ../Cellar/gcc/8.1.0
lrwxr-xr-x 23 fhiyo  4 Apr 15:57 gcc@5 -> ../Cellar/gcc@5/5.5.0_2
lrwxr-xr-x 21 fhiyo 19 May 18:08 gcc@7 -> ../Cellar/gcc@7/7.3.0
```

とりあえずこれ書き換えるしかないっぽいのでそうする。

```sh
$ mv gcc gcc@8
$ ln -s gcc@7 gcc
```

これで動くようにはなったけど嫌だなぁ

#### ちなみに

Issueが4日前に上がっていた。[LightGBM and gcc 8 in MacOS: `Library not loaded: /usr/local/opt/gcc/lib/gcc/7/libgomp.1.dylib` · Issue #1369 · Microsoft/LightGBM](https://github.com/Microsoft/LightGBM/issues/1369)

PyPIの現時点の最新である2.1.1はgcc-7になっているようだが、8を使うように変更したプルリクは既に上がっているとのこと。問題が解決するまでそんなには時間かからなさそう？


#### 追記 (2018/05/21)

```sh
$ R
dyld: Library not loaded: /usr/local/opt/gcc/lib/gcc/8/libgfortran.5.dylib
  Referenced from: /usr/local/Cellar/r/3.5.0_1/lib/libR.dylib
  Reason: image not found
zsh: abort      R
```

やっぱり他に影響があった。
