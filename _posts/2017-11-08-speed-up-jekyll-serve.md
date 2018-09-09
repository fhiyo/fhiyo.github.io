---
layout: article
date: 2017-11-08 16:39
title: "jekyll serveを速くする"
tags: [ruby, jekyll]
comments: true
published: true
---

このブログのジェネレーターはJekyllを使っているが，書いた記事をプレビューするために`jekyll serve`コマンドを使うことでローカル上にサーバーが構築され，`localhost:4000`でアクセスすることができる．  
しかし，記事の数が増えるに従って`jekyll serve`でサーバーが立ち上がるまでにかかる時間がどんどん遅くなり，1分くらいまで遅くなってしまったので調べて何とかした．

[ここ](http://okzk.org/blog/how-to-speed-up-jekyll-lsi/)を参考にした．方法は以下．

1. `brew install gsl`で`gsl`をインストール
1. `Gemfile`に以下の記述を追加して`bundle install`でJekyllで`gsl`を使うための準備
```ruby
# For gsl
gem "sass"
gem "bourbon"
gem "neat"
gem "rb-gsl"
```

そもそも遅かったのは原因は，関連記事を表示させるのに`lsi: true`をしていたためっぽい．

この処理をしたところ，約1分かかっていたサーバーの立ち上げが3秒程度まで高速化された．

## gslというのは何か？
[GNU Scientific Library - Wikipedia](https://ja.wikipedia.org/wiki/GNU_Scientific_Library)より

>
GNU Scientific Library (GSL) は、C言語で記述された科学技術計算関数のライブラリである。オープンソースであり、GNU General Public Licenseのもとで配布されている。 このプロジェクトは1996年にロスアラモス国立研究所のDr. M. GalassiとDr. J. Theilerの着想に始まり、計算物理の専門家集団（Dr G. Jungman、Dr B. Gough、Dr J. Davies、R. Priedhorsky、Dr M. Booth、Dr F. Rossi、Dr D. Eddelbuettelら）を中心に作成された。

数値計算を高速に行えるライブラリのようだ．これで記事同士の距離みたいなのを計算して，関連度が高い記事を表示させているのかな？
