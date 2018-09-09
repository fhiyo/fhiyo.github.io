---
layout: article
date: 2018-04-16 02:06
title: "Juliaを使ってKaggleのhouse pricesの問題を試してみた"
tags: [machine-learning, Kaggle]
comments: true
published: true
---

今回は[Julia](https://julialang.org/)という比較的新しめの言語でKaggleをやってみることにした．
まだKaggleも大してやれていないのになぜ今まで触ったこと無い言語を試してみたかというと，社内のハッカソンでやることになったから．

Juliaは実行時コンパイルすることでPythonより高速に計算できることを売りにしている言語のようだ (だが今回の簡単なKernelではその良さは享受できず．．．)．どんな言語かはこちらとかを見ればわかると思う→[Python使いをJuliaに引き込むサンプル集](http://www.mwsoft.jp/programming/julia/python_to_julia.html)

題材にしたKaggleのコンペは[House Prices: Advanced Regression Techniques](https://www.kaggle.com/c/house-prices-advanced-regression-techniques)で，アイオワ州のエイムズにある家についての様々な情報 (特徴量の数が79個) から家の値段を予測する，というもの．要するに家の値段について回帰をしろ，という問題．

試してみた結果は[こちらのGist](https://gist.github.com/fhiyo/d514d0f83efe22bdd262a8593d9b6efa)にKernelを公開している．本当はKaggle上で公開したかったのだが，現時点ではJuliaのKernel公開用の環境が整っていないみたい (昔のversionのKernelは投稿できたようだ)．

Juliaのライブラリはまだあまり成熟していないのか，versionの依存関係によって動かなくなる関数があったりwarningが出たりすることが多く (dependency hellというやつ？)，中々苦戦した． (まだ見れていないが，[Pkg3](https://github.com/JuliaLang/Pkg3.jl)というパッケージ管理システムを使えば問題が解決されるのかも？)
またJuliaは動的型付けながらも，コンパイルされる関係からか型にうるさい印象だった．
Pythonではよしなにやってくれていた処理が全然通らず，NAという欠損値を表す特別な値が型チェックですぐに引っかかったり，GridSearchCVがなぜか全然動かなかったりでハッカソンの一日では前処理をまともに行うことができずじまいだった．

言語自体にはLispのようにマクロが使えたり，関数型言語としての側面もあるので慣れれば楽しそうではあるものの，データ分析でぱっと使うにはライブラリが整っていないためまだまだ厳しいという印象．ネット上のサンプルが0.4とか昔のversionのものが多くてコピペで動作確認ができないし，helpで呼び出せるリファレンスもまだ充実していない感じなので勉強する環境が十分に整っていないのでは，と思う．まあ一日ちょっと触ってみただけなので，もう少し触っていればやりやすくなるのかもしれないが．

KaggleのKernelを見てやり方を勉強する方法はとてもよいと思うので，ぜひJuliaのKernelも早く公開できるようになってほしい．今の状態でコード書くのは自分で試行錯誤しなければならない部分が多くて無駄に時間を食ってしまった．
