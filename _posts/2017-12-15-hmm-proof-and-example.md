---
layout: article
date: 2017-12-15 17:19
title: "HMMに使われているアルゴリズムの詳細とサンプル"
tags: [machine-learning, hmm]
comments: true
published: true
---

HMMに使われているアルゴリズム (Baum-Welch, Forward-Backward, Viterbi) の理解を深めるために証明を行ったので，ここにリンクを載せておく．~~ただし一部証明ができない部分があった．~~

(~~Baum-Welchの$$\Sigma_k$$の期待値を最大化するところの行列の微分が分からなかった．~~最終的な答えがわかっているので辿り着けているところが大いにあるので，式変形がどっかで間違っている可能性も十分ある．)

18/01/16追記．$$\Sigma_k$$の期待値を最大化する式の導出ができた．PRMLのAppendixに行列の微分についての公式があり，それを参考にすることでできた．


[hmmアルゴリズムメモ](https://github.com/fhiyo/hmm-algo-memo/blob/master/main.pdf)

また，アルゴリズムの実装ではないけど使い方として練習したものはこちら．  
初期確率と遷移確率，出力確率を与えたHMMのモデルを構築し，観測変数と潜在変数をサンプリングして出力している．  
また，観測変数を適当に与えてやり，その潜在変数を推定する操作を後半で行っている．  
[Hidden Markov models sample](https://gist.github.com/fhiyo/e0373db6f4f8c92d4fddd7467117f8de)

やっつけだし[hmmlearn](http://hmmlearn.readthedocs.io/en/latest/auto_examples/plot_hmm_sampling.html#sphx-glr-auto-examples-plot-hmm-sampling-py) の例とほとんど同じだけどとりあえず載せておこう．
