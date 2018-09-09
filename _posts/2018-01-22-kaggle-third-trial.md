---
layout: article
date: 2018-01-22 01:28
title: "KaggleのTitanic問題を試してみた3"
tags: [machine-learning, Kaggle]
comments: true
published: true
---

[前回](/blog/2018/01/kaggle-second-trial/)の続編．色んな機械学習のアルゴリズムを試して，良さそうなやつ数個をVotingClassifierでアンサンブル学習して精度を出している．  
結論から言うと今回のスコアは0.72で，前回の0.71よりもちょっとだけ良くなってる (誤差の範囲か？)

やったことはいつもの通り[gistのjupyter notebook](https://gist.github.com/fhiyo/58efcfd245c238e00d0cd6eace911893)に書いた．

アルゴリズムを色々なものを試してみたものの，データの前処理の方は特に触っていないのでそんなに精度が上がることは期待していなかったが，ほとんど変わらないとはちょっと意外だった．やはりまだほとんどのデータの属性を捨てているのでアルゴリズムをいくら頑張っても効果は薄いんだろうな．

使ってみたアルゴリズムは以下の通り．

- SVC
- Decision Tree
- AdaBoost
- Random Forest
- Extra Trees
- Gradient Boosting
- Multiple layer perceprton (neural network)
- KNN
- Logistic regression
- Linear Discriminant Analysis
- xgboost

最終的にVotingClassifierに使用したアルゴリズムは，

- Random Forest
- Adaboost
- ExtraTrees
- Gradient boosting
- SVC
- (xgboost: これをVotingClassifierに入れるとモデルの学習に時間がかかって結果が出なかったので，今回は除外した．)

である．アルゴリズムの選定はパラメータチューニングしていない状態でcross validation scoreを算出して，比較して高いものを選んでいる (と思う．この辺りは[参考にしたKernel](https://www.kaggle.com/yassineghouzam/titanic-top-4-with-ensemble-modeling)のやり方をそのままにしてるので，勘違いかも．)


次は捨てていた属性を学習に使えるように前処理をもう少し頑張ってみよう．次で精度が向上することを期待．
