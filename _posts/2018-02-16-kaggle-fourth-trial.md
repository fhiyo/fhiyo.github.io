---
layout: article
date: 2018-02-16 18:24
title: "KaggleのTitanic問題を試してみた4"
tags: [machine-learning, Kaggle]
comments: true
published: true
---

[前回](/blog/2018/01/kaggle-third-trial/)の続編．feature engineeringした．スコアは0.78になった(前回は0.72)．

結果を載せたjupyter notebookは[こちら](https://gist.github.com/fhiyo/618c58615e664442075d396e32dcc39b)．

前回までは"Ticket", "SibSp", "Parch", "Fare", "Cabin", "Embarked"の属性は扱うのが面倒だったので削除していたが，今回はそれぞれの属性と"Survived"との相関を可視化しつつ処理の方法を決めていった．また，データの分布を見て偏りが大きい場合はlogを取って正規分布に近づける処理も行った．結果として前回よりも0.06精度が向上しているので，やったかいはあったといえる．

Titanicは飽きてきたので，今後は別の問題にチャレンジするか機械学習のアルゴリズムの理解を行っていくことにする．
