---
layout: article
date: 2017-11-05 21:45
title: "haskellの中置記法について"
tags: [haskell, programming]
comments: true
published: true
---

## 疑問

Haskellのコードを見ていて，このように書いてあるものを見つけたのだが，これがなぜ上手く動作するのかがわからなかった．

```haskell
> (`elem` [1,2,3]) 3
True
```

```3 `elem` [1,2,3]```や```elem 3 [1,2,3]```の形ならわかるし，```(3 `elem`) [1,2,3]```でも部分適応をしてるのはわかる．しかし，問題のコードがこれでいける理由が謎だった．

結果から見ると，どうも```(`elem` [1,2,3])```は```(\ x -> x `elem` [1,2,3])```と等価なようだが．．．？

## 答え

Stackoverflowに同じような疑問を聞いてる人がいて，その回答を見てわかった．

リンク: [functional programming - Why do both map (^2) xs and map (2^) xs work as expected in Haskell? - Stack Overflow](https://stackoverflow.com/questions/4703576/why-do-both-map-2-xs-and-map-2-xs-work-as-expected-in-haskell/4703600)

>
The Haskell grammar has special support for construct like this, called "operator sections". If you have any infix operator, like say #$%, then the following notation is supported:
```haskell
(#$%)   = \x y -> x #$% y
(#$% y) = \x   -> x #$% y
(x #$%) = \y   -> x #$% y
```
So you are expecting some mathematical consistency to break this, and if Haskell were a miniscule language like Forth, I would be inclined to agree with your intuition. The reason it works is basically "because they wrote it to work like that".

この答えを見る限り，やはり ```(`elem` [1,2,3])```は```(\ x -> x `elem` [1,2,3])```と等価であるらしい．[wiki.haskell.orgのページ](https://wiki.haskell.org/Section_of_an_infix_operator)を見てもやはり同じことが書いてあった．

>
In Haskell there is a special syntax for partial application on infix operators. Essentially, you only give one of the arguments to the infix operator, and it represents a function which intuitively takes an argument and puts it on the "missing" side of the infix operator.
>
```haskell
(2^) (left section) is equivalent to (^) 2, or more verbosely \x -> 2 ^ x
(^2) (right section) is equivalent to flip (^) 2, or more verbosely \x -> x ^ 2
```

つまり，中置記法の関数に関してはこのような感じで定義がされているということが予想できる．

```haskell
(infix_op e) = \ x -> x infix_op e
(e infix_op) = \ x -> e infix_op x
```

上記のような理解で問題なさそうだ．自分にとってはいまいちわかりづらい見た目をしていたので，しっかり定義を覚えておこう．

## 参考
- [functional programming - Why do both map (^2) xs and map (2^) xs work as expected in Haskell? - Stack Overflow](https://stackoverflow.com/questions/4703576/why-do-both-map-2-xs-and-map-2-xs-work-as-expected-in-haskell/4703600)
- [Section of an infix operator - HaskellWiki](https://wiki.haskell.org/Section_of_an_infix_operator)
