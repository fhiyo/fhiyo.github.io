---
layout: article
date: 2017-10-20 11:45
title: "`show $ head [1, 2]`は良くて`show . head [1, 2]`がエラーになる理由"
tags: [haskell, difference]
comments: true
published: true
---

Haskellでプログラミングの問題を解いていて，`.`だとエラーになるけど`$`に変えてみるとコンパイルが通る，という状況に何度か出くわした．理由を調査したのでメモする．

### 現象

```haskell
> show . head [1, 2]

<interactive>:15:1: error:
    • Non type-variable argument in the constraint: Num (a1 -> a)
      (Use FlexibleContexts to permit this)
    • When checking the inferred type
        it :: forall a a1. (Num (a1 -> a), Show a) => a1 -> String
> show $ head [1, 2]
"1"
```

現象が再現される最小の構成にしたところ，どうも演算子の結合規則辺りの話かと予測がついた．

### 調査

まず，なぜ`show . head [1,2]`のコンパイルが通ると思ったのか，を考えてみる．  
`head`はリストモナドを引数に取ってそのリストの先頭の要素を返す函数であり，  
`show`は`Show`の型クラスのインスタンスである引数を取って，String型にして返す函数である．  
この2つを`.`で合成して，`show . head`という合成函数を作り，引数としてリストモナドを取ることでString型の返却値を取ることができるだろう，と考えたのが理由である．

この考え方は，暗黙的に次のように演算の順番を考えていることになる．

```haskell
-- 間違った解釈
> show . head [1,2]
-- => (show . head) [1,2]
-- => show 1
-- => "1"
```

一見良さそうにも見えるが，実際の挙動を見てみると以下のように解釈されているようである．

```haskell
> show . head [1,2]
-- => show . (head [1,2])
-- => show . 1              -- ERROR!
```

`.`は第二引数 (面倒なのでそう呼ぶ) は函数を指定するので，型が異なりエラーが吐かれる．  
実際に`head [1,2]`に括弧をつけて実行すると，問題のコンパイルエラーと同様のログが吐かれることもわかる．

```haskell
> show . (head [1,2])

<interactive>:18:1: error:
    • Non type-variable argument in the constraint: Num (a1 -> a)
      (Use FlexibleContexts to permit this)
    • When checking the inferred type
        it :: forall a a1. (Num (a1 -> a), Show a) => a1 -> String
```

というわけで，問題は`.`演算子と函数適用の優先順位がどちらが高いか，という問題に帰着できそうだ．`.`よりも函数適用の方が優先順位が高ければOK.

[Learn You a Haskell for Great Good!](http://learnyouahaskell.com/chapters)を調べたところ，それらしい記述があった．

[Higher order functions](http://learnyouahaskell.com/higher-order-functions)
> Putting a space between two things is simply function application. The space is sort of like an operator and it has the highest precedence.

というわけで，**函数適用は最も優先順位の高い演算子としてみることができる**ことがわかったので，コンパイルエラーが吐かれる理由もわかった．

`.`演算子を使いたいなら，`(show . head) [1,2]`のように括弧を使うか，`(.) show head [1,2]`のように`.`演算子を函数をして使うことで優先順位の高い函数適用を使えばOK．

一方，`$`の方は，

```haskell
> :type ($)
($) :: (a -> b) -> a -> b
```

なので第二引数は函数の形ではなく，

```haskell
> show $ head [1,2]
-- => show $ (head [1,2])
-- => show $ 1            -- ここが許される！
-- => "1"
```

となり，予想と同じ結果が返ってくる．



### 感想

中置記法は混乱することがあるので注意が必要ということがわかった．  
atcoderの問題解くときに，よく

```haskell
main = getContents >>= print . map fooFunc . map words . lines
```

みたいな感じで問題解いてるんだから調査しなくても気づきたかった感はある (`map fooFunc . map words`が`map (fooFunc . map) words`のようには解釈されていないので)．


### 参考
- [Learn You a Haskell for Great Good!](http://learnyouahaskell.com/chapters)
