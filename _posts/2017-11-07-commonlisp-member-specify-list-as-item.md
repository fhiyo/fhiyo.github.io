---
layout: article
date: 2017-11-07 14:11
title: "CommonLispのmemberでitemとしてリストを指定する"
tags: [common-lisp, programming]
comments: true
published: true
---

Common Lispの`member`はリスト内にitemがあるかどうかを検索するための関数だが，itemがリストのときに上手くいかずにハマっていた．

```lisp
> (member '(2 3) '(() (1 2 3 4) (2 3) (2)))
NIL  ; '(2 3)は検索対象内にあるのにNILが返ってくる
```

結論としては，`member`のキーワード引数の`test`に`equal`関数を指定することで上手くいった．

```lisp
> (member '(2 3) '(() (1 2 3 4) (2 3) (2)) :test #'equal)
((2 3) (2))  ; 見つけた要素から後ろの部分のリストが返る
```

`equal`は比較関数だが．異なるリスト同士の比較でも値が同じなら真を返す．，`eq`だとS式として同一でないと偽を返すのでダメ．

```lisp
> (eq '(1 2) '(1 2))
NIL
```

余談 (これを書いているときに知った) だが，単なる数値の比較でも値が大きい場合は別のS式としてみなされるようだ．

```lisp
> (eq 100000000000000 100000000000000)
T
> (eq 1000000000000000 1000000000000000)
NIL
```

これ，Pythonとかでも同じ現象が起きる (`a=1000; id(a)==id(1000)`がFalseになる) ので，気をつけていたほうがいいかもしれない．

最後に`member`関数の定義を確認する． (参考: [CLHS: Function MEMBER, MEMBER-IF, MEMBER-IF-NOT](http://clhs.lisp.se/Body/f_mem_m.htm))

```
Syntax:

member item list &key key test test-not => tail

item---an object.
list---a proper list.
predicate---a designator for a function of one argument that returns a generalized boolean.
test---a designator for a function of two arguments that returns a generalized boolean.
test-not---a designator for a function of two arguments that returns a generalized boolean.
key---a designator for a function of one argument, or nil.
tail---a list.
```

リファレンスを見ても`test`のデフォルトの関数がわからなかったのだが，`eq`がデフォルトなのだろうか？
