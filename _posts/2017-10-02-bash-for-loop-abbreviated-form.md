---
layout: article
date: 2017-10-02 13:13
title: "bashのforループのin句は省略できる"
tags: [bash, command]
comments: true
published: true
---

bashで配列内にある要素が存在しているか確認するための関数を書こうと思って，stackexchangeとかにあるかなと思って見てみたところ，見慣れない`for`の使い方があったのでメモする．

```sh
containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}
```

上の関数の使い方は以下の通り．

```sh
$ array=("something to search for" "a string" "test2000")
$ containsElement "a string" "${array[@]}"
$ echo $?
0
$ containsElement "blaha" "${array[@]}"
$ echo $?
1
```

出典元: [Check if a Bash array contains a value](https://stackoverflow.com/questions/3685970/check-if-a-bash-array-contains-a-value#answer-8574392)

この`for e; do hogehoge; done`の`for`の使い方が初めてだった．`in`のstatementが省略されている場合，どういう挙動になるのだろうか．

`for`文中で`in`が省略されている場合，ループ変数には関数の引数が順に格納され，引数を全て格納し終わったらループを抜ける，という実装になっているようだ．

つまり，以下の書き方と同じ．

```sh
for e in "$@"; do hogehoge; done
```

一応確認をしてみる．

```
#!/usr/bin/env bash

loop_test () {
  for e; do
    echo ${e}
  done
}

loop_test foo bar baz
```

これを`sample.sh`の名前で保存し，デバッグモードで実行してみる．

```sh
$ bash -x sample.sh
+ loop_test foo bar baz
+ for e in "$@"
+ echo foo
foo
+ for e in "$@"
+ echo bar
bar
+ for e in "$@"
+ echo baz
baz
```

デバッグしたときの命令文は`for e in "$@"`となっていた．というわけで確かに`for e`は`for e in "$@"`と同じであることがわかった．自分でこのように書く機会は少ないかもしれないが，他の人のコード見るときにわからなくなるからしっかり覚えておこう．

注: 検証環境はmacOS 10.12.16のbash 4.4.12と3.2.57である．
