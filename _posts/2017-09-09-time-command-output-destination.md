---
layout: article
date: 2017-09-09 12:21
title: "timeコマンドの出力先はどこなのか？"
tags: [bash, command, shell]
comments: true
published: true
---

`time`コマンドの出力先が標準出力だったか標準エラー出力だったかを確認しようと思って，コマンドを叩いてみた．

```sh
$ time echo foo 1>/dev/null

real	0m0.000s
user	0m0.000s
sys	0m0.000s
$ time echo foo 2>/dev/null
foo

real	0m0.000s
user	0m0.000s
sys	0m0.000s
```

標準出力も標準エラー出力も出力先ではない？？じゃあ出力されているのはなんなんだ．．`man time`で確認してみる．

>      ...When command finishes, time writes a message to  standard
       error giving timing statistics about this program run.

標準エラー出力に書き込むと書いてある．ではなぜ`2>/dev/null`で出力先が`/dev/null`に切り替わらないのか？

### 理由

`time`の文法は以下の通り．

```sh
time [options] command [arguments...]
```

`time echo foo 2>/dev/null`の場合，
- options: なし
- command: `echo`
- arguments: `foo 2>/dev/null`
となる．つまり，`time`は`echo foo 2>dev/null`の実行時間を計測して標準エラー出力に書き込んでいた，ということ．

では，`echo foo`までを`time`コマンドの引数にするにはどうすればよいかというと，code blockかsubshellを使って`time`コマンドに対する引数の終わりを明示してやればよいようだ．

```sh
$ time { echo foo; } 2>/dev/null  # code blockの場合
foo
$ time ( echo foo ) 2>/dev/null  # subshellの場合
foo
```

ちなみにcode blockを使うときは，最後のセミコロンを忘れないこと．セミコロンがないと文法エラーになる．
