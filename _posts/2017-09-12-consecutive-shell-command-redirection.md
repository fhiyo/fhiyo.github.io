---
layout: article
date: 2017-09-12 22:47
title: "shell commandの連続リダイレクト"
tags: [shell, bash, redirect]
comments: true
published: true
---

```sh
$ echo foo 1>&2 1>file.txt
```

と

```sh
$ echo foo 1>file.txt 1>&2
```

は異なる結果になる．

```sh
$ echo foo 1>&2 1>file.txt
# 何も出力されない
$ cat file.txt
foo
$ echo foo 1>file.txt 1>&2
foo
$ cat file.txt
# 何も出力されない
```

これはなぜかというと，リダイレクトの解釈の順番は

1. **左から右へ行われ**，
2. ファイルディスクリプタへのリダイレクトは**現在のファイルディスクリプタの参照先へリダイレクトされる**

からである．

`echo foo 1>&2 1>file.txt`は，  
`1>&2`でファイルディスクリプタ1へのストリームをファイルディスクリプタの2が**現在**参照している参照先である標準エラー出力にリダイレクトし，  
`1>file.txt`で現在標準エラー出力にリダイレクトされているファイルディスクリプタの1へのストリームを`file.txt`にリダイレクトしている．  
よって`echo`の結果としては端末には何も出力されず，`file.txt`には`foo`の文字が出力されている．

一方，`echo foo 1>file.txt 1>&2`は，  
`1>file.txt`でファイルディスクリプタの1へのストリームを`file.txt`にリダイレクトし，  
`1>&2`で現在`file.txt`にリダイレクトされているファイルディスクリプタの1へのストリームをファイルディスクリプタの2が**現在**参照している参照先である標準エラー出力にリダイレクトしている．  
よって`echo`の結果としては端末に`foo`の文字が出力され (これは標準エラー出力として表示されている)，`file.txt`には何も出力されていない．

この挙動を理解すると，標準出力と標準エラー出力を入れ替えることができる．

```sh
$ { { echo stdout; echo stderr 1>&2; } 3>&1 1>&2 2>&3; } 1>1.txt 2>2.txt;
$ cat 1.txt
stderr
$ cat 2.txt
stdout
```

ファイルディスクリプタ1が標準エラー出力，ファイルディスクリプタ2が標準出力にストリームを流していることが確認できる．  
ただ，この標準出力と標準エラー出力の入れ替え，探せば解説しているところは結構出てくるんだけど，この入れ替えをする使い所が書いてあるサイトが見つからなかった．．もし知っている方がいらっしゃれば教えて欲しい．

別の例としては，`/dev/null`に出力全部を捨てるときもリダイレクトの仕組みを理解していないとはまる場合がある．

```sh
$ { echo stdout; echo stderr 1>&2; } 2>&1 >/dev/null
stderr
$ { echo stdout; echo stderr 1>&2; } >/dev/null 2>&1
# 何も出力されない
```

上の`{ echo stdout; echo stderr 1>&2; } 2>&1 >/dev/null`はファイルディスクリプタ2が標準出力へストリームを流し，ファイルディスクリプタ1が`/dev/null`へストリームを流すので標準エラー出力が出力されたままになってしまっている．  
下の`{ echo stdout; echo stderr 1>&2; } >/dev/null 2>&1`はファイルディスクリプタ1が`/dev/null`へリダイレクトされた状態でファイルディスクリプタ2が**現在の**ファイルディスクリプタ1の参照先である`/dev/null`へストリームを流すので全て出力を捨てることができる．

`/dev/null`へ出力を捨てるときはこういう書き方もできる．

```sh
$ { echo stdout; echo stderr 1>&2; } >/dev/null 2>/dev/null
```

これは，同じファイルに標準出力と標準エラー出力の両方をリダイレクトする場合には上手くいかない．

```sh
$ { echo -e "output\noutput" && echo -e "error" 1>&2; } 1>file.txt 2>file.txt
$ cat file.txt
error

output
```

これは別々のopen file descriptionが独自のファイルポジションを持っており，同期していないために後の`echo`文字列が前の`echo`文字列を書き換えてしまうかららしい．  
なので，標準出力，標準エラー出力の両方を同じファイルにリダイレクトしたい場合は
```sh
{ echo -e "output\noutput" && echo -e "error" 1>&2; } 1>file.txt 2>&1
```
か，または
```sh
{ echo -e "output\noutput" && echo -e "error" 1>&2; } &>file.txt
```
のように記述すること (下の書き方はbashのv4以降のみ対応)．この記法の場合はopen file descriptionは1つしか作成されず，ファイルポジションはディスクリプタ間で共有されるので上書きされることはない．

詳細はstackexchangeのここの質問を参照: [Why is the behavior of command 1>file.txt 2>file.txt different from command 1>file.txt 2>&1?](https://unix.stackexchange.com/questions/391456/why-is-the-behavior-of-command-1file-txt-2file-txt-different-from-command-1)  
↑実はこの`{ echo -e "output\noutput" && echo -e "error" 1>&2; } 1>file.txt 2>file.txt`の挙動の理由が不明だったので自分がした質問．open file descriptionとファイルディスクリプタの関係を教えてくれたので非常に良かった．

## 参考
- [What does the ampersand indicate in this bash command 1>&2](https://stackoverflow.com/questions/2341023/what-does-the-ampersand-indicate-in-this-bash-command-12)
- [Linuxのファイルディスクリプタをハックする](http://qiita.com/tajima_taso/items/7e696944cb521958e63d)
