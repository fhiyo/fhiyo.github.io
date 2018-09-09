---
layout: article
date: 2017-09-21 01:54
title: "execコマンドについて"
tags: [bash, command]
comments: true
published: true
---

stackoverflowとかで時々，`exec 3>&1`という記述を見かける．何してる？


`man bash`の`exec`の項

```sh
       exec [-cl] [-a name] [command [arguments]]
              If  command  is  specified,  it  replaces  the  shell.  No new
              process is created.  The arguments  become  the  arguments  to
              command.   If  the  -l  option is supplied, the shell places a
              dash at the beginning of the zeroth argument  passed  to  com‐
              mand.   This is what login(1) does.  The -c option causes com‐
              mand to be executed with an empty environment.  If -a is  sup‐
              plied,  the  shell  passes  name as the zeroth argument to the
              executed command.  If command cannot be executed for some rea‐
              son,  a  non-interactive  shell exits, unless the shell option
              execfail is enabled, in which case  it  returns  failure.   An
              interactive  shell  returns failure if the file cannot be exe‐
              cuted.  If command is not  specified,  any  redirections  take
              effect  in  the current shell, and the return status is 0.  If
              there is a redirection error, the return status is 1.
```

最後の2行見る限り，`exec`の引数にcommandが指定されていない場合，リダイレクトが現在のシェルに適応されるようだ．つまり，`exec 3>&1`は現在のシェルでファイルディスクリプタの3を標準出力へ (1が標準出力を指しているなら) リダイレクトするように設定している，ということか．

`exec`の引数にcommandを指定するとどうなるのだろう？

```sh
$ echo $SHLVL
1
$ bash
$ echo $SHLVL
2
$ exec ls
foo.txt bar.txt
$ echo $SHLVL
1
```

`exec`の終了のタイミングで現在のシェルのプロセスが終了してしまった．これが`man`の説明にあった，"If  command  is  specified,  it  replaces  the shell."の指定したコマンドがシェルを置き換える，言っている意味なのだろうか？現在動作しているシェルのプロセスを`exec`で指定したコマンドである`ls`に持っていったため，`ls`の終了と同時にシェルが終了する，のだろうか．．

ということは，`exec fish`とすれば`bash`が`fish`に置き換わる？

```sh
$ echo $$
21185
$ echo $SHLVL
1
$ exec fish
Welcome to fish, the friendly interactive shell
Type help for instructions on how to use fish
vagrant@localhost ~> echo %self  # fishでは`echo $$`の代わりにこう書くらしい
21185
vagrant@localhost ~> echo $SHLVL
1
```

プロセスidが一致している．`fish`シェルも子プロセスとして起動しているわけではない．  
確かに`bash`のプロセスが`fish`のプロセスに置き換わっていることがわかった．  
(プロセスidが一致していることがプロセスを置き換えている証明になるのか怪しいけど)


関係ないけど，`fish`使いやすそうだな．．今まで`zsh`で俺は一生やっていくんだからいいんだと思って目を背けてたけどちょっと見てみたほうがいいかもしれない．
