---
layout: article
date: 2017-10-30 23:17
title: "Macのghci上でCtrl+yを入力するとsuspendされる"
tags: [haskell, ghci]
comments: true
published: true
---

MacのアプリケーションではCtrl-Yはヤンクに割り当てられていることが多いが，Haskellの対話環境であるghciを使っているときに，Ctrl-Yがおかしな動作をすることに気づいた．

```sh
$ ghci
Configuring GHCi with the following packages:
GHCi, version 8.0.2: http://www.haskell.org/ghc/  :? for help
Loaded GHCi configuration from /private/var/folders/5q/q2ghsy6s29z88zqbtthn1x5w0000gn/T/ghci57891/ghci-script
Prelude> -- ここでCtrl-Yを押す
Prelude> zsh: suspended (signal)  stack ghci --
$
```

Ctrl-Yを押すと，`zsh: suspend (signal)  stack ghci --`の文字列が出力されてghciがサスペンドしてしまう．

## 解決策

StackExchangeに書いてあった．[Haskell: GHCi treats Ctrl-Y like Ctrl-Z](https://stackoverflow.com/questions/46290504/haskell-ghci-treats-ctrl-y-like-ctrl-z)によると，要するに

```sh
stty dsusp undef
```

を`.zshrc`に設定すればOKのようだ.

`dsusp`が何かとか，なんでこれを設定すればOKなのかを少し調べてみる．

まず，`stty -a`で現在のterminalの設定を調べられるので，それを見てみる．

```sh
$ stty -a
speed 38400 baud; 40 rows; 158 columns;
lflags: icanon isig iexten echo echoe echok echoke -echonl echoctl
	-echoprt -altwerase -noflsh -tostop -flusho pendin -nokerninfo
	-extproc
iflags: -istrip icrnl -inlcr -igncr -ixon -ixoff ixany imaxbel iutf8
	-ignbrk brkint -inpck -ignpar -parmrk
oflags: opost onlcr -oxtabs -onocr -onlret
cflags: cread cs8 -parenb -parodd hupcl -clocal -cstopb -crtscts -dsrflow
	-dtrflow -mdmbuf
cchars: discard = ^O; dsusp = ^Y; eof = ^D; eol = <undef>;
	eol2 = <undef>; erase = ^?; intr = ^C; kill = ^U; lnext = ^V;
	min = 1; quit = ^\; reprint = ^R; start = ^Q; status = ^T;
	stop = ^S; susp = ^Z; time = 0; werase = ^W;
```

`dsusp = ^Y;`とあるのでCtrl-Yに`dsusp`なるものが割り当てられているのがわかる．

では`dsusp`は何かというと，ジョブコントロールのシグナルであるSIGTSTPをジョブに送出するようだ．  
で，`SIGTSTP`がジョブに送出されるとジョブはsuspendする．
ジョブをサスペンドするといえばメジャーなのはCtrl-Zを入力することだと思うけど，Ctrl-Zだと`SIGSTOP`シグナルを送出する．`SIGTSTP`の`SIGSTOP`の違いは，`SIGSTOP`が信号をすぐに送出するのに対して`SIGTSTP`はプログラムが入力としてシグナルを読み込んだときに送出されるらしい．要するに`SIGTSTP`の場合だとちょっと待ってからサスペンドされる，ということか？

ghciはどうもその`dsusp`の機能が有効になっているらしく，そのせいでCtrl-Yを押すとサスペンドされてしまったらしい．

よって`dsusp`をCtrl-Yに割り当てられている状態を解除することで，Ctrl-Yを入力してもシグナルが送出されなくなった．またアプリ側で上書きされていたヤンクの機能も使えるようになり大分使いやすくなった．


## 参考
- [pryがCtrl-Yでサスペンドしてしまう](http://neocat.hatenablog.com/category/mac)
- [17.4.9.2 Characters that Cause Signals](https://www.gnu.org/software/libc/manual/html_node/Signal-Characters.html)
