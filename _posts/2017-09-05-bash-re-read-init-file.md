---
layout: article
date: 2017-09-05 00:59
title: "bashで.inputrcの設定を読み込むには"
tags: [bash, readline]
comments: true
published: true
---

仮想環境を構築した後，そこでしばらく作業する予定のときは`.inputrc`ファイルに大文字小文字に依らず補完できる設定を書いて作業をするようにしている (`set completion-ignore-case on`を記述) のだけど，その設定を反映させるのにどうすればいいかやるたびに忘れるのでメモ．

```sh
bind -f ~/.inputrc
# もしくは，\C-x\C-r (re-read-init-file)でもOK.
```

`bind`はReadline (bashのコマンドラインを制御するライブラリ) で有効になっているkeyやfunctionを表示したり設定したりするshell builtin command．`-f`オプションで次に続く引数のファイルから設定を読み込む．

デフォルトでは`\C-x\C-r`がre-read-init-fileにbindされている．どのkeyがどの関数にbindされているかは`bind -p`で確認することができる．

```sh
$ bind -p | grep "\\\\C-x\\\\C-r"
"\C-x\C-r": re-read-init-file
```

-p: Display readline function names and bindings in such a way that they can be re-read.
