---
layout: article
date: 2017-08-31 00:39
title: "zshでタブ補完のときにエラーが出る"
tags: [zsh, shell]
comments: true
published: true
---

zshでタブ補完をしようとすると，

```sh
(eval):setopt:3: no such option: NO_warnnestedvar
```

のようにエラーが出る．補完自体はできるのだが，補完するたびに出力されるので困る．

```sh
$ /bin/zsh --version
zsh 5.2 (x86_64-apple-darwin16.0)
$ /usr/local/bin/zsh --version
zsh 5.4.1 (x86_64-apple-darwin16.7.0)
```

`/usr/local/bin/zsh`を起動させたところ，タブ補完でエラーメッセージは出なかった．tab-completionで認識してるzshとログインしているzshにバージョンに食い違いがある，ということだろうか？  
試しにログインシェルを`/usr/local/bin/zsh`に変更しようとしたが，エラーが出て失敗．

```sh
$ chsh -s /usr/local/bin/zsh
Changing shell for user.
Password for user:
chsh: /usr/local/bin/zsh: non-standard shell
```

[OS X refuses to setting fish as default shell(installed via Homebrew)](https://github.com/fish-shell/fish-shell/issues/989)を参考に，上記のエラーを解決．ログインシェルとして認識されていないのが問題だったらしく，`/usr/local/bin/zsh`を`/etc/shells`に追加してやれば上手くいった．

```sh
$ sudo vim /etc/shells  # /usr/local/bin/zsh をリストに追加
$ cat /etc/shells
# List of acceptable shells for chpass(1).
# Ftpd will not allow users to connect who are not using
# one of these shells.

/bin/bash
/bin/csh
/bin/ksh
/bin/sh
/bin/tcsh
/bin/zsh
/usr/local/bin/zsh
$ chsh -s /usr/local/bin/zsh  # シェルを再起動

```

で問題は解決した．
