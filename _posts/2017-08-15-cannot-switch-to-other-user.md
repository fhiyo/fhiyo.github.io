---
layout: article
date: 2017-08-15 23:39
title: "suで他のユーザーに切り替えができない"
tags: [shell, security]
comments: true
published: true
---

EC2インスタンスで作業を行っていたとき，他のユーザーでの操作を確認したい場面があったので

```sh
$ sudo su user1
```

と入力したところ，

This account is currently not available.

と言われてユーザーの切り替えに失敗した．とりあえず出力されたメッセージでググってみる．  
["This account is currently not available" ... Please help](https://www.centos.org/forums/viewtopic.php?t=56978)がヒットした．`/etc/passwd`の中身を変更して`/sbin/nologin`にしたせいでrootアカウントに入れなくなった？らしい．対象ユーザーがrootかどうかの違いはあるものの，自分と似た現象が起こっている．`/etc/passwd`の中身を確認する．

```sh
$ cat /etc/passwd
root:x:0:0:root:/root:/bin/bash
bin:x:1:1:bin:/bin:/sbin/nologin
...
user1:x:1000:1000:bin/bin:/sbin/nologin
...
```

まず，`/etc/passwd`の見方をよく知らなかったので調べる．  
[Understanding /etc/passwd File Format](https://www.cyberciti.biz/faq/understanding-etcpasswd-file-format/)によると，`/etc/passwd`は以下のフォーマットになっているようだ．  

```
USERNAME:PASSWORD:UID:GID:USER_ID_INFO:HOME_DIRECTORY:COMMAND/SHELL
```

各属性がコロン区切りで置かれている．
- USERNAME: ユーザーがログインするときに使われるユーザー名．1-32文字の長さであることが条件．
- PASSWORD: ここが`x`となっている場合は`/etc/shadow`で暗号化されているパスワードである，という意味らしい．`/etc/passwd`ファイル自体はデフォルトで誰でも読める (permissionが644) のでここにパスワードが書かれないのは当然．
- UID: ユーザーid．0はroot, 1-99はあらかじめ定義されているアカウント．100-999は管理用にsystem側が予約しているようだ．
- GID: グループid．`/etc/group`で一覧を確認できる．
- USER_ID_INFO: ユーザーに関するコメント．
- HOME_DIRECTORY: ユーザーがログインしたときのホームディレクトリ．指定なしだと`/`になる．
- COMMAND/SHELL: ログイン時に実行されるコマンド/シェル．


確かに，user1でログインを行うと`/sbin/nologin`が実行されるようになっているようだ．
以下のコマンドを実行してログイン時に実行されるコマンドを書き換える．

```sh
$ sudo usermod -s /bin/bash user1
```

`usermod -s`でユーザーのログインシェルを変更できる．もう一度`/etc/passwd`を確認すると，

```sh
$ cat /etc/passwd
root:x:0:0:root:/root:/bin/bash
bin:x:1:1:bin:/bin:/sbin/nologin
...
user1:x:1000:1000:bin/bin:/bin/bash
...
```

書き換わっていることを確認．この状態で`sudo su user1`をすると，正常にuser1としてbashを起動させることができた．

おそらく他の人がsambaユーザーを作成するときにsshでログインできないように設定をしたのを自分が知らずにその人でログインしようとして失敗した，というところなのだろう．．．

ちなみに，`nologin`実行時に出力されるメッセージは`/etc/nologin.txt`を追加してそこに好きなメッセージを入れることで編集できるようだ．
