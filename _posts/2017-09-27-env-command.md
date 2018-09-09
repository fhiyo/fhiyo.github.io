---
layout: article
date: 2017-09-27 22:53
title: "envコマンドについて"
tags: [bash, command]
comments: true
published: true
---

pythonスクリプトの1行目に`#!/usr/bin/env python`とかを書くが，`env`コマンドが何をやっているのかいまいちわかっていなかった．のでメモする．

まず，`env`コマンドをそのまま実行すると，環境変数の一覧が見れる．

```sh
$ env
XDG_SESSION_ID=214
HOSTNAME=localhost.localdomain
SELINUX_ROLE_REQUESTED=
TERM=xterm
...(中略)...
OLDPWD=/home
```

また，`env`の後に引数としてコマンドを指定すると，`${PATH}`を頭から順に見ていき，指定したコマンドを実行したところでそいつを実行する．(ようだ)  
そのため，`#!/usr/bin/env python`が頭に書いてあるスクリプトは，shell上で`python`の引数としてではなく直接実行 (`./foo.py`のように) されたとき，ユーザーの`${PATH}`環境変数を頭から見て最初に見つかったパスの`python`からスクリプトが実行されるようになっている．  
これは環境によって`/usr/bin/python`だったり`/usr/local/bin/python`だったりと`python`の居場所が違うことがあるので，どこに`python`がいてもパスさえ通っていればちゃんとスクリプトが動作するようにしているのが`#!/usr/bin/env python`の意味になるようだ．

`env COMMAND`の挙動の実験をしてみる．

```sh
$ cat foo.sh
#!/bin/bash

set -u

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "This function\'s location is ${DIR}."
```

この`foo.sh`を`/home/vagrant`と`/home/vagrant/tmp`の両方に置いた．これらが呼ばれるとそのスクリプトが存在しているフルパスを返す．

この状態で以下の実験を行った．

```sh
$ echo ${PATH}
/sbin:/bin:/usr/sbin:/usr/bin
$ # パスを追加
$ export PATH=/home/vagrant:/home/vagrant/tmp:${PATH}
$ echo ${PATH}
/home/vagrant:/home/vagrant/tmp:/sbin:/bin:/usr/sbin:/usr/bin
$ env foo.sh
This function\'s location is /home/vagrant.
$ # 追加したパスをワンライナーで削除 (関数化したほうがいいかも)
$ P=/home/vagrant; export PATH=$(echo ${PATH} | awk -v RS=: -v ORS=: '$0 != "'${P}'"' | sed 's/:$//')
$ P=/home/vagrant/tmp; export PATH=$(echo ${PATH} | awk -v RS=: -v ORS=: '$0 != "'${P}'"' | sed 's/:$//')
$ echo ${PATH}
/sbin:/bin:/usr/sbin:/usr/bin
$ # 今度はパスの順番を入れ替えて再度export
$ export PATH=/home/vagrant/tmp:/home/vagrant:${PATH}
$ env foo.sh
This function\'s location is /home/vagrant/tmp.
```

確かに`${PATH}`の前にあるパスから先に読まれている．  
本当はこういうの，カーネルのソースコードとか読んでどうなってるのか見てみたいんだけど，どこから探して見ればいいのかあまり分かってないので仕方なくこうやって実験．．．

`man env`や`info coreutils 'env invocation'`でドキュメント見てみたんだが，"`env COMMAND`は`${PATH}`を前から見て最初に見つかったコマンドを実行します"という旨の記述が見当たらなかった．



## 参考
- [What type of path in shebang is more preferable?](https://askubuntu.com/questions/88257/what-type-of-path-in-shebang-is-more-preferable)
