---
layout: article
date: 2017-09-08 00:04
title: "declare -rとreadonlyの関係について"
tags: [bash, command]
comments: true
published: true
---

bashのshellscriptを書くときに，変数を書き換え不可にして定義をしたいときがよくある．

```sh
#!/bin/bash

FILEPATH=$0
...
FILEPATH=another_value
```

この場合，`${FILEPATH}`はこのスクリプトファイルへのパスなのだから後で別の値を代入できるようにしたくない．今までこのようなときは`readonly`のbuiltin commandを使って実現していたが，`declare -r`という定義の方法でも同じことができるらしい．何か違いがあるのかを調べてみた．

結論から言うと，変数のスコープに違いが生じる．`declare -r`はローカルスコープ, `readonly`はグローバルスコープになる．

```sh
#!/bin/bash

function func () {
  foo=10
  local bar=20
  declare -r baz=30
  readonly qux=40
  declare -rg quux=50
}

func

echo ${foo}
echo ${bar}
echo ${baz}
echo ${qux}
echo ${quux}
```

Output:  
```sh
$ ./example.sh
10


40
50
```

foo: 修飾子なし．関数内で定義された変数は関数外でも参照可  
bar: `local`修飾子あり．`local`で修飾された変数は，定義された関数内とその子の関数内でしか可視でなくなる．  
baz: `declare -r`で修飾された変数は書き換え不可となる．スコープはローカルになる．  
qux: `readonly`で修飾された変数は書き換え不可となる．スコープはグローバルとなる．  
quux: `declare -rg`で修飾された変数は書き換え不可となる．スコープはグローバルとなる．  

`readonly`と`declare -rg`は同じと見てよさそう．また，`man`を読んでいて`declare`の説明のところに`typeset`というbuiltin commandがあったが，`bash`においては完全に`declare`と同じものでobsoluteなコマンドらしい．

ちなみに，Mac環境のbashでは`declare -g`のオプションは存在していないようだった．  

```sh
$ cat ex.sh
#!/bin/bash

func () {
  declare -gr foo=30
}

func

echo "foo: " ${foo}
$ ./ex.sh
./ex.sh: line 4: declare: -g: invalid option
declare: usage: declare [-afFirtx] [-p] [name[=value] ...]
foo:
$ sw_vers
ProductName:	Mac OS X
ProductVersion:	10.12.6
BuildVersion:	16G29
$ bash --version
GNU bash, version 4.4.12(1)-release (x86_64-apple-darwin16.3.0)
Copyright (C) 2016 Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>

This is free software; you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
```

元々調査した環境は  
CentOS 7.3 (GNU bash, version 4.2.46(1)-release (x86_64-redhat-linux-gnu))



おまけ: `local`で修飾された変数は孫の関数内でも参照できるのか試したもの  
```sh
#!/bin/bash

function f () {
  local hoge=1000
  echo ${hoge}
  function g () {
    echo ${hoge}
    function h () {
      echo ${hoge}
    }
    h
  }
  g
}
f
```

Output:  
```sh
$ ./example.sh
1000
1000
1000
```

#### 参考
[In bash script, what's the different between declare and a normal variable?](https://unix.stackexchange.com/questions/254367/in-bash-script-whats-the-different-between-declare-and-a-normal-variable)

[what is difference in `declare -r` and `readonly` in bash?](https://stackoverflow.com/questions/30362831/what-is-difference-in-declare-r-and-readonly-in-bash)
