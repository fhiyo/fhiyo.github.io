---
layout: article
date: 2017-08-31 14:18
title: "treeコマンドで文字化けを防ぐ (+他のオプション)"
tags: [shell, command]
comments: true
published: true
---

`tree`コマンドを使ってると、日本語のファイル名が文字化けしてしまい困っていた。

```sh
$ ls
file.txt         ファイル.txt
$ tree .
.
├── file.txt
└── �\203\225�\202��\202��\203�.txt

0 directories, 2 files

```

`tree -N`のオプションを使えば、マルチバイト文字もそのまま表示してくれるので文字化けは起きないようだ．

- -N: Print non-printable characters as is instead of as escaped octal numbers.

```sh
$ tree -N
.
├── file.txt
└── ファイル.txt

0 directories, 2 files
```

ついでなので他のオプションも見てみよう。 (version: 1.7.0)

- -l: Follows  symbolic  links if they point to directories, as if they were directories. Symbolic links that will result in recursion are avoided when detected.

```sh
$ ls -1
fuga.txt
hoge.txt
tmp
$ ls -1 ../tmp3
piyo.txt
$ tree
.
├── fuga.txt
├── hoge.txt
└── tmp -> ../tmp3

1 directory, 2 files
$ tree -l
.
├── fuga.txt
├── hoge.txt
└── tmp -> ../tmp3
    └── piyo.txt

1 directory, 3 files
```

最後に返すディレクトリとファイルの数も`-l`オプションの有無によって変わるようだ．  
(`--noreport`オプションで最後のディレクトリとファイル数の出力は返らなくなる)

- -f: Prints the full path prefix for each file.

```sh
$ tree -Nf
.
├── ./file.txt
├── ./tmp
│   └── ./tmp/tmp -> ..
└── ./ファイル.txt

2 directories, 2 files
```

- -D: Print the date of the last modification time or if -c is used, the last status change time for the file listed.

```sh
$ tree -DNc
.
├── [Aug 31 23:09]  ファイル.txt
├── [Aug 31 23:10]  file.txt
└── [Aug 31 23:16]  tmp
    └── [Aug 31 23:16]  tmp -> ..

2 directories, 2 files
```

`-D`でファイルの更新時刻を表示，`-c`で更新時刻順にソート．

- -I pattern: Do not list those files that match the wild-card pattern.

```sh
$ tree -N
.
├── file.txt
├── tmp
│   └── tmp -> ..
└── ファイル.txt

2 directories, 2 files
$ tree -IN "file*"
.
├── tmp
│   └── tmp -> ..
└── ファイル.txt

2 directories, 1 file
$ tree -IN "file*|tmp"
.
└── ファイル.txt

0 directories, 1 file
```

除外したいパターンが複数ある場合は，`|`を使って区切って指定すればよい．

- -P pattern: List only those files that match the wild-card pattern.

```sh
$ tree -N
.
├── file.txt
├── tmp
│   └── tmp -> ..
└── ファイル.txt

2 directories, 2 files
$ tree -NP file.txt
.
├── file.txt
└── tmp

1 directory, 1 file
$ tree -NP file.txt --prune
.
└── file.txt

0 directories, 1 file
```

ディレクトリはパターンマッチングの結果に関係なくそのまま表示されるようだ．  
該当しないディレクトリを表示したくない場合は，`--prune`オプションを使えばよい．  
ただし現行バージョン (1.7.0) では`--prune`オプション使用時は`tree`の結果がメモリに全部載ってから出力されるようなので，結果が巨大な場合は出力が遅くなるなどの現象が起こるようである．

> The --prune and --du options cause tree to accumulate the entire tree in memory before emitting it. For large directory trees this can cause a significant delay in output and the use of large amounts of memory.

- -J: Turn on JSON output. Outputs the directory tree as an JSON formatted array.
- -X: Turn on XML output. Outputs the directory tree as an XML formatted file.

```sh
$ tree -J
[{"type":"directory","name": ".","contents":[
    {"type":"file","name":"file.txt"},
    {"type":"directory","name":"tmp","contents":[
      {"type":"link","name":"tmp","target":"..","contents":[]}
    ]},
    {"type":"file","name":"ファイル.txt"}
  ]},
  {"type":"report","directories":2,"files":2}
]
$ tree -X
<?xml version="1.0" encoding="UTF-8"?>
<tree>
  <directory name=".">
    <file name="file.txt"></file>
    <directory name="tmp">
      <link name="tmp" target=".."></link>
    </directory>
    <file name="ファイル.txt"></file>
  </directory>
  <report>
    <directories>2</directories>
    <files>2</files>
  </report>
</tree>
```

`-J`, `-X`のオプションのときは`-N`を指定しなくてもマルチバイト文字が文字化けしなかった．
