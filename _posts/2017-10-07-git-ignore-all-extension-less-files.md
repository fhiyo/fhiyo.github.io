---
layout: article
date: 2017-10-07 11:43
title: ".gitignoreに拡張子なしのファイルを登録する"
tags: [git, gitignore]
comments: true
published: true
---

この前gccでコンパイルした実行ファイルを無視するために，拡張子がないファイルだけをgitの管理対象から外したいと思ってやり方を調べたのでメモする．

[How do I add files without dots in them (all extension-less files) to the gitignore file?](https://stackoverflow.com/questions/19023550/how-do-i-add-files-without-dots-in-them-all-extension-less-files-to-the-gitign) にあるStackoverflowの回答によると，以下のように`.gitignore`などのファイルに記述すれば (or 変数に値を入れれば) 実現できると書いてあった．

```
*
!*/
!*.*
```

ただし，上の記述は先頭に加えないと動作しない (後勝ちなのでこのルールで前のルールが全て無効になってしまう)．  
また，2行目は`!/**/`と書いてもOK.  

このルールは，

- 1行目の`*`で全てのファイル，ディレクトリを管理対象から外し，
- 2行目の`!*/`で全てのディレクトリを管理対象にre-includeして，
- 3行目の`!*.*`で全ての拡張子付きのファイルを管理対象にre-includeする

という意味になっている．2行目がなければディレクトリ以下のファイルをre-includeすることができない．理由はgitのdocumentationにある．

>     The pattern dir/ excludes a directory named dir and (implicitly) everything under it.
With dir/, Git will never look at anything under dir, and thus will never apply any of the “un-exclude” patterns to anything under dir.
      The pattern dir/* says nothing about dir itself; it just excludes everything under dir. With dir/*, Git will process the direct contents of dir, giving other patterns a chance to “un-exclude” some bit of the content (!dir/sub/).

あるディレクトリがgitの管理対象から外れてしまうと，その下の全てのディレクトリ及びファイルはパフォーマンスの理由からどうあっても再び管理対象にすることはできないらしい．そのため`!*.*`で拡張子付きのファイルを管理対象にし直す前に，`!*/`や`!/**/`で全てのディレクトリを管理対象にし直す必要がある．


余談だが，gitで管理しているリポジトリ内にあるファイルをgitの管理対象から外す方法は`.gitignore`にパターンを記述する以外にもあり，全部で3つあるようだ．

1. `.gitignore`にパターンを記述する  
1. `${GIT_DIR}/info/exclude`にパターンを記述する  
1. `core.excludesFile` (`.gitconfig`内の変数)にパターンを記述する  

1はそのリポジトリをクローンする他の開発者にも同じように無視してほしいパターンを，  
2は自分の環境のリポジトリ内でのみ無視したいパターン (他の環境ではそのパターンはいらない場合) を，  
3は自分のユーザー環境全体で常に無視したいパターンを  

それぞれ記述するという住み分けになっており，上のものほど優先順位が高い．

参考:
- [How do I add files without dots in them (all extension-less files) to the gitignore file?](https://stackoverflow.com/questions/19023550/how-do-i-add-files-without-dots-in-them-all-extension-less-files-to-the-gitign)
- [Git docs -- gitignore PATTERN FORMAT](https://git-scm.com/docs/gitignore#_pattern_format)
- [.gitignore exclude folder but include specific subfolder](https://stackoverflow.com/questions/5533050/gitignore-exclude-folder-but-include-specific-subfolder/20652768#20652768)
