---
layout: article
date: 2017-11-03 10:12
title: "Macのファイル検索CUIのmdfindについて"
tags: [mac, mdfind]
comments: true
published: true
---

ファイル検索をするとき，`find`を使うのは一般的だと思うが，`locate`の方がデータベースを先に構築してそこにクエリをかけるから速い，と聞いてからはそちらを使うようにしていた．  
しかし`locate`コマンドを使っていると，なぜかDropbox上のディレクトリが検索できなくて困っていた．これの原因を調べて見ると，以下の質問に答えが書いてあった．  
[Locate command can't find anything inside Documents folder on mac](https://stackoverflow.com/questions/15887431/locate-command-cant-find-anything-inside-documents-folder-on-mac)

> It's because your Documents-folder isn't world readable, which is a good thing, specially on shared systems.
>
> The BUGS section of the locate(1) man-page explains it:
>
> The locate database is typically built by user ''nobody'' and the locate.updatedb(8) utility skips directories which are not readable for user ''nobody'', group ''nobody'', or world. For example, if your HOME directory is not world-readable, none of your files are in the database.
Try running ls -ld ~/Documents and you'll see the permissions. Wikipedia have an article on Unix permissions if you are unfamiliar with these.

`locate`で使用しているデータベースである`/var/db/locate.database`のownerが`nobody`であるため，`nobody`がread権限を持ってないディレクトリは検索することができないとのこと．解決方法は質問の答えにもあるように，GNU版の`locate`をhomebrewでインストールするか，`mdfind`を使え，とある．`mdfind`の方はmacに標準で入っているやつだから，今回はそちらを使うことにしよう．

## mdfindについて
`mdfind`はmacのSpotlight検索のCUI版らしい．Spotlightも`locate`と同じで検索用のデータベースを持っていて，それを参照して検索結果を表示するのでファイルシステム全体を検索しても非常に速く結果が返ってくる．また，daemonによって更新は即座に反映されるようだ．

使い方としては，`mdfind -name foo`で`foo`の文字列が入っているファイルやディレクトリをrootから検索ができる．

```sh
$ mdfind -name foo
/Users/fhiyo/Project/example/server/system/composer/vendor/symfony/class-loader/Symfony/Component/ClassLoader/Tests/Fixtures/Namespaced/Foo.php
/Users/fhiyo/Project/example/server/application/view/layouts/admin/footer.html
/Users/fhiyo/Project/example/html/img/_common/dummy/dummy_footer.png
/Users/fhiyo/Project/example/html/_scss/footer/_footer.scss
/Users/fhiyo/Project/example/html/_scss/footer
/Users/fhiyo/Project/example/design/www/img/_common/dummy/dummy_footer.png
/Users/fhiyo/Project/example/design/www/_scss/footer/_footer.scss
/Users/fhiyo/Project/example/design/www/_scss/footer
```

のような感じ．`-name`オプションを無くせば全文検索 (ファイル名だけでなく，中身も見て検索) をしてくれるらしい．また，特定のパス下に検索対象に絞りたい場合は，`-onlyin <パス>`のオプションでOK.

しかし，しばらく使っていると`mdfind`にも問題があることが発覚．EvernoteやGoogle Driveとかの一部のアプリケーションが検索結果に反映されない．この現象はSpotlight検索でも見られていて，不便だな，と思っていたけど大きな問題にはならなかったので放置していた問題だった．

これについて調べて見ると，[Qiitaの記事](https://qiita.com/delphinus/items/438046d2bbeb3e63f8fa)に解決策が書いてあった．Spotlight検索がされないのはそのファイルやディレクトリがシンボリックリンクだからで，Spotlightはシンボリックリンクを検索対象にしない仕様になっているのだそう．[^1]シンボリックリンクではなくfinderのエイリアスというものを使えば検索対象にしてくれるらしい(bashとかのコマンドのエイリアスとはまた別物)．

エイリアスの作り方は，

```sh
$ osascript -e 'tell application "Finder" to make alias file to POSIX file "/file/to/make/link/from" at POSIX file "/folder/where/to/make/link"'
```

のようにコマンドを打てばよい．("/file/to/make/link/from"と"/folder/where/to/make/link"にはそれぞれエイリアスを作成したいファイルと作成先のパスを入力する)

問題も解決したので，後は`mdfind`を自分で使いやすくシェルコマンドのエイリアスを作成しておこう．

```sh
l () {
  # Usage: l [search-start-point] pattern
  if [ $# -eq 1 ]; then
    mdfind -onlyin $(pwd) -name $1
  else
    mdfind -onlyin $2 -name $1
  fi
}
```

とりあえずこんな感じで登録しておく．`l <pattern>`で現在のディレクトリ以下を，`l <file-path> <pattern>`で任意のディレクトリ以下を検索する．

ファイル検索関係は長いこと不便なまま放置してしまっていたので，ようやく楽になると思うと嬉しい．


## 参考
- [Locate command can't find anything inside Documents folder on mac](https://stackoverflow.com/questions/15887431/locate-command-cant-find-anything-inside-documents-folder-on-mac)
- [mdfind – macOSのCLIでもSpotlightを使ってファイルを高速全文検索する ｜ Developers.IO](https://dev.classmethod.jp/etc/spotlight-via-terminal/)
- [Mac/Linuxの「locate」コマンドで高速ファイル検索｜「find」コマンドとの違いから「mdfind」の紹介まで - QiitaQiita](https://qiita.com/kenju/items/862d59d29bb033eb0a8c)
- [brew linkapps したアプリが Spotlight でヒットしないとお嘆きのあなたに - QiitaQiita](https://qiita.com/delphinus/items/438046d2bbeb3e63f8fa)
- [osx - How do I create a Macintosh Finder Alias From the Command Line? - Stack Overflow](https://stackoverflow.com/questions/7072208/how-do-i-create-a-macintosh-finder-alias-from-the-command-line/10067437#10067437)

[^1]: 自分の場合は昔homebrewのディレクトリをSpotlightの検索対象に外していたからこんな現象が起こったのであって，本来rootから全てのディレクトリを検索するならシンボリックリンクは検索しないのが普通な実装である気がする．まあSpotlight検索対象を制限している人はこういう回避策もあるよ，ということで．
