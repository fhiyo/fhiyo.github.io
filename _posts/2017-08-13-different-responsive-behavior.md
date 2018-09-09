---
layout: article
date: 2017-08-13 12:32
title: "本番環境のレスポンシブ対応の挙動がローカルと異なる"
tags: [css, responsive-web-design, jekyll, github-pages]
comments: true
published: true
---

自分が運営しているサイトのレスポンシブ対応の挙動が，ローカルで見たときと本番環境で見たときで異なっていることに気づいた．

![responsive-different-behavior](/assets/images/different-responsive-behavior/different-responsive-behavior.png)  
左がローカルで動かしたときのサイトのトップページ，右が本番環境でのサイトのトップページ．

圧倒的に本番の方が見にくい．これは何とかしないと困る．挙動が異なるということは，ローカルでは読み込まれていたcssが本番環境では読み込まれていない，ということが起こっているのだと思われる．

ググってみたところ，[ローカルとサーバーでCSSの挙動が異なるときに考えられる２つの原因](http://www.sho-dw.com/web/41/)というサイトがヒットした．ここにはOSの違いやブラウザのレンダリングモードの違いによって，ファイルパスやcssセレクタの大文字小文字の違いの解釈が異なることにより読み込まれない場合があると書かれている．
自分のローカル環境はWindowsではなくMacであり，本番環境 (多分Linux) とのファイルパスの解釈の違いはないはず．また，レンダリングモードも自分はChromeで確認してるのでこの場合問題ない．

もう少しローカルとgithub pages上のサイトを比較してみる．デベロッパーツールを使って適応されているcssを一つ一つ確認．すると，`/assets/css/style.css`に書いていたcssの内容が適用されていないことがわかった．なぜかと思い，デベロッパーツール上でSourcesタブから`/assets/css/style.css`を確認．すると，自分で定義していたはずのstylesheetが違う内容になっていた．中身を見ると，`normalize.css v4.1.1 | MIT License | github.com/necolas/normalize.css`と書いてある．

どうもnormalize.cssはブラウザ側の環境によらず，同じようにサイトを見せてくれる機能を提供するためのstylesheetのようだ．しかしデフォルトのstyle.cssはなぜ削除されてしまうのだろうか？

[GitHub Pages - CSS is incorrect on live site, fine on local host](https://stackoverflow.com/questions/44557181/github-pages-css-is-incorrect-on-live-site-fine-on-local-host)  
[Jekyll site works locally but not on Github Pages](https://stackoverflow.com/questions/42450554/jekyll-site-works-locally-but-not-on-github-pages?rq=1)  
同じような現象で困っている人もいた．Stackoverflowで聞こうかと思ったけど，同じ内容の質問はダメなので[問い合わせフォーム](https://github.com/contact)からgithubのサポートに聞いてみることにする．

送ってもらった回答は[Customizing your Jekyll theme's CSS](https://help.github.com/articles/customizing-css-and-html-in-your-jekyll-theme/#customizing-your-jekyll-themes-css)を見てみてくれとのこと．リンクを見ると，cssを自分でいじるときはcssではなくscssのファイルを使い，publishのときにcssファイルを生成する形にするようだ．うーん，もっとwebサイト構築の方法勉強しなきゃいけないなぁ．．．


というわけで，`/assets/css/style.css`を`/assets/css/style.scss`にrenameして再びリポジトリをpushしてみる．

![fixed-responsive-iphone](/assets/images/different-responsive-behavior/fixed-responsive-iphone.png)  
直った！スタイルシートも自分が定義したものになっていることを確認．

これで直ったということは，github pagesではデフォルトでstyle.scssがstyle.cssを生成していて，その際に自分で作成したstyle.cssは上書きされてしまうのだろうか？自分でstyle.scssを定義してやればデフォルトのstyle.scssは動かないので想定通りの挙動になるのか．とりあえずは直ったが，workaroundな解決法だからどこかのタイミングでサイトを改修してやらなきゃいけないかもなぁ．


<!--site-->
