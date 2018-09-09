---
layout: article
date: 2017-11-11 17:24
title: "リンクコピーを簡単にできるChrome拡張を作ってみた"
tags: [chrome-extension, programming]
comments: true
published: true
---

Google Chromeの拡張機能であるChrome Extensionを作ったのでメモを残しておく．Chrome Extensionの作成方法が知りたい場合は，[Chromeのオリジナル拡張機能を開発しよう（ソースコードあり） 株式会社LIG](https://liginc.co.jp/web/tool/browser/163575)や[Chrome拡張の開発方法まとめ その1：概念編 - QiitaQiita](https://qiita.com/edit-mode/items/26d7a22233ecdf48fed8)のリンクが参考になったので見ておくといいかも．

## 何を作ったか

ショートカットキーを押すと，今見ているページのリンクをコピーする．markdown形式と単純にurlをシングルクォーテーションかダブルクォーテーションで括ったものかを選択できるようにしている．

作ったものへのリンク: [fhiyo/simple-link-clipper: Clip web page's link](https://github.com/fhiyo/simple-link-clipper)

## 動機

今回自分が作ったものは最初いくらでもあるだろうと思って探したが，なぜかバッジをクリックしなきゃいけないものが多く，マウス使いたくない人間にとっては痒い所に手が届いてない感じだった．やりたいことは非常にシンプルなので，練習がてら作ってみるか，と思い作成をした．

## 上手くいかなかったところ

最初はvimiumのように表示しているページによって活性/非活性を変えられる機能を持たせようと思ったけど，どうも簡単にはいかないらしい．作りたい機能としてはそこまで優先度が高いものでは無かったので，今回は断念．

また，ポップアップのページからショートカットキーを変更できるようにしたかったんだけど，それも上手くいかなかったので断念した．`chrome://extensions/configureCommands`への遷移ができればいいんだけど，セキュリティの関係かポップアップページにハイパーリンクをつけても移動してくれなかった．．．

## 作ってみて

メモ書きの作業効率が地味に結構改善されている．markdown形式のコピーもそうだし，昔からシングルクォートでurlを囲った文字列が欲しいときが時々あるので，それを一発でできるようになったのもありがたい．簡単な機能だけど割りと気に入っているし，今後他に拡張が欲しくなったときにすぐに開発できるようになったのは大きい．

## 参考
- [Chromeのオリジナル拡張機能を開発しよう（ソースコードあり） 株式会社LIG](https://liginc.co.jp/web/tool/browser/163575)
- [Chrome拡張の開発方法まとめ その1：概念編 - QiitaQiita](https://qiita.com/edit-mode/items/26d7a22233ecdf48fed8)
- [Getting Started: Building a Chrome Extension - Google Chrome](https://developer.chrome.com/extensions/getstarted)
