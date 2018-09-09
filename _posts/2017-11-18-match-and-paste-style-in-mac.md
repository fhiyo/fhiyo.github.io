---
layout: article
date: 2017-11-18 19:50
title: "貼り付け先の書式に合わせてペーストする(macで)"
tags: [mac]
comments: true
published: true
---

evernoteとかにメモを取るとき，webページの文章をペーストする際に書式がコピー元のものと同じになってしまい，面倒なときがある．

<img src="/assets/images/match-and-paste-style-in-mac/before.png" alt="before" style="width: 400px;"/>


調べたところ，macの場合は**"Option + Shift + Command + V"**で書式を合わせてペーストしてくれるようだ．

<img src="/assets/images/match-and-paste-style-in-mac/after.png" alt="after" style="width: 400px;"/>

しかしこのコマンドでは指がつるのでデフォルトの"Command + V"の挙動をこちらに変更する．方法は以下の通り．

1. system preferencesのApplication shortcutを選択
1. ShortcutsタブのApp Shortcutsを選択し，+ボタンをクリック
1. "Paste and match style"の文字列を入力し，Keyboard Shortcutの欄で"Command + V"を入力する．
1. Addを押してOK.

<img src="/assets/images/match-and-paste-style-in-mac/config.png" alt="config" style="width: 400px;"/>

もし書式も一緒にペーストしたければ上のリボンからEdit→Pasteを選択すればよい．もし書式も一緒にペーストしたい機会が増えればまたショートカットを別に割り当てればよい．


## 参考

[Take Control of Your Mac's Clipboard](https://computers.tutsplus.com/tutorials/take-control-of-your-macs-clipboard--mac-30472)
