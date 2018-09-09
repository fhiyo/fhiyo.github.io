---
layout: article
date: 2017-10-23 17:24
title: "vimでsnippetのインデント深くするときはハードタブじゃないとダメ"
tags: [vim, neosnippet]
comments: true
published: true
---


## 問題

vimの[neosnippet](https://github.com/Shougo/neosnippet.vim)を使っているときに，インデントが上手く動作しない場合があった．

```vim
snippet     initial
options     head
    #include <iostream>
    int main(int argc, char const* argv[])
    {
      ${1}
      ${2}

      ${3}
      return 0;
    }
```

のように書いて，snippetを起動させると以下のように，

<img src="/assets/images/neosnippet-indent-must-be-hard-tab/vim-cpp-snippet.png" alt="cpp-snip-wrong-result" style="width: 400px;"/>

一つ目のplaceholderより下のインデントでおかしくなる部分が出てしまう．

## 解決策

ググってみたところ，docsの[ここ](https://github.com/Shougo/neosnippet.vim/blob/master/doc/neosnippet.txt#L699)を見ろという内容のissue[^1]を発見した．

> If you use hard-tab for indentation inside a snippet file, neosnippet will use
'shiftwidth' instead of Vim indent plugin. This feature is useful while some
languages' indent files do not work very well (e.g.: PHP, Python).

ハードタブのインデントをスニペットの登録側ですると，呼び出しの際には対応するfiletypeのshiftwidthの値分だけインデントされるように変換されて展開されるようだ．  
というわけで以下のように先程のスニペット内のスペースをタブに変換して，

```vim
snippet     initial
options     head
	#include <iostream>
	int main(int argc, char const* argv[])
	{
		${1}
		${2}

		${3}
		return 0;
	}
```

再度snippetを起動させる．

<img src="/assets/images/neosnippet-indent-must-be-hard-tab/vim-cpp-snippet-2.png" alt="cpp-snip-correct-result" style="width: 400px;"/>

今度は正しくインデントされて展開された！

今後の対策としては，`ftplugin/snippet.vim`に以下の記述を書いてsnippetファイルを編集するときだけハードタブが挿入されるように変更する[^2]．

```vim
setl noexpandtab
```


## 感想

vimを使い始めて1年半くらいになるが，まだ時々Emacsが恋しくなる．Emacs時代は毎日エディタの設定をしていたから最終的にかなり使いやすくなったが，小指痛めてvimに移ってからあまりいじらなくなったからなぁ．．もう少し頑張って軽くVim scriptを触れるくらいにはならなければ．

[^1]: リンク: [Improper indentation #370](https://github.com/Shougo/neosnippet.vim/issues/370)
[^2]: ここを参考にした: [neosnippet のスニペットファイルを書く上での注意点](http://rhysd.hatenablog.com/entry/2012/11/11/222150)
