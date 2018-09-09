---
layout: article
date: 2017-10-25 18:15
title: "vim-surroundで文字列を括弧で囲むときに空白をつけない方法"
tags: [vim, code-reading]
comments: true
published: true
---

普段[vim-surround](https://github.com/tpope/vim-surround)を使っているが，括弧をつけるときに選択した文字列の間にスペースが入ってしまうのを解決する方法があったのでメモする．


## 先に結論

vim-surroundをインストールした環境では，選択モードに入ってSを押した後，
- bで`(`
- Bで`{`
- rで`[`
- aで`<`

をスペースなしで選択部分を囲むことができる．

他にも，
- tで任意のタグ (<でもOK)
- sでスペース (スペースでもOK)
- pで空行
- lでLatexのEnvironments (`\begin{foo} ~ \end{foo}`のやつ)

で囲むことができる．

括弧を削除するときは，対応する括弧の中に入って  
`ds(` (括弧の種類は消したいものに対応) と打てばよい．


## 現象

vim内でテキストを選択モードで範囲選択した後，`S`でvim-surroundのコマンド受け付けの状態に入り，`(`を押すと範囲選択した部分を括弧で囲める．  
しかし，囲んだ部分の間にスペースが入ってしまい困っていた．

```txt
foo bar baz
```

のようなテキストがあるとき，barを選択した状態で`S(`と打つと，

```txt
foo ( bar ) baz
```

のようになる．この`bar`と括弧の間にあるスペースを何とかしたい．

## 解決法

ドキュメントを見てもいいが，今回はVimscriptを読む練習ということでソースコード見てみた．

Visual modeで`S`を押してvim-surroundの入力待ち状態にする処理は，

```vimscript
  xmap S   <Plug>VSurround
```

の部分で行っているようだ．調べると`xmap`は`Visual mode`でのコマンドをmapするらしい．  (参考: [Vim documentation: map  1.3 MAPPING AND MODES](http://vimdoc.sourceforge.net/htmldoc/map.html))  
`VSurround`の関数を辿っていくと，

```vimscript
function! s:wrap(string,char,type,removed,special)
```

で問題の文字列を囲む処理を行っているようだ．関数の中を見てみると，`stridx(pairs, newchar)`という関数呼び出しの結果を`idx`という変数に格納している部分がある．`newchar`はユーザーが入力した文字のようである．`pairs`は`"b()B{}r[]a<>"`という文字列が格納されている．

ここで，`stridx()`のhelpをvim内のnormal modeで`:help stridx`と入力して見てみると，

```vimscript
stridx({haystack}, {needle} [, {start}])		*stridx()*
		The result is a Number, which gives the byte index in
		{haystack} of the first occurrence of the String {needle}.
		If {start} is specified, the search starts at index {start}.
		This can be used to find a second match: >
			:let colon1 = stridx(line, ":")
			:let colon2 = stridx(line, ":", colon1 + 1)
		The search is done case-sensitive.
		For pattern searches use |match()|.
		-1 is returned if the {needle} does not occur in {haystack}.
		See also |strridx()|.
		Examples: >
		  :echo stridx("An Example", "Example")	     3
		  :echo stridx("Starting point", "Start")    0
		  :echo stridx("Starting point", "start")   -1
```

とある．最後の例だけ見れば雰囲気はつかめて，{haystack}内に{needle}がある場合，{haystack}内で一番最初に{needle}が出現する開始位置を返し，{needle}が無ければ-1を返す関数のようだ．`wrap()`関数のもう少し先を見ると，

```vimscript
  elseif idx >= 0
    let spc = (idx % 3) == 1 ? " " : ""
    let idx = idx / 3 * 3
    let before = strpart(pairs,idx+1,1) . spc
    let after  = spc . strpart(pairs,idx+2,1)
```

というelse文がある．`spc`には`stridx`の{needle}が`(`,`{`,`[`,`<`のいずれかの場合に空白が入る．この部分が自分が問題にしていた原因の箇所だろう．  
`before`と`after`にはユーザーが選択した括弧 (`(`，`{`，`[`，`<`のいずれか) で囲むための文字列が入る．

ここまで見ると，`pairs = "b()B{}r[]a<>"`と`spc`の処理から，`b`，`B`，`r`，`a`を入力すれば`spc`に空白が入らないため，上手くできるのではないかと予測がつく．

実際に，`foo bar baz`の文字列でbarを選択して`Sb`と入力したところ，`foo (bar) baz`と空白をつけずに括弧を挿入できた．

`SB`で`[`，`Sr`で`{`，`Sa`で`<`を使って文字列を囲める．  
ついでなので他のコマンドも少し調べて"先に結論"の項に書いておいた．


## 感想

括弧の前後に空白が入らないようなオプションがあるかと思って探したけど，空白入れない用のコマンドを別に作っていたのは意外だった．  
ドキュメント読むだけじゃなく実際にコードを追うのも中々楽しい．
