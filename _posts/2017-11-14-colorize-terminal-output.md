---
layout: article
date: 2017-11-14 12:34
title: "コマンドの出力に色をつける"
tags: [command, linux]
comments: true
published: true
---

`cat`や`less`など，terminalで何か文字列を出力する際には色がついてないと読むのがキツイ．  
今までも一部のコマンドには色をつけて表示されるようにしていたが，出力系コマンドで色をつけて表示してくれるようなコマンドを調べたのでメモする．  
調査したのはMacOS上でだが，見る限りどれもLinuxでもインストールできるっぽい．


## ls→exa

`exa`: [ogham/exa: Replacement for 'ls' written in Rust.](https://github.com/ogham/exa)

`exa`はRustで書かれた`ls`の代わりになるツール．

`brew install exa`でインストール可能．


<img src="/assets/images/colorize-terminal-output/exa-sample.png" alt="exa-sample" style="width: 600px;"/>  
↑綺麗に表示してくれる

普通の`ls`と違って`-T`でtree表示ができたり，`git`のrepositoryはbranchの表示もされる．

## cat→ccat

`ccat`: [jingweno/ccat: Colorizing cat](https://github.com/jingweno/ccat)

`brew install ccat`でインストール．

<img src="/assets/images/colorize-terminal-output/ccat.png" alt="ccat" style="width: 500px;"/>

## less→source-highlight

`brew install sourch-highlight`でインストール．

```sh
export LESS='-R'
export LESSOPEN='| /usr/local/bin/src-hilite-lesspipe.sh  %s'
```

を`.zshrc`などに設定後，読み込みすればOK．


## diff→colordiff

`brew install colordiff`でインストール．


## pingなど→grcを使って色をつける

`brew install grc`でインストール．

```sh
[[ -s "/usr/local/etc/grc.zsh" ]] && source /usr/local/etc/grc.zsh
```

を.zshrcに追記すればOK．bashやfishなら`grc.zsh`を`grc.bashrc`や`grc.fish`に変えてそれぞれの設定ファイルに書けばOK.

<img src="/assets/images/colorize-terminal-output/color-ping.png" alt="color-ping" style="width: 500px;"/>

色をつけることができるコマンドは以下．

```
cc           ls
configure    make
cvs          mount
df           mtr
diff         netstat
dig          ping
gcc          ping6
gmake        ps
ifconfig     traceroute
last         traceroute6
ldap         wdiff
```

## メモ

自分は設定を以下のようにした．

```sh
####################
#     Colorize     #
####################

# read grc setting
[[ -s "/usr/local/etc/grc.zsh" ]] && source /usr/local/etc/grc.zsh

# ls
alias ls="exa"
alias ll="exa -la"
alias la="exa -l"
alias l1="exa -1"
alias lll="exa -abghHliS --git"

# less
export LESS='-R'
export LESSOPEN='| /usr/local/bin/src-hilite-lesspipe.sh  %s'

# cat
alias cat="ccat"

# diff
alias diff='colordiff -u'

### End Colorize ###
```

aliasをしているが，元々の`ls`とか使いたい場合は`\ls`のようにバックスラッシュをコマンドの前につければaliasは無効化されるので当面はこの設定でやってみる．`exa`やら`colordiff`やらがインストールされてない場合の処理は省略しているが，つけたい場合は

```sh
if [[ -x `which ccat` ]]; then
  alias cat='ccat'
fi
```

のように記述すればよい．


`grc`の`diff`の設定と`colordiff`の設定がバッティングする．`.zshrc`の設定読み込みの順番で制御してもいいが，変更に弱いのであまりやりたくない．  
ので，`/usr/local/etc/grc.zsh`を直接書き換えて`diff`の部分をコメントアウトすることにした．気になる人は`grc.zsh`をコピーしてユーザーローカルなところに置いてからそいつを変更して読み込む感じにすれば良さそう．

`cat`と`less`は`pygments`使って色つける方法もあるようだ．どちらが良さそうかは未検証．`pygments`は`pip install pygments`でインストールして使えば良さそう．細かい設定方法は
- [catやlessをシンタックスハイライト+行番号表示してコードを見やすくする - QiitaQiita](https://qiita.com/zaru/items/4a7811ac21f74a13481c)
- [cat, less コマンドの表示を Syntax Highlight させる - xykのブログ](http://xyk.hatenablog.com/entry/2014/12/24/161507)

この辺りのブログで紹介されている．
