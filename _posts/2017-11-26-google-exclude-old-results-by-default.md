---
layout: article
date: 2017-11-26 10:38
title: "googleの検索で古い情報をフィルタリングする"
tags: [google]
comments: true
published: true
---

googleで検索するときに，情報が古いものに当たることが最近増えてしまったので，デフォルトの検索で古い情報を除外するようにする．今回はGoogle Chrome (Mac) で行ったが，他のブラウザでもデフォルトの検索を設定できれば変更可能なはず．

## 手順

1. Chromeの環境設定を開き，"Manage search engines"を開く  
<img src="/assets/images/google-exclude-old-results-by-default/search-engine.png" alt="search-engine" style="width: 600px;"/>
2. "Default search engines"の右下にある"Add"の文字をクリックする．  
<img src="/assets/images/google-exclude-old-results-by-default/add.png" alt="add" style="width: 600px;"/>
3. "Add search engine"のポップアップが出るので，以下のように入力する． (Search engineとKeywordは任意)  
<img src="/assets/images/google-exclude-old-results-by-default/edit-search-engine.png" alt="edit-search-engine" style="width: 600px;"/>  
4. "Add"をクリックして"Other search engines"の項に追加する．  
5. 追加したカラムの右にある縦3点のマークをクリックして，"Make default"をクリックする．

やったことは要するにデフォルトのサーチエンジンのクエリを設定して，1年前までに更新されているページのみを検索の対象にする，というもの．


## 解説

`https://www.google.com/search?q=%s&tbs=qdr:y&tbo=1`で一年前までに更新されているページのみを表示する．  
`tbs=qdr:y`が該当の部分．検索する更新日時を変更したい場合は，以下のようなクエリパラメータに変えてやればよい．

- 指定なし: `tbs=qdr:a`
- 1秒前までに更新: `tbs=qdr:s`
- 1分前までに更新: `tbs=qdr:n`
- 10分前までに更新: `tbs=qdr:n10`
- 1時間前までに更新: `tbs=qdr:h`
- 10時間前までに更新: `tbs=qdr:h10`
- 1日前までに更新: `tbs=qdr:d`
- 1週間前までに更新: `tbs=qdr:w`
- 1ヶ月前までに更新: `tbs=qdr:m`
- 1年前までに更新: `tbs=qdr:y`
- 1984年3月2日から1987年6月5日の間に更新: `tbs=cdr:1,cd_min:3/2/1984,cd_max:6/5/1987`
- 更新日時でソート: `tbs=sbd:1`
- 関連性が高い順でソート: `tbs=sbd:0`


`tbo=1`は検索のツールバーを表示させるクエリオプションのようだが，`tbs`の指定をするとパラメータ無くても勝手にパラメータつけられているようだ．自分が試したところ`tbo=1`の有無で結果に違いが現れなかった．


## 参考

- [URL parameters of the Google search engine results page - SEO Hero](https://seoheronews.com/url-google)
- [Google Search URL Request Parameters DETECTED](https://stenevang.wordpress.com/2013/02/22/google-advanced-power-search-url-request-parameters/)
- [Setting Google to Return Recent Search Results by Default - cameronreilly.com](http://cameronreilly.com/setting-google-to-return-recent-search-results-by-default/)
- [【小ワザ】Googleの検索結果にデフォルトでフィルターをかける方法](https://conpath.net/columns/filter-google-search-results-by-default/)
