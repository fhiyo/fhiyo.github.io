---
layout: article
date: 2018-02-25 19:37
title: "Jekyllでjupyter notebookの画像を使う方法"
tags: [python, jekyll, jupyter-notebook]
comments: true
published: true
---

最近jupyter notebookを使うことが多くなったが，notebookで作ったグラフをブログに載せるのが面倒で，いい方法はないかと探っていた．  
[こちらのブログ](https://briancaffey.github.io/2016/03/14/ipynb-with-jekyll.html)にずばりな方法があったのでここでもやり方をメモしておく．

方法:
1. `jupyter notebook`を編集する
1. 編集したnotebookを以下のコマンドでmarkdown形式に変換する
    `jupyter nbconvert --to markdown jupyter-notebook-with-jekyll.ipynb jupyter-notebook-with-jekyll_files`  
    (jupyter-notebook-with-jekyll_filesは画像の出力されるディレクトリ名．省略可．)
1. markdownファイルがカレントディレクトリに，グラフのpng画像が`jupyter-notebook-with-jekyll_files`のディレクトリの下にそれぞれ出力されるので，それを記事に貼り付ける

これでOK．出力されたmarkdownファイルを記事に貼り付ける場合は，画像を各自適切な場所へ配置してその場所へのパスに置換すること．

```md
![png](jupyter-notebook-with-jekyll_files/jupyter-notebook-with-jekyll_2_1.png)
```

↑最初はこのようなパスになっているので，画像を再配置した場所へのパスに置換する．

上の方法で作成したjupyter notebookを貼り付けた結果がこちら:

```python
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set(style='white', context='notebook', palette='deep')

%matplotlib inline
```


```python
from sklearn import datasets
iris = datasets.load_iris()
iris_df = pd.DataFrame(iris.data, columns=iris.feature_names)
iris_df['target'] = pd.Series(iris.target)
iris_df.head()
```

    <class 'sklearn.utils.Bunch'>





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sepal length (cm)</th>
      <th>sepal width (cm)</th>
      <th>petal length (cm)</th>
      <th>petal width (cm)</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>5.1</td>
      <td>3.5</td>
      <td>1.4</td>
      <td>0.2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4.9</td>
      <td>3.0</td>
      <td>1.4</td>
      <td>0.2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>4.7</td>
      <td>3.2</td>
      <td>1.3</td>
      <td>0.2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4.6</td>
      <td>3.1</td>
      <td>1.5</td>
      <td>0.2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5.0</td>
      <td>3.6</td>
      <td>1.4</td>
      <td>0.2</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>




```python
sns.heatmap(iris_df.corr(), annot=True, fmt='.2f', cmap='coolwarm')
```




    <matplotlib.axes._subplots.AxesSubplot at 0x11043f978>




![png](/assets/images/jupyter-notebook-with-jekyll/jupyter-notebook-with-jekyll_2_1.png)
