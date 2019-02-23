---
layout: article
date: 2019-02-23 14:34
title: "凸関数の和は凸関数"
tags: [math]
comments: true
published: true
---

凸関数の和は凸関数になるが、それの証明を確認したので記録しておく。

以下では簡単化のために集合$X$は有限次元ユークリッド空間$R^m$の部分集合であるとする。

まず用語の定義から。

### 凸結合


$X$の2つの点$\boldsymbol{x}_0$, $\boldsymbol{x}_1$をとり、$0 \le \lambda \le 1, \lambda \in \mathbb{R}$である$\lambda$を一つ定める。このとき、

$$
\boldsymbol{x} = \lambda \boldsymbol{x}_0 + (1 - \lambda) \boldsymbol{x}_1
$$

で表される点の集合$\boldsymbol{x}$を、$\boldsymbol{x}_0$, $\boldsymbol{x}_1$の凸結合という。

===

この凸結合は2点を結ぶ線分のことです (3点以上に対しても凸結合は定義できる)。

(念の為補足するとこれは、

$$
\begin{align*}
\boldsymbol{x} &= \boldsymbol{x}_1 + \lambda (\boldsymbol{x}_1 - \boldsymbol{x}_0)
\end{align*}
$$

となることから、$\boldsymbol{x}$は$\boldsymbol{x}_1$を通り$\boldsymbol{x}_1 - \boldsymbol{x}_0$の方向に
$|\boldsymbol{x}_1 - \boldsymbol{x}_0|$の$\lambda$倍だけ動く点の集合だから2点を結ぶ線分になる。)


### 凸集合

集合$X$がそこに含まれる任意の2点から作られる凸結合が$X$に含まれるとき、集合$X$は凸集合であるという。


### 凸関数

凸集合$X$を定義域とする関数$f$が、任意の2点 $\boldsymbol{x}_0, \boldsymbol{x}_1 \in X$, $0 \le \lambda \le 1$である任意の実数$\lambda$に対して、

$$
f(\lambda \boldsymbol{x}_0 + (1 - \lambda)\boldsymbol{x}_1) \le \lambda f(\boldsymbol{x}_0) + (1 - \lambda) f(\boldsymbol{x}_1)
$$

が成立するとき、$f$は凸関数とよばれる。


### 凸関数の和は凸関数

$n$を正の整数として、$f_1, \dots, f_n$が凸関数であるとき、$\sum_{i=1}^{n} f_i$は凸関数になる。

[証明]  
$\boldsymbol{x}_0, \boldsymbol{x}_1 \in X$, $0 \le \lambda \le 1$とすると、

$$
\begin{align*}
(\sum_{i=1}^{n} f_i)(\lambda \boldsymbol{x}_0 + (1 - \lambda)\boldsymbol{x}_1) &= \sum_{i=1}^{n} f_i(\lambda \boldsymbol{x}_0 + (1 - \lambda)\boldsymbol{x}_1)\\
&\le \sum_{i=1}^{n} \lambda f_i(\boldsymbol{x}_0) + (1 - \lambda) f_i(\boldsymbol{x}_1)\\
&= \lambda \sum_{i=1}^{n} f_i(\boldsymbol{x}_0) + (1 - \lambda) \sum_{i=1}^{n} f_i(\boldsymbol{x}_1)\\
&= \lambda (\sum_{i=1}^{n} f_i)(\boldsymbol{x}_0) + (1 - \lambda) (\sum_{i=1}^{n} f_i)(\boldsymbol{x}_1)\\
\end{align*}
$$


## Reference

`http://web.econ.keio.ac.jp/staff/ito/pdf97/me97fuku.pdf`
