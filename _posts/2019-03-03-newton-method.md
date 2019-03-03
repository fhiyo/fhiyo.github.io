---
layout: article
date: 2019-03-03 17:41
title: "Newton method"
tags: [math]
comments: true
published: true
---

Newton法について勉強したので記録する<sup><a href="#1">1</a></sup>。なお、この記事ではベクトルをボールド体にせず、スカラと同じように表記する。

## 語句の定義

### 狭義凸関数 (strictly convex function)

凸集合$X$を定義域とする関数$f$が、任意の2点 $x_0, x_1 \in X$, $0 \le \lambda \le 1$である任意の実数$\lambda$に対して、

$$
f(\lambda x_0 + (1 - \lambda)x_1) < \lambda f(x_0) + (1 - \lambda) f(x_1)
$$

が成立するとき、$f$は狭義凸関数とよばれる。

(凸集合の定義は[前回の記事](/2019/02/23/sum-of-convex-function-is-convex-function.html)を参照のこと。)

つまり狭義凸関数は凸関数の定義から等号を除いたバージョンということである。


## Newton法がやっていること

Newton法は、狭義凸関数 (strictly convex function) $f$に対して、step wiseに$f$を最小化するパラメータを求めるアルゴリズムである。

なお、以下の説明で出てくる$k$は非負整数であり、$f$を最小化するためにパラメータを更新する回数のことである。

最小化したい$f$は$C^n (n \ge 2)$級の狭義凸関数であるとする。$f(\theta)$について、2次のテーラー展開で近似した$\theta_k$周りの関数を$f_{\text{quad}}(\theta_k)$とすると、 ($\theta \in \mathbb{R}^n, n \in \mathbb{N}$)

$$
\begin{align}
f_{\text{quad}} &= f_k + (\nabla f(\theta_k))^{T}(\theta - \theta_k) + \frac{1}{2}(\theta - \theta_k)^TH_k(\theta - \theta_k) \tag{1}\\
\end{align}
$$

ただし、$H_k$は$f_k$についてのヘッセ行列であり、

$$
H_k = \begin{bmatrix}
\frac{\partial}{\partial \theta_1 \partial \theta_1}&\dots&\frac{\partial}{\partial \theta_1 \partial \theta_n}\\
\vdots&\ddots&\vdots\\
\frac{\partial}{\partial \theta_n \partial \theta_1}&\dots&\frac{\partial}{\partial \theta_n \partial \theta_n}\\
\end{bmatrix} f(\theta_k) \tag{2}
$$

である。

ここで、以下のように記号を導入する。

$A = \frac{1}{2}H_k$

$b = \nabla f(\theta_k) - H_k\theta_k$

$c = f_k - (\nabla f(\theta_k))^{T}\theta_k + \frac{1}{2}\theta_k^TH_k\theta_k$

このようにおくと、式$(1)$は

$$
f_{\text{quad}} = \theta^TA\theta + b^T\theta + c
$$

のように変形できる。これは$f_{\text{quad}}$を$\theta$について整理した形である。

さて、ここで$f_{\text{quad}}$の極値を求める。

$$
\begin{align}
\nabla f_{\text{quad}} &= \frac{\partial}{\partial \theta}\theta^TA\theta + \frac{\partial}{\partial \theta}b^T\theta + \frac{\partial}{\partial \theta}c = 0\\
&\Rightarrow \frac{\partial}{\partial \theta}\theta^TA\theta + \frac{\partial}{\partial \theta}b^T\theta = 0\\
&\Rightarrow (A + A^T)\theta + \frac{\partial}{\partial \theta}b^T\theta = 0\\
&\Rightarrow (A + A^T)\theta + b = 0\\
\end{align}
$$

ところで、式$(2)$と、$f$が$C^n (n \ge 2)$級の関数であることから$A$の各成分の微分は順序交換が可能である。よって$A$は対称行列であることがわかるので、

$$
2A\theta + b = 0\\
$$

となる。よって、__$A$が正則であるならば__、

$$
\begin{align}
\tilde \theta &= \frac{1}{2}A^{-1}b\\
&= H_k^{-1}(\nabla f(\theta_k) - H_k\theta_k)\\
&= \theta_k - H_k^{-1}\nabla f(\theta_k)\\
\end{align}
$$

上式を満たす$\tilde \theta$は元の関数$f$を二次近似した関数$f_{\text{quad}}$が最小値となるパラメータ (minimizer)なので、$\tilde \theta$は$\theta_k$よりも$f$の最小値に近いことが期待できる。

よって、$k$番目のstep sizeを$\eta_k \in \mathbb{R}$とすると、<sup><a href="#2">2</a></sup>

$$\theta_{k+1} \leftarrow \theta_k - \eta_k H_k^{-1}\nabla f(\theta_k)$$

のように$\theta$を更新することを$\theta$の値が収束するまで繰り返せば、$\theta$がminimizerに近い値になっているであろう、というアルゴリズムの考え方である。

...

ところで、$A$、つまりヘッセ行列が正則である場合はminimizerを求めることができたが、$f$が狭義凸関数であるときにヘッセ行列が正則になる証明はしていない。これは次回以降の記事に回そうと思う (本記事執筆時点で証明はできているはず)。

### Reference

1. Machine Learning: A Probabilistic Perspective, Kevin P. Murphy, 2012
2. [偏微分の順序交換の十分条件とその証明 高校数学の美しい物語](https://mathtrain.jp/henbibunexchange)

---

<span id="1" style="font-size:x-small">1: Newton法、
$$
x _ { n + 1 } = x _ { n } - \frac { f \left( x _ { n } \right) } { f ^ { \prime } \left( x _ { n } \right) }
$$
のように値を更新する、という説明が書いてある記事がちらほら見られるけど、これはどのように違うものなのだろうか？違う手法に同じ名前がついてしまっているのか、同じものだけど流派が違うよ、くらいのものなのか。。？
</span>

<span id="2" style="font-size:x-small">2: なんでstep sizeが出てくるのかよく分かっていない。。そのまま$\tilde \theta$に更新するんじゃダメなんだろうか？</span>


<!-- 「半正定値」ってことは、固有値が0のやつがn個中何個かある状態を許すということだよな？ -->
<!-- ということはその状態は正則じゃないということだ。。。 -->
<!--  -->
<!-- $f_{\text{quad}}$を計算するときにhessian matrixの逆行列を使ったが、こいつが存在しない場合でも計算できる、ということか？ -->
<!--  -->
<!-- 狭義凸関数じゃない場合に逆行列が存在しない場合があるな。 -->
<!-- quasi-newton methodsを使うとかになるのか。 -->
<!--  -->
<!-- 別の、一次の微分しか使わないバージョンのnewton methodsは何なんだろうか？？ -->
<!-- どう違うの？？ -->
<!--  -->
<!--  -->
<!-- 微分の順序を交換できることについて言及しないとhessian matrixが対称行列であることは言えない。 -->

