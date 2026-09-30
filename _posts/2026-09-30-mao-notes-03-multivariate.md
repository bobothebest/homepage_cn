---
layout: post
title: "茆书复习笔记 03｜多维分布、协方差与条件期望"
author: 张梓源
tags:
- 茆书复习笔记
- 概率论
- 多维分布
- 条件期望
- 茆诗松
- 考研
- 读书笔记
date: 2026-09-30 10:03 -0400
toc: true
---
> **导读**
>
> 多维随机变量这一章，是笔记里标注“要会证”最密集的部分：切比雪夫不等式、施瓦茨不等式、相关系数为 ±1 的充要条件、协方差矩阵非负定、重期望公式……
>
> 好消息是，这些证明的思路其实很统一，基本都是“**构造一个非负的量，然后看它什么时候等于 0**”。读的时候可以留意这条主线。

---

## 1. 常用多维分布

| 分布 | 页码 | 要点 |
|---|---|---|
| 多项分布 $$M(n;p_1,\dots,p_r)$$ | P145 | $$P(X_1=n_1,\dots,X_r=n_r)=\dfrac{n!}{n_1!\cdots n_r!}p_1^{n_1}\cdots p_r^{n_r}$$，$$\sum n_i=n$$ |
| 多维超几何分布 | P146 | $$\dfrac{\binom{N_1}{n_1}\cdots\binom{N_r}{n_r}}{\binom Nn}$$ |
| 多维均匀分布 | P147 | 在区域 $$D$$ 上 $$p(x)=\frac1{S_D}$$ |
| ☆ 二元正态 $$N(\mu_1,\mu_2,\sigma_1^2,\sigma_2^2,\rho)$$ | P148 | 见下 |
| 二维指数分布 | P152 | $$F(x,y)=1-e^{-x}-e^{-y}+e^{-x-y-\lambda xy}$$，$$x,y>0$$ |
| n 元正态分布 | P189 | $$p(\boldsymbol x)=\dfrac{1}{(2\pi)^{n/2}\vert \boldsymbol B\vert ^{1/2}}\exp\left\{-\tfrac12(\boldsymbol x-\boldsymbol a)^T\boldsymbol B^{-1}(\boldsymbol x-\boldsymbol a)\right\}$$ |

**二元正态的密度**

$$p(x,y)=\frac{1}{2\pi\sigma_1\sigma_2\sqrt{1-\rho^2}}\exp\left\{-\frac{1}{2(1-\rho^2)}\left[\frac{(x-\mu_1)^2}{\sigma_1^2}-2\rho\frac{(x-\mu_1)(y-\mu_2)}{\sigma_1\sigma_2}+\frac{(y-\mu_2)^2}{\sigma_2^2}\right]\right\}$$

---

## 2. 多维分布与边际分布的关系（P135）

- 多项分布的一维边际分布仍为二项分布：$$X_i\sim b(n,p_i)$$（会记）。
- 二维正态的边际分布为一维正态：$$X\sim N(\mu_1,\sigma_1^2)$$，$$Y\sim N(\mu_2,\sigma_2^2)$$（P150，会推）。
- **具有相同边际分布的联合分布可以不同（P156）**：例如二维指数分布中，不管 $$\lambda$$ 取什么值，两个边际分布都是 $$Exp(1)$$，可见边际分布不能决定联合分布。同理，二维正态的边际分布与 $$\rho$$ 无关。

<details markdown="1"><summary>二维正态边际分布的推导要点</summary>

对 $$y$$ 配方：令 $$u=\frac{x-\mu_1}{\sigma_1}$$，$$v=\frac{y-\mu_2}{\sigma_2}$$，指数部分可写为

$$-\frac{u^2}{2}-\frac{(v-\rho u)^2}{2(1-\rho^2)}.$$

对 $$y$$ 积分时，后一项是均值 $$\rho u$$、方差 $$1-\rho^2$$ 的正态核，积分后只剩一个常数，于是

$$p_X(x)=\frac{1}{\sqrt{2\pi}\sigma_1}e^{-\frac{(x-\mu_1)^2}{2\sigma_1^2}}.$$

</details>

---

## 3. 独立性

- 多维随机变量相互独立时，其中一部分随机变量与另一部分随机变量也相互独立（指导书 P149）。
- 多维随机变量相互独立时，联合分布可由边际分布唯一确定（边际分布之积）。
- 独立性可以按定义判断（$$p(x,y)=p_X(x)p_Y(y)$$），也可以从实际背景判断。
- 对离散随机变量，考试以**会应用**为主：熟悉书中各种典型的判别方式即可。

**连续场合的卷积公式（P167）**：$$X,Y$$ 独立时，$$Z=X+Y$$ 的密度为

$$p_Z(z)=\int_{-\infty}^{\infty}p_X(z-y)\,p_Y(y)\,dy=\int_{-\infty}^{\infty}p_X(x)\,p_Y(z-x)\,dx.$$

**$$\max$$ 与 $$\min$$ 的恒等式**（习题 3.3.12）

$$\min\{a,b\}=\tfrac12\big[a+b-\vert b-a\vert \big],\qquad \max\{a,b\}=\tfrac12\big[a+b+\vert b-a\vert \big].$$

应用（习题 3.4.14、3.4.18）：$$X,Y\overset{iid}{\sim}N(0,1)$$ 时 $$E\max\{X,Y\}=\frac1{\sqrt\pi}$$；$$X,Y\overset{iid}{\sim}N(a,\sigma^2)$$ 时 $$E\max\{X,Y\}=a+\frac{\sigma}{\sqrt\pi}$$。

<details markdown="1"><summary>推导</summary>

$$X-Y\sim N(0,2)$$，而 $$E\vert Z\vert =\sigma_Z\sqrt{2/\pi}$$，所以 $$E\vert X-Y\vert =\sqrt2\cdot\sqrt{2/\pi}=\frac{2}{\sqrt\pi}$$。  
于是 $$E\max\{X,Y\}=\frac12\left[E(X+Y)+E\vert X-Y\vert \right]=\frac1{\sqrt\pi}$$。一般情形令 $$X=a+\sigma X'$$ 即可。

</details>

---

## 4. 期望、方差与协方差

**切比雪夫不等式（P90，要求会证明）**

$$P\big(\vert X-E(X)\vert \ge\varepsilon\big)\le\frac{\mathrm{Var}(X)}{\varepsilon^2},\qquad\forall\varepsilon>0.$$

<details markdown="1"><summary>证明（连续型）</summary>

$$P(\vert X-EX\vert \ge\varepsilon)=\int_{\vert x-EX\vert \ge\varepsilon}p(x)\,dx\le\int_{\vert x-EX\vert \ge\varepsilon}\frac{(x-EX)^2}{\varepsilon^2}p(x)\,dx\le\frac1{\varepsilon^2}\int_{-\infty}^{\infty}(x-EX)^2p(x)\,dx=\frac{\mathrm{Var}(X)}{\varepsilon^2}.$$

更一般地（习题 2.3.12）：若 $$g(x)$$ 非负不减且 $$E g(X)$$ 存在，则 $$P(X>\varepsilon)\le\dfrac{E\,g(X)}{g(\varepsilon)}$$。

</details>

**定理 2.3.2：$$\mathrm{Var}(X)=0\iff P(X=a)=1$$（要会证）**

<details markdown="1"><summary>证明</summary>

($$\Leftarrow$$) $$X$$ 几乎处处等于常数 $$a$$，所以 $$E(X)=a$$，$$\mathrm{Var}(X)=E(X-a)^2=0$$。

($$\Rightarrow$$) 记 $$a=E(X)$$。由切比雪夫不等式，对每个 $$n$$ 有 $$P\left(\vert X-a\vert \ge\frac1n\right)\le n^2\,\mathrm{Var}(X)=0$$，所以

$$P(X\ne a)=P\left(\bigcup_{n=1}^{\infty}\left\{\vert X-a\vert \ge\tfrac1n\right\}\right)\le\sum_{n=1}^{\infty}P\left(\vert X-a\vert \ge\tfrac1n\right)=0.$$

</details>

**均值与中位数的最优性**（习题 2.3.9）
- $$E(X-EX)^2\le E(X-c)^2$$：平方偏差在均值处最小（P264，定理 5.3.2 是其样本版本）。
- $$E\vert X-m\vert \le E\vert X-c\vert$$：绝对偏差在中位数 $$m$$ 处最小。

**用分布函数求期望**（习题 2.2.20、2.2.21）
- 连续随机变量：$$E(X)=\int_0^{\infty}[1-F(x)]\,dx-\int_{-\infty}^{0}F(x)\,dx$$
- 非负连续随机变量：$$E(X)=\int_0^{\infty}P(X>x)\,dx$$，$$E(X^n)=\int_0^{\infty}nx^{n-1}P(X>x)\,dx$$

**☆ 施瓦茨不等式（P183，会证）**

$$\big[\mathrm{Cov}(X,Y)\big]^2\le\mathrm{Var}(X)\,\mathrm{Var}(Y).$$

**性质 3.4.12（P183，会证）**：$$\mathrm{Corr}(X,Y)=\pm1$$ 的充要条件是 $$X$$ 与 $$Y$$ 几乎处处有线性关系，即存在 $$a\ne0$$ 与 $$b$$，使得 $$P(Y=aX+b)=1$$。

<details markdown="1"><summary>证明</summary>

记 $$X^*=X-EX$$，$$Y^*=Y-EY$$。对任意实数 $$t$$，

$$g(t)=E(tX^*+Y^*)^2=t^2\mathrm{Var}(X)+2t\,\mathrm{Cov}(X,Y)+\mathrm{Var}(Y)\ge0.$$

这个关于 $$t$$ 的二次函数非负，所以判别式 $$\le0$$：$$4[\mathrm{Cov}(X,Y)]^2-4\mathrm{Var}(X)\mathrm{Var}(Y)\le0$$，即施瓦茨不等式。

$$\vert \mathrm{Corr}\vert =1$$ 等价于判别式 $$=0$$，即存在 $$t_0$$ 使 $$g(t_0)=E(t_0X^*+Y^*)^2=0$$，也就是 $$\mathrm{Var}(t_0X+Y)=0$$。由定理 2.3.2，$$t_0X+Y$$ 几乎处处为常数，即 $$Y=aX+b$$，其中 $$a=-t_0$$。$$a\ne0$$，否则 $$Y$$ 为常数，相关系数没有定义。反之，若 $$Y=aX+b$$，直接计算可得 $$\mathrm{Corr}=\mathrm{sgn}(a)=\pm1$$。

</details>

**$$N(\mu_1,\mu_2,\sigma_1^2,\sigma_2^2,\rho)$$ 的相关系数就是 $$\rho$$（P182，会证）**。所以对二维正态来说，**不相关与独立等价**。

**定理 3.4.2（P188，会证）**：$$n$$ 维随机向量的协方差矩阵 $$\mathrm{Cov}(\boldsymbol X)=\big(\mathrm{Cov}(X_i,X_j)\big)_{n\times n}$$ 是对称的非负定矩阵。

<details markdown="1"><summary>证明</summary>

对称性显然。对任意 $$\boldsymbol c=(c_1,\dots,c_n)^T$$，

$$\boldsymbol c^T\mathrm{Cov}(\boldsymbol X)\boldsymbol c=\sum_{i,j}c_ic_j\mathrm{Cov}(X_i,X_j)=\mathrm{Var}\left(\sum_i c_iX_i\right)\ge0.$$

</details>

**和的方差**：

$$\mathrm{Var}\left(\sum_{i=1}^nX_i\right)=\sum_{i=1}^n\mathrm{Var}(X_i)+2\sum_{i<j}\mathrm{Cov}(X_i,X_j).$$

若只有相邻项相关（$$\vert i-j\vert \ge2$$ 时协方差为 0），则简化为 $$\sum\mathrm{Var}(X_i)+2\sum_{i=1}^{n-1}\mathrm{Cov}(X_i,X_{i+1})$$（习题 4.3.4，挺常用）。

**配对模型**（抽礼物问题）：$$n$$ 个人各自随机抽回一件礼物，$$S_n$$ 为抽到自己礼物的人数，则 $$E(S_n)=1$$，$$\mathrm{Var}(S_n)=1$$，与 $$n$$ 无关。

<details markdown="1"><summary>推导</summary>

令 $$X_i=1$$ 表示第 $$i$$ 人抽到自己的礼物，则 $$P(X_i=1)=\frac1n$$，$$P(X_iX_j=1)=\frac1{n(n-1)}$$。  
$$E(S_n)=n\cdot\frac1n=1$$；  
$$\mathrm{Var}(S_n)=n\cdot\frac1n\left(1-\frac1n\right)+n(n-1)\left[\frac{1}{n(n-1)}-\frac1{n^2}\right]=1-\frac1n+\frac1n=1$$。

</details>

---

## 5. 条件分布与条件期望

**条件分布函数的推导（P197）**：连续场合下，$$P(X\le x\mid Y=y)$$ 定义为 $$\lim_{h\to0^+}P(X\le x\mid y\le Y\le y+h)$$，由此得到条件密度

$$p(x\mid y)=\frac{p(x,y)}{p_Y(y)}.$$

**二维正态的条件分布仍为正态（P198，会证）**

$$X\mid Y=y\ \sim\ N\left(\mu_1+\rho\frac{\sigma_1}{\sigma_2}(y-\mu_2),\ \sigma_1^2(1-\rho^2)\right).$$

条件均值 $$g_1(y)=E(X\mid Y=y)=\mu_1+\rho\frac{\sigma_1}{\sigma_2}(y-\mu_2)$$ 是 $$y$$ 的线性函数，也就是回归直线；条件方差不依赖于 $$y$$，且比 $$\sigma_1^2$$ 小（知道 $$Y$$ 之后，对 $$X$$ 的不确定性减少了）。（指导书 P204）

**条件期望**：$$E(X\mid Y=y)$$ 是 $$y$$ 的函数；把 $$y$$ 换成 $$Y$$，得到的 $$E(X\mid Y)$$ 是一个随机变量。

**☆ 定理 3.5.1 重期望公式（P202，会证）**

$$E(X)=E\big(E(X\mid Y)\big).$$

<details markdown="1"><summary>证明（离散型）</summary>

$$E\big(E(X\mid Y)\big)=\sum_jE(X\mid Y=y_j)P(Y=y_j)=\sum_j\sum_ix_iP(X=x_i\mid Y=y_j)P(Y=y_j)=\sum_ix_i\sum_jP(X=x_i,Y=y_j)=\sum_ix_iP(X=x_i)=E(X).$$

</details>

直观理解：先分组求平均，再对各组的平均按组的权重求平均，结果等于总平均。

**随机个随机变量之和的期望（P204，会证）**：设 $$X_1,X_2,\dots$$ 独立同分布，$$N$$ 是取正整数值的随机变量且与 $$\{X_i\}$$ 独立，则

$$E\left(\sum_{i=1}^NX_i\right)=E(X_1)\,E(N).$$

<details markdown="1"><summary>证明</summary>

$$E\left(\sum_{i=1}^NX_i\right)=E\left[E\left(\sum_{i=1}^NX_i\,\Big\vert \,N\right)\right]=\sum_nP(N=n)\,E\left(\sum_{i=1}^nX_i\right)=\sum_nP(N=n)\,nE(X_1)=E(X_1)E(N).$$

</details>

---

**本系列目录**

1. [01｜概念简答索引]({% post_url 2026-09-30-mao-notes-01-concepts %})
2. [02｜常用分布]({% post_url 2026-09-30-mao-notes-02-distributions %})
3. **03｜多维分布、协方差与条件期望**（本篇）
4. [04｜特征函数、大数定律与中心极限定理]({% post_url 2026-09-30-mao-notes-04-limit-theorems %})
5. [05｜抽样分布、次序统计量与充分统计量]({% post_url 2026-09-30-mao-notes-05-sampling %})
6. [06｜参数估计]({% post_url 2026-09-30-mao-notes-06-estimation %})
7. [07｜区间估计、假设检验与习题结论]({% post_url 2026-09-30-mao-notes-07-intervals-tests %})

[← 上一篇：常用分布]({% post_url 2026-09-30-mao-notes-02-distributions %})　·　[下一篇：特征函数、大数定律与中心极限定理 →]({% post_url 2026-09-30-mao-notes-04-limit-theorems %})

*本系列是我复习茆诗松、程依明、濮晓龙《概率论与数理统计教程》时整理的笔记，共 7 篇。文中 P××× 为茆书页码，“指导书”为配套的学习指导与习题解答。如有错漏，欢迎指出。*
