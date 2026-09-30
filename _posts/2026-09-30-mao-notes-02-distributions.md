---
layout: post
title: "茆书复习笔记 02｜常用分布：特征、无记忆性、可加性与分布间关系"
author: 张梓源
tags:
- 茆书复习笔记
- 概率论
- 常用分布
- 茆诗松
- 考研
- 读书笔记
date: 2026-09-30 10:02 -0400
toc: true
---
> **导读**
>
> 常用分布是整本概率论的地基：后面的特征函数、抽样分布、参数估计，处处都要用到它们的期望、方差和相互关系。
>
> 这一篇的要求是“**会推 + 记住**”。建议先看第 5 节的关系网：只要记住 $$Ga(1,\lambda)=Exp(\lambda)$$、$$Ga(\frac n2,\frac12)=\chi^2(n)$$、$$Be(1,1)=U(0,1)$$ 这几个等式，很多分布就能互相推出来，不用死记。每个分布都要能写出分布列或密度、期望和方差，标 ☆ 的是重点。

---

## 1. 17 个常用分布一览

记号：$$q=1-p$$。

| # | 分布 | 分布列 / 密度 | $$E(X)$$ | $$\mathrm{Var}(X)$$ |
|---|---|---|---|---|
| 1 | 二项 $$b(n,p)$$ | $$\binom nk p^kq^{n-k},\ k=0,\dots,n$$ | $$np$$ | $$npq$$ |
| 2 | 二点/0-1/伯努利 $$b(1,p)$$ | $$p^xq^{1-x},\ x=0,1$$ | $$p$$ | $$pq$$ |
| 3 | 泊松 $$P(\lambda)$$ | $$\frac{\lambda^k}{k!}e^{-\lambda},\ k=0,1,\dots$$ | $$\lambda$$ | $$\lambda$$ |
| 4 | 超几何 $$h(n,N,M)$$ | $$\dfrac{\binom Mk\binom{N-M}{n-k}}{\binom Nn}$$ | $$n\frac MN$$ | $$\dfrac{nM(N-M)(N-n)}{N^2(N-1)}$$ |
| 5 | 几何 $$Ge(p)$$ | $$q^{k-1}p,\ k=1,2,\dots$$ | $$\frac1p$$ | $$\frac{q}{p^2}$$ |
| 6 | 负二项/帕斯卡 $$Nb(r,p)$$ | $$\binom{k-1}{r-1}p^rq^{k-r},\ k=r,r+1,\dots$$ | $$\frac rp$$ | $$\frac{rq}{p^2}$$ |
| 7 | 正态 $$N(\mu,\sigma^2)$$ | $$\frac1{\sqrt{2\pi}\sigma}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$$ | $$\mu$$ | $$\sigma^2$$ |
| 8 | 均匀 $$U(a,b)$$ | $$\frac1{b-a},\ a<x<b$$ | $$\frac{a+b}2$$ | $$\frac{(b-a)^2}{12}$$ |
| 9 | 指数 $$Exp(\lambda)$$ | $$\lambda e^{-\lambda x},\ x>0$$ | $$\frac1\lambda$$ | $$\frac1{\lambda^2}$$ |
| 10 | ☆ 伽玛 $$Ga(\alpha,\lambda)$$ | $$\frac{\lambda^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\lambda x},\ x>0$$ | $$\frac\alpha\lambda$$ | $$\frac\alpha{\lambda^2}$$ |
| 11 | ☆ 贝塔 $$Be(a,b)$$ | $$\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}x^{a-1}(1-x)^{b-1},\ 0<x<1$$ | $$\frac a{a+b}$$ | $$\frac{ab}{(a+b)^2(a+b+1)}$$ |
| 12 | 对数正态 $$LN(\mu,\sigma^2)$$ | $$\frac1{x\sigma\sqrt{2\pi}}e^{-\frac{(\ln x-\mu)^2}{2\sigma^2}},\ x>0$$ | $$e^{\mu+\sigma^2/2}$$ | $$(e^{\sigma^2}-1)e^{2\mu+\sigma^2}$$ |
| 13 | 柯西 $$Cau(\mu,\lambda)$$ | $$\frac1\pi\cdot\frac{\lambda}{\lambda^2+(x-\mu)^2}$$ | 不存在 | 不存在 |
| 14 | 卡方 $$\chi^2(n)$$（P283） | 即 $$Ga(\frac n2,\frac12)$$ | $$n$$ | $$2n$$ |
| 15 | $$F(m,n)$$（P286） | 见第 5 篇 | $$\frac n{n-2}\ (n>2)$$ | $$\frac{2n^2(m+n-2)}{m(n-2)^2(n-4)}\ (n>4)$$ |
| 16 | $$t(n)$$（P288） | 见第 5 篇 | $$0\ (n>1)$$ | $$\frac n{n-2}\ (n>2)$$ |
| 17 | 韦布尔 | $$F(x)=1-e^{-(x/\eta)^m},\ x>0$$ | $$\eta\Gamma(1+\frac1m)$$ | $$\eta^2\left[\Gamma(1+\frac2m)-\Gamma^2(1+\frac1m)\right]$$ |

补充结论：
- **对数正态**（习题 2.7.3、2.7.11）：中位数为 $$e^\mu$$，$$p$$ 分位数为 $$x_p=\exp\{\mu+\sigma u_p\}$$（因为 $$\ln X\sim N(\mu,\sigma^2)$$，而取对数是单调变换，分位数跟着一起变）。
- **伽玛**（习题 2.7.4）：$$E(X^k)=\dfrac{\Gamma(\alpha+k)}{\Gamma(\alpha)\lambda^k}$$。
- **指数**（习题 2.7.5）：$$E(X^k)=\dfrac{k!}{\lambda^k}$$，变异系数 $$C_v=1$$，偏度 $$\beta_s=2$$，峰度 $$\beta_k=6$$。

<details markdown="1"><summary>推导示例：对数正态的期望与方差</summary>

令 $$Y=\ln X\sim N(\mu,\sigma^2)$$，则 $$E(X^k)=E(e^{kY})$$ 就是正态的矩母函数在 $$k$$ 处的值：

$$E(e^{kY})=e^{k\mu+\frac12k^2\sigma^2}.$$

取 $$k=1$$ 得 $$E(X)=e^{\mu+\sigma^2/2}$$；取 $$k=2$$ 得 $$E(X^2)=e^{2\mu+2\sigma^2}$$，于是

$$\mathrm{Var}(X)=e^{2\mu+2\sigma^2}-e^{2\mu+\sigma^2}=(e^{\sigma^2}-1)e^{2\mu+\sigma^2}.$$

</details>

---

## 2. 伽玛函数与贝塔函数

**伽玛函数（P115）**：$$\Gamma(\alpha)=\int_0^{\infty}x^{\alpha-1}e^{-x}\,dx$$，$$\alpha>0$$。
- $$\Gamma(1)=1$$，$$\Gamma(\tfrac12)=\sqrt\pi$$
- $$\Gamma(\alpha+1)=\alpha\Gamma(\alpha)$$（分部积分），所以 $$\Gamma(n+1)=n!$$

**伽玛分布的参数含义**：$$\alpha$$ 是形状参数，$$\lambda$$ 是尺度参数（准确说 $$1/\lambda$$ 才是尺度）。$$\alpha\le1$$ 时密度单调递减，$$\alpha>1$$ 时密度先增后减、呈单峰。

**贝塔函数（P117）**：$$B(a,b)=\int_0^1x^{a-1}(1-x)^{b-1}\,dx$$，$$a,b>0$$。
- 对称性：$$B(a,b)=B(b,a)$$
- 与伽玛函数的关系：$$B(a,b)=\dfrac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}$$

**贝塔分布的形状**：$$a=b$$ 时关于 $$\frac12$$ 对称；$$a=b=1$$ 就是 $$U(0,1)$$；$$a<b$$ 时右偏，$$a>b$$ 时左偏。

**三个常用积分**：

$$\int_{-\infty}^{\infty}e^{-u^2/2}\,du=\sqrt{2\pi},\qquad \int_{-\infty}^{\infty}e^{-x^2}\,dx=\sqrt\pi,\qquad \int_0^{\infty}x^2e^{-x^2}\,dx=\frac{\sqrt\pi}4.$$

---

## 3. 无记忆性（都要会证）

**几何分布（P103）**：$$X\sim Ge(p)$$，则对任意正整数 $$m,n$$，

$$P(X>m+n\mid X>m)=P(X>n).$$

**指数分布（P114）**：$$X\sim Exp(\lambda)$$，则对任意 $$s,t>0$$，

$$P(X>s+t\mid X>s)=P(X>t).$$

<details markdown="1"><summary>证明</summary>

几何分布：$$P(X>n)=\sum_{k=n+1}^{\infty}q^{k-1}p=q^n$$，所以

$$P(X>m+n\mid X>m)=\frac{P(X>m+n)}{P(X>m)}=\frac{q^{m+n}}{q^m}=q^n=P(X>n).$$

指数分布：$$P(X>t)=e^{-\lambda t}$$，所以

$$P(X>s+t\mid X>s)=\frac{e^{-\lambda(s+t)}}{e^{-\lambda s}}=e^{-\lambda t}=P(X>t).$$

</details>

直观理解：“已经等了 $$s$$ 分钟”这件事，不会改变“还要再等多久”的分布。几何分布是唯一具有无记忆性的离散分布，指数分布是唯一具有无记忆性的连续分布。

---

## 4. 可加性（都要会推导）

下面各式中的随机变量**相互独立**：

| 分布 | 可加性 | 页码 |
|---|---|---|
| 二项 | $$b(n,p)+b(m,p)=b(n+m,p)$$（$$p$$ 相同） | P163 |
| 泊松 | $$P(\lambda_1)+P(\lambda_2)=P(\lambda_1+\lambda_2)$$ | P163 |
| 正态 | $$N(\mu_1,\sigma_1^2)+N(\mu_2,\sigma_2^2)=N(\mu_1+\mu_2,\sigma_1^2+\sigma_2^2)$$ | P167 |
| 伽玛 | $$Ga(\alpha_1,\lambda)+Ga(\alpha_2,\lambda)=Ga(\alpha_1+\alpha_2,\lambda)$$（$$\lambda$$ 相同） | P168 |
| 卡方 | $$\chi^2(n)+\chi^2(m)=\chi^2(n+m)$$ | P169 |

> **注意：指数分布没有可加性。** $$Exp(\lambda)+Exp(\lambda)=Ga(2,\lambda)$$，已经不是指数分布了。

<details markdown="1"><summary>推导方法（以泊松、伽玛为例）</summary>

**方法一：卷积公式。** 泊松：

$$P(Z=k)=\sum_{i=0}^{k}\frac{\lambda_1^i}{i!}e^{-\lambda_1}\frac{\lambda_2^{k-i}}{(k-i)!}e^{-\lambda_2}=\frac{e^{-(\lambda_1+\lambda_2)}}{k!}\sum_{i=0}^k\binom ki\lambda_1^i\lambda_2^{k-i}=\frac{(\lambda_1+\lambda_2)^k}{k!}e^{-(\lambda_1+\lambda_2)}.$$

**方法二：特征函数**（更快）。独立和的特征函数等于特征函数之积：
- 泊松：$$e^{\lambda_1(e^{it}-1)}\cdot e^{\lambda_2(e^{it}-1)}=e^{(\lambda_1+\lambda_2)(e^{it}-1)}$$
- 伽玛：$$(1-\frac{it}\lambda)^{-\alpha_1}(1-\frac{it}\lambda)^{-\alpha_2}=(1-\frac{it}\lambda)^{-(\alpha_1+\alpha_2)}$$

再由特征函数的唯一性定理得到结论。

</details>

---

## 5. 分布间关系（会写最好）

| 关系 | 页码 |
|---|---|
| $$r$$ 个独立同分布的 $$Ge(p)$$ 之和服从 $$Nb(r,p)$$ | P104 |
| 泊松过程中，单位时间事件数 $$\sim P(\lambda)$$ 时，首次发生的等待时间 $$\sim Exp(\lambda)$$ | P114 |
| $$Ga(1,\lambda)=Exp(\lambda)$$ | P116 |
| $$Ga(\tfrac n2,\tfrac12)=\chi^2(n)$$ | P116 |
| 第 $$n$$ 个事件发生的时间 $$\sim Ga(n,\lambda)$$，与泊松分布的关系如下 | P116 |
| $$Be(1,1)=U(0,1)$$ | P118 |
| $$X\sim Ga(\alpha,\lambda)\Rightarrow kX\sim Ga(\alpha,\lambda/k)$$，$$k>0$$ | P125 |
| $$X\sim N(0,1)\Rightarrow X^2\sim\chi^2(1)$$ | P127 |
| $$t(1)=Cau(0,1)$$ | |
| $$X\sim t(n)\Rightarrow X^2\sim F(1,n)$$ | |

**泊松与指数、伽玛的关系**：设 $$N(t)$$ 为 $$[0,t]$$ 内事件发生的次数，$$N(t)\sim P(\lambda t)$$。
- 首次发生时间 $$T_1$$：$$P(T_1>t)=P(N(t)=0)=e^{-\lambda t}$$，所以 $$T_1\sim Exp(\lambda)$$。
- 第 $$n$$ 次发生时间 $$T_n$$：$$P(T_n>t)=P(N(t)<n)=\sum_{k=0}^{n-1}\frac{(\lambda t)^k}{k!}e^{-\lambda t}$$，所以 $$T_n\sim Ga(n,\lambda)$$。

**$$Y=F_X(X)\sim U(0,1)$$（P126，会证）**：设 $$X$$ 的分布函数 $$F_X$$ 连续且严格单增，则 $$Y=F_X(X)\sim U(0,1)$$。它是逆变换法生成随机数的依据。

<details markdown="1"><summary>证明</summary>

当 $$0<y<1$$ 时，

$$F_Y(y)=P(F_X(X)\le y)=P\big(X\le F_X^{-1}(y)\big)=F_X\big(F_X^{-1}(y)\big)=y;$$

当 $$y\le0$$ 时 $$F_Y(y)=0$$，当 $$y\ge1$$ 时 $$F_Y(y)=1$$。所以 $$Y\sim U(0,1)$$。  
反过来，若 $$U\sim U(0,1)$$，则 $$X=F^{-1}(U)$$ 的分布函数就是 $$F$$。

</details>

**正态的线性变换仍为正态（P124）**：$$X\sim N(\mu,\sigma^2)\Rightarrow aX+b\sim N(a\mu+b,a^2\sigma^2)$$，$$a\ne0$$。

---

## 6. 分布形状与近似

**形状**
- 二项分布：随着 $$p$$ 增大，峰逐渐右移（P95）。$$p<0.5$$ 时右偏，$$p=0.5$$ 时对称，$$p>0.5$$ 时左偏。
- 泊松分布：随着 $$\lambda$$ 增大，逐渐趋于对称（由右偏变为对称）。

**近似**

| 近似 | 页码 | 条件与用法 |
|---|---|---|
| 二项 → 泊松（泊松定理，定理 2.4.1） | P98 | $$n$$ 大、$$p$$ 小、$$np=\lambda$$ 适中时，$$\binom nkp^k(1-p)^{n-k}\approx\frac{\lambda^k}{k!}e^{-\lambda}$$ |
| 超几何 → 二项 | P102 | $$N$$ 很大、$$n$$ 相对 $$N$$ 很小时，不放回近似于有放回，$$h(n,N,M)\approx b(n,\frac MN)$$（知道结论即可） |
| 二项 → 正态（棣莫弗-拉普拉斯） | P243 | $$n$$ 大时，$$P(a\le X\le b)\approx\Phi\left(\frac{b+0.5-np}{\sqrt{npq}}\right)-\Phi\left(\frac{a-0.5-np}{\sqrt{npq}}\right)$$（修正项 0.5 可提高精度） |

<details markdown="1"><summary>泊松定理的证明（要会）</summary>

**定理**：若 $$\lim_{n\to\infty}np_n=\lambda>0$$，则对任意固定的 $$k$$，$$\lim_{n\to\infty}\binom nkp_n^k(1-p_n)^{n-k}=\frac{\lambda^k}{k!}e^{-\lambda}$$。

记 $$\lambda_n=np_n$$，则

$$\binom nkp_n^k(1-p_n)^{n-k}=\frac{\lambda_n^k}{k!}\cdot\underbrace{\frac{n(n-1)\cdots(n-k+1)}{n^k}}_{\to1}\cdot\underbrace{\left(1-\frac{\lambda_n}n\right)^{n}}_{\to e^{-\lambda}}\cdot\underbrace{\left(1-\frac{\lambda_n}n\right)^{-k}}_{\to1}\to\frac{\lambda^k}{k!}e^{-\lambda}.$$

</details>

---

## 7. 矩与偏度、峰度

**前四阶中心矩与原点矩的关系（P120）**（$$\mu_k$$ 为原点矩，$$\nu_k$$ 为中心矩）

$$
\begin{aligned}
\nu_2&=\mu_2-\mu_1^2\\
\nu_3&=\mu_3-3\mu_2\mu_1+2\mu_1^3\\
\nu_4&=\mu_4-4\mu_3\mu_1+6\mu_2\mu_1^2-3\mu_1^4
\end{aligned}
$$

**正态分布的矩（P120，P130 例 2.7.1）**：$$X\sim N(0,\sigma^2)$$ 时，$$k$$ 为奇数则 $$\mu_k=0$$，$$k$$ 为偶数则 $$\mu_k=\sigma^k(k-1)!!$$。前四阶为 $$0,\ \sigma^2,\ 0,\ 3\sigma^4$$。

**几种常见分布的偏度与峰度（P137）**

| 分布 | $$E(X)$$ | $$\mathrm{Var}(X)$$ | 偏度 $$\beta_s$$ | 峰度 $$\beta_k$$ |
|---|---|---|---|---|
| $$U(a,b)$$ | $$\frac{a+b}{2}$$ | $$\frac{(b-a)^2}{12}$$ | 0 | −1.2 |
| $$N(\mu,\sigma^2)$$ | $$\mu$$ | $$\sigma^2$$ | 0 | 0 |
| $$Exp(\lambda)$$ | $$\frac1\lambda$$ | $$\frac1{\lambda^2}$$ | 2 | 6 |
| $$Ga(\alpha,\lambda)$$ | $$\frac\alpha\lambda$$ | $$\frac{\alpha}{\lambda^2}$$ | $$\frac{2}{\sqrt\alpha}$$ | $$\frac6\alpha$$ |

伽玛分布随 $$\alpha\to\infty$$ 时偏度和峰度都趋于 0，也就是越来越接近正态。

**正态分布的 3σ 原则（P111）**

$$P(\vert X-\mu\vert <k\sigma)=\begin{cases}0.6826,&k=1\\0.9545,&k=2\\0.9973,&k=3\end{cases}$$

**应用：过程能力指数** $$C_p=\dfrac{\text{上规格限}-\text{下规格限}}{6\sigma}$$，用来判断生产过程是否稳定；$$C_p\ge1.33$$ 时认为生产过程正常。（$$C_{pk}$$ 进一步考虑了均值偏离中心的情况。）

---

## 8. 其他相关结论

- 一个连续分布的密度函数不唯一（P74）：改变有限个点处的值不影响积分。
- 方差存在时期望必存在（P87），因为 $$\vert x\vert \le1+x^2$$。
- 有界随机变量的期望与方差总存在：若 $$X\in[a,b]$$，则 $$a\le E(X)\le b$$，$$\mathrm{Var}(X)\le\left(\frac{b-a}2\right)^2$$。

---

**本系列目录**

1. [01｜概念简答索引]({% post_url 2026-09-30-mao-notes-01-concepts %})
2. **02｜常用分布**（本篇）
3. [03｜多维分布、协方差与条件期望]({% post_url 2026-09-30-mao-notes-03-multivariate %})
4. [04｜特征函数、大数定律与中心极限定理]({% post_url 2026-09-30-mao-notes-04-limit-theorems %})
5. [05｜抽样分布、次序统计量与充分统计量]({% post_url 2026-09-30-mao-notes-05-sampling %})
6. [06｜参数估计]({% post_url 2026-09-30-mao-notes-06-estimation %})
7. [07｜区间估计、假设检验与习题结论]({% post_url 2026-09-30-mao-notes-07-intervals-tests %})

[← 上一篇：概念简答索引]({% post_url 2026-09-30-mao-notes-01-concepts %})　·　[下一篇：多维分布、协方差与条件期望 →]({% post_url 2026-09-30-mao-notes-03-multivariate %})

*本系列是我复习茆诗松、程依明、濮晓龙《概率论与数理统计教程》时整理的笔记，共 7 篇。文中 P××× 为茆书页码，“指导书”为配套的学习指导与习题解答。如有错漏，欢迎指出。*
