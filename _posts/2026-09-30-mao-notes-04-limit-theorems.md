---
layout: post
title: "茆书复习笔记 04｜特征函数、大数定律与中心极限定理"
author: 张梓源
tags:
- 茆书复习笔记
- 概率论
- 特征函数
- 大数定律
- 中心极限定理
- 茆诗松
- 考研
- 读书笔记
date: 2026-09-30 10:04 -0400
toc: true
---
> **导读**
>
> 大数定律和中心极限定理，最容易混淆的不是结论，而是**条件**：哪个要求独立同分布？哪个只要两两不相关？哪个连方差都不需要？
>
> 这一篇先整理特征函数这个“工具”，再把几个极限定理的条件放进表格里对比，最后给出用特征函数证明辛钦大数定律和林德伯格-莱维中心极限定理的主线。

---

## 1. 特征函数

**定义（P216）**

$$\varphi(t)=E(e^{itX})=\begin{cases}\sum_k e^{itx_k}P(X=x_k),&\text{离散}\\[4pt]\int_{-\infty}^{\infty}e^{itx}p(x)\,dx,&\text{连续}\end{cases}$$

记住欧拉公式 $$e^{i\theta}=\cos\theta+i\sin\theta$$。因为 $$\vert e^{itX}\vert =1$$，所以任何随机变量的特征函数都存在（矩母函数就不一定）。

**性质**
1. $$\vert \varphi(t)\vert \le\varphi(0)=1$$，$$\varphi(-t)=\overline{\varphi(t)}$$
2. 若 $$Y=aX+b$$，则 $$\varphi_Y(t)=e^{ibt}\varphi_X(at)$$
3. 若 $$X$$ 与 $$Y$$ 独立，则 $$\varphi_{X+Y}(t)=\varphi_X(t)\varphi_Y(t)$$
4. 若 $$E(X^k)$$ 存在，则 $$\varphi^{(k)}(0)=i^kE(X^k)$$
5. **逆转公式**：若 $$x_1<x_2$$ 是 $$F$$ 的连续点，则

   $$F(x_2)-F(x_1)=\lim_{T\to\infty}\frac1{2\pi}\int_{-T}^{T}\frac{e^{-itx_1}-e^{-itx_2}}{it}\varphi(t)\,dt$$

6. **唯一性定理**：特征函数与分布函数一一对应。若 $$\int\vert \varphi(t)\vert \,dt<\infty$$，则 $$X$$ 为连续型，且

   $$p(x)=\frac1{2\pi}\int_{-\infty}^{\infty}e^{-itx}\varphi(t)\,dt$$

**狄利克雷积分（P223，证明逆转公式时要用）**

$$D(a)=\frac1\pi\int_0^{\infty}\frac{\sin at}{t}\,dt=\begin{cases}\tfrac12,&a>0\\0,&a=0\\-\tfrac12,&a<0\end{cases}$$

**常用分布的特征函数（P220，指导书 P227）**

| 分布 | $$\varphi(t)$$ |
|---|---|
| 单点 $$P(X=a)=1$$ | $$e^{ita}$$ |
| 0-1 分布 $$b(1,p)$$ | $$pe^{it}+q$$ |
| 二项 $$b(n,p)$$ | $$(pe^{it}+q)^n$$ |
| 泊松 $$P(\lambda)$$ | $$e^{\lambda(e^{it}-1)}$$ |
| 几何 $$Ge(p)$$ | $$\dfrac{pe^{it}}{1-qe^{it}}$$ |
| 负二项 $$Nb(r,p)$$ | $$\left(\dfrac{pe^{it}}{1-qe^{it}}\right)^r$$ |
| 正态 $$N(\mu,\sigma^2)$$ | $$\exp\left\{i\mu t-\frac{\sigma^2t^2}2\right\}$$ |
| 标准正态 $$N(0,1)$$ | $$e^{-t^2/2}$$ |
| 柯西 $$Cau(\mu,\lambda)$$ | $$\exp\{i\mu t-\lambda\lvert t\rvert\}$$ |
| 均匀 $$U(a,b)$$ | $$\dfrac{e^{itb}-e^{ita}}{it(b-a)}$$ |
| 均匀 $$U(-a,a)$$ | $$\dfrac{\sin at}{at}$$ |
| 指数 $$Exp(\lambda)$$ | $$\left(1-\frac{it}{\lambda}\right)^{-1}$$ |
| 伽玛 $$Ga(\alpha,\lambda)$$ | $$\left(1-\frac{it}{\lambda}\right)^{-\alpha}$$ |
| 卡方 $$\chi^2(n)$$ | $$(1-2it)^{-n/2}$$ |
| 贝塔 $$Be(a,b)$$ | $$\dfrac{\Gamma(a+b)}{\Gamma(a)}\sum_{j=0}^{\infty}\dfrac{\Gamma(a+j)\,(it)^j}{\Gamma(a+b+j)\,\Gamma(j+1)}$$ |

<details markdown="1"><summary>推导示例：泊松与指数</summary>

泊松：$$\varphi(t)=\sum_{k=0}^{\infty}e^{itk}\frac{\lambda^k}{k!}e^{-\lambda}=e^{-\lambda}\sum_{k=0}^{\infty}\frac{(\lambda e^{it})^k}{k!}=e^{-\lambda}e^{\lambda e^{it}}$$。

指数：$$\varphi(t)=\int_0^{\infty}e^{itx}\lambda e^{-\lambda x}\,dx=\frac{\lambda}{\lambda-it}=\left(1-\frac{it}\lambda\right)^{-1}$$。

记忆技巧：伽玛是 $$\alpha$$ 个“指数”相加，所以是指数的 $$\alpha$$ 次方；卡方是 $$Ga(\frac n2,\frac12)$$，代入即得。

</details>

---

## 2. 收敛性之间的关系

**定理 4.1.1（P209）**：若 $$X_n\xrightarrow{P}a$$，$$Y_n\xrightarrow{P}b$$，则
1. $$X_n\pm Y_n\xrightarrow{P}a\pm b$$
2. $$X_nY_n\xrightarrow{P}ab$$
3. $$X_n/Y_n\xrightarrow{P}a/b$$（$$b\ne0$$）

更一般地，若 $$g$$ 在 $$a$$ 处连续，则 $$g(X_n)\xrightarrow{P}g(a)$$。

**定理 4.1.2（P212）**：$$X_n\xrightarrow{P}X\ \Rightarrow\ X_n\xrightarrow{L}X$$（反之一般不成立）。

**定理 4.1.3**：若 $$c$$ 为常数，则 $$X_n\xrightarrow{P}c\iff X_n\xrightarrow{L}c$$。

**连续性定理**：$$F_n\to F$$（弱收敛）$$\iff\varphi_n(t)\to\varphi(t)$$ 对每个 $$t$$ 成立。这是用特征函数证明极限定理的关键工具。

---

## 3. 大数定律

**统一形式**：

$$\lim_{n\to\infty}P\left(\left\vert \frac1n\sum_{i=1}^nX_i-\frac1n\sum_{i=1}^nE(X_i)\right\vert <\varepsilon\right)=1.$$

| 名称 | 页码 | 条件 | 掌握要求 |
|---|---|---|---|
| 伯努利大数定律 | P230 | $$S_n\sim b(n,p)$$，结论为 $$\frac{S_n}n\xrightarrow{P}p$$ | 会证 |
| 切比雪夫大数定律 | P233 | $$\{X_i\}$$ 两两不相关，方差存在且有共同上界 | |
| 马尔可夫大数定律 | P234 | $$\frac1{n^2}\mathrm{Var}\left(\sum_{i=1}^nX_i\right)\to0$$（马尔可夫条件） | |
| 辛钦大数定律 | P235 | $$\{X_i\}$$ 独立同分布，$$E(X_i)$$ 存在（不要求方差存在） | 会证 |

直观理解：伯努利大数定律说明“频率稳定于概率”，辛钦大数定律说明“样本均值稳定于总体均值”，这是矩估计和蒙特卡罗方法的依据。

<details markdown="1"><summary>证明：伯努利（用切比雪夫不等式）</summary>

$$E\left(\frac{S_n}{n}\right)=p$$，$$\mathrm{Var}\left(\frac{S_n}{n}\right)=\frac{p(1-p)}{n}$$，所以

$$P\left(\left\vert \frac{S_n}n-p\right\vert \ge\varepsilon\right)\le\frac{p(1-p)}{n\varepsilon^2}\le\frac{1}{4n\varepsilon^2}\to0.$$

切比雪夫、马尔可夫大数定律的证明完全相同，只是方差的上界换一下。

</details>

<details markdown="1"><summary>证明：辛钦（用特征函数）</summary>

设 $$E(X_i)=a$$，$$\varphi(t)$$ 为 $$X_i$$ 的特征函数。因为 $$\varphi'(0)=ia$$，所以 $$\varphi(t)=1+iat+o(t)$$。  
$$Y_n=\frac1n\sum X_i$$ 的特征函数为

$$\varphi_{Y_n}(t)=\left[\varphi\left(\frac tn\right)\right]^n=\left[1+\frac{iat}{n}+o\left(\frac1n\right)\right]^n\to e^{iat}.$$

$$e^{iat}$$ 是单点分布 $$P(X=a)=1$$ 的特征函数，所以 $$Y_n\xrightarrow{L}a$$，再由定理 4.1.3 得 $$Y_n\xrightarrow{P}a$$。

</details>

**习题 4.1.18**：$$\{X_n\}$$ 独立同分布，$$E(X_n)=0$$，$$\mathrm{Var}(X_n)=\sigma^2$$，则 $$\frac1n\sum_{k=1}^nX_k^2\xrightarrow{P}\sigma^2$$。（对 $$X_k^2$$ 用辛钦大数定律，$$E(X_k^2)=\sigma^2$$。）

---

## 4. 中心极限定理

| 名称 | 页码 | 条件 | 结论 |
|---|---|---|---|
| 林德伯格-莱维（会证） | P240 | $$\{X_n\}$$ 独立同分布，$$E(X_i)=\mu$$，$$\mathrm{Var}(X_i)=\sigma^2>0$$ 存在 | $$\dfrac{\sum X_i-n\mu}{\sigma\sqrt n}\xrightarrow{L}N(0,1)$$ |
| 棣莫弗-拉普拉斯 | P242 | $$S_n\sim b(n,p)$$（伯努利） | $$\dfrac{S_n-np}{\sqrt{npq}}\xrightarrow{L}N(0,1)$$ |

<details markdown="1"><summary>证明：林德伯格-莱维（特征函数法）</summary>

令 $$Y_i=\frac{X_i-\mu}{\sigma}$$，其特征函数 $$\varphi(t)=1-\frac{t^2}{2}+o(t^2)$$（因为 $$E Y_i=0$$，$$E Y_i^2=1$$）。  
$$Z_n=\frac{1}{\sqrt n}\sum Y_i$$ 的特征函数为

$$\left[\varphi\left(\frac t{\sqrt n}\right)\right]^n=\left[1-\frac{t^2}{2n}+o\left(\frac1n\right)\right]^n\to e^{-t^2/2},$$

这正是 $$N(0,1)$$ 的特征函数，由连续性定理即得结论。

</details>

**应用：二项分布的正态近似（P243）**

$$P(a\le S_n\le b)\approx\Phi\left(\frac{b+0.5-np}{\sqrt{npq}}\right)-\Phi\left(\frac{a-0.5-np}{\sqrt{npq}}\right)$$

（加减 0.5 是连续性修正。）

直观理解：不管单个 $$X_i$$ 长什么样，只要独立同分布且方差有限，大量叠加之后就会“磨成”钟形。这也是正态分布在自然界中无处不在的原因。

---

**本系列目录**

1. [01｜概念简答索引]({% post_url 2026-09-30-mao-notes-01-concepts %})
2. [02｜常用分布]({% post_url 2026-09-30-mao-notes-02-distributions %})
3. [03｜多维分布、协方差与条件期望]({% post_url 2026-09-30-mao-notes-03-multivariate %})
4. **04｜特征函数、大数定律与中心极限定理**（本篇）
5. [05｜抽样分布、次序统计量与充分统计量]({% post_url 2026-09-30-mao-notes-05-sampling %})
6. [06｜参数估计]({% post_url 2026-09-30-mao-notes-06-estimation %})
7. [07｜区间估计、假设检验与习题结论]({% post_url 2026-09-30-mao-notes-07-intervals-tests %})

[← 上一篇：多维分布、协方差与条件期望]({% post_url 2026-09-30-mao-notes-03-multivariate %})　·　[下一篇：抽样分布、次序统计量与充分统计量 →]({% post_url 2026-09-30-mao-notes-05-sampling %})

*本系列是我复习茆诗松、程依明、濮晓龙《概率论与数理统计教程》时整理的笔记，共 7 篇。文中 P××× 为茆书页码，“指导书”为配套的学习指导与习题解答。如有错漏，欢迎指出。*
