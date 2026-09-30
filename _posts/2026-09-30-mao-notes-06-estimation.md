---
layout: post
title: "茆书复习笔记 06｜参数估计：无偏、相合、MLE、Fisher 信息与 C-R 下界"
author: 张梓源
tags:
- 茆书复习笔记
- 数理统计
- 参数估计
- 最大似然估计
- UMVUE
- 茆诗松
- 考研
- 读书笔记
date: 2026-09-30 10:06 -0400
toc: true
math: true
---
> **导读**
>
> 参数估计这一章的性质特别容易混：无偏性有没有不变性？矩估计唯一吗？MLE 呢？
>
> 这一篇先用一张表把这些“唯一 / 不变”的结论对照清楚，再整理相合性、UMVUE、充分性原则和 Cramér-Rao 不等式这几个需要会证的定理，最后用一个共轭先验的例子讲贝叶斯估计。

---

## 1. 两种点估计方法

**矩估计（替换原理）**：用样本矩替换总体矩，用样本矩的函数替换总体矩的同一函数。  
例：$$E(X)=\mu\Rightarrow\hat\mu=\bar x$$；$$\mathrm{Var}(X)=\sigma^2\Rightarrow\hat\sigma^2=s_n^2$$。

**最大似然估计**：$$L(\theta)=\prod_{i=1}^np(x_i;\theta)$$，通常解对数似然方程 $$\frac{\partial\ln L}{\partial\theta}=0$$。
- 正态总体（P316）：$$\hat\mu=\bar x$$，$$\hat\sigma^2=s_n^2=\frac1n\sum(x_i-\bar x)^2$$
- 支撑依赖参数时不能求导，要直接看 $$L$$ 的单调性。例：$$U(0,\theta)$$ 时 $$L(\theta)=\theta^{-n}I_{\{\theta\ge x_{(n)}\}}$$，在 $$\theta\ge x_{(n)}$$ 上递减，所以 $$\hat\theta=x_{(n)}$$。

---

## 2. 性质对照（哪些“唯一”、哪些“不变”）

| 结论 | 页码 | 说明 / 例子 |
|---|---|---|
| **无偏性不具有不变性** | P303 | $$s^2$$ 是 $$\sigma^2$$ 的无偏估计，但 $$s$$ 不是 $$\sigma$$ 的无偏估计（由 Jensen 不等式，$$E s<\sigma$$） |
| 矩估计：$$g(\theta)$$ 的矩估计可取 $$g(\hat\theta)$$ | P309 | 替换原理的直接推论 |
| **矩估计不唯一** | P309 | 泊松总体中 $$\bar x$$ 和 $$s_n^2$$ 都是 $$\lambda$$ 的矩估计（因为 $$E X=\mathrm{Var}X=\lambda$$） |
| **相合估计不唯一** | P310 | 若 $$\hat\theta_n$$ 相合，则 $$\hat\theta_n+\frac1n$$ 也相合 |
| **相合估计具有不变性**（定理 6.2.2，会证） | P311 | $$\hat\theta_n\xrightarrow{P}\theta$$，$$g$$ 连续 $$\Rightarrow g(\hat\theta_n)\xrightarrow{P}g(\theta)$$ |
| **矩估计一般都具有相合性** | P312 | 辛钦大数定律 + 定理 6.2.2 |
| **最大似然估计具有不变性** | P316 | $$\hat\theta$$ 是 MLE $$\Rightarrow g(\hat\theta)$$ 是 $$g(\theta)$$ 的 MLE。例：$$\sigma$$ 的 MLE 是 $$s_n$$ |
| **MLE 可能不止一个**（习题 6.3.2） | | $$U(\theta-\frac12,\theta+\frac12)$$ 中，$$[x_{(n)}-\frac12,\ x_{(1)}+\frac12]$$ 内的任一点都是 MLE |

---

## 3. 相合性

**定理 6.2.1（P310，会证）**：若 $$\lim_{n\to\infty}E(\hat\theta_n)=\theta$$，且 $$\lim_{n\to\infty}\mathrm{Var}(\hat\theta_n)=0$$，则 $$\hat\theta_n$$ 是 $$\theta$$ 的相合估计。

<details markdown="1"><summary>证明</summary>

由马尔可夫（切比雪夫型）不等式，

$$P(\vert \hat\theta_n-\theta\vert \ge\varepsilon)\le\frac{E(\hat\theta_n-\theta)^2}{\varepsilon^2}=\frac{\mathrm{Var}(\hat\theta_n)+\big(E\hat\theta_n-\theta\big)^2}{\varepsilon^2}\to0.$$

</details>

**定理 6.2.2（P311，会证）**：若 $$\hat\theta_{n1},\dots,\hat\theta_{nk}$$ 分别是 $$\theta_1,\dots,\theta_k$$ 的相合估计，$$\eta=g(\theta_1,\dots,\theta_k)$$ 是连续函数，则 $$\hat\eta_n=g(\hat\theta_{n1},\dots,\hat\theta_{nk})$$ 是 $$\eta$$ 的相合估计。

<details markdown="1"><summary>证明要点（一维）</summary>

由 $$g$$ 在 $$\theta$$ 处连续，对任意 $$\varepsilon>0$$，存在 $$\delta>0$$，使得 $$\vert \hat\theta_n-\theta\vert <\delta$$ 时 $$\vert g(\hat\theta_n)-g(\theta)\vert <\varepsilon$$。于是

$$P\big(\vert g(\hat\theta_n)-g(\theta)\vert \ge\varepsilon\big)\le P(\vert \hat\theta_n-\theta\vert \ge\delta)\to0.$$

</details>

**MLE 的渐近正态性（指导书 P328）**：在一定正则条件下，

$$\hat\theta_{MLE}\ \dot\sim\ N\left(\theta,\ \frac{1}{nI(\theta)}\right),$$

其中 $$n$$ 为样本容量，$$I(\theta)$$ 为费希尔信息量。所以 MLE 同时具有**相合性**和**渐近正态性**，而且渐近方差达到了 C-R 下界（渐近有效）。

---

## 4. 均方误差与有效性

**均方误差（P323）**

$$\mathrm{MSE}(\hat\theta)=E(\hat\theta-\theta)^2=\mathrm{Var}(\hat\theta)+\big(E\hat\theta-\theta\big)^2.$$

无偏估计的 MSE 就是方差。有偏估计的 MSE 可能反而更小。例：正态总体中，$$\sigma^2$$ 的估计 $$\frac1{n+1}\sum(x_i-\bar x)^2$$ 的 MSE 比 $$s^2$$ 更小。

---

## 5. UMVUE

**定理 6.4.1（P325，会证）**：设 $$\hat\theta$$ 是 $$\theta$$ 的无偏估计，$$\mathrm{Var}(\hat\theta)<\infty$$。则 $$\hat\theta$$ 是 $$\theta$$ 的 UMVUE，当且仅当对任意满足 $$E(\varphi)=0$$、$$\mathrm{Var}(\varphi)<\infty$$ 的统计量 $$\varphi$$（即“零的无偏估计”），都有

$$\mathrm{Cov}_\theta(\hat\theta,\varphi)=0,\qquad\forall\theta\in\Theta.$$

一句话：UMVUE 必与任一零的无偏估计不相关。

<details markdown="1"><summary>证明</summary>

**充分性**：设 $$\tilde\theta$$ 是任一无偏估计，令 $$\varphi=\tilde\theta-\hat\theta$$，则 $$E\varphi=0$$，于是

$$\mathrm{Var}(\tilde\theta)=\mathrm{Var}(\hat\theta+\varphi)=\mathrm{Var}(\hat\theta)+\mathrm{Var}(\varphi)+2\underbrace{\mathrm{Cov}(\hat\theta,\varphi)}_{=0}\ge\mathrm{Var}(\hat\theta).$$

**必要性**（反证）：若对某个 $$\theta_0$$ 有 $$\mathrm{Cov}(\hat\theta,\varphi)\ne0$$，考虑无偏估计 $$\hat\theta+\lambda\varphi$$，

$$\mathrm{Var}(\hat\theta+\lambda\varphi)=\mathrm{Var}(\hat\theta)+2\lambda\mathrm{Cov}(\hat\theta,\varphi)+\lambda^2\mathrm{Var}(\varphi).$$

取 $$\lambda=-\mathrm{Cov}/\mathrm{Var}(\varphi)$$，方差严格小于 $$\mathrm{Var}(\hat\theta)$$，与 UMVUE 矛盾。

</details>

**定理 6.4.2 充分性原则（P327，会证）**：设 $$T$$ 是充分统计量，$$\hat\theta$$ 是 $$\theta$$ 的无偏估计，令 $$\tilde\theta=E(\hat\theta\mid T)$$，则 $$\tilde\theta$$ 也是 $$\theta$$ 的无偏估计，并且 $$\mathrm{Var}(\tilde\theta)\le\mathrm{Var}(\hat\theta)$$。

<details markdown="1"><summary>证明</summary>

因为 $$T$$ 充分，条件分布与 $$\theta$$ 无关，所以 $$E(\hat\theta\mid T)$$ 不含 $$\theta$$，是一个统计量。
- 无偏性：由重期望公式，$$E\tilde\theta=E\big(E(\hat\theta\mid T)\big)=E\hat\theta=\theta$$。
- 方差：由条件方差公式，$$\mathrm{Var}(\hat\theta)=\mathrm{Var}\big(E(\hat\theta\mid T)\big)+E\big(\mathrm{Var}(\hat\theta\mid T)\big)\ge\mathrm{Var}(\tilde\theta)$$。

</details>

直观理解：好的估计只应依赖于充分统计量，对充分统计量取条件期望只会让估计更好（至少不会更差）。所以找 UMVUE，只需在充分统计量的函数中找。

---

## 6. 费希尔信息量与 C-R 不等式

**费希尔信息量（P329）**

$$I(\theta)=E\left[\frac{\partial}{\partial\theta}\ln p(x;\theta)\right]^2=-E\left[\frac{\partial^2}{\partial\theta^2}\ln p(x;\theta)\right].$$

第二个等号见习题 6.4.5，计算时通常用二阶导更方便。$$I(\theta)$$ 越大，一个样本中关于 $$\theta$$ 的信息越多。

**☆ 定理 6.4.3 Cramér-Rao 不等式（P329，会证）**：在正则条件下，若 $$T$$ 是 $$g(\theta)$$ 的无偏估计，则

$$\mathrm{Var}(T)\ge\frac{[g'(\theta)]^2}{nI(\theta)}.$$

特别地，$$\theta$$ 的任一无偏估计的方差 $$\ge\dfrac1{nI(\theta)}$$。达到下界的无偏估计称为**有效估计**。

<details markdown="1"><summary>证明要点</summary>

记得分函数 $$S=\sum_{i=1}^n\frac{\partial}{\partial\theta}\ln p(x_i;\theta)$$，则 $$E(S)=0$$，$$\mathrm{Var}(S)=nI(\theta)$$。  
对 $$E_\theta(T)=g(\theta)$$ 两边关于 $$\theta$$ 求导（正则条件保证求导与积分可交换）：

$$g'(\theta)=\int T(\boldsymbol x)\frac{\partial}{\partial\theta}p(\boldsymbol x;\theta)\,d\boldsymbol x=E(TS)=\mathrm{Cov}(T,S).$$

由施瓦茨不等式，$$[g'(\theta)]^2=[\mathrm{Cov}(T,S)]^2\le\mathrm{Var}(T)\,\mathrm{Var}(S)=\mathrm{Var}(T)\cdot nI(\theta)$$。

</details>

**常见分布的费希尔信息量（指导书 P341）**

| 分布 | 费希尔信息量 |
|---|---|
| $$b(1,p)$$ | $$I(p)=\dfrac1{p(1-p)}$$ |
| $$P(\lambda)$$ | $$I(\lambda)=\dfrac1\lambda$$ |
| $$Exp(\lambda)$$，密度 $$\lambda e^{-\lambda x}$$ | $$I(\lambda)=\dfrac1{\lambda^2}$$ |
| $$N(\mu,1)$$ | $$I(\mu)=1$$ |
| $$N(0,\sigma^2)$$ | $$I(\sigma^2)=\dfrac1{2\sigma^4}$$ |
| $$N(\mu,\sigma^2)$$ | $$I(\mu,\sigma^2)=\begin{pmatrix}\frac1{\sigma^2}&0\\0&\frac1{2\sigma^4}\end{pmatrix}$$ |

例：泊松总体中 $$\mathrm{Var}(\bar x)=\frac\lambda n=\frac{1}{nI(\lambda)}$$，所以 $$\bar x$$ 是 $$\lambda$$ 的有效估计。

---

## 7. 贝叶斯估计

**后验分布公式**

$$\pi(\theta\mid x_1,\dots,x_n)=\frac{h(x_1,\dots,x_n,\theta)}{m(x_1,\dots,x_n)}=\frac{p(x_1,\dots,x_n\mid\theta)\,\pi(\theta)}{\int_\Theta p(x_1,\dots,x_n\mid\theta)\,\pi(\theta)\,d\theta},$$

其中 $$\pi(\theta)$$ 为先验分布，$$m$$ 为样本的边际分布。一般用后验均值 $$E(\theta\mid\boldsymbol x)$$ 作为 $$\theta$$ 的贝叶斯估计。

**共轭先验示例**：$$x_i\sim b(1,p)$$，先验 $$p\sim Be(a,b)$$，则后验为 $$Be\left(a+\sum x_i,\ b+n-\sum x_i\right)$$，贝叶斯估计为

$$\hat p_B=\frac{a+\sum x_i}{a+b+n}.$$

直观理解：先验相当于“事先看过 $$a+b$$ 次试验，其中成功了 $$a$$ 次”，再和真实数据合在一起算频率。

---

**本系列目录**

1. [01｜概念简答索引]({{ site.baseurl }}{% post_url 2026-09-30-mao-notes-01-concepts %})
2. [02｜常用分布]({{ site.baseurl }}{% post_url 2026-09-30-mao-notes-02-distributions %})
3. [03｜多维分布、协方差与条件期望]({{ site.baseurl }}{% post_url 2026-09-30-mao-notes-03-multivariate %})
4. [04｜特征函数、大数定律与中心极限定理]({{ site.baseurl }}{% post_url 2026-09-30-mao-notes-04-limit-theorems %})
5. [05｜抽样分布、次序统计量与充分统计量]({{ site.baseurl }}{% post_url 2026-09-30-mao-notes-05-sampling %})
6. **06｜参数估计**（本篇）
7. [07｜区间估计、假设检验与习题结论]({{ site.baseurl }}{% post_url 2026-09-30-mao-notes-07-intervals-tests %})

[← 上一篇：抽样分布、次序统计量与充分统计量]({{ site.baseurl }}{% post_url 2026-09-30-mao-notes-05-sampling %})　·　[下一篇：区间估计、假设检验与习题结论 →]({{ site.baseurl }}{% post_url 2026-09-30-mao-notes-07-intervals-tests %})

*本系列是我复习茆诗松、程依明、濮晓龙《概率论与数理统计教程》时整理的笔记，共 7 篇。文中 P××× 为茆书页码，“指导书”为配套的学习指导与习题解答。如有错漏，欢迎指出。*
