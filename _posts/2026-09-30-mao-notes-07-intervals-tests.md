---
layout: post
title: "茆书复习笔记 07｜区间估计、假设检验与习题结论汇总"
author: 张梓源
tags:
- 茆书复习笔记
- 数理统计
- 置信区间
- 假设检验
- 茆诗松
- 考研
- 读书笔记
date: 2026-09-30 10:07 -0400
toc: true
---
> **导读**
>
> 系列最后一篇，把区间估计和假设检验全部整理成表：单正态总体、两正态总体的各种情形，大样本方法，以及 u、t、χ²、F 检验。
>
> 记这些表格有个窍门：**区间和检验用的是同一个枢轴量**。先想清楚“在什么条件下，哪个统计量服从什么分布”，区间和拒绝域都能自己推出来。文末附上从习题中推出的常用结论。

---

> 记号：$$u_p$$、$$t_p$$、$$\chi^2_p$$、$$F_p$$ 都是**下侧** $$p$$ 分位数。置信水平为 $$1-\alpha$$。

---

## 1. 置信区间与枢轴量（P343～P352，指导书 P369）

**枢轴量法**：构造一个只含待估参数、分布已知且与参数无关的函数 $$G$$（枢轴量），由 $$P(a\le G\le b)=1-\alpha$$ 反解出参数的区间。

### 单个正态总体 $$N(\mu,\sigma^2)$$

| 待估 | 条件 | 枢轴量 | $$1-\alpha$$ 置信区间 |
|---|---|---|---|
| $$\mu$$ | $$\sigma$$ 已知 | $$\dfrac{\bar x-\mu}{\sigma/\sqrt n}\sim N(0,1)$$ | $$\bar x\pm u_{1-\alpha/2}\dfrac{\sigma}{\sqrt n}$$ |
| $$\mu$$ | $$\sigma$$ 未知 | $$\dfrac{\sqrt n(\bar x-\mu)}{s}\sim t(n-1)$$ | $$\bar x\pm t_{1-\alpha/2}(n-1)\dfrac{s}{\sqrt n}$$ |
| $$\sigma^2$$ | $$\mu$$ 未知 | $$\dfrac{(n-1)s^2}{\sigma^2}\sim\chi^2(n-1)$$ | $$\left[\dfrac{(n-1)s^2}{\chi^2_{1-\alpha/2}(n-1)},\ \dfrac{(n-1)s^2}{\chi^2_{\alpha/2}(n-1)}\right]$$ |

### 两个正态总体 $$X\sim N(\mu_1,\sigma_1^2)$$（样本量 $$m$$），$$Y\sim N(\mu_2,\sigma_2^2)$$（样本量 $$n$$）

**$$\mu_1-\mu_2$$ 的置信区间**

| 情形 | 区间 |
|---|---|
| ① $$\sigma_1^2,\sigma_2^2$$ 已知 | $$\bar x-\bar y\pm u_{1-\alpha/2}\sqrt{\dfrac{\sigma_1^2}{m}+\dfrac{\sigma_2^2}{n}}$$ |
| ② $$\sigma_1^2=\sigma_2^2$$ 未知 | $$\bar x-\bar y\pm\sqrt{\dfrac{m+n}{mn}}\,s_w\,t_{1-\alpha/2}(m+n-2)$$，$$s_w^2=\dfrac{(m-1)s_x^2+(n-1)s_y^2}{m+n-2}$$ |
| ③ $$\sigma_1^2/\sigma_2^2=\theta$$ 已知 | $$\bar x-\bar y\pm\sqrt{\dfrac{m\theta+n}{mn}}\,s_t\,t_{1-\alpha/2}(m+n-2)$$，$$s_t^2=\dfrac{(m-1)s_x^2+(n-1)s_y^2/\theta}{m+n-2}$$ |
| ④ $$m,n$$ 都很大（近似） | $$\bar x-\bar y\pm u_{1-\alpha/2}\sqrt{\dfrac{s_x^2}{m}+\dfrac{s_y^2}{n}}$$ |
| ⑤ 一般场合，小样本（近似 t，Welch） | $$\bar x-\bar y\pm s_0\,t_{1-\alpha/2}(l)$$，$$s_0^2=\dfrac{s_x^2}m+\dfrac{s_y^2}n$$，$$l=\dfrac{s_0^4}{\frac{s_x^4}{m^2(m-1)}+\frac{s_y^4}{n^2(n-1)}}$$（取最接近的整数） |

**$$\sigma_1^2/\sigma_2^2$$ 的置信区间**（$$\mu_1,\mu_2$$ 未知）

$$F=\frac{s_x^2/\sigma_1^2}{s_y^2/\sigma_2^2}\sim F(m-1,n-1),\qquad\left[\frac{s_x^2}{s_y^2}\cdot\frac{1}{F_{1-\alpha/2}(m-1,n-1)},\ \frac{s_x^2}{s_y^2}\cdot\frac{1}{F_{\alpha/2}(m-1,n-1)}\right].$$

### 大样本置信区间（会推）

$$X\sim b(1,p)$$，由中心极限定理 $$\dfrac{\bar x-p}{\sqrt{p(1-p)/n}}\ \dot\sim\ N(0,1)$$，将分母中的 $$p$$ 用 $$\bar x$$ 代替，得到

$$\bar x\pm u_{1-\alpha/2}\sqrt{\frac{\bar x(1-\bar x)}{n}}.$$

<details markdown="1"><summary>更精确的做法（不代入 \(\bar x\)，直接解二次不等式）</summary>

由 $$\left\vert \bar x-p\right\vert \le u_{1-\alpha/2}\sqrt{p(1-p)/n}$$，两边平方后得到关于 $$p$$ 的二次不等式：

$$\left(1+\frac{u^2}n\right)p^2-\left(2\bar x+\frac{u^2}n\right)p+\bar x^2\le0,\qquad u=u_{1-\alpha/2}.$$

两根之间就是 $$p$$ 的置信区间。$$n$$ 很大时，它与上面的简化区间几乎相同。

</details>

---

## 2. 假设检验汇总

**假设检验的基本步骤（P357）**：见第 1 篇。  
**p 值**：在 $$H_0$$ 成立的条件下，得到“和当前观测值一样极端或更极端”结果的概率。p 值 $$\le\alpha$$ 时拒绝 $$H_0$$。

### 单个正态总体均值 $$\mu$$ 的检验（P369）

| 检验 | $$H_0$$ | $$H_1$$ | 检验统计量 | 拒绝域 $$W$$ | p 值 |
|---|---|---|---|---|---|
| u 检验（$$\sigma$$ 已知） | $$\mu\le\mu_0$$ | $$\mu>\mu_0$$ | $$u=\dfrac{\bar x-\mu_0}{\sigma/\sqrt n}$$ | $$u\ge u_{1-\alpha}$$ | $$1-\Phi(u_0)$$ |
| | $$\mu\ge\mu_0$$ | $$\mu<\mu_0$$ | | $$u\le u_\alpha$$ | $$\Phi(u_0)$$ |
| | $$\mu=\mu_0$$ | $$\mu\ne\mu_0$$ | | $$\lvert u\rvert\ge u_{1-\alpha/2}$$ | $$2[1-\Phi(\lvert u_0\rvert)]$$ |
| t 检验（$$\sigma$$ 未知） | $$\mu\le\mu_0$$ | $$\mu>\mu_0$$ | $$t=\dfrac{\bar x-\mu_0}{s/\sqrt n}$$ | $$t\ge t_{1-\alpha}(n-1)$$ | $$P(T\ge t_0)$$ |
| | $$\mu\ge\mu_0$$ | $$\mu<\mu_0$$ | | $$t\le t_\alpha(n-1)$$ | $$P(T\le t_0)$$ |
| | $$\mu=\mu_0$$ | $$\mu\ne\mu_0$$ | | $$\lvert t\rvert\ge t_{1-\alpha/2}(n-1)$$ | $$2P(T\ge\lvert t_0\rvert)$$ |

（$$u_0$$、$$t_0$$ 为统计量的观测值，$$T\sim t(n-1)$$。）

### 两个正态总体均值差的检验（P373）

原假设都写成 $$H_0:\mu_1-\mu_2=0$$（或 $$\le0$$、$$\ge0$$），拒绝域的方向同上表。

| 检验 | 条件 | 检验统计量 | $$H_0$$ 下的分布 |
|---|---|---|---|
| 两样本 u 检验 | $$\sigma_1,\sigma_2$$ 已知 | $$u=\dfrac{\bar x-\bar y}{\sqrt{\sigma_1^2/m+\sigma_2^2/n}}$$ | $$N(0,1)$$ |
| 两样本 t 检验 | $$\sigma_1=\sigma_2=\sigma$$ 未知 | $$t=\dfrac{\bar x-\bar y}{s_w\sqrt{\frac1m+\frac1n}}$$ | $$t(m+n-2)$$ |
| 大样本 u 检验 | $$m,n$$ 充分大 | $$u=\dfrac{\bar x-\bar y}{\sqrt{s_x^2/m+s_y^2/n}}$$ | 近似 $$N(0,1)$$ |
| 近似 t 检验 | $$m,n$$ 不很大，方差未知且不等 | $$t=\dfrac{\bar x-\bar y}{\sqrt{s_x^2/m+s_y^2/n}}$$ | 近似 $$t(l)$$，$$l$$ 同上 |

### 方差的检验

| 检验 | $$H_0$$ | 检验统计量 | 拒绝域（双侧 / 右侧 / 左侧） |
|---|---|---|---|
| $$\chi^2$$ 检验 | $$\sigma^2=\sigma_0^2$$ | $$\chi^2=\dfrac{(n-1)s^2}{\sigma_0^2}\sim\chi^2(n-1)$$ | $$\chi^2\le\chi^2_{\alpha/2}$$ 或 $$\chi^2\ge\chi^2_{1-\alpha/2}$$ ／ $$\chi^2\ge\chi^2_{1-\alpha}$$ ／ $$\chi^2\le\chi^2_\alpha$$ |
| F 检验 | $$\sigma_1^2=\sigma_2^2$$ | $$F=\dfrac{s_x^2}{s_y^2}\sim F(m-1,n-1)$$ | $$F\le F_{\alpha/2}$$ 或 $$F\ge F_{1-\alpha/2}$$ ／ $$F\ge F_{1-\alpha}$$ ／ $$F\le F_\alpha$$ |

**区间估计与检验的对偶关系**：$$\mu=\mu_0$$ 的双侧检验在水平 $$\alpha$$ 下被接受，当且仅当 $$\mu_0$$ 落在 $$1-\alpha$$ 置信区间内。

**势函数（P360）**：例如对 u 检验 $$H_0:\mu\le\mu_0$$ vs $$H_1:\mu>\mu_0$$，

$$g(\mu)=P_\mu(u\ge u_{1-\alpha})=1-\Phi\left(u_{1-\alpha}-\frac{\mu-\mu_0}{\sigma/\sqrt n}\right),$$

它是 $$\mu$$ 的增函数；在 $$\mu=\mu_0$$ 处等于 $$\alpha$$。

---

## 3. 习题推得的结论汇总

| 习题 | 结论 |
|---|---|
| 1.2.1(3) | $$\binom n0+\binom n1+\cdots+\binom nn=2^n$$ |
| 1.2.1(4) | $$\binom n1+2\binom n2+\cdots+n\binom nn=n2^{n-1}$$ |
| 1.2.1(5) | 范德蒙恒等式：$$\sum_{k=0}^{n}\binom ak\binom b{n-k}=\binom{a+b}n$$，$$n\le\min\{a,b\}$$ |
| 1.3.2 | 概率不能反推事件：$$P(AB)=0\not\Rightarrow AB=\varnothing$$（例如连续型随机变量取单点的概率为 0，但这个事件并非不可能事件） |
| 2.2.20 | 连续型：$$E(X)=\int_0^{\infty}[1-F(x)]dx-\int_{-\infty}^0F(x)dx$$ |
| 2.2.21 | 非负连续型：$$E(X)=\int_0^{\infty}P(X>x)dx$$，$$E(X^n)=\int_0^{\infty}nx^{n-1}P(X>x)dx$$ |
| 2.3.8 | 熟练应用上面两条；$$\int_0^{\infty}x^2e^{-x^2}dx=\frac{\sqrt\pi}4$$ |
| 2.3.9 | $$E(X-EX)^2\le E(X-c)^2$$（均值处最小）；$$E\lvert X-m\rvert\le E\lvert X-c\rvert$$（中位数处最小） |
| 2.3.12 | $$g$$ 非负不减，$$E g(X)$$ 存在 $$\Rightarrow P(X>\varepsilon)\le\frac{E g(X)}{g(\varepsilon)}$$ |
| — | 有界随机变量 $$X\in[a,b]$$：$$a\le EX\le b$$，$$\mathrm{Var}X\le\left(\frac{b-a}2\right)^2$$ |
| 2.7.3(3) | $$LN(\mu,\sigma^2)$$ 的中位数为 $$e^\mu$$ |
| 2.7.4 | $$Ga(\alpha,\lambda)$$：$$E(X^k)=\frac{\Gamma(\alpha+k)}{\Gamma(\alpha)\lambda^k}$$ |
| 2.7.5 | $$Exp(\lambda)$$：$$E(X^k)=\frac{k!}{\lambda^k}$$，$$C_v=1$$，$$\beta_s=2$$，$$\beta_k=6$$ |
| 2.7.11 | $$LN(\mu,\sigma^2)$$ 的 $$p$$ 分位数为 $$e^{\mu+\sigma u_p}$$ |
| 3.3.12 | $$\min\{a,b\}=\frac12[a+b-\lvert b-a\rvert]$$，$$\max\{a,b\}=\frac12[a+b+\lvert b-a\rvert]$$ |
| 3.4.14 | $$X,Y\overset{iid}{\sim}N(0,1)$$：$$E\max\{X,Y\}=\frac1{\sqrt\pi}$$ |
| 3.4.18 | $$X,Y\overset{iid}{\sim}N(a,\sigma^2)$$：$$E\max\{X,Y\}=a+\frac\sigma{\sqrt\pi}$$ |
| 配对模型 | 抽礼物问题：$$E(S_n)=1$$，$$\mathrm{Var}(S_n)=1$$ |
| 4.1.18 | iid，$$EX=0$$，$$\mathrm{Var}X=\sigma^2$$ $$\Rightarrow\frac1n\sum X_k^2\xrightarrow{P}\sigma^2$$ |
| 4.3.4 | 只有相邻项相关时：$$\mathrm{Var}(\sum X_i)=\sum\mathrm{Var}(X_i)+2\sum_{i=1}^{n-1}\mathrm{Cov}(X_i,X_{i+1})$$ |
| 5.3.4 | $$s_{n+1}^2=\frac{n-1}ns_n^2+\frac1{n+1}(x_{n+1}-\bar x_n)^2$$ |
| 5.3.10 | $$(\sum x_i)^2=\sum x_i^2+2\sum_{i<j}x_ix_j$$ |
| 5.3.22–23 | 定理 5.3.6（两个次序统计量的联合密度）的应用 |
| 5.3.28 | $$\eta_i=F(x_{(i)})$$ 是 $$U(0,1)$$ 的次序统计量；$$E\eta_i=\frac i{n+1}$$，$$\mathrm{Var}\eta_i=\frac{i(n+1-i)}{(n+1)^2(n+2)}$$ |
| 5.3.34 | $$U(0,\theta)$$ 的全体次序统计量的联合密度为 $$\frac{n!}{\theta^n}$$，$$0<u_1\le\dots\le u_n\le\theta$$ |
| 6.3.2 | MLE 可能不止一个 |
| 6.4.5 | $$I(\theta)=-E\left[\frac{\partial^2}{\partial\theta^2}\ln p(x;\theta)\right]$$ |

---

**本系列目录**

1. [01｜概念简答索引]({% post_url 2026-09-30-mao-notes-01-concepts %})
2. [02｜常用分布]({% post_url 2026-09-30-mao-notes-02-distributions %})
3. [03｜多维分布、协方差与条件期望]({% post_url 2026-09-30-mao-notes-03-multivariate %})
4. [04｜特征函数、大数定律与中心极限定理]({% post_url 2026-09-30-mao-notes-04-limit-theorems %})
5. [05｜抽样分布、次序统计量与充分统计量]({% post_url 2026-09-30-mao-notes-05-sampling %})
6. [06｜参数估计]({% post_url 2026-09-30-mao-notes-06-estimation %})
7. **07｜区间估计、假设检验与习题结论**（本篇）

**系列完结**，感谢一路看到这里。　·　[← 上一篇：参数估计]({% post_url 2026-09-30-mao-notes-06-estimation %})

*本系列是我复习茆诗松、程依明、濮晓龙《概率论与数理统计教程》时整理的笔记，共 7 篇。文中 P××× 为茆书页码，“指导书”为配套的学习指导与习题解答。如有错漏，欢迎指出。*
