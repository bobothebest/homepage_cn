---
layout: post
title: "随机过程笔记 03｜随机积分"
author: 张梓源
tags:
- 随机过程笔记
- 学习笔记
date: 2026-10-21 10:03 -0400
toc: true
math: true
---
> **导读**
>
> 随机积分其实就是**累计的交易盈亏（P/L）**：持有 $$\Delta$$ 股，价格变动 $$dW$$，加起来。这一篇分三步构造 $$\int_0^t\Delta_s\,dW_s$$——离散时间、"买入持有"的简单策略、用逼近得到的一般策略，并说明同样的五条性质每一步都保留下来：
>
> **线性 · 适应 · 鞅 · Itô 等距 $$\mathrm{Var}[I_t]=E\int_0^t\Delta_s^2ds$$ · $$[I,I]_t=\int_0^t\Delta_s^2ds$$。**

---

## 1. 离散时间：积分就是交易盈亏

设 $$M=\{M_n\}$$ 是鞅（想成股价）。对交易策略 $$\Delta=\{\Delta_n\}$$ 定义

$$I_0=0,\qquad I_n=\sum_{j=1}^n\underbrace{\Delta_j}_{\text{持股数}}\ \underbrace{dM_j}_{\text{价格变动}},\qquad dM_j=M_j-M_{j-1}.$$

$$\Delta_j\,dM_j$$ 是 $$[j-1,j]$$ 上的盈亏，$$I_n$$ 是 $$[0,n]$$ 上的累计盈亏。

**可料（predictable）策略。** 我们在 $$[j-1,j]$$ 上持有 $$\Delta_j$$ 股，而投资时不能看到未来，所以 $$\Delta_j$$ 必须在 $$j-1$$ 时刻就已知：$$\Delta_j$$ 是 $$\mathcal F_{j-1}$$ 可测的。这样的 $$\Delta=\{\Delta_n\}_{n\ge1}$$ 叫**可料**的（$$\Delta_0$$ 没有定义，因为 $$\mathcal F_{-1}$$ 不存在）。

### $$I=\sum\Delta\,dM$$ 的性质（$$\Delta$$ 可料，$$M$$ 是鞅）

1. **对 $$\Delta$$ 线性：** $$I(\alpha\Delta+\beta\psi)=\alpha I(\Delta)+\beta I(\psi)$$——拆开求和即可。
2. **适应：** $$\Delta_1,\dots,\Delta_n$$ 在 $$n-1$$ 时刻已知，$$M_0,\dots,M_n$$ 在 $$n$$ 时刻已知，所以 $$I_n$$ 是 $$\mathcal F_n$$ 可测的。
3. **鞅：** $$I_{n+1}-I_n=\Delta_{n+1}(M_{n+1}-M_n)$$，$$\Delta_{n+1}$$ 是 $$\mathcal F_n$$ 可测的，所以

   $$E[I_{n+1}-I_n\mid\mathcal F_n]=\Delta_{n+1}\,E[M_{n+1}-M_n\mid\mathcal F_n]=0.$$

4. **二次变差：** $$[I,I]_n=\sum_j\Delta_j^2(dM_j)^2=\sum_j\Delta_j^2\,d[M,M]_j$$。
5. **方差：** $$I$$ 是鞅且 $$I_0=0$$，由第 02 篇结论 1，

   $$\mathrm{Var}[I_n]=E\big[[I,I]_n\big]=E\Big[\sum_{j=1}^n\Delta_j^2\,d[M,M]_j\Big].$$

对随机游走，$$d[M,M]_j=1$$，所以 $$[I,I]_n=\sum_j\Delta_j^2$$，$$\mathrm{Var}[I_n]=E\sum_j\Delta_j^2$$。

"鞅"这条的直观：**只要下注前不知道硬币结果，怎么调整下注大小都赢不了公平的赌局。**

## 2. 连续时间：为什么直接定义行不通

我们想要 $$I_t=\int_0^t\Delta_s\,dW_s$$（$$\Delta$$：被积函数，$$W$$：积分器；盈亏的解释不变）。自然的想法是黎曼和：

$$\int_0^t\Delta_s\,dW_s\approx\sum_{j=1}^n\Delta_{t_{j-1}}\big(W_{t_j}-W_{t_{j-1}}\big).$$

如果 $$W$$ 光滑，中值定理会把它变成 $$\int_0^t\Delta_s\dot W_s\,ds$$。但 $$W$$ 处处不可微（第 02 篇），$$\dot W_s$$ 不存在，所以**逐条路径**取极限在通常意义下不存在。我们需要一种新的极限。

## 3. 简单（买入持有）被积函数

先回到"差不多是离散"的情形。固定 $$0=t_0<t_1<\dots<t_n$$ 和 $$\mathcal F_{t_j}$$ 可测的 $$\Delta_j$$，令

$$\Delta_t=\Delta_0\mathbf 1_{\{t=0\}}+\sum_{j=1}^n\Delta_{j-1}\mathbf 1_{\{t\in(t_{j-1},t_j]\}}:$$

在 $$(t_{j-1},t_j]$$ 上持有 $$\Delta_{j-1}$$ 股。对这种**简单被积函数**，积分可以直接写出来：

$$I_t=\sum_{j=1}^n\underbrace{\Delta_{j-1}}_{\text{头寸}}\ \underbrace{\big(W_{t\wedge t_j}-W_{t\wedge t_{j-1}}\big)}_{\text{价格变动}},\qquad t\wedge s=\min\{t,s\}.$$

分段来看：$$t_{j-1}\le t\le t_j$$ 时

$$I_t=\underbrace{\sum_{i=1}^{j-1}\Delta_{i-1}\big(W_{t_i}-W_{t_{i-1}}\big)}_{[0,t_{j-1}]\text{ 上的盈亏}}+\underbrace{\Delta_{j-1}\big(W_t-W_{t_{j-1}}\big)}_{[t_{j-1},t]\text{ 上的盈亏}},$$

$$t\ge t_n$$ 之后不再变化——我们已经停止交易了。

### 简单 $$\Delta$$ 的性质

- **线性**、**适应**：和离散时间一样。
- **路径连续**（新性质）：每一项 $$t\mapsto\Delta_{j-1}(W_{t\wedge t_j}-W_{t\wedge t_{j-1}})$$ 都连续，因为 $$W$$ 连续；有限个连续过程之和仍连续。
- **鞅。** 取 $$t_{j-1}<s<t_j<t$$（相邻两段）：

  $$\begin{aligned}
  E[I_t-I_s\mid\mathcal F_s]&=E[\Delta_{j-1}(W_{t_j}-W_s)\mid\mathcal F_s]+E\big[E[\Delta_j(W_t-W_{t_j})\mid\mathcal F_{t_j}]\mid\mathcal F_s\big]\\
  &=\Delta_{j-1}\cdot0+E[\Delta_j\cdot0\mid\mathcal F_s]=0.
  \end{aligned}$$

  关键在于头寸是在每段**开始时**选定的。

- **Itô 等距：** $$\mathrm{Var}[I_t]=E\Big[\int_0^t\Delta_s^2\,ds\Big]$$。
- **二次变差：** $$[I,I]_t=\int_0^t\Delta_s^2\,ds$$。

<details markdown="1"><summary>证明：简单被积函数的 Itô 等距</summary>

$$I_0=0$$ 且 $$I$$ 是鞅，所以 $$\mathrm{Var}[I_t]=E[I_t^2]$$。重新编号使 $$t=t_j$$，$$I_t=\sum_{i=1}^j\Delta_{t_{i-1}}(W_{t_i}-W_{t_{i-1}})$$。交叉项和第 02 篇一样消失（对较晚的左端点取条件），对角项为

$$E\big[\Delta_{t_{i-1}}^2(W_{t_i}-W_{t_{i-1}})^2\big]=E\big[\Delta_{t_{i-1}}^2\,E[(W_{t_i}-W_{t_{i-1}})^2\mid\mathcal F_{t_{i-1}}]\big]=E\big[\Delta_{t_{i-1}}^2(t_i-t_{i-1})\big].$$

$$\Delta$$ 在每段 $$(t_{i-1},t_i]$$ 上是常数，所以

$$\mathrm{Var}[I_t]=E\Big[\sum_{i=1}^j\Delta_{t_{i-1}}^2(t_i-t_{i-1})\Big]=E\Big[\int_0^t\Delta_s^2\,ds\Big].$$

</details>

<details markdown="1"><summary>证明：简单被积函数的二次变差</summary>

把一段 $$[t_{i-1},t_i]$$ 细分为 $$s_0<s_1<\dots<s_n$$。段内头寸不变，$$I_{s_k}-I_{s_{k-1}}=\Delta_{t_{i-1}}(W_{s_k}-W_{s_{k-1}})$$，所以

$$[I,I]^\Pi_{[t_{i-1},t_i]}=\Delta_{t_{i-1}}^2\,[W,W]^\Pi_{[t_{i-1},t_i]}\ \longrightarrow\ \Delta_{t_{i-1}}^2(t_i-t_{i-1}).$$

对各段求和，$$[I,I]_t=\sum_i\Delta_{t_{i-1}}^2(t_i-t_{i-1})=\int_0^t\Delta_s^2\,ds$$。

</details>

注意后两条的区别：$$[I,I]_t$$ 是逐条路径的**随机**量；等距是它的**期望**。

## 4. 一般被积函数：先逼近，再取 $$L^2$$ 极限

固定 $$T>0$$。记 $$\mathcal L$$ 为满足以下条件的被积函数 $$\Delta$$ 的全体：

- 适应；
- 平方可积：$$E\Big[\int_0^T\Delta_s^2\,ds\Big]<\infty$$。

对 $$\Delta\in\mathcal L$$ 和步长趋于 0 的划分 $$\Pi^n$$，把 $$\Delta$$ 冻结在左端点：

$$\Delta^n_t=\Delta_0\mathbf 1_{\{t=0\}}+\sum_{j=1}^n\Delta_{t_{j-1}}\mathbf 1_{\{t\in(t_{j-1},t_j]\}},\qquad I^n_t=I(\Delta^n)_t .$$

$$\Delta^n$$ 是 $$\Delta$$ "只在网格时刻调仓"的简单版本，$$I^n$$ 有明确定义。固定 $$\omega$$ 时 $$\lim I^n_t(\omega)$$ 不一定存在；我们改为取**均方意义下**的极限。

**第 1 步——逼近是好的。** 课上的一个关键结论：对 $$\Delta\in\mathcal L$$，

$$\lim_{n\to\infty}E\Big[\int_0^T(\Delta^n_s-\Delta_s)^2\,ds\Big]=0,$$

记作 $$\Delta^n\to\Delta$$ 于 $$L^2([0,T]\times\Omega)$$。（直观：调仓网格越来越细，冻结后的策略在时间和情景上平均来说就越来越贴近真实策略。）

**第 2 步——所以逼近序列是柯西列。** 由 $$(a+b)^2\le2a^2+2b^2$$，

$$(\Delta^n_u-\Delta^m_u)^2\le2(\Delta^n_u-\Delta_u)^2+2(\Delta_u-\Delta^m_u)^2\ \Longrightarrow\ \lim_{n,m\to\infty}E\int_0^T(\Delta^n_s-\Delta^m_s)^2ds=0.$$

**第 3 步——等距把这一点传给积分。** 由线性，$$I^n_t-I^m_t=I(\Delta^n-\Delta^m)_t$$，而 $$\Delta^n-\Delta^m$$ 是简单的，所以

$$\mathrm{Var}\big[I^n_t-I^m_t\big]=E\Big[\int_0^t(\Delta^n_u-\Delta^m_u)^2du\Big]\le E\Big[\int_0^T(\Delta^n_u-\Delta^m_u)^2du\Big]\xrightarrow{n,m\to\infty}0.$$

**第 4 步——把积分定义为极限。** $$I^n$$ 在均方意义下是柯西列，所以（作为过程）收敛到某个 $$I$$，定义

$$\Big(\int_0^t\Delta_s\,dW_s\Big)(\omega):=I_t(\omega),\qquad t\le T.$$

**等距为什么是主角。** 它说"策略 $$\mapsto$$ 盈亏"这个映射保持距离：两个策略在 $$E\int(\cdot)^2$$ 意义下接近，盈亏在方差意义下就接近。于是"逼近策略"自动给出"逼近积分"——即使没有任何一条路径按经典意义收敛。

## 5. 一般被积函数的性质

对 $$\Delta\in\mathcal L$$ 全部保留：

- **线性：** $$I(\alpha\Delta+\beta\psi)=\alpha I(\Delta)+\beta I(\psi)$$；
- **路径连续**，$$I_0=0$$；
- **适应**；
- **鞅**；
- **Itô 等距：** $$\mathrm{Var}[I_t]=E\Big[\int_0^t\Delta_s^2\,ds\Big]$$；
- **二次变差：** 逐条路径 $$[I,I]_t=\int_0^t\Delta_s^2\,ds$$。

用微分记号：$$dI_t=\Delta_t\,dW_t$$，$$d[I,I]_t=\Delta_t^2\,dt=(\Delta_t\,dW_t)^2$$——与 $$(dW_t)^2=dt$$ 一致。

---

**本系列目录**

1. [01｜概率空间、信息流、条件期望与鞅]({{ site.baseurl }}{% post_url 2026-10-21-mf790-01-foundations %})
2. [02｜随机游走、二次变差与布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-02-brownian-motion %})
3. **03｜随机积分**（本篇）
4. [04｜Itô 公式、Itô 过程与几何布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-04-ito-formula %})
5. [05｜Black–Scholes：复制、偏微分方程与希腊字母]({{ site.baseurl }}{% post_url 2026-10-21-mf790-05-black-scholes %})
6. [06｜风险中性定价、测度变换与 Girsanov 定理]({{ site.baseurl }}{% post_url 2026-10-21-mf790-06-risk-neutral-pricing %})
7. [07｜随机微分方程、Feynman–Kac 与停时]({{ site.baseurl }}{% post_url 2026-10-21-mf790-07-sde-feynman-kac %})

[← 上一篇：随机游走、二次变差与布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-02-brownian-motion %})　·　[下一篇：Itô 公式、Itô 过程与几何布朗运动 →]({{ site.baseurl }}{% post_url 2026-10-21-mf790-04-ito-formula %})

*本系列是我在波士顿大学 Questrom 商学院修读 QST MF 790 Stochastic Calculus（随机分析）时整理的课程笔记，共 7 篇。内容按课程顺序整理，表述、直观解释和几处勘误是我自己加的。文中 $$N(\mu,\sigma^2)$$ 的第二个参数一律为方差。如有错漏，欢迎指出。本篇英文版：[English version](https://bobothebest.github.io/homepage_en/2026/10/21/mf790-03-stochastic-integration/)。*
