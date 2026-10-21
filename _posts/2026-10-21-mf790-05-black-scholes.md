---
layout: post
title: "随机过程笔记 05｜Black–Scholes：复制、偏微分方程与希腊字母"
author: 张梓源
tags:
- 随机过程笔记
- BU MF 790
- 随机分析
- Black-Scholes
- 期权定价
- 希腊字母
- 隐含波动率
- 读书笔记
date: 2026-10-21 10:05 -0400
toc: true
math: true
---
> **导读**
>
> 期权的价格就是**复制它所需要的本金**。这一篇搭好 Black–Scholes 市场，写出自融资交易者的财富，再把它和期权价格过程逐项对齐（$$dW$$ 项和 $$dt$$ 项）。对齐 $$dW$$ 项得到 **Delta 对冲**，对齐 $$dt$$ 项得到 **Black–Scholes 偏微分方程**。解方程得到 Black–Scholes 公式，对公式求导得到希腊字母。
>
> 一个值得早点注意的意外：股票的期望收益率 $$\mu$$ 从价格里**消失**了。

---

## 1. 市场

固定 $$(\Omega,\mathcal F,\mathbb F,\mathbb P)$$ 和布朗运动 $$W$$。两种资产可交易：

1. **股票**（风险资产），$$S\sim GBM(\mu,\sigma)$$：$$\dfrac{dS_t}{S_t}=\mu\,dt+\sigma\,dW_t$$，$$S_0>0$$。$$S_t(\omega)$$ 是情景 $$\omega$$ 下 $$t$$ 时刻一股的价格。
2. **货币市场**（无风险资产），$$\dfrac{dB_t}{B_t}=r\,dt$$，$$B_0=1$$：0 时刻的 1 元到 $$t$$ 时刻变成 $$e^{rt}$$ 元，$$r>0$$ 为常数。

## 2. 自融资策略与财富

**交易策略** $$\Delta$$ 是适应过程；$$\Delta_t(\omega)$$ 是 $$(t,\omega)$$ 时持有的股数。初始资金为 $$x$$。**自融资**：不在股票里的钱全部放在货币市场，0 时刻之后没有钱进出。

财富 $$X=X^{x,\Delta}$$ 拆成

$$X_t=\underbrace{\Delta_tS_t}_{\text{股票}}+\underbrace{X_t-\Delta_tS_t}_{\text{货币市场}},$$

货币市场部分相当于持有 $$(X_t-\Delta_tS_t)/B_t$$ "股" $$B$$。

**财富动态。** 在 $$[t,t+dt]$$ 开始时定好头寸：

$$X_{t+dt}-X_t=\Delta_t(S_{t+dt}-S_t)+\frac{X_t-\Delta_tS_t}{B_t}(B_{t+dt}-B_t)\ \xrightarrow{dt\to0}\ dX_t=\Delta_t\,dS_t+(X_t-\Delta_tS_t)\frac{dB_t}{B_t}.$$

代入 $$S$$ 和 $$B$$ 的动态：

$$dX_t=\big(rX_t+(\mu-r)\Delta_tS_t\big)dt+\sigma\Delta_tS_t\,dW_t,$$

$$X$$ 是 Itô 过程。读法：全部资金赚 $$r$$，股票部分再多赚超额收益 $$\mu-r$$，再加上噪声。

## 3. "价格"指的是什么

我们要给欧式**看涨期权**定价：在到期日 $$T$$ 以执行价 $$K$$ 买入 $$S$$ 的权利，收益 $$(S_T-K)^+=\max(S_T-K,0)$$。

**无套利价格 = 复制成本。** 如果从 $$x$$ 出发、用某个 $$\Delta$$ 自融资交易，最后以概率 1 得到 $$X^{x,\Delta}_T=(S_T-K)^+$$，那么看涨期权的价格必须是 $$x$$。任何别的价格都能让你卖贵的、买便宜的，无风险地赚差价。

**起点相同 + 动态相同 = 同一个过程。** 若 $$dY_t=\Theta^Y_tdt+\Delta^Y_tdW_t$$，$$dZ_t=\Theta^Z_tdt+\Delta^Z_tdW_t$$，且 $$Y_0=Z_0$$，$$\Theta^Y\equiv\Theta^Z$$，$$\Delta^Y\equiv\Delta^Z$$，则对所有 $$t\le T$$ 以概率 1 有 $$Y_t=Z_t$$。这就是一句话版的复制论证：对齐起点，再对齐 $$dt$$ 项和 $$dW$$ 项。

## 4. 把财富和期权价格对齐

令 $$c(t,s)$$ 为 $$S_t=s$$ 时 $$t$$ 时刻的看涨期权价格。已知 $$c(T,s)=(s-K)^+$$，要求 $$t<T$$ 时的 $$c$$。用**贴现**后的过程（代数更简单）：$$e^{-rt}X_t$$ 和 $$e^{-rt}c(t,S_t)$$。

**贴现财富。** 用第 04 篇的贴现公式，

$$d(e^{-rt}X_t)=e^{-rt}(\mu-r)\Delta_tS_t\,dt+e^{-rt}\sigma\Delta_tS_t\,dW_t .$$

**贴现价格。** 假设 $$c$$ 光滑（严格说收益函数在 $$s=K$$ 处不可微，但 $$t<T$$ 时 $$c$$ 是光滑的）。由 $$d[S,S]_t=\sigma^2S_t^2dt$$，记 $$(\cdot)=(t,S_t)$$：

$$d\big(e^{-rt}c(t,S_t)\big)=e^{-rt}\Big(-rc+c_t+c_s\,\mu S_t+\tfrac12c_{ss}\,\sigma^2S_t^2\Big)(\cdot)\,dt+e^{-rt}c_s(\cdot)\,\sigma S_t\,dW_t .$$

令 $$x=c(0,S_0)$$ 并对齐：

- **$$dW$$ 项：** $$e^{-rt}\Delta_t\sigma S_t=e^{-rt}c_s(t,S_t)\sigma S_t$$，所以

  $$\boxed{\ \Delta_t=c_s(t,S_t)\ }\quad\text{——Delta 对冲。}$$

- **$$dt$$ 项：** $$\Delta_t(\mu-r)S_t=-rc+c_t+c_s\mu S_t+\frac12c_{ss}\sigma^2S_t^2$$。代入 $$\Delta_t=c_s$$，两边的 $$c_s\mu S_t$$ **抵消**，剩下

  $$0=-rc(t,S_t)+c_t(t,S_t)+rS_t\,c_s(t,S_t)+\tfrac12\sigma^2S_t^2\,c_{ss}(t,S_t).$$

只要对所有 $$t<T$$、$$s>0$$ 都成立就够了，这就是 **Black–Scholes 偏微分方程**：

$$\begin{cases}0=-rc(t,s)+c_t(t,s)+rs\,c_s(t,s)+\tfrac12\sigma^2s^2\,c_{ss}(t,s),&t<T,\ s>0,\\[4pt] c(T,s)=(s-K)^+,&s>0.\end{cases}$$

**$$\mu$$ 不见了。** 对冲把股票的风险完全抵消掉了，所以股票预期涨多快和价格无关。

## 5. 解方程

两步：（1）假设解 $$c$$ 存在，求出它的表达式；（2）验证这个表达式确实满足方程。

**第 1 步。** 令 $$\widetilde S\sim GBM(r,\sigma)$$，在某个测度 $$\widetilde{\mathbb P}$$ 下、布朗运动为 $$\widetilde W$$——注意**换了测度**，漂移是 **$$r$$ 而不是 $$\mu$$**。由 Itô，

$$d\big(e^{-rt}c(t,\widetilde S_t)\big)=e^{-rt}\big(\underbrace{-rc+c_t+r\widetilde S_tc_s+\tfrac12\sigma^2\widetilde S_t^2c_{ss}}_{=0\ \text{（由方程）}}\big)dt+e^{-rt}c_s\,\sigma\widetilde S_t\,d\widetilde W_t,$$

所以贴现价格是 **$$\widetilde{\mathbb P}$$ 下的鞅**。结合终值条件，

$$e^{-rt}c(t,\widetilde S_t)=E^{\widetilde{\mathbb P}}\big[e^{-rT}c(T,\widetilde S_T)\mid\widetilde{\mathcal F}_t\big]\quad\Longrightarrow\quad c(t,\widetilde S_t)=E^{\widetilde{\mathbb P}}\big[e^{-r(T-t)}(\widetilde S_T-K)^+\mid\widetilde{\mathcal F}_t\big].$$

用 $$\widetilde S_t$$ 表示 $$\widetilde S_T$$：

$$\widetilde S_T=\widetilde S_t\exp\Big\{\big(r-\tfrac12\sigma^2\big)(T-t)+\sigma(\widetilde W_T-\widetilde W_t)\Big\}.$$

$$\widetilde S_t$$ 在 $$t$$ 时已知，$$\widetilde W_T-\widetilde W_t\perp\!\!\!\perp\widetilde{\mathcal F}_t$$，所以可以把 $$\widetilde S_t$$ 固定为 $$s$$。由 $$\widetilde W_T-\widetilde W_t\sim N(0,T-t)$$：

$$c(t,s)=\frac{e^{-r(T-t)}}{\sqrt{2\pi(T-t)}}\int_{-\infty}^{\infty}\Big(s\,e^{(r-\frac12\sigma^2)(T-t)+\sigma z}-K\Big)^+e^{-\frac{z^2}{2(T-t)}}\,dz .$$

算出这个积分就是 **Black–Scholes 公式**：

$$\boxed{\ c(t,s)=s\,N\big(d_+(T-t,s)\big)-Ke^{-r(T-t)}N\big(d_-(T-t,s)\big)\ }$$

$$d_\pm(\tau,s)=\frac1{\sigma\sqrt\tau}\Big(\log\frac sK+\big(r\pm\tfrac12\sigma^2\big)\tau\Big),\qquad N(y)=\frac1{\sqrt{2\pi}}\int_{-\infty}^ye^{-z^2/2}\,dz .$$

公式里只出现**剩余期限** $$\tau=T-t$$，不出现绝对时间。

**第 2 步。** 直接验证这个 $$c$$ 满足方程且 $$c(T,s)=(s-K)^+$$（是一道作业题；关键一步就是下面第 7 节验证的恒等式）。

## 6. 希腊字母

把 $$c$$ 看成 $$c(t,s;T,K,r,\sigma)$$。对这些输入的敏感度就是**希腊字母**。记 $$\tau=T-t$$，$$d_\pm=d_\pm(\tau,s)$$。

| 希腊字母 | 定义 | 公式 | 符号 |
|---|---|---|---|
| Delta | $$c_s$$ | $$N(d_+)$$ | $$>0$$ |
| Gamma | $$c_{ss}$$ | $$\dfrac{N'(d_+)}{\sigma s\sqrt\tau}$$ | $$>0$$ |
| Theta | $$c_t$$ | $$-rKe^{-r\tau}N(d_-)-\dfrac{\sigma s}{2\sqrt\tau}N'(d_+)$$ | $$<0$$ |
| Vega | $$c_\sigma$$ | $$s\sqrt\tau\,N'(d_+)$$ | $$>0$$ |

**Delta。** 股价越高，看涨期权越值钱。复制策略是 $$\Delta_t=c_s(t,S_t)$$：只看 $$S$$ 的一阶变化（假装 $$dt=0$$），$$dX_t=c_s\,dS_t=dc(t,S_t)$$，对冲和期权同步变动——我们是 **Delta 中性**的。

**Gamma。** $$c_{ss}>0$$：看涨期权关于股价是**凸**的。看到 $$S$$ 的二阶变化，

$$dX_t=c_s\,dS_t,\qquad dc(t,S_t)=c_s\,dS_t+\tfrac12c_{ss}\,d[S,S]_t .$$

*"做多凸性"的直观。* 持有看涨期权、卖空 $$c_s$$ 股。Delta 对冲是凸价格曲线的一条切线；股价无论往哪边跳，曲线都在切线上方，所以对冲后的组合从任意方向的变动里都赚到 $$\frac12c_{ss}(dS)^2\ge0$$。Gamma 在执行价附近最大，越临近到期越集中。

**Theta。** $$c_t<0$$：其他条件不变时，随着时间流逝期权贬值。这是做多 Gamma 的代价：把所有二阶变化都算上，

$$dX_t=c_s\,dS_t+(c-c_sS_t)\,r\,dt,\qquad dc(t,S_t)=c_t\,dt+c_s\,dS_t+\tfrac12c_{ss}\sigma^2S_t^2\,dt,$$

而方程说这两者相等——**完美对冲**。Theta 恰好支付了 Gamma 的收益加上融资成本。

**Vega。** $$c_\sigma>0$$：波动率越大，价格变动越大；上涨对看涨期权有利，下跌最多跌到 0。Vega 为正 = **做多波动率**。

## 7. 验证：对冲是完美的

由 $$c=SN(d_+)-Ke^{-r\tau}N(d_-)$$，$$c_s=N(d_+)$$，$$c_{ss}=\frac{N'(d_+)}{\sigma S\sqrt\tau}$$ 以及表中的 $$c_t$$，需要验证 $$(c-c_sS)\,r=c_t+\frac12c_{ss}\sigma^2S^2$$：

$$\text{左边}=\big[SN(d_+)-Ke^{-r\tau}N(d_-)-SN(d_+)\big]r=-rKe^{-r\tau}N(d_-),$$

$$\text{右边}=-rKe^{-r\tau}N(d_-)-\frac{\sigma S}{2\sqrt\tau}N'(d_+)+\frac12\cdot\frac{\sigma S}{\sqrt\tau}N'(d_+)=-rKe^{-r\tau}N(d_-).\ \checkmark$$

Theta 项和 Gamma 项正好抵消。

## 8. 隐含波动率

在 Black–Scholes 下，收益率的条件均值和方差不随时间变化：

$$E\Big[\frac{dS_t}{S_t}\,\Big\vert\,\mathcal F_t\Big]=\mu\,dt,\qquad\mathrm{Var}\Big[\frac{dS_t}{S_t}\,\Big\vert\,\mathcal F_t\Big]=\sigma^2\,d[W,W]_t=\sigma^2\,dt .$$

现实中波动率不是常数。估计它有两种办法：

- 用历史收益率算**历史波动率**（向后看）；
- 用期权价格反推**隐含波动率**（向前看）。

对到期日 $$T$$、执行价 $$K$$ 的欧式看涨期权，记 $$C_{BS}(\sigma,K,T)$$ 为 Black–Scholes 价格，$$C_m(K,T)$$ 为市场价格。**隐含波动率** $$\sigma_{imp}=\sigma_{imp}(K,T)$$ 由下式定义：

$$C_{BS}(\sigma_{imp},K,T)=C_m(K,T).$$

因为 Vega $$>0$$，$$C_{BS}$$ 关于 $$\sigma$$ 单调递增，所以解唯一。期权通常**用隐含波动率报价**，而不是直接报价格。如果世界完全符合 Black–Scholes，对所有 $$K,T$$ 都有 $$\sigma_{imp}(K,T)=\sigma$$。实际数据中 $$\sigma_{imp}$$ 随执行价和期限变化——这就是**波动率微笑 / 偏斜**，是市场在告诉我们这个模型只是近似。

---

**本系列目录**

1. [01｜概率空间、信息流、条件期望与鞅]({{ site.baseurl }}{% post_url 2026-10-21-mf790-01-foundations %})
2. [02｜随机游走、二次变差与布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-02-brownian-motion %})
3. [03｜随机积分]({{ site.baseurl }}{% post_url 2026-10-21-mf790-03-stochastic-integration %})
4. [04｜Itô 公式、Itô 过程与几何布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-04-ito-formula %})
5. **05｜Black–Scholes：复制、偏微分方程与希腊字母**（本篇）
6. [06｜风险中性定价、测度变换与 Girsanov 定理]({{ site.baseurl }}{% post_url 2026-10-21-mf790-06-risk-neutral-pricing %})
7. [07｜随机微分方程、Feynman–Kac 与停时]({{ site.baseurl }}{% post_url 2026-10-21-mf790-07-sde-feynman-kac %})

[← 上一篇：Itô 公式、Itô 过程与几何布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-04-ito-formula %})　·　[下一篇：风险中性定价、测度变换与 Girsanov 定理 →]({{ site.baseurl }}{% post_url 2026-10-21-mf790-06-risk-neutral-pricing %})

*本系列是我在波士顿大学 Questrom 商学院修读 QST MF 790 Stochastic Calculus（随机分析）时整理的课程笔记，共 7 篇。内容按课程顺序整理，表述、直观解释和几处勘误是我自己加的。文中 $$N(\mu,\sigma^2)$$ 的第二个参数一律为方差。如有错漏，欢迎指出。本篇英文版：[English version](https://bobothebest.github.io/homepage_en/2026/10/21/mf790-05-black-scholes/)。*
