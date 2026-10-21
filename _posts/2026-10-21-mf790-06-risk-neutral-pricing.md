---
layout: post
title: "随机过程笔记 06｜风险中性定价、测度变换与 Girsanov 定理"
author: 张梓源
tags:
- 随机过程笔记
- BU MF 790
- 随机分析
- 风险中性定价
- 测度变换
- Girsanov 定理
- 读书笔记
date: 2026-10-21 10:06 -0400
toc: true
math: true
---
> **导读**
>
> 第 05 篇给看涨期权定了价。这一篇用一个公式给**任意**收益定价：
>
> $$V_t=E^{\mathbb Q}\big[e^{-r(T-t)}V_T\mid\mathcal F_t\big]$$。
>
> 价格是**净现值——但要在另一个概率测度 $$\mathbb Q$$ 下算**，在 $$\mathbb Q$$ 下股票只赚无风险利率。然后回答几个自然的问题：测度变换是什么？为什么 $$\mathbb P$$ 和 $$\mathbb Q$$ 对什么*可能*发生意见一致，对*多大可能*意见不同？Girsanov 定理又是怎么通过平移布朗运动把 $$\mu$$ 变成 $$r$$ 的？

---

## 1. 从看涨期权到任意未定权益

**回顾 Black–Scholes 设定。** $$\mathbb P$$ 下 $$S\sim GBM(\mu,\sigma)$$；货币市场利率为常数 $$r$$；自融资策略 $$\Delta$$。贴现财富：

$$d(e^{-rt}X_t)=e^{-rt}\big((\mu-r)\Delta_tS_t\,dt+\sigma\Delta_tS_t\,dW_t\big).$$

对看涨期权我们得到：令 $$\widetilde{\mathbb P}$$ 下 $$\widetilde S\sim GBM(r,\sigma)$$，

$$c(t,s)=E^{\widetilde{\mathbb P}}\big[e^{-r(T-t)}(\widetilde S_T-K)^+\mid\widetilde{\mathcal F}_t,\ \widetilde S_t=s\big],$$

初始资金 $$x=c(0,S_0)$$，策略 $$\Delta_t=c_s(t,S_t)$$，终值财富 $$X_T=(S_T-K)^+$$，由无套利，$$c(0,S_0)$$ 就是价格。现在推广。

**未定权益。** 在 $$T$$ 时支付 $$V_T$$，$$V_T$$ 可以是任何足够可积的 $$\mathcal F_T$$ 可测随机变量：

| 权益 | 收益 $$V_T$$ |
|---|---|
| 看涨期权 | $$(S_T-K)^+$$ |
| 看跌期权 | $$(K-S_T)^+$$ |
| 亚式看涨 | $$\Big(\frac1T\int_0^TS_t\,dt-K\Big)^+$$ |
| 回望期权 | $$\max_{t\le T}S_t-S_T$$ |
| 向上敲出看涨 | $$(S_T-K)^+\,\mathbf 1_{\{\max_{t\le T}S_t\le L\}}$$ |

## 2. 两个假设与定价公式

1. **风险中性测度。** $$(\Omega,\mathcal F)$$ 上存在测度 $$\mathbb Q$$，使 $$\mathbb Q$$ 下 $$S\sim GBM(r,\sigma)$$。等价地，存在 $$\mathbb Q$$-布朗运动 $$W^{\mathbb Q}$$ 使

   $$\frac{dS_t}{S_t}=\underbrace{\mu\,dt+\sigma\,dW_t}_{\mathbb P\text{ 下}}=\underbrace{r\,dt+\sigma\,dW^{\mathbb Q}_t}_{\mathbb Q\text{ 下}}.$$

2. **可复制。** 存在 $$x$$ 和 $$\Delta$$ 使 $$\mathbb P[X^{x,\Delta}_T=V_T]=\mathbb Q[X^{x,\Delta}_T=V_T]=1$$。

**贴现财富是 $$\mathbb Q$$-鞅。** 在 $$\mathbb Q$$ 下漂移项消失：

$$d(e^{-rt}X_t)=e^{-rt}\sigma\Delta_tS_t\,dW^{\mathbb Q}_t .$$

所以 $$\mathbb Q$$ 也叫**鞅测度**。

**在 $$[t,T]$$ 上复制。** 无套利论证从任何起始时间都成立：$$t$$ 时持有 $$X_t$$（$$t$$ 时已知），在 $$[t,T]$$ 上自融资交易，最后得到 $$V_T$$，那么 $$t$$ 时的价格必须是 $$X_t$$。

**确定 $$V_t$$。** 由可复制，$$e^{-rT}X_T=e^{-rT}V_T$$。由鞅性质，

$$e^{-rt}X_t=E^{\mathbb Q}\big[e^{-rT}X_T\mid\mathcal F_t\big]=E^{\mathbb Q}\big[e^{-rT}V_T\mid\mathcal F_t\big].$$

又因为 $$V_T=X_T=X_t+\int_t^TdX_u$$，$$t$$ 时手里有 $$X_t$$ 恰好足以交付 $$V_T$$，所以

$$\boxed{\ V_t=X_t=E^{\mathbb Q}\big[e^{-r(T-t)}V_T\mid\mathcal F_t\big]\ }$$

**怎么读。** 推理链是：（i）复制组合*就是*这份权益，所以权益价格等于组合价值；（ii）$$\mathbb Q$$ 下每个贴现自融资组合都是公平赌局，所以今天的贴现价值等于终值贴现后的期望；（iii）终值就是收益。价格是净现值——**不是在真实测度下**，而是在风险资产和无风险资产期望收益相同的测度下。

这两个假设引出五个问题：$$\mathbb Q$$ 怎么得到？$$\mathbb P$$ 和 $$\mathbb Q$$ 什么关系？怎么在两者之间切换？为什么 $$\mathbb Q$$ 下 $$S\sim GBM(r,\sigma)$$？$$W^{\mathbb Q}$$ 从哪来？本篇剩下的部分逐一回答。

## 3. 换测度会改变什么

**不受测度影响的**（只涉及集合和映射）：

1. $$\Omega$$，结果的集合；
2. $$\mathcal F$$，信息（$$\Omega$$ 的子集族）；
3. $$\mathbb F$$，信息流；
4. 随机变量 $$X:\Omega\to\mathbb R$$——只是一个映射；
5. 过程 $$\{X_t\}$$；
6. 适应性：$$\{\omega:X_t(\omega)\in E\}\in\mathcal F_t$$。

**受测度影响的**（涉及可能性）：

1. 概率 $$P(A)=E^P[\mathbf 1_A]$$；
2. 期望 $$E\big[\sum_i\alpha_i\mathbf 1_{A_i}\big]=\sum_i\alpha_iP(A_i)$$；
3. 矩母函数 $$\alpha\mapsto E^P[e^{\alpha X}]$$；
4. 分布 $$P(X\le t)$$；
5. 鞅性质 $$E[X_t\mid\mathcal F_s]=X_s$$；
6. 独立性 $$P(A\cap B)=P(A)P(B)$$。

所以股价路径还是那些路径，变的只是我们给它们的**权重**——漂移也随之改变。

## 4. 测度变换的构造

给定 $$(\Omega,\mathcal F,\mathbb P)$$，任取随机变量 $$Z$$，$$\mathbb P[Z>0]=1$$（记作 $$Z>0$$）且 $$E^{\mathbb P}[Z]=1$$。定义

$$\mathbb Q(A):=E^{\mathbb P}[Z\,\mathbf 1_A],\qquad A\in\mathcal F .$$

$$\mathbb Q$$ 是概率测度：$$\mathbb Q(A)\ge0$$；$$\mathbb Q(\Omega)=E^{\mathbb P}[Z]=1$$；可列可加，比如 $$A,B$$ 不交时 $$\mathbb Q(A\cup B)=E^{\mathbb P}[(\mathbf 1_A+\mathbf 1_B)Z]=\mathbb Q(A)+\mathbb Q(B)$$。

**例：公平硬币变不公平硬币。** 抛两次，$$\mathbb P$$ 下公平且独立。令

$$Z(HH)=\tfrac49,\quad Z(HT)=Z(TH)=\tfrac89,\quad Z(TT)=\tfrac{16}9,\qquad E^{\mathbb P}Z=\tfrac14\Big(\tfrac49+\tfrac89+\tfrac89+\tfrac{16}9\Big)=1 .$$

则 $$\mathbb Q(HH)=\frac14\cdot\frac49=\frac19$$，$$\mathbb Q(HT)=\mathbb Q(TH)=\frac29$$，$$\mathbb Q(TT)=\frac49$$。$$\mathbb Q$$ 下两次抛掷仍然独立（$$\frac19=\frac13\cdot\frac13$$ 等），但硬币**不公平**了：$$\mathbb Q(H)=\frac13$$。

**等价。** 因为 $$Z>0$$：$$\mathbb P(A)>0\Rightarrow\mathbb Q(A)=E^{\mathbb P}[Z\mathbf 1_A]>0$$，反之亦然。所以 $$\mathbb P(A)>0\iff\mathbb Q(A)>0$$，称 $$\mathbb P,\mathbb Q$$ **等价**。等价不是相等：它们对什么能发生、什么不能发生意见一致，对可能性大小不必一致（$$\mathbb P(A)=0.01$$、$$\mathbb Q(A)=0.99$$ 完全可以）。

**Radon–Nikodym 导数。** $$\dfrac{d\mathbb Q}{d\mathbb P}:=Z$$，也叫**似然比**。在有限空间上

$$\frac{d\mathbb Q}{d\mathbb P}(\omega)=\frac{\mathbb Q(\{\omega\})}{\mathbb P(\{\omega\})}=Z(\omega).$$

## 5. 平移正态分布的均值

$$\mathbb P$$ 下 $$X\sim N(0,\sigma^2)$$，$$\theta\in\mathbb R$$，$$Z=e^{-\theta X-\frac12\sigma^2\theta^2}$$。则 $$Z>0$$，$$E[Z]=e^{-\frac12\sigma^2\theta^2}E[e^{-\theta X}]=1$$。在 $$\mathbb Q(A)=E^{\mathbb P}[Z\mathbf 1_A]$$ 下：

$$E^{\mathbb Q}[e^{\alpha X}]=E^{\mathbb P}\big[e^{(\alpha-\theta)X}\big]e^{-\frac12\sigma^2\theta^2}=e^{\frac12\sigma^2(\alpha-\theta)^2-\frac12\sigma^2\theta^2}=e^{-\sigma^2\theta\alpha+\frac12\sigma^2\alpha^2},$$

所以 $$X\overset{\mathbb Q}{\sim}N(-\sigma^2\theta,\ \sigma^2)$$，即 $$X+\sigma^2\theta\overset{\mathbb Q}{\sim}N(0,\sigma^2)$$。**乘一个指数权重，均值平移，方差不变。**

**一般情形。** 要把 $$\mathbb P$$ 下的 $$N(\mu_P,\sigma^2)$$ 变成 $$\mathbb Q$$ 下的 $$N(\mu_Q,\sigma^2)$$，取 $$\frac{d\mathbb Q}{d\mathbb P}=\frac{e^{-\theta X}}{E^{\mathbb P}[e^{-\theta X}]}=e^{-\theta X-\frac12\sigma^2\theta^2+\theta\mu_P}$$。则

$$E^{\mathbb Q}[e^{\alpha X}]=\exp\Big\{\alpha(\mu_P-\theta\sigma^2)+\tfrac12\alpha^2\sigma^2\Big\},\qquad\text{所以需要}\quad \theta=\frac{\mu_P-\mu_Q}{\sigma^2}.$$

## 6. 用到布朗运动和股价上

固定 $$T$$；$$\mathbb P$$ 下 $$W_T\sim N(0,T)$$。取 $$\theta=\frac{\mu-r}{\sigma}$$，

$$Z=e^{-\theta W_T-\frac12\theta^2T}.$$

由第 5 节（$$\sigma^2=T$$，$$\mu_P=0$$），$$\mathbb Q$$ 下 $$W^{\mathbb Q}_T:=W_T+\theta T\sim N(0,T)$$。因为 $$\mu-\sigma\theta=r$$，

$$S_T=S_0\exp\Big\{\big(\mu-\tfrac12\sigma^2\big)T+\sigma W_T\Big\}=S_0\exp\Big\{\big(r-\tfrac12\sigma^2\big)T+\sigma W^{\mathbb Q}_T\Big\}.$$

所以**在固定的时刻 $$T$$**，$$\mathbb Q$$ 下 $$S_T$$ 看起来像 $$GBM(r,\sigma)$$。要让 $$\mathbb Q$$ 下 $$S\sim GBM(r,\sigma)$$ **对所有** $$t\le T$$ 成立，还需要更多：$$t\mapsto W^{\mathbb Q}_t=W_t+\theta t$$ 是不是 $$[0,T]$$ 上的 $$\mathbb Q$$-布朗运动？等价地问

$$\frac{dS_t}{S_t}=\mu\,dt+\sigma\,dW_t\overset{?}{=}r\,dt+\sigma\,dW^{\mathbb Q}_t,\qquad W^{\mathbb Q}_t=W_t+\frac{\mu-r}{\sigma}t ?$$

核心就是验证 $$W^{\mathbb Q}$$ 是布朗运动。（注意 $$Z_t=E^{\mathbb P}[Z_T\mid\mathcal F_t]=e^{-\theta W_t-\frac12\theta^2t}$$ 是 $$\mathbb P$$-鞅——第 04 篇的指数鞅。）

## 7. Girsanov 定理

**简单版。** 若 $$W$$ 是 $$\mathbb P$$-布朗运动，$$Z_T=e^{-\theta W_T-\frac12\theta^2T}$$，$$\mathbb Q(A)=E^{\mathbb P}[\mathbf 1_AZ_T]$$，则

$$t\mapsto W^{\mathbb Q}_t:=W_t+\theta t\ \text{ 是 }[0,T]\text{ 上的 }\mathbb Q\text{-布朗运动}.$$

**一般版。** 固定 $$(\widetilde\Omega,\widetilde{\mathcal F},\widetilde{\mathbb F},\widetilde{\mathbb P})$$ 和 $$\widetilde{\mathbb P}$$-布朗运动 $$\widetilde W$$，以及适应的 $$\widetilde\theta$$。令

$$\widetilde Z_t=\exp\Big(-\int_0^t\widetilde\theta_u\,d\widetilde W_u-\frac12\int_0^t\widetilde\theta_u^2\,du\Big),\quad t\le T,\qquad \frac{d\widetilde{\mathbb Q}}{d\widetilde{\mathbb P}}=\widetilde Z_T\ \text{于 }\widetilde{\mathcal F}_T .$$

则 $$\widetilde W^{\widetilde{\mathbb Q}}_t=\widetilde W_t+\int_0^t\widetilde\theta_u\,du$$（$$t\le T$$）是 $$\widetilde{\mathbb Q}$$-布朗运动。

直观：用 $$Z$$ 给路径重新加权，让往下漂的路径更可能、往上漂的路径更不可能，调整的幅度正好把"布朗运动加漂移 $$\theta$$"变回无漂移的布朗运动。

### Lévy 刻画

直接从定义证明某个过程是布朗运动很麻烦。Lévy 刻画是捷径：**若 $$\widetilde M$$ 是路径连续的 $$\widetilde{\mathbb Q}$$-鞅，$$\widetilde M_0=0$$ 且 $$[\widetilde M,\widetilde M]_t=t$$，则 $$\widetilde M$$ 是 $$\widetilde{\mathbb Q}$$-布朗运动。**

所以证明 Girsanov 只需要：（1）$$[\widetilde W^{\widetilde{\mathbb Q}},\widetilde W^{\widetilde{\mathbb Q}}]_t=t$$；（2）$$\widetilde W^{\widetilde{\mathbb Q}}$$ 是 $$\widetilde{\mathbb Q}$$-鞅。

<details markdown="1"><summary>Girsanov 的证明（θ 为常数）</summary>

**第 1 步。** 二次变差不依赖测度，漂移也没有贡献：

$$d[\widetilde W^{\widetilde{\mathbb Q}},\widetilde W^{\widetilde{\mathbb Q}}]_t=(d\widetilde W_t+\widetilde\theta\,dt)^2=(d\widetilde W_t)^2=dt .$$

**第 2 步。** 对 $$s<t$$ 和 $$A_s\in\widetilde{\mathcal F}_s$$，换到 $$\widetilde{\mathbb P}$$ 下，用 Tower 和鞅 $$\widetilde Z_t=E^{\widetilde{\mathbb P}}[\widetilde Z_T\mid\widetilde{\mathcal F}_t]$$：

$$E^{\widetilde{\mathbb Q}}[\mathbf 1_{A_s}\widetilde W^{\widetilde{\mathbb Q}}_t]=E^{\widetilde{\mathbb P}}[\mathbf 1_{A_s}\widetilde W^{\widetilde{\mathbb Q}}_t\widetilde Z_T]=E^{\widetilde{\mathbb P}}[\mathbf 1_{A_s}\widetilde W^{\widetilde{\mathbb Q}}_t\widetilde Z_t]=E^{\widetilde{\mathbb P}}\big[\mathbf 1_{A_s}f(t,\widetilde W_t)\big],$$

其中 $$f(t,w)=(w+\widetilde\theta t)\,e^{-\widetilde\theta w-\frac12\widetilde\theta^2t}$$。直接求导：

$$f_t=e^{-\widetilde\theta w-\frac12\widetilde\theta^2t}\,\widetilde\theta\big(1-\tfrac12\widetilde\theta w-\tfrac12\widetilde\theta^2t\big),\quad f_w=e^{-\widetilde\theta w-\frac12\widetilde\theta^2t}\big(1-\widetilde\theta w-\widetilde\theta^2t\big),\quad f_{ww}=e^{-\widetilde\theta w-\frac12\widetilde\theta^2t}\big(-2\widetilde\theta+\widetilde\theta^2w+\widetilde\theta^3t\big),$$

所以 $$f_t+\frac12f_{ww}=0$$，由 Itô $$df(t,\widetilde W_t)=f_w\,d\widetilde W_t$$：$$f(t,\widetilde W_t)$$ 是 $$\widetilde{\mathbb P}$$-鞅。于是

$$E^{\widetilde{\mathbb Q}}[\mathbf 1_{A_s}\widetilde W^{\widetilde{\mathbb Q}}_t]=E^{\widetilde{\mathbb P}}[\mathbf 1_{A_s}f(s,\widetilde W_s)]=E^{\widetilde{\mathbb P}}\big[\mathbf 1_{A_s}(\widetilde W_s+\widetilde\theta s)e^{-\widetilde\theta\widetilde W_s-\frac12\widetilde\theta^2s}\big].$$

对 $$\widetilde W^{\widetilde{\mathbb Q}}_s$$ 做同样的 Tower/TOWK 计算，得到完全相同的右边，所以

$$E^{\widetilde{\mathbb Q}}[\mathbf 1_{A_s}\widetilde W^{\widetilde{\mathbb Q}}_t]=E^{\widetilde{\mathbb Q}}[\mathbf 1_{A_s}\widetilde W^{\widetilde{\mathbb Q}}_s]\quad\text{对所有 }A_s\in\widetilde{\mathcal F}_s .$$

$$\widetilde W^{\widetilde{\mathbb Q}}_s$$ 是 $$\widetilde{\mathcal F}_s$$ 可测的，由部分平均得 $$E^{\widetilde{\mathbb Q}}[\widetilde W^{\widetilde{\mathbb Q}}_t\mid\widetilde{\mathcal F}_s]=\widetilde W^{\widetilde{\mathbb Q}}_s$$。$$\widetilde\theta_t$$ 随时间变化时思路相同，只是记账多一些。

</details>

### $$\mathbb Q$$ 下的价格动态

由 Girsanov，$$W^{\mathbb Q}_t=W_t+\theta t$$ 是 $$[0,T]$$ 上的 $$\mathbb Q$$-布朗运动；又 $$\mu-\sigma\theta=r$$，

$$dS_t=\mu S_t\,dt+\sigma S_t\big(dW^{\mathbb Q}_t-\theta\,dt\big)=rS_t\,dt+\sigma S_t\,dW^{\mathbb Q}_t\iff S\overset{\mathbb Q}{\sim}GBM(r,\sigma).$$

$$\theta=\dfrac{\mu-r}{\sigma}$$ 是**风险的市场价格**：每单位波动率的超额收益，也就是**夏普比率**。这就回答了第 2 节关于假设 1 的五个问题。

## 8. 假设 2：找对冲策略

先说好消息：$$\frac{d\mathbb Q}{d\mathbb P}\big\vert_{\mathcal F_T}=Z_T>0$$，所以 $$\mathbb P$$ 和 $$\mathbb Q$$ 等价，**只需证明 $$\mathbb Q[X_T=V_T]=1$$**。

两个事实：

1. $$e^{-rt}X_t$$ 是 $$\mathbb Q$$-鞅：$$d(e^{-rt}X_t)=e^{-rt}\Delta_tS_t\sigma\,dW^{\mathbb Q}_t$$。
2. $$M^V_t:=E^{\mathbb Q}[e^{-rT}V_T\mid\mathcal F_t]$$ 是 $$[0,T]$$ 上的 $$\mathbb Q$$-鞅，由 Tower：$$E^{\mathbb Q}[M^V_t\mid\mathcal F_s]=E^{\mathbb Q}[e^{-rT}V_T\mid\mathcal F_s]=M^V_s$$。

我们想让这两个鞅成为同一个过程。**假设** $$M^V$$ 能写成随机积分，$$dM^V_t=\Gamma_t\,dW^{\mathbb Q}_t$$，$$\Gamma$$ 适应。那么取

$$\Delta_t=\frac{\Gamma_t}{e^{-rt}\sigma S_t},\qquad t\le T,$$

就有 $$d(e^{-rt}X_t)=e^{-rt}\Delta_tS_t\sigma\,dW^{\mathbb Q}_t=\Gamma_t\,dW^{\mathbb Q}_t=dM^V_t$$。两者都从 $$x=M^V_0=E^{\mathbb Q}[e^{-rT}V_T]$$ 出发，"起点相同 + 动态相同"给出对所有 $$t$$ 都有 $$e^{-rt}X_t=M^V_t$$，到 $$T$$ 时 $$X_T=V_T$$。

剩下的问题——*是不是每个 $$\mathbb Q$$-鞅都能写成 $$\int\Gamma\,dW^{\mathbb Q}$$？*——由**鞅表示定理**回答（信息流由布朗运动生成时，答案是肯定的）。换句话说：在布朗运动驱动的世界里，每一个公平赌局都是一种股票交易策略，所以每一份权益都能被对冲——市场是**完备**的。

## 9. 二次协变差与分部积分

对连续的 $$X,Y$$ 和 $$[0,t]$$ 的划分 $$\Pi$$，

$$[X,Y]^\Pi_t=\sum_{j}\big(X_{t_j}-X_{t_{j-1}}\big)\big(Y_{t_j}-Y_{t_{j-1}}\big),\qquad [X,Y]_t=\lim_{\Vert \Pi\Vert \to0}[X,Y]^\Pi_t .$$

对同一个 $$W$$ 驱动的 Itô 过程 $$dX_t=\Theta_tdt+\Delta_tdW_t$$，$$dY_t=\phi_tdt+\psi_tdW_t$$，启发式 $$d[X,Y]_t=dX_t\,dY_t$$ 给出

$$d[X,Y]_t=\Delta_t\psi_t\,dt,\qquad [X,Y]_t=\int_0^t\Delta_s\psi_s\,ds,$$

$$\Delta\equiv0$$ 或 $$\psi\equiv0$$ 时 $$[X,Y]\equiv0$$。

<details markdown="1"><summary>用极化恒等式证明</summary>

对任意数有 $$4(a-b)(c-d)=\big((a+c)-(b+d)\big)^2-\big((a-c)-(b-d)\big)^2$$，所以

$$[X,Y]=\tfrac14\big([X+Y,X+Y]-[X-Y,X-Y]\big).$$

$$X\pm Y$$ 是 Itô 过程，$$dW$$ 系数为 $$\Delta_t\pm\psi_t$$，所以

$$d[X,Y]_t=\tfrac14\big((\Delta_t+\psi_t)^2-(\Delta_t-\psi_t)^2\big)dt=\Delta_t\psi_t\,dt .$$

</details>

**随机分部积分。**

$$d(X_tY_t)=X_t\,dY_t+Y_t\,dX_t+d[X,Y]_t=\big(X_t\phi_t+Y_t\Theta_t+\Delta_t\psi_t\big)dt+\big(X_t\psi_t+Y_t\Delta_t\big)dW_t,$$

仍是 Itô 过程。对光滑（没有 $$dW$$ 项）的 $$f,g$$，括号项为零，回到普通的乘积法则 $$d(f_tg_t)=f_t\,dg_t+g_t\,df_t$$。

---

**本系列目录**

1. [01｜概率空间、信息流、条件期望与鞅]({{ site.baseurl }}{% post_url 2026-10-21-mf790-01-foundations %})
2. [02｜随机游走、二次变差与布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-02-brownian-motion %})
3. [03｜随机积分]({{ site.baseurl }}{% post_url 2026-10-21-mf790-03-stochastic-integration %})
4. [04｜Itô 公式、Itô 过程与几何布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-04-ito-formula %})
5. [05｜Black–Scholes：复制、偏微分方程与希腊字母]({{ site.baseurl }}{% post_url 2026-10-21-mf790-05-black-scholes %})
6. **06｜风险中性定价、测度变换与 Girsanov 定理**（本篇）
7. [07｜随机微分方程、Feynman–Kac 与停时]({{ site.baseurl }}{% post_url 2026-10-21-mf790-07-sde-feynman-kac %})

[← 上一篇：Black–Scholes：复制、偏微分方程与希腊字母]({{ site.baseurl }}{% post_url 2026-10-21-mf790-05-black-scholes %})　·　[下一篇：随机微分方程、Feynman–Kac 与停时 →]({{ site.baseurl }}{% post_url 2026-10-21-mf790-07-sde-feynman-kac %})

*本系列是我在波士顿大学 Questrom 商学院修读 QST MF 790 Stochastic Calculus（随机分析）时整理的课程笔记，共 7 篇。内容按课程顺序整理，表述、直观解释和几处勘误是我自己加的。文中 $$N(\mu,\sigma^2)$$ 的第二个参数一律为方差。如有错漏，欢迎指出。本篇英文版：[English version](https://bobothebest.github.io/homepage_en/2026/10/21/mf790-06-risk-neutral-pricing/)。*
