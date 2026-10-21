---
layout: post
title: "随机过程笔记 07｜随机微分方程、Feynman–Kac 与停时"
author: 张梓源
tags:
- 随机过程笔记
- BU MF 790
- 随机分析
- 随机微分方程
- Feynman-Kac
- 停时
- 读书笔记
date: 2026-10-21 10:07 -0400
toc: true
math: true
---
> **导读**
>
> 最后一篇把视角拉远。股价、利率、波动率都被建模成**随机微分方程**的解。我们显式解出 Ornstein–Uhlenbeck 方程，认识 CIR 过程，并说明为什么所有 SDE 的解都是**马尔可夫**的。然后 **Feynman–Kac** 把整门课首尾接上：收益的条件期望和一个偏微分方程是同一个东西——Black–Scholes 只是其中一例。最后讲派息股票、停时和可选抽样定理。

---

## 1. 随机微分方程

固定 $$(\Omega,\mathcal F,\mathbb F,\mathbb P)$$、布朗运动 $$W$$，$$\mathbb F=\mathbb F^W$$。**SDE** 是如下形式的方程：

$$dX_t=\mu(t,X_t)\,dt+\sigma(t,X_t)\,dW_t,\qquad X_0=x,$$

$$\mu:[0,\infty)\times\mathbb R\to\mathbb R$$ 是**漂移**函数，$$\sigma:[0,\infty)\times\mathbb R\to\mathbb R$$ 是**扩散**（离散）函数。和 Itô 过程不同，这里的系数依赖未知的 $$X$$ 本身。

**解**是满足下式的过程 $$X$$：对所有 $$t\ge0$$

$$\mathbb P\Big[X_t=x+\int_0^t\mu(u,X_u)\,du+\int_0^t\sigma(u,X_u)\,dW_u\Big]=1,$$

且所有积分都有定义。

**GBM 就是一个 SDE。** 取 $$\mu(t,x)=\mu x$$，$$\sigma(t,x)=\sigma x$$：$$dX_t=\mu X_tdt+\sigma X_tdW_t$$，解为 $$X_t=x\,e^{(\mu-\frac12\sigma^2)t+\sigma W_t}$$——正是 Black–Scholes 里的股票。

**存在唯一性。** 不是每个 SDE 都有解，有显式解的更少。在温和的条件下可以证明解存在且唯一；关键是 $$\mu$$ 和 $$\sigma$$ 在空间变量上满足 **Lipschitz 条件**：

$$\vert f(x)-f(y)\vert\le K\vert x-y\vert,$$

和常微分方程一样。这门课里我们总是假设解存在。

## 2. Ornstein–Uhlenbeck 过程

$$dX_t=\kappa(\theta-X_t)\,dt+\sigma\,dW_t,\qquad \kappa>0 .$$

$$\theta$$：长期均值；$$\kappa$$：均值回复速度；$$\sigma$$：波动率。$$X$$ 高于 $$\theta$$ 时漂移把它往下拉，反之往上拉——像一根带噪声的弹簧。

**求解。** （1）假设解存在，求出形式；（2）验证。乘以积分因子 $$e^{\kappa t}$$ 再用 Itô（和确定性因子的乘积法则）：

$$d(e^{\kappa t}X_t)=e^{\kappa t}\big(\kappa X_t\,dt+\kappa(\theta-X_t)\,dt+\sigma\,dW_t\big)=e^{\kappa t}\big(\kappa\theta\,dt+\sigma\,dW_t\big).$$

积分：$$e^{\kappa t}X_t=x+\theta(e^{\kappa t}-1)+\int_0^te^{\kappa u}\sigma\,dW_u$$，所以

$$X_t=\theta+e^{-\kappa t}\Big((x-\theta)+\int_0^te^{\kappa u}\sigma\,dW_u\Big).$$

验证：求微分，$$dX_t=-\kappa e^{-\kappa t}\big((x-\theta)+\int_0^te^{\kappa u}\sigma dW_u\big)dt+\sigma\,dW_t=\kappa(\theta-X_t)dt+\sigma dW_t$$。✓

**分布。** 被积函数 $$e^{\kappa u}\sigma$$ 是确定性的，由第 04 篇，随机积分服从正态，方差由 Itô 等距给出：

$$e^{-\kappa t}\int_0^te^{\kappa u}\sigma\,dW_u\sim N\big(0,\sigma_t^2\big),\qquad \sigma_t^2=\sigma^2e^{-2\kappa t}\int_0^te^{2\kappa u}\,du=\frac{\sigma^2}{2\kappa}\big(1-e^{-2\kappa t}\big).$$

因此

$$X_t\sim N\Big(\theta+e^{-\kappa t}(x-\theta),\ \frac{\sigma^2}{2\kappa}\big(1-e^{-2\kappa t}\big)\Big)\ \xrightarrow{t\to\infty}\ N\Big(\theta,\ \frac{\sigma^2}{2\kappa}\Big).$$

起点以指数速度被遗忘；$$\theta$$ 是**长期均值**，$$\sigma^2/(2\kappa)$$ 是**长期方差**。

**时变系数。** 若 $$\kappa(t),\theta(t),\sigma(t)$$ 非随机（可以随 $$t$$ 变，但不随 $$\omega$$ 变），$$dX_t=\kappa(t)\big(\theta(t)-X_t\big)dt+\sigma(t)dW_t$$ 的解为

$$X_t=e^{-\int_0^t\kappa(u)du}\Big(x+\int_0^te^{\int_0^u\kappa(v)dv}\kappa(u)\theta(u)\,du+\int_0^te^{\int_0^u\kappa(v)dv}\sigma(u)\,dW_u\Big).$$

## 3. Cox–Ingersoll–Ross 过程

$$dX_t=\kappa(\theta-X_t)\,dt+\xi\sqrt{X_t}\,dW_t,\qquad X_0=x>0 .$$

均值回复和 OU 一样，但瞬时条件方差是 $$\xi^2X_t\,dt$$：$$X$$ 越接近 0，噪声越小。解必须保持为正，$$\sqrt{X_t}$$ 才有意义；当 $$\kappa\theta\ge\frac12\xi^2$$ 时解严格为正（Feller 条件）。没有显式解，但解存在，$$X_t$$ 的分布可以求出来。它是常用的利率模型和方差模型。

## 4. 马尔可夫性

$$dX_t=\mu(t,X_t)dt+\sigma(t,X_t)dW_t$$ 的解是**马尔可夫**过程：在时刻 $$t$$ 想评估过程的未来，知道 $$X_t$$ 就够了——$$\mathcal F_t$$ 里 $$X_t$$ 以外的信息都无关。实际用起来就是：

$$E[g(X_u)\mid\mathcal F_t]=E[g(X_u)\mid X_t],\qquad E\Big[e^{-\int_t^Tg(X_u)du}\,\Big\vert\,\mathcal F_t\Big]=E\Big[e^{-\int_t^Tg(X_u)du}\,\Big\vert\,X_t\Big].$$

为什么：在时刻 $$t$$，下一步 $$dX_t$$ 由 $$X_t$$（通过 $$\mu$$ 和 $$\sigma$$）和 $$dW_t$$ 决定，而 $$dW_t$$ 与 $$\mathcal F_t$$ 独立。分别展开 $$E[g(X_t+dX_t)\mid\mathcal F_t]$$ 和 $$E[g(X_t+dX_t)\mid X_t]$$，一步一步都相同。

## 5. Feynman–Kac

设 $$dX_t=\mu(t,X_t)dt+\sigma(t,X_t)dW_t$$，$$g$$ 为收益函数，$$r(t,x)$$ 为贴现率。那么

$$f(t,x):=E\Big[e^{-\int_t^Tr(u,X_u)du}\,g(X_T)\,\Big\vert\,X_t=x\Big],\qquad t\le T,\ x\in\mathbb R,$$

与下面的偏微分方程一一对应：

$$\begin{cases}0=\Big(f_t+\mu f_x+\tfrac12\sigma^2f_{xx}-rf\Big)(t,x),&t<T,\ x\in\mathbb R,\\[4pt] f(T,x)=g(x).\end{cases}$$

不贴现时，$$u(t,x)=E[g(X_T)\mid X_t=x]$$ 满足 $$0=u_t+\mu u_x+\frac12\sigma^2u_{xx}$$，$$u(T,x)=g(x)$$。

**为什么成立（一句话）。** 由马尔可夫性，$$e^{-\int_0^tr\,du}f(t,X_t)=E\big[e^{-\int_0^Tr\,du}g(X_T)\mid\mathcal F_t\big]$$ 是鞅。对它用 Itô：$$dt$$ 系数是 $$e^{-\int_0^tr}\big(f_t+\mu f_x+\frac12\sigma^2f_{xx}-rf\big)$$，而鞅没有 $$dt$$ 项。

**Black–Scholes 是特例。** 取 $$\mathbb Q$$ 下的 $$X=S$$，即 $$\mu(t,s)=rs$$，$$\sigma(t,s)=\sigma s$$，$$r$$ 为常数，$$g(s)=(s-K)^+$$。偏微分方程就是第 05 篇的 Black–Scholes 方程，期望就是第 06 篇的风险中性公式。整门课的两半——偏微分方程和期望——是同一件事。

## 6. 派息股票

若股票以速率 $$a_t$$ 连续派息，股价会因派息而下降：

$$\frac{dS_t}{S_t}=\mu_t\,dt+\sigma_t\,dW_t-a_t\,dt .$$

**财富。** 在 $$[t,t+dt]$$ 上持有 $$\Delta_t$$ 股，还会收到股息 $$\Delta_tS_ta_t\,dt$$：

$$X_{t+dt}=\Delta_tS_{t+dt}+\Delta_tS_ta_t\,dt+(X_t-\Delta_tS_t)\frac{B_{t+dt}}{B_t},$$

$$\begin{aligned}
dX_t&=\Delta_t\,dS_t+\Delta_tS_ta_t\,dt+(X_t-\Delta_tS_t)r_t\,dt\\
&=\Delta_tS_t\big((\mu_t-a_t)dt+\sigma_tdW_t\big)+\Delta_tS_ta_t\,dt+(X_t-\Delta_tS_t)r_t\,dt\\
&=r_tX_t\,dt+\Delta_t\sigma_tS_t\big(dW_t+\theta_t\,dt\big),\qquad\theta_t=\frac{\mu_t-r_t}{\sigma_t}.
\end{aligned}$$

**财富动态和不派息时一样**（收到的股息正好抵消股价的下跌），但**股价动态变了**。所以欧式期权的价格会变，因为股价不同了。在 $$\mathbb Q$$ 下，

$$dX_t=r_tX_t\,dt+\Delta_tS_t\sigma_t\,dW^{\mathbb Q}_t,\qquad V_t=E^{\mathbb Q}\Big[e^{-\int_t^Tr_udu}V_T\,\Big\vert\,\mathcal F_t\Big],$$

而 $$\mathbb Q$$ 下股票的漂移变成 $$r_t-a_t$$。实际计算时，把收益写成 $$S_T=S_t\cdot\frac{S_T}{S_t}$$ 的函数，固定 $$S_t=s$$ 再算——和 Black–Scholes 的推导完全一样。

## 7. 停时

随机时间 $$\tau$$ 是**停时**，如果

$$\{\omega:\tau(\omega)\le t\}\in\mathcal F_t\quad\text{对所有 }t,$$

即在每个时刻，根据当前信息就能判断是否已经停下。例：非随机时间 $$\tau\equiv t_0$$，$$t_0\le t$$ 时 $$\{\tau\le t\}=\Omega$$，否则为 $$\varnothing$$。

**首中时（hitting time）。** $$\tau(\omega)=\inf\{t:X_t(\omega)=\lambda\}$$，$$X$$ 第一次碰到水平 $$\lambda$$ 的时间。对路径连续的适应过程 $$X$$，首中时是停时。（反过来，"$$T$$ 之前最后一次碰到 $$\lambda$$ 的时间"不是停时：要知道它必须看未来。）

**停时处的值与停止过程。** $$X_\tau$$ 是随机变量 $$\omega\mapsto X_{\tau(\omega)}(\omega)$$——比如 $$\tau(\omega)=1$$ 时是 $$X_1(\omega)$$，$$\tau(\omega)=3$$ 时是 $$X_3(\omega)$$。**停止过程** $$X^\tau$$ 为

$$X^\tau_t(\omega):=X_{t\wedge\tau(\omega)}(\omega):$$

到 $$\tau$$ 之前跟着 $$X$$ 走，之后保持不动。

## 8. 可选抽样定理

**定理（OST）。** 设 $$M$$ 是鞅，$$\alpha\le\tau$$ 是两个**有界**停时。则

$$E[M_\tau\mid\mathcal F_\alpha]=M_\alpha .$$

对有界停时 $$\tau$$ 的**推论**：

1. 停止过程 $$t\mapsto M_{t\wedge\tau}$$ 是鞅；
2. $$E[M_\tau]=M_0$$；
3. 若 $$\tau\le T$$，权益在随机时刻 $$\tau$$ 的价格为 $$V_\tau=E^{\mathbb Q}\big[e^{-r(T-\tau)}V_T\mid\mathcal F_\tau\big]$$——就是第 06 篇的定价公式把 $$t$$ 换成 $$\tau$$。

直观：只要决定"什么时候看"的规则不偷看未来，公平的赌局就仍然公平。有界性很重要："赢了就走"也是停时，但它是无界的，用它就可能 $$E[M_\tau]>M_0$$。

---

**本系列目录**

1. [01｜概率空间、信息流、条件期望与鞅]({{ site.baseurl }}{% post_url 2026-10-21-mf790-01-foundations %})
2. [02｜随机游走、二次变差与布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-02-brownian-motion %})
3. [03｜随机积分]({{ site.baseurl }}{% post_url 2026-10-21-mf790-03-stochastic-integration %})
4. [04｜Itô 公式、Itô 过程与几何布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-04-ito-formula %})
5. [05｜Black–Scholes：复制、偏微分方程与希腊字母]({{ site.baseurl }}{% post_url 2026-10-21-mf790-05-black-scholes %})
6. [06｜风险中性定价、测度变换与 Girsanov 定理]({{ site.baseurl }}{% post_url 2026-10-21-mf790-06-risk-neutral-pricing %})
7. **07｜随机微分方程、Feynman–Kac 与停时**（本篇）

[← 上一篇：风险中性定价、测度变换与 Girsanov 定理]({{ site.baseurl }}{% post_url 2026-10-21-mf790-06-risk-neutral-pricing %})

*本系列是我在波士顿大学 Questrom 商学院修读 QST MF 790 Stochastic Calculus（随机分析）时整理的课程笔记，共 7 篇。内容按课程顺序整理，表述、直观解释和几处勘误是我自己加的。文中 $$N(\mu,\sigma^2)$$ 的第二个参数一律为方差。如有错漏，欢迎指出。本篇英文版：[English version](https://bobothebest.github.io/homepage_en/2026/10/21/mf790-07-sde-feynman-kac/)。*
