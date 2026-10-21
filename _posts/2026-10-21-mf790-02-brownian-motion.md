---
layout: post
title: "随机过程笔记 02｜随机游走、二次变差与布朗运动"
author: 张梓源
tags:
- 随机过程笔记
- 学习笔记
date: 2026-10-21 10:02 -0400
toc: true
math: true
---
> **导读**
>
> 布朗运动最好理解的方式，是把它看成抛硬币游戏的极限。这一篇从随机游走出发，引入让随机分析和普通微积分分道扬镳的那个量——**二次变差**，然后说明布朗运动继承了这一切：它是鞅，$$W_t^2-t$$ 是鞅，并且 $$[W,W]_t=t$$。
>
> 要记住的结论：**$$(dW_t)^2=dt$$。** 光滑曲线的二次变差是 0，布朗运动的路径不是，Itô 公式多出来的那一项就来自这里。

---

## 1. 随机游走

取无穷次公平抛硬币的空间，$$\omega=(\omega_1,\omega_2,\dots)$$，各次独立。定义

$$Z_j(\omega)=\begin{cases}+1,&\omega_j=H\\-1,&\omega_j=T\end{cases},\qquad M_0=0,\quad M_n=\sum_{j=1}^nZ_j\ (n\ge1).$$

$$Z_j$$ 独立同分布，$$EZ_j=0$$，$$\mathrm{Var}(Z_j)=EZ_j^2=1$$。$$M$$ 就是**随机游走**，$$Z_n=M_n-M_{n-1}$$ 是它的增量。

**性质。**

1. **$$M$$ 是鞅。** $$\mathcal F_n$$ = 前 $$n$$ 次抛掷的信息；$$Z_{n+1}$$ 只依赖第 $$n+1$$ 次，所以 $$Z_{n+1}\perp\!\!\!\perp\mathcal F_n$$：

   $$E[M_{n+1}\mid\mathcal F_n]=E[M_n+Z_{n+1}\mid\mathcal F_n]=M_n+E[Z_{n+1}]=M_n.$$

2. **独立增量。** $$0\le m<n\le p<q$$ 时，$$M_n-M_m=\sum_{j=m+1}^nZ_j$$ 与 $$M_q-M_p=\sum_{j=p+1}^qZ_j$$ 用的是不相交的抛掷，所以独立。
3. **增量的均值和方差。** $$E[M_n-M_m]=0$$，且

   $$\mathrm{Var}(M_n-M_m)=\sum_{i,j=m+1}^nE[Z_iZ_j]=n-m,\qquad E[Z_iZ_j]=\begin{cases}1,&i=j\\0,&i\ne j.\end{cases}$$

## 2. 随机游走的平方

$$M_n^2$$ **不是**鞅——它会往上漂：

$$E[M_{n+1}^2\mid\mathcal F_n]=E[(M_n+Z_{n+1})^2\mid\mathcal F_n]=M_n^2+2M_nE[Z_{n+1}]+E[Z_{n+1}^2]=M_n^2+1.$$

所以 $$M^2$$ 是下鞅；减掉漂移就修好了：$$Y_n:=M_n^2-n$$ 满足 $$E[Y_{n+1}\mid\mathcal F_n]=M_n^2+1-(n+1)=Y_n$$，是**鞅**。

## 3. 二次变差（离散时间）

过程 $$X=\{X_n\}$$ 的**二次变差**是一个过程：

$$[X,X]_0=0,\qquad [X,X]_n=\sum_{j=1}^n\big(X_j-X_{j-1}\big)^2 .$$

它把增量的平方加起来，衡量路径**有多粗糙**。（定义本身不要求 $$X$$ 是鞅。）

- 随机游走：每个增量都是 $$\pm1$$，所以对**每个** $$\omega$$ 都有 $$[M,M]_n=n$$——硬币不公平也一样。这是特例；一般来说 $$[X,X]_n$$ 同时依赖 $$n$$ 和 $$\omega$$。
- $$Y_n=M_n^2-n$$：用 $$M_j=M_{j-1}+Z_j$$ 和 $$Z_j^2=1$$，

  $$[Y,Y]_n=\sum_{j=1}^n\big(Z_j^2+2Z_jM_{j-1}-1\big)^2=\sum_{j=1}^n4M_{j-1}^2,$$

  确实随 $$\omega$$ 变化。

### 鞅与二次变差的两个结论

**结论 1。** 若 $$X$$ 是鞅且 $$X_0=0$$，则 $$\mathrm{Var}(X_n)=E\big[[X,X]_n\big]$$。

<details markdown="1"><summary>证明</summary>

因为 $$EX_n=EX_0=0$$，$$\mathrm{Var}(X_n)=EX_n^2$$。把 $$X_n=\sum_j(X_j-X_{j-1})$$ 平方展开：

$$X_n^2=\sum_{i=1}^n(X_i-X_{i-1})^2+\sum_{i\ne j}(X_i-X_{i-1})(X_j-X_{j-1}).$$

$$i<j$$ 时，对 $$\mathcal F_{j-1}$$ 取条件（Tower：对第二大的时间取条件）：

$$E[(X_i-X_{i-1})(X_j-X_{j-1})]=E\big[(X_i-X_{i-1})\,E[X_j-X_{j-1}\mid\mathcal F_{j-1}]\big]=E[(X_i-X_{i-1})\cdot0]=0.$$

交叉项消失，$$EX_n^2=E[X,X]_n$$。对随机游走就是 $$\mathrm{Var}(M_n)=n=E[M,M]_n$$。

</details>

**结论 2（二次变差是补偿项）。** 若 $$X$$ 是鞅且 $$EX_n^2<\infty$$，则 $$\{X_n^2-[X,X]_n\}$$ 是鞅。

<details markdown="1"><summary>证明</summary>

令 $$Y_n=X_n^2-[X,X]_n$$，则

$$Y_n-Y_{n-1}=X_n^2-X_{n-1}^2-(X_n-X_{n-1})^2=2X_{n-1}(X_n-X_{n-1}),$$

而 $$E[2X_{n-1}(X_n-X_{n-1})\mid\mathcal F_{n-1}]=2X_{n-1}\cdot0=0$$。

</details>

$$[X,X]_n$$ 只增不减，所以 $$X_n^2$$ 本身是下鞅；必须**补偿**（减掉）$$[X,X]$$ 才是鞅，这就是 $$[X,X]$$ 又叫补偿项（compensator）的原因。前面的 $$M_n^2-n$$ 就是 $$[M,M]_n=n$$ 的特例。

## 4. 从随机游走到布朗运动

随机游走和布朗运动的共同点：独立增量，增量均值为 0、方差等于经过的时间。不同点：$$M$$ 是离散时间、整数跳；$$W$$ 是连续时间、取任意实数。

**缩放随机游走。** 时间加速 $$n$$ 倍，跳幅缩到 $$1/\sqrt n$$：

$$W^{(n)}_t:=\frac1{\sqrt n}M_{\lfloor nt\rfloor}.$$

- **独立增量：** 从 $$M$$ 继承。
- **均值和方差对得上：** 取 $$ns,nt$$ 为整数，

  $$E[W^{(n)}_t-W^{(n)}_s]=0,\qquad \mathrm{Var}[W^{(n)}_t-W^{(n)}_s]=\frac1n\,\mathrm{Var}[M_{nt}-M_{ns}]=\frac{nt-ns}{n}=t-s.$$

- **极限是正态（中心极限定理）：**

  $$W^{(n)}_t-W^{(n)}_s=\underbrace{\sqrt{\frac{\lfloor nt\rfloor-\lfloor ns\rfloor}{n}}}_{\to\sqrt{t-s}}\cdot\underbrace{\frac1{\sqrt{\lfloor nt\rfloor-\lfloor ns\rfloor}}\sum_{j=\lfloor ns\rfloor+1}^{\lfloor nt\rfloor}Z_j}_{\to N(0,1)\text{（CLT）}}\ \Longrightarrow\ N(0,t-s).$$

所以 $$W^{(n)}\to W$$：**布朗运动是存在的**，它就是从很远的地方看一场抛硬币游戏。（注意 CLT 里的"$$n$$"是这段时间里的抛掷次数 $$\lfloor nt\rfloor-\lfloor ns\rfloor$$。）

## 5. 布朗运动的性质

回忆第 01 篇的定义：$$W_0=0$$；路径连续；$$W_t-W_s\sim N(0,t-s)$$；独立增量。若 $$W$$ 还适应于 $$\mathbb F$$ 且 $$W_t-W_s\perp\!\!\!\perp\mathcal F_s$$，就说它是**关于 $$\mathbb F$$ 的**布朗运动。

**1. 鞅。** 第 01 篇已证。连续时间下上鞅、下鞅以及 $$EX_t$$ 的单调性和离散时间完全一样。

**2. 协方差。** $$s<t$$ 时，减一项加一项凑出独立增量：

$$\mathrm{Cov}(W_s,W_t)=E[W_sW_t]=E[W_s(W_t-W_s)]+E[W_s^2]=0+s=\min(s,t).$$

所以 $$(W_s,W_t)$$ 联合正态：

$$\begin{pmatrix}W_s\\W_t\end{pmatrix}\sim N\left(\begin{pmatrix}0\\0\end{pmatrix},\begin{pmatrix}s&s\\s&t\end{pmatrix}\right).$$

**3. 条件分布与马尔可夫性。** 用 $$Y\sim N(\mu,\sigma^2)$$ 时 $$E[e^{\alpha Y}]=e^{\alpha\mu+\frac12\alpha^2\sigma^2}$$：

$$E[e^{\alpha W_t}\mid\mathcal F_s]=e^{\alpha W_s}\,E[e^{\alpha(W_t-W_s)}]=e^{\alpha W_s+\frac12\alpha^2(t-s)}\ \Longrightarrow\ W_t\mid\mathcal F_s\sim N(W_s,\ t-s).$$

给定到 $$s$$ 为止的一切，$$W_t$$ 的分布只依赖 $$W_s$$——这就是**马尔可夫性**。**转移密度**为

$$p(\tau,y,x)=\frac1{\sqrt{2\pi\tau}}\,e^{-\frac{(x-y)^2}{2\tau}},\qquad E[h(W_{s+\tau})\mid W_s=y]=\int_{\mathbb R}h(x)\,p(\tau,y,x)\,dx,$$

读作"用时间 $$\tau$$ 从 $$y$$ 走到 $$x$$ 附近"。

**4. $$W_t^2-t$$ 是鞅。** 用 $$W_t^2=(W_t-W_s)^2+2(W_t-W_s)W_s+W_s^2$$：

$$E[W_t^2-t\mid\mathcal F_s]=(t-s)+0+W_s^2-t=W_s^2-s.$$

这是 $$M_n^2-n$$ 的连续时间版本。

## 6. 连续时间的二次变差

固定 $$t$$ 和划分 $$\Pi=\{0=t_0<t_1<\dots<t_n=t\}$$，步长 $$\Vert \Pi\Vert =\max_j(t_j-t_{j-1})$$。定义

$$[X,X]^\Pi_t=\sum_{j=1}^n\big(X_{t_j}-X_{t_{j-1}}\big)^2,\qquad [X,X]_t=\lim_{\Vert \Pi\Vert \to0}[X,X]^\Pi_t .$$

区间上的 $$[X,X]_{[s,t]}$$ 同样定义，且二次变差是累加的：$$[X,X]_t=[X,X]_s+[X,X]_{[s,t]}$$。

**定理。** 以概率 1，$$[W,W]_t=t$$，$$[W,W]_{[s,t]}=t-s$$。

*直观。* 缩放随机游走 $$[W^{(n)},W^{(n)}]_t=\frac1n\lfloor nt\rfloor\in[t-\frac1n,t]$$，而 $$W^{(n)}\to W$$。

<details markdown="1"><summary>证明梗概</summary>

取均匀划分 $$t_j=jt/n$$，则

$$[W,W]^\Pi_t=\frac tn\sum_{j=1}^n\left(\frac{W_{t_j}-W_{t_{j-1}}}{\sqrt{t_j-t_{j-1}}}\right)^2=t\cdot\frac1n\sum_{j=1}^nZ_j^2,$$

$$Z_j$$ 独立同分布 $$N(0,1)$$。由强大数定律 $$\frac1n\sum Z_j^2\to E[Z^2]=1$$，所以 $$[W,W]^\Pi_t\to t$$。一般划分下，证明 $$E\big[[W,W]^\Pi_t\big]=t$$ 且 $$\mathrm{Var}\big[[W,W]^\Pi_t\big]\le K\Vert \Pi\Vert $$（$$K$$ 只依赖 $$t$$），即可推出收敛到 $$t$$。

</details>

结论 2 也推广到连续时间：若 $$X$$ 是连续鞅且 $$EX_t^2<\infty$$，则 $$X_t^2-[X,X]_t$$ 是鞅（证明要用第 04 篇的 Itô 公式）。

## 7. 布朗运动的路径为什么必须粗糙

路径的**一阶变差**是 $$FV(X)_t=\lim_{\Vert \Pi\Vert \to0}\sum_j\vert X_{t_j}-X_{t_{j-1}}\vert$$。对光滑路径，由中值定理

$$FV(X)_t=\lim\sum_j(t_j-t_{j-1})\vert X'(\xi_j)\vert=\int_0^t\vert X'_s\vert\,ds<\infty.$$

**一阶变差有限 ⟹ 二次变差为 0**（路径连续时）：

$$\sum_j(X_{t_j}-X_{t_{j-1}})^2\le\underbrace{\max_j\vert X_{t_j}-X_{t_{j-1}}\vert}_{\to0\ \text{（连续）}}\cdot\underbrace{\sum_j\vert X_{t_j}-X_{t_{j-1}}\vert}_{\to FV(X)_t<\infty}\ \to0.$$

而 $$[W,W]_t=t>0$$，所以布朗运动的路径一阶变差不可能有限；事实上它以概率 1 **处处不可微**。

这种粗糙不是缺点：如果价格路径是光滑的，局部就是可预测的，可预测就意味着套利。

**大家都在用的启发式写法。** 对 $$[W,W]_t=t$$ 求微分得 $$d[W,W]_t=dt$$，人们写成

$$(dW_t)^2=dt,$$

想的是 $$[W,W]_{[t,t+dt]}\approx\sum(W_{t_j}-W_{t_{j-1}})^2$$。开方说"$$dW_t=\sqrt{dt}$$"没有意义，但乘法规则 $$(dW_t)^2=dt$$ 完全正确，第 04 篇整篇都靠它。

---

**本系列目录**

1. [01｜概率空间、信息流、条件期望与鞅]({{ site.baseurl }}{% post_url 2026-10-21-mf790-01-foundations %})
2. **02｜随机游走、二次变差与布朗运动**（本篇）
3. [03｜随机积分]({{ site.baseurl }}{% post_url 2026-10-21-mf790-03-stochastic-integration %})
4. [04｜Itô 公式、Itô 过程与几何布朗运动]({{ site.baseurl }}{% post_url 2026-10-21-mf790-04-ito-formula %})
5. [05｜Black–Scholes：复制、偏微分方程与希腊字母]({{ site.baseurl }}{% post_url 2026-10-21-mf790-05-black-scholes %})
6. [06｜风险中性定价、测度变换与 Girsanov 定理]({{ site.baseurl }}{% post_url 2026-10-21-mf790-06-risk-neutral-pricing %})
7. [07｜随机微分方程、Feynman–Kac 与停时]({{ site.baseurl }}{% post_url 2026-10-21-mf790-07-sde-feynman-kac %})

[← 上一篇：概率空间、信息流、条件期望与鞅]({{ site.baseurl }}{% post_url 2026-10-21-mf790-01-foundations %})　·　[下一篇：随机积分 →]({{ site.baseurl }}{% post_url 2026-10-21-mf790-03-stochastic-integration %})

*本系列是我在波士顿大学 Questrom 商学院修读 QST MF 790 Stochastic Calculus（随机分析）时整理的课程笔记，共 7 篇。内容按课程顺序整理，表述、直观解释和几处勘误是我自己加的。文中 $$N(\mu,\sigma^2)$$ 的第二个参数一律为方差。如有错漏，欢迎指出。本篇英文版：[English version](https://bobothebest.github.io/homepage_en/2026/10/21/mf790-02-brownian-motion/)。*
