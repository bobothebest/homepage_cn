---
layout: post
title: "茆书复习笔记 05｜抽样分布、次序统计量与充分统计量"
author: 张梓源
tags:
- 茆书复习笔记
- 数理统计
- 抽样分布
- 次序统计量
- 充分统计量
- 茆诗松
- 考研
- 读书笔记
date: 2026-09-30 10:05 -0400
toc: true
---
> **导读**
>
> 从这一篇开始进入数理统计。样本均值 $$\bar x$$、样本方差 $$s^2$$ 和 $$\chi^2$$、t、F 三大分布，是后面所有区间估计和假设检验的“原材料”。
>
> 其中最核心的一个定理是：**正态样本下，$$\bar x$$ 与 $$s^2$$ 相互独立**（定理 5.4.1）。后面几乎所有枢轴量都建立在它上面，所以一定要会证。

---

## 1. 样本均值与样本方差

**样本方差（P267）**

$$s^2=\frac1{n-1}\sum_{i=1}^n(x_i-\bar x)^2,\qquad s_n^2=\frac1n\sum_{i=1}^n(x_i-\bar x)^2.$$

计算用的恒等式（记住，方便计算）：

$$\sum_{i=1}^n(x_i-\bar x)^2=\sum x_i^2-\frac{\left(\sum x_i\right)^2}{n}=\sum x_i^2-n\bar x^2.$$

配凑证明时常用（习题 5.3.10）：$$\left(\sum_{i=1}^nx_i\right)^2=\sum_{i=1}^nx_i^2+2\sum_{i<j}x_ix_j$$。

**递推公式（习题 5.3.4）**：记 $$\bar x_n=\frac1n\sum_{i=1}^nx_i$$，$$s_n^2=\frac1{n-1}\sum_{i=1}^n(x_i-\bar x_n)^2$$，则

$$s_{n+1}^2=\frac{n-1}{n}s_n^2+\frac1{n+1}(x_{n+1}-\bar x_n)^2.$$

**样本均值的平方偏差最小（P264，定理 5.3.2）**：$$\sum(x_i-\bar x)^2\le\sum(x_i-c)^2$$，对任意常数 $$c$$ 成立。

**定理 5.3.3（P265）：$$\bar x$$ 的分布**
- 总体为 $$N(\mu,\sigma^2)$$ 时：$$\bar x\sim N(\mu,\sigma^2/n)$$（精确成立）
- 总体任意，但方差 $$\sigma^2$$ 存在时：$$n$$ 较大时 $$\bar x\ \dot\sim\ N(\mu,\sigma^2/n)$$（由中心极限定理）

**定理 5.3.4（P268，会证）**：若总体二阶矩存在，则

$$E(\bar x)=\mu,\qquad\mathrm{Var}(\bar x)=\frac{\sigma^2}n,\qquad E(s^2)=\sigma^2.$$

<details markdown="1"><summary>证明</summary>

前两个直接由期望和方差的性质得到。对于 $$s^2$$：

$$E\left[\sum(x_i-\bar x)^2\right]=\sum E(x_i^2)-nE(\bar x^2)=n(\sigma^2+\mu^2)-n\left(\frac{\sigma^2}n+\mu^2\right)=(n-1)\sigma^2.$$

所以 $$E(s^2)=\sigma^2$$。这也是 $$s^2$$ 分母取 $$n-1$$ 的原因：用 $$\bar x$$ 代替 $$\mu$$ 会“用掉”一个自由度。

</details>

---

## 2. 次序统计量

**单个次序统计量的分布（P273，背记）**：总体密度为 $$p(x)$$、分布函数为 $$F(x)$$ 时，

$$p_k(x)=\frac{n!}{(k-1)!\,(n-k)!}\,[F(x)]^{k-1}\,[1-F(x)]^{n-k}\,p(x).$$

特别地：

$$p_{(1)}(x)=n[1-F(x)]^{n-1}p(x),\qquad p_{(n)}(x)=n[F(x)]^{n-1}p(x).$$

> **注意**：密度公式只适用于连续总体，对离散分布不成立。离散时应从分布函数出发：$$F_{(n)}(x)=[F(x)]^n$$，$$F_{(1)}(x)=1-[1-F(x)]^n$$。

<details markdown="1"><summary>记忆方法：“三段分类”</summary>

$$x_{(k)}$$ 落在 $$x$$ 附近，需要：$$k-1$$ 个样本落在 $$x$$ 左边（概率 $$F$$），$$1$$ 个落在 $$x$$ 附近（概率 $$p(x)dx$$），$$n-k$$ 个落在右边（概率 $$1-F$$）。系数就是把 $$n$$ 个样本分成三组的多项式系数 $$\frac{n!}{(k-1)!\,1!\,(n-k)!}$$。

</details>

**多个次序统计量的分布（P274；定理 5.3.6 的应用见习题 5.3.22、5.3.23）**：$$(x_{(i)},x_{(j)})$$，$$i<j$$，联合密度为

$$p_{ij}(y,z)=\frac{n!}{(i-1)!\,(j-i-1)!\,(n-j)!}[F(y)]^{i-1}[F(z)-F(y)]^{j-i-1}[1-F(z)]^{n-j}p(y)p(z),\quad y\le z.$$

全部次序统计量的联合密度为 $$n!\prod_{i=1}^np(y_i)$$，$$y_1<\dots<y_n$$。  
例（习题 5.3.34）：$$X\sim U(0,\theta)$$ 时，$$(x_{(1)},\dots,x_{(n)})$$ 的联合密度为 $$\frac{n!}{\theta^n}$$，$$0<u_1\le\dots\le u_n\le\theta$$。

**均匀分布的次序统计量**
- 例 5.3.8（P274）：总体 $$U(0,1)$$ 时，$$x_{(k)}\sim Be(k,\ n-k+1)$$
- 例 5.3.9（P275）：总体 $$U(0,1)$$ 时，极差 $$x_{(n)}-x_{(1)}\sim Be(n-1,\ 2)$$

**☆ 习题 5.3.28（把任意连续总体化为均匀分布）**：设总体分布函数 $$F$$ 连续，令 $$\eta_i=F(x_{(i)})$$，则
1. $$\eta_1\le\dots\le\eta_n$$ 恰是来自 $$U(0,1)$$ 的次序统计量（因为 $$F(X)\sim U(0,1)$$ 且 $$F$$ 单调）；
2. $$E(\eta_i)=\frac{i}{n+1}$$，$$\mathrm{Var}(\eta_i)=\frac{i(n+1-i)}{(n+1)^2(n+2)}$$（由 $$Be(i,n-i+1)$$ 的期望和方差得到）；
3. $$\eta_i,\eta_j$$（$$i<j$$）的协方差矩阵为

$$\begin{pmatrix}\frac{a_1(1-a_1)}{n+2}&\frac{a_1(1-a_2)}{n+2}\\[2pt]\frac{a_1(1-a_2)}{n+2}&\frac{a_2(1-a_2)}{n+2}\end{pmatrix},\qquad a_1=\frac i{n+1},\ a_2=\frac j{n+1}.$$

**样本 $$p$$ 分位数的渐近分布（定理 5.3.7，P276）**：设总体密度 $$p(x)$$ 在 $$x_p$$ 处连续且 $$p(x_p)>0$$，则

$$m_p\ \dot\sim\ N\left(x_p,\ \frac{p(1-p)}{n\,p^2(x_p)}\right).$$

特别地，样本中位数 $$m_{0.5}\ \dot\sim\ N\left(x_{0.5},\ \frac{1}{4n\,p^2(x_{0.5})}\right)$$。

**中位数具有稳健性（P277）**：少数异常值对中位数几乎没有影响，对均值的影响却可能很大。

---

## 3. 三大抽样分布

| 分布 | 构造 | 密度 | 页码 |
|---|---|---|---|
| $$\chi^2(n)$$ | $$X_1,\dots,X_n\overset{iid}{\sim}N(0,1)$$，$$\sum X_i^2$$ | $$\frac{1}{2^{n/2}\Gamma(n/2)}x^{\frac n2-1}e^{-x/2},\ x>0$$ | P283 |
| $$F(m,n)$$ | $$\dfrac{X/m}{Y/n}$$，$$X\sim\chi^2(m)$$，$$Y\sim\chi^2(n)$$ 独立 | $$\frac{\Gamma(\frac{m+n}2)}{\Gamma(\frac m2)\Gamma(\frac n2)}\left(\frac mn\right)^{\frac m2}x^{\frac m2-1}\left(1+\frac mnx\right)^{-\frac{m+n}2},\ x>0$$ | P286 |
| $$t(n)$$ | $$\dfrac{X}{\sqrt{Y/n}}$$，$$X\sim N(0,1)$$，$$Y\sim\chi^2(n)$$ 独立 | $$\frac{\Gamma(\frac{n+1}2)}{\sqrt{n\pi}\,\Gamma(\frac n2)}\left(1+\frac{x^2}n\right)^{-\frac{n+1}2}$$ | P288 |

- **F 分布的分位数（P287）**：$$F_\alpha(n,m)=\dfrac{1}{F_{1-\alpha}(m,n)}$$（因为 $$F\sim F(m,n)$$ 时 $$\frac1F\sim F(n,m)$$）
- $$t(1)$$ 就是标准柯西分布；$$n\to\infty$$ 时 $$t(n)\to N(0,1)$$；若 $$T\sim t(n)$$，则 $$T^2\sim F(1,n)$$
- $$\chi^2(n)$$：期望 $$n$$，方差 $$2n$$；$$F(m,n)$$、$$t(n)$$ 的期望和方差见第 2 篇

<details markdown="1"><summary>F 分布密度的导出（会导）</summary>

先求 $$Z=X/Y$$ 的密度。由商的密度公式

$$p_Z(z)=\int_0^{\infty}p_X(zy)\,p_Y(y)\,y\,dy=\frac{z^{\frac m2-1}}{2^{\frac{m+n}2}\Gamma(\frac m2)\Gamma(\frac n2)}\int_0^{\infty}y^{\frac{m+n}2-1}e^{-\frac{(1+z)y}{2}}\,dy.$$

后面的积分用伽玛函数算出，等于 $$\Gamma(\frac{m+n}2)\left(\frac{2}{1+z}\right)^{\frac{m+n}2}$$，于是

$$p_Z(z)=\frac{\Gamma(\frac{m+n}2)}{\Gamma(\frac m2)\Gamma(\frac n2)}z^{\frac m2-1}(1+z)^{-\frac{m+n}2}.$$

最后 $$F=\frac nmZ$$，做线性变换即得。$$t$$ 分布的推导类似：先求 $$\sqrt{Y/n}$$ 的密度，再用商的密度公式。

</details>

**☆ 定理 5.4.1（P284，会证）**：设 $$x_1,\dots,x_n$$ 来自 $$N(\mu,\sigma^2)$$，则
1. $$\bar x$$ 与 $$s^2$$ 相互独立；
2. $$\bar x\sim N(\mu,\sigma^2/n)$$；
3. $$\dfrac{(n-1)s^2}{\sigma^2}\sim\chi^2(n-1)$$。

<details markdown="1"><summary>证明要点（正交变换）</summary>

不妨设 $$\mu=0$$。取正交矩阵 $$A$$，第一行为 $$\left(\frac1{\sqrt n},\dots,\frac1{\sqrt n}\right)$$，令 $$\boldsymbol y=A\boldsymbol x$$。
- 正交变换保持独立同正态：$$y_1,\dots,y_n\overset{iid}{\sim}N(0,\sigma^2)$$
- $$y_1=\sqrt n\,\bar x$$
- $$\sum y_i^2=\sum x_i^2$$，所以 $$(n-1)s^2=\sum x_i^2-n\bar x^2=\sum_{i=2}^ny_i^2$$

于是 $$\bar x$$ 只依赖于 $$y_1$$，$$s^2$$ 只依赖于 $$y_2,\dots,y_n$$，二者独立；并且 $$\frac{(n-1)s^2}{\sigma^2}=\sum_{i=2}^n\left(\frac{y_i}\sigma\right)^2\sim\chi^2(n-1)$$。

</details>

**推论（构造枢轴量的基础）**
- $$\dfrac{\sqrt n(\bar x-\mu)}{s}\sim t(n-1)$$
- 两个独立正态样本：$$\dfrac{s_x^2/\sigma_1^2}{s_y^2/\sigma_2^2}\sim F(m-1,n-1)$$
- 若 $$\sigma_1=\sigma_2$$：$$\dfrac{(\bar x-\bar y)-(\mu_1-\mu_2)}{s_w\sqrt{\frac1m+\frac1n}}\sim t(m+n-2)$$，其中 $$s_w^2=\dfrac{(m-1)s_x^2+(n-1)s_y^2}{m+n-2}$$

---

## 4. 充分统计量

**定理 5.5.1 因子分解定理（P297）**：$$T(\boldsymbol x)$$ 是 $$\theta$$ 的充分统计量，当且仅当存在函数 $$g(t;\theta)$$ 和 $$h(\boldsymbol x)$$，使得

$$p(\boldsymbol x;\theta)=g\big(T(\boldsymbol x);\theta\big)\,h(\boldsymbol x).$$

其中 $$h$$ 不含 $$\theta$$，$$g$$ 只通过 $$T$$ 依赖于样本。

**常见分布的充分统计量（指导书 P300）**

| 分布 | 参数 | 充分统计量 |
|---|---|---|
| $$b(1,p)$$ | $$p$$ | $$\sum x_i$$ |
| $$P(\lambda)$$ | $$\lambda$$ | $$\sum x_i$$ |
| $$Ge(\theta)$$ | $$\theta$$ | $$\sum x_i$$ |
| $$Exp(\lambda)$$ | $$\lambda$$ | $$\sum x_i$$ |
| $$U(0,\theta)$$ | $$\theta$$ | $$x_{(n)}$$ |
| $$U(\theta_1,\theta_2)$$ | $$\theta_1,\theta_2$$ | $$(x_{(1)},x_{(n)})$$ |
| $$U(\theta,2\theta)$$ | $$\theta$$ | $$(x_{(1)},x_{(n)})$$ |
| $$N(\mu,\sigma^2)$$ | $$\mu,\sigma^2$$ | $$\left(\bar x,\ \sum(x_i-\bar x)^2\right)$$ |
| 幂分布 $$\theta x^{\theta-1},\ 0<x<1$$ | $$\theta$$ | $$\prod x_i$$ 或 $$\sum\ln x_i$$ |
| 双参数指数 $$\frac1\theta e^{-\frac{x-\mu}\theta},\ x>\mu$$ | $$\mu,\theta$$ | $$\left(x_{(1)},\ \sum x_i\right)$$ |
| $$Ga(\alpha,\lambda)$$ | $$\alpha,\lambda$$ | $$\left(\sum x_i,\ \prod x_i\right)$$ |
| $$LN(\mu,\sigma^2)$$ | $$\mu,\sigma^2$$ | $$\left(\sum\ln x_i,\ \sum(\ln x_i)^2\right)$$ |
| $$Be(a,b)$$ | $$a,b$$ | $$\left(\sum\ln x_i,\ \sum\ln(1-x_i)\right)$$ |

- 充分统计量的一一对应变换仍是充分统计量（指导书 P300），例如 $$\sum x_i$$ 与 $$\bar x$$ 等价。
- 看规律：指数族分布的充分统计量就是密度指数上与参数相乘的那些样本函数之和；支撑依赖参数（如均匀分布）时，充分统计量是最值。

<details markdown="1"><summary>示例：\(U(0,\theta)\)</summary>

$$p(\boldsymbol x;\theta)=\prod_{i=1}^n\frac1\theta I_{\{0<x_i<\theta\}}=\underbrace{\theta^{-n}I_{\{x_{(n)}<\theta\}}}_{g(x_{(n)};\theta)}\cdot\underbrace{I_{\{x_{(1)}>0\}}}_{h(\boldsymbol x)},$$

所以 $$x_{(n)}$$ 是充分统计量。

</details>

---

**本系列目录**

1. [01｜概念简答索引]({% post_url 2026-09-30-mao-notes-01-concepts %})
2. [02｜常用分布]({% post_url 2026-09-30-mao-notes-02-distributions %})
3. [03｜多维分布、协方差与条件期望]({% post_url 2026-09-30-mao-notes-03-multivariate %})
4. [04｜特征函数、大数定律与中心极限定理]({% post_url 2026-09-30-mao-notes-04-limit-theorems %})
5. **05｜抽样分布、次序统计量与充分统计量**（本篇）
6. [06｜参数估计]({% post_url 2026-09-30-mao-notes-06-estimation %})
7. [07｜区间估计、假设检验与习题结论]({% post_url 2026-09-30-mao-notes-07-intervals-tests %})

[← 上一篇：特征函数、大数定律与中心极限定理]({% post_url 2026-09-30-mao-notes-04-limit-theorems %})　·　[下一篇：参数估计 →]({% post_url 2026-09-30-mao-notes-06-estimation %})

*本系列是我复习茆诗松、程依明、濮晓龙《概率论与数理统计教程》时整理的笔记，共 7 篇。文中 P××× 为茆书页码，“指导书”为配套的学习指导与习题解答。如有错漏，欢迎指出。*
