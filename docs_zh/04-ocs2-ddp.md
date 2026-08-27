# 04 · ocs2_ddp：SLQ 与 iLQR 的完整推导

这是 OCS2 的招牌模块。本章从最优控制的第一性原理出发，
把连续时间 SLQ 与离散时间 iLQR 的**每一步推导写全**，
并逐一对应到代码。

**符号表**（代码变量 ↔ 数学符号）：

| 代码 | 数学 | 含义 |
|---|---|---|
| `Sm`, `Sv`, `s` | $S_m$, $S_v$, $s$ | 值函数 $V=\frac12\delta x^\top S_m\delta x+S_v^\top\delta x+s$ |
| `Am`, `Bm` (`dynamics.dfdx/dfdu`) | $A$, $B$ | 动力学 Jacobian |
| `Hv` (`dynamicsBias`) | $h$ | 动力学常数项（投影后产生） |
| `Qm`,`Qv`,`q` (`cost.dfdxx/dfdx/f`) | $Q$, $q_v$, $q$ | 代价对状态的二阶/一阶/零阶 |
| `Rm`,`Rv`,`Pm` (`cost.dfduu/dfdu/dfdux`) | $R$, $r$, $P$ | 代价对输入的项 |
| `Hm` | $H$ | 哈密顿量对输入的 Hessian |
| `Gm`, `Gv` | $G_m$, $G_v$ | 哈密顿量的 $\partial^2/\partial u\partial x$ 与 $\partial/\partial u$ |
| `Km`, `Lv` | $K$, $\ell$ | 反馈增益、前馈项 |
| `Dm` (`stateInputEqConstraint.dfdu`) | $D$ | 状态-输入约束对输入的 Jacobian |
| `Cm` (`.dfdx`), `Ev` (`.f`) | $C$, $e$ | 同上，对状态 / 常数项 |
| `Qu` (`constraintNullProjector_`) | $P_u$ | 零空间投影矩阵 |
| `DmDagger` (`constraintRangeProjector_`) | $D^\dagger$ | 加权伪逆 |

---

## 4.1 问题设定与 Bellman 方程

### 4.1.1 连续时间最优控制问题

$$
\min_{u(\cdot)}\ J = \phi(x(t_f)) + \int_{t_0}^{t_f} l(t,x,u)\,\mathrm{d}t
\qquad\text{s.t.}\qquad \dot x = f(t,x,u),\ x(t_0)=x_0
$$

定义**值函数**（cost-to-go）：

$$
V(t,x) = \min_{u(\cdot)}\left\{\phi(x(t_f))+\int_t^{t_f} l(\tau,x,u)\,\mathrm{d}\tau\right\}
$$

### 4.1.2 HJB 方程

在 $[t, t+\mathrm{d}t]$ 上用动态规划原理：

$$
V(t,x)=\min_u\Big\{l(t,x,u)\,\mathrm{d}t + V\big(t+\mathrm{d}t,\ x+f\,\mathrm{d}t\big)\Big\}
$$

对右端做泰勒展开：

$$
V(t+\mathrm{d}t, x+f\mathrm{d}t) = V(t,x) + \frac{\partial V}{\partial t}\mathrm{d}t
 + \frac{\partial V}{\partial x}^{\!\top} f\,\mathrm{d}t + O(\mathrm{d}t^2)
$$

代入并消去 $V(t,x)$、除以 $\mathrm{d}t$：

$$
\boxed{\;-\frac{\partial V}{\partial t}=\min_u\underbrace{\left\{l(t,x,u)+\frac{\partial V}{\partial x}^{\!\top}f(t,x,u)\right\}}_{\text{哈密顿量 }\mathcal{H}}\;}
$$

边界条件 $V(t_f,x)=\phi(x)$。这是 **Hamilton–Jacobi–Bellman 方程**。

**HJB 的困难**：对一般非线性系统这是无穷维 PDE，不可解。
DDP 的思路是：**沿标称轨迹做二阶泰勒展开，把 HJB 局部化为 Riccati ODE**。

### 4.1.3 二次值函数假设

设标称轨迹 $(\bar x(t),\bar u(t))$，偏差 $\delta x = x-\bar x$、$\delta u = u-\bar u$。
假设值函数在标称轨迹附近是二次的：

$$
\boxed{\;V(t,\bar x+\delta x)\;\approx\;\tfrac12\,\delta x^{\!\top}S_m(t)\,\delta x + S_v(t)^{\!\top}\delta x + s(t)\;}
$$

于是
$$
\frac{\partial V}{\partial x}=S_m\,\delta x + S_v,\qquad
\frac{\partial V}{\partial t}=\tfrac12\delta x^{\!\top}\dot S_m\delta x+\dot S_v^{\!\top}\delta x+\dot s
 \;-\;\big(S_m\delta x+S_v\big)^{\!\top}\dot{\bar x}
$$

（最后一项来自 $\delta x$ 依赖 $t$：$\frac{\mathrm{d}}{\mathrm{d}t}\delta x = \dot x - \dot{\bar x}$；
在标称轨迹上 $\delta x=0$，标准 DDP 推导中把这一项吸收进定义，
下面按 OCS2 代码的约定直接给出结果。）

---

## 4.2 ⭐ 连续时间 Riccati 方程（SLQ）的完整推导

### 4.2.1 展开哈密顿量

哈密顿量

$$
\mathcal{H}(t,x,u)=l(t,x,u)+\frac{\partial V}{\partial x}^{\!\top}f(t,x,u)
$$

在 $(\bar x,\bar u)$ 处做二阶展开。记 LQ 近似量：

$$
l \approx q + q_v^{\!\top}\delta x + r^{\!\top}\delta u
 + \tfrac12\delta x^{\!\top}Q\,\delta x + \delta u^{\!\top}P\,\delta x + \tfrac12\delta u^{\!\top}R\,\delta u
$$
$$
f \approx h + A\,\delta x + B\,\delta u
$$

（$h$ 是 `dynamicsBias`；在未投影的原始问题里 $h=0$，因为
$f(\bar x,\bar u)=\dot{\bar x}$ 已被吸收；投影后会产生非零 $h$，见 §4.5）

代入 $\frac{\partial V}{\partial x}=S_m\delta x+S_v$：

$$
\begin{aligned}
\mathcal{H} &\approx q + q_v^{\!\top}\delta x + r^{\!\top}\delta u
 + \tfrac12\delta x^{\!\top}Q\delta x + \delta u^{\!\top}P\delta x + \tfrac12\delta u^{\!\top}R\delta u\\
&\quad + (S_m\delta x+S_v)^{\!\top}\big(h + A\delta x + B\delta u\big)
\end{aligned}
$$

**按幂次收集**：

| 项 | 系数 |
|---|---|
| 常数 | $q + S_v^{\!\top}h$ |
| $\delta x$ | $q_v + A^{\!\top}S_v + S_m h$ |
| $\delta u$ | $\underbrace{r + B^{\!\top}S_v}_{=:G_v}$ |
| $\frac12\delta x^\top(\cdot)\delta x$ | $Q + S_m A + A^{\!\top}S_m$ |
| $\delta u^\top(\cdot)\delta x$ | $\underbrace{P + B^{\!\top}S_m}_{=:G_m}$ |
| $\frac12\delta u^\top(\cdot)\delta u$ | $\underbrace{R}_{=:H}$ |

**这三个量在代码里就叫 `projectedGv_`、`projectedGm_`、`projectedRm_`**
（`ContinuousTimeRiccatiEquations.cpp:207-210`）：

```cpp
creCache.projectedGm_.noalias() += creCache.projectedBm_.transpose() * Sm;  // P + BᵀSm
creCache.projectedGv_.noalias() += creCache.projectedBm_.transpose() * Sv;  // r + BᵀSv
```

（`projectedGm_`/`projectedGv_` 初值已从 `cost.dfdux`/`cost.dfdu` 插值得到，
所以 `+=` 后就是完整的 $G_m$、$G_v$。）

### 4.2.2 对输入最小化

$$
\frac{\partial\mathcal{H}}{\partial\delta u}=G_v + G_m\,\delta x + H\,\delta u = 0
$$

若 $H\succ0$，最优输入偏差为

$$
\boxed{\;\delta u^\star = -H^{-1}\big(G_m\,\delta x + G_v\big) =: K\,\delta x + \ell\;}
$$

其中

$$
K = -H^{-1}G_m,\qquad \ell = -H^{-1}G_v
$$

**⭐ OCS2 的重要设计**：代码里 $K$、$\ell$ 存的是**没有乘 $H^{-1}$ 的版本**：

```cpp
creCache.projectedKm_ = -(creCache.projectedGm_ + creCache.projectedKm_);   // = -(Gm + ΔGm)
creCache.projectedLv_ = -(creCache.projectedGv_ + creCache.projectedLv_);   // = -(Gv + ΔGv)
```
（`ContinuousTimeRiccatiEquations.cpp:213-215`）

为什么？因为 **§4.4 的约束投影已经把 $H$ 变成了单位阵**！
`computeProjections()` 构造的零空间投影 $P_u$ 满足 $P_u^\top H P_u = I$
（`GaussNewtonDDP.cpp:775-780` 显式验证了这一点）。
所以在投影坐标下 $H^{-1}=I$，$K=-G_m$、$\ell=-G_v$。

`ΔGm`、`ΔGv` 是搜索策略贡献的修正项（Levenberg–Marquardt 用，见 §4.7）。

### 4.2.3 代回 HJB 得到 Riccati 方程

把 $\delta u^\star$ 代回 $\mathcal{H}$。利用 $H=I$（投影后）：

$$
\begin{aligned}
\mathcal{H}^\star &= \big(q+S_v^{\!\top}h\big)
 + \big(q_v+A^{\!\top}S_v+S_m h\big)^{\!\top}\delta x
 + \tfrac12\delta x^{\!\top}\big(Q+S_mA+A^{\!\top}S_m\big)\delta x\\
&\quad + G_v^{\!\top}\delta u^\star + \delta u^{\star\top}G_m\delta x
 + \tfrac12\delta u^{\star\top}H\,\delta u^\star
\end{aligned}
$$

代入 $\delta u^\star = K\delta x+\ell$（$K=-G_m$、$\ell=-G_v$）：

- **$\delta u^\star$ 的线性项贡献**：
  $G_v^\top(K\delta x+\ell) = G_v^\top K\delta x + G_v^\top\ell$
- **交叉项贡献**：$(K\delta x+\ell)^\top G_m\delta x = \delta x^\top K^\top G_m\delta x + \ell^\top G_m\delta x$
- **二次项贡献**：$\frac12(K\delta x+\ell)^\top H(K\delta x+\ell)
  = \frac12\delta x^\top K^\top HK\delta x + \ell^\top HK\delta x + \frac12\ell^\top H\ell$

HJB 要求 $-\dot V = \mathcal{H}^\star$，即逐阶匹配：

$$
\boxed{
\begin{aligned}
-\dot S_m &= Q + S_m A + A^{\!\top}S_m + K^{\!\top}G_m + G_m^{\!\top}K + K^{\!\top}HK\\
-\dot S_v &= q_v + A^{\!\top}S_v + S_m h + G_m^{\!\top}\ell + K^{\!\top}G_v + K^{\!\top}H\ell\\
-\dot s   &= q + S_v^{\!\top}h + \ell^{\!\top}G_v + \tfrac12\,\ell^{\!\top}H\ell
\end{aligned}}
$$

终端条件：$S_m(t_f)=\phi_{xx}$、$S_v(t_f)=\phi_x$、$s(t_f)=\phi$。

**这就是 `computeFlowMapSLQ()` 的 `reducedFormRiccati_ == false` 分支**
（`ContinuousTimeRiccatiEquations.cpp:243-248, 265-272, 286-289`）：

```cpp
// Sm
dSm += creCache.deltaQm_ + SmTrans_projectedAm_ + SmTrans_projectedAm_.transpose();
dSm += projectedKm_T_projectedGm_ + projectedKm_T_projectedGm_.transpose();  // KᵀGm + GmᵀK
dSm.noalias() += projectedKm_.transpose() * projectedRm_projectedKm_;        // KᵀHK
// Sv
dSv.noalias() += Sm.transpose() * projectedHv_;      // Sm·h
dSv.noalias() += projectedAm_.transpose() * Sv;      // Aᵀ·Sv
dSv.noalias() += projectedGm_.transpose() * projectedLv_;   // Gmᵀ·ℓ
dSv.noalias() += projectedKm_.transpose() * projectedGv_;   // Kᵀ·Gv
dSv.noalias() += projectedRm_projectedKm_.transpose() * projectedLv_;  // KᵀHℓ
// s
ds += projectedHv_.dot(Sv);                          // hᵀSv
ds += projectedLv_.dot(projectedGv_);                // ℓᵀGv
ds += 0.5 * projectedLv_.dot(projectedRm_projectedLv_);   // ½ℓᵀHℓ
```

（`dSm`、`dSv`、`ds` 初值是从 `projectedModelData` 插值的 $Q$、$q_v$、$q$。）

### 4.2.4 ⭐ Reduced form：为什么可以少算一半

`reducedFormRiccati_` 为 `true` 时（`preComputeRiccatiTerms_` 设置，默认开），
代码只算：

```cpp
dSm += projectedKm_T_projectedGm_;                       // 只有 KᵀGm，没有 GmᵀK、KᵀHK
dSv.noalias() += projectedGm_.transpose() * projectedLv_; // 只有 Gmᵀℓ
ds += 0.5 * projectedLv_.dot(projectedGv_);              // 只有 ½ℓᵀGv
```

**推导**：在**精确最优**时 $K=-H^{-1}G_m$、$\ell=-H^{-1}G_v$，于是

$$
K^{\!\top}HK = G_m^{\!\top}H^{-1}HH^{-1}G_m = G_m^{\!\top}H^{-1}G_m = -G_m^{\!\top}K = -K^{\!\top}G_m
$$

代入 $S_m$ 方程的后三项：

$$
K^{\!\top}G_m + G_m^{\!\top}K + K^{\!\top}HK
= K^{\!\top}G_m + G_m^{\!\top}K - K^{\!\top}G_m
= G_m^{\!\top}K = K^{\!\top}G_m
$$

（最后一步用了 $K^\top G_m = -G_m^\top H^{-1}G_m$ 对称）。**三项塌缩成一项**。

同理 $S_v$ 方程：

$$
G_m^{\!\top}\ell + K^{\!\top}G_v + K^{\!\top}H\ell
= G_m^{\!\top}\ell + K^{\!\top}G_v - K^{\!\top}G_v = G_m^{\!\top}\ell
$$

$s$ 方程：

$$
\ell^{\!\top}G_v + \tfrac12\ell^{\!\top}H\ell = \ell^{\!\top}G_v - \tfrac12\ell^{\!\top}G_v
=\tfrac12\,\ell^{\!\top}G_v
$$

**完全对上代码**。

**计算复杂度对比**（代码注释里明确写出，`:234-236`）：

| 形式 | $S_m$ 的复杂度 |
|---|---|
| reduced | $O(n_x^3) + 2\,O(n_x^2 n_p)$ |
| full | $O(n_x^3) + 3\,O(n_x^2 n_p) + O(n_x n_p^2)$ |

（$n_p$ 是投影后的输入维度）reduced form 还省掉了对 $R$ 的插值。

**代价**：reduced form 只在 $K,\ell$ 严格取最优值时等价。
当 Levenberg–Marquardt 引入 $\Delta G_m\ne0$ 时 $K\ne-H^{-1}G_m$，
此时 reduced form 引入误差。这就是为什么它是可配置的。

### 4.2.5 值函数的状态一致性

`getValueFunctionImpl`（`GaussNewtonDDP.cpp`）返回的值函数是**相对标称轨迹**的：

$$
V(t,x)=\tfrac12(x-\bar x)^{\!\top}S_m(x-\bar x)+S_v^{\!\top}(x-\bar x)+s
$$

要转成绝对坐标 $\frac12 x^\top S_m x + \tilde S_v^\top x + \tilde s$，需

$$
\tilde S_v = S_v - S_m\bar x
$$

SQP 里也有同样的修正（`SqpSolver.cpp:315`）：

```cpp
valueFunction_[i].dfdx.noalias() -= valueFunction_[i].dfdxx * x[i];
```

---

## 4.3 ⭐ 离散时间 Riccati 方程（iLQR）

**文件**：`ocs2_ddp/src/riccati_equations/DiscreteTimeRiccatiEquations.cpp:65`

### 4.3.1 设定

离散动力学 $\delta x_{k+1}=A_k\delta x_k+B_k\delta u_k+h_k$，
阶段代价 $l_k$ 的二次近似同前。Bellman 递推：

$$
V_k(\delta x) = \min_{\delta u}\Big\{l_k(\delta x,\delta u) + V_{k+1}\big(A\delta x+B\delta u+h\big)\Big\}
$$

代入 $V_{k+1}(\xi)=\frac12\xi^\top S_m'\xi + S_v'^\top\xi + s'$（撇号表示 $k+1$）：

$$
\begin{aligned}
\mathcal{Q}(\delta x,\delta u)&= q+q_v^{\!\top}\delta x+r^{\!\top}\delta u
 +\tfrac12\delta x^{\!\top}Q\delta x+\delta u^{\!\top}P\delta x+\tfrac12\delta u^{\!\top}R\delta u\\
&\quad+\tfrac12\big(A\delta x+B\delta u+h\big)^{\!\top}S_m'\big(A\delta x+B\delta u+h\big)\\
&\quad+S_v'^{\!\top}\big(A\delta x+B\delta u+h\big)+s'
\end{aligned}
$$

### 4.3.2 收集各项

展开二次型：

$$
\tfrac12(A\delta x+B\delta u+h)^{\!\top}S_m'(\cdots)
=\tfrac12\delta x^{\!\top}A^{\!\top}S_m'A\delta x
+\delta u^{\!\top}B^{\!\top}S_m'A\delta x
+\tfrac12\delta u^{\!\top}B^{\!\top}S_m'B\delta u
$$
$$
\quad+\,h^{\!\top}S_m'A\delta x + h^{\!\top}S_m'B\delta u + \tfrac12 h^{\!\top}S_m'h
$$

于是：

$$
\boxed{
\begin{aligned}
H &:= R + B^{\!\top}S_m'B  &&\text{（输入 Hessian）}\\
G_m &:= P + B^{\!\top}S_m'A &&\text{（交叉项）}\\
G_v &:= r + B^{\!\top}\big(S_v' + S_m'h\big) &&\text{（输入梯度）}
\end{aligned}}
$$

对照代码（`:71-82`）：

```cpp
dreCache.Sm_projectedHv_ = SmNext * dynamicsBias;                       // Sm'·h
dreCache.Sm_projectedAm_ = SmNext * dynamics.dfdx;                      // Sm'·A
dreCache.Sm_projectedBm_ = SmNext * dynamics.dfdu;                      // Sm'·B
dreCache.Sv_plus_Sm_projectedHv_ = SvNext + Sm_projectedHv_;            // Sv' + Sm'h

projectedGm_ = cost.dfdux + dynamics.dfdu.transpose() * Sm_projectedAm_;   // P + BᵀSm'A ✓
projectedGv_ = cost.dfdu  + dynamics.dfdu.transpose() * Sv_plus_Sm_projectedHv_; // r + Bᵀ(Sv'+Sm'h) ✓
projectedHm_ = cost.dfduu + Sm_projectedBm_.transpose() * dynamics.dfdu;   // R + BᵀSm'B ✓
```

### 4.3.3 最小化与回代

$$
\delta u^\star = -H^{-1}(G_m\delta x + G_v) = K\delta x+\ell
$$

（同样地，投影后 $H=I$，故代码里 `projectedKm = -Gm - ΔGm`、`projectedLv = -Gv - ΔGv`，`:85-87`）

回代得到 $V_k$ 的系数：

$$
\boxed{
\begin{aligned}
S_m &= Q + A^{\!\top}S_m'A + K^{\!\top}G_m + G_m^{\!\top}K + K^{\!\top}HK\\
S_v &= q_v + A^{\!\top}\big(S_v'+S_m'h\big) + G_m^{\!\top}\ell + K^{\!\top}G_v + K^{\!\top}H\ell\\
s   &= s' + q + h^{\!\top}\big(S_v'+S_m'h\big) - \tfrac12 h^{\!\top}S_m'h
       + \ell^{\!\top}G_v + \tfrac12\ell^{\!\top}H\ell
\end{aligned}}
$$

对照代码（`:103-146`）：

```cpp
Sm = cost.dfdxx + deltaQm_;
Sm.noalias() += Sm_projectedAm_.transpose() * dynamics.dfdx;      // AᵀSm'A
Sm += projectedKm_T_projectedGm_ + projectedKm_T_projectedGm_.transpose();
Sm.noalias() += projectedKm.transpose() * projectedHm_projectedKm_;

Sv = cost.dfdx;
Sv.noalias() += dynamics.dfdx.transpose() * Sv_plus_Sm_projectedHv_;
Sv.noalias() += projectedGm_.transpose() * projectedLv;
Sv.noalias() += projectedKm.transpose() * projectedGv_;
Sv.noalias() += projectedHm_projectedKm_.transpose() * projectedLv;

s  = sNext + cost.f;
s += dynamicsBias.dot(Sv_plus_Sm_projectedHv_);   // hᵀ(Sv' + Sm'h)
s -= 0.5 * dynamicsBias.dot(Sm_projectedHv_);     // −½hᵀSm'h
```

**注意 $s$ 里的 $h^\top(S_v'+S_m'h)-\frac12 h^\top S_m'h = h^\top S_v' + \frac12 h^\top S_m'h$**，
代码用两步写是为了复用已算好的中间量 `Sv_plus_Sm_projectedHv_`，避免额外的矩阵-向量乘法。

**Reduced form** 的塌缩推导与 §4.2.4 完全一致。

### 4.3.4 连续 vs 离散：一个自洽性检查

令离散量取 $A=I+\Delta t\,A_c$、$B=\Delta t\,B_c$、$Q=\Delta t\,Q_c$ 等（$\Delta t\to0$）：

$$
S_m^{(k)} = \Delta t\,Q_c + (I+\Delta t A_c)^{\!\top}S_m^{(k+1)}(I+\Delta t A_c)+\cdots
$$

保留 $O(\Delta t)$：

$$
S_m^{(k)}-S_m^{(k+1)} = \Delta t\big(Q_c + A_c^{\!\top}S_m + S_mA_c + \cdots\big)
$$

而 $S_m^{(k)}-S_m^{(k+1)}\approx -\Delta t\,\dot S_m$，故 $-\dot S_m = Q_c+A_c^\top S_m+S_mA_c+\cdots$，
与 §4.2.3 一致 ✓。

---

## 4.4 ⭐ 状态-输入等式约束的投影

这是 OCS2 处理硬约束的核心技巧，也是 SLQ 与经典 DDP 的最大区别。

### 4.4.1 问题

在每个时刻有约束

$$
C\,\delta x + D\,\delta u + e = 0,\qquad D\in\mathbb{R}^{n_c\times n_u}
$$

我们希望在**约束流形内**做无约束优化。§3.5.3 已给出通解形式，
但那里用的是标准 QR（$Q^\top Q=I$）。DDP 需要一个**更强的性质**：
投影后哈密顿量 Hessian 变成单位阵。

### 4.4.2 加权投影的构造

**文件**：`ocs2_core/src/misc/LinearAlgebra.cpp:129` 的 `computeConstraintProjection`

设 $H\succ0$（哈密顿量 Hessian）。第一步做 $H^{-1}$ 的 **UUT 分解**
（`computeInverseMatrixUUT`，`:119`）：

$$
H = L L^{\!\top}\ (\text{Cholesky})\ \Longrightarrow\
H^{-1}=L^{-\top}L^{-1} = U U^{\!\top},\quad U:=L^{-\top}
$$

代码里 `AmInvUmUmT` 存的就是 $U$（对 $L^\top$ 做原地上三角求解得到 $L^{-\top}$）。

第二步对 $U^{\!\top}D^{\!\top}$ 做 Householder QR（`:135`）：

$$
U^{\!\top}D^{\!\top}=\begin{bmatrix}\mathcal{Q}_c & \mathcal{Q}_u\end{bmatrix}
\begin{bmatrix}\mathcal{R}_c\\0\end{bmatrix}
$$

其中 $\mathcal{Q}=[\mathcal{Q}_c\ \mathcal{Q}_u]$ 正交，$\mathcal{R}_c$ 上三角。

**第三步定义两个投影**：

$$
\boxed{
\begin{aligned}
P_u &:= U\,\mathcal{Q}_u && \text{（零空间投影，`constraintNullProjector_`）}\\
D^\dagger &:= U\,\mathcal{Q}_c\,\mathcal{R}_c^{-\top} && \text{（加权伪逆，`constraintRangeProjector_`）}
\end{aligned}}
$$

对照代码（`:151-154`）：

```cpp
DmDaggerTRmDmDaggerUUT.setIdentity(numConstraints, numConstraints);
QRof_RmInvUmUmTT_DmT_Rc.triangularView<Eigen::Upper>().solveInPlace(DmDaggerTRmDmDaggerUUT); // = Rc^{-1}
DmDagger.noalias() = RmInvUmUmT * (QRof_..._Qc * DmDaggerTRmDmDaggerUUT.transpose());  // U Qc Rc^{-T}
RmInvConstrainedUUT.noalias() = RmInvUmUmT * QRof_..._Qu;                              // U Qu
```

### 4.4.3 三条关键性质的验证

**性质 1：$P_u$ 张成 $D$ 的零空间**

$$
D\,P_u = D\,U\,\mathcal{Q}_u = \big(U^{\!\top}D^{\!\top}\big)^{\!\top}\mathcal{Q}_u
= \big(\mathcal{Q}_c\mathcal{R}_c\big)^{\!\top}\mathcal{Q}_u
= \mathcal{R}_c^{\!\top}\underbrace{\mathcal{Q}_c^{\!\top}\mathcal{Q}_u}_{=0}=0\quad\checkmark
$$

**性质 2：$P_u^{\!\top}HP_u = I$** ← 这是 §4.2.2 简化的依据

$$
P_u^{\!\top}HP_u=\mathcal{Q}_u^{\!\top}U^{\!\top}HU\,\mathcal{Q}_u
$$

注意 $U=L^{-\top}$，$H=LL^\top$，故

$$
U^{\!\top}HU = L^{-1}\,LL^{\!\top}\,L^{-\top}=I
$$

于是 $P_u^\top HP_u = \mathcal{Q}_u^\top\mathcal{Q}_u = I\quad\checkmark$

**代码里显式验证了这一点**（`GaussNewtonDDP.cpp:774-781`）：

```cpp
if (ddpSettings_.checkNumericalStability_) {
  matrix_t HmProjected = constraintNullProjector.transpose() * Hm * constraintNullProjector;
  if (!HmProjected.isApprox(matrix_t::Identity(nullSpaceDim, nullSpaceDim), 1e-6))
    throw std::runtime_error("HmProjected should be identity!");
}
```

**性质 3：$D\,D^\dagger = I$**

$$
D\,D^\dagger = D\,U\mathcal{Q}_c\mathcal{R}_c^{-\top}
= \big(\mathcal{Q}_c\mathcal{R}_c\big)^{\!\top}\mathcal{Q}_c\mathcal{R}_c^{-\top}
= \mathcal{R}_c^{\!\top}\mathcal{Q}_c^{\!\top}\mathcal{Q}_c\mathcal{R}_c^{-\top}
= \mathcal{R}_c^{\!\top}\mathcal{R}_c^{-\top}=I\quad\checkmark
$$

**$D^\dagger$ 是什么伪逆？** 它是最小化 $\|\delta u\|_H^2 = \delta u^\top H\delta u$
意义下的最小范数解，即**加权最小二乘伪逆**：

$$
D^\dagger = H^{-1}D^{\!\top}\big(DH^{-1}D^{\!\top}\big)^{-1}
$$

验证：$H^{-1}D^\top = UU^\top D^\top = U(\mathcal{Q}_c\mathcal{R}_c)$，
而 $DH^{-1}D^\top = (\mathcal{Q}_c\mathcal{R}_c)^\top(\mathcal{Q}_c\mathcal{R}_c)=\mathcal{R}_c^\top\mathcal{R}_c$，
故

$$
H^{-1}D^{\!\top}(DH^{-1}D^{\!\top})^{-1}
= U\mathcal{Q}_c\mathcal{R}_c\,\mathcal{R}_c^{-1}\mathcal{R}_c^{-\top}
= U\mathcal{Q}_c\mathcal{R}_c^{-\top} = D^\dagger\quad\checkmark
$$

**为什么用加权伪逆而不是普通伪逆？** 因为我们要最小化的是代价，
$H$ 是代价对输入的曲率。用 $H$ 加权意味着"在代价意义下最省力地满足约束"，
这与最优性一致。

### 4.4.4 数值鲁棒性

`setTriangularMinimumEigenvalues`（`LinearAlgebra.cpp:38`）在 QR 之后
把 $\mathcal{R}_c$ 的对角元的绝对值下限设为 `minEigenValue`：

```cpp
if (eigenValue < 0.0) eigenValue = std::min(-minEigenValue, eigenValue);
else                  eigenValue = std::max( minEigenValue, eigenValue);
```

这防止约束近乎冗余（$D$ 接近行亏秩）时 $\mathcal{R}_c^{-1}$ 爆炸。
**注意它保号**——直接置 $\pm\epsilon$ 而不是取绝对值，保持分解的符号结构。

`Dm.rows()==0`（无约束）时走快速路径：$P_u=U$（即 $H^{-1}$ 的 UUT 因子），
$D^\dagger$ 为空（`GaussNewtonDDP.cpp:762-765`）。
此时 $P_u^\top HP_u = U^\top HU = I$ 依然成立 ✓

---

## 4.5 ⭐ `projectLQ`：把 LQ 问题投到零空间

**文件**：`ocs2_ddp/src/DDP_HelperFunctions.cpp:143`

有了投影，输入替换为

$$
\delta u = \underbrace{P_u}_{\texttt{Qu}}\,\delta\tilde u
\;\underbrace{-\,D^\dagger C}_{\texttt{Px}=-\texttt{CmProjected}}\,\delta x
\;\underbrace{-\,D^\dagger e}_{\texttt{u0}=-\texttt{EvProjected}}
$$

代码把 $D^\dagger C$、$D^\dagger e$ 存进 `projectedModelData.stateInputEqConstraint`
的 `dfdx`、`f` 字段（复用了这个已经不再需要的结构体）：

```cpp
projectedModelData.stateInputEqConstraint.f    = constraintRangeProjector * modelData...f;    // D†e
projectedModelData.stateInputEqConstraint.dfdx = constraintRangeProjector * modelData...dfdx; // D†C
const auto& Pu = constraintNullProjector;
const matrix_t Px = -projectedModelData.stateInputEqConstraint.dfdx;   // -D†C
const matrix_t u0 = -projectedModelData.stateInputEqConstraint.f;      // -D†e
```

然后对动力学与代价施加 §3.3.3 的变量替换：

```cpp
changeOfInputVariables(projectedModelData.dynamics, Pu, Px, u0);
projectedModelData.dynamicsBias += modelData.dynamics.dfdu * u0;    // ← h 的来源！
changeOfInputVariables(projectedModelData.cost, Pu, Px, u0);
```

**⭐ `dynamicsBias`（数学里的 $h$）从哪来**：
原动力学 $\delta\dot x = A\delta x+B\delta u$，代入 $\delta u = P_u\delta\tilde u+P_x\delta x+u_0$：

$$
\delta\dot x = \underbrace{(A+BP_x)}_{\tilde A}\delta x + \underbrace{BP_u}_{\tilde B}\delta\tilde u
+ \underbrace{B\,u_0}_{h}
$$

**这就是为什么 SLQ 的 Riccati 方程里有 $h$ 项，而教科书上的 LQR 没有**——
$h = -B D^\dagger e$ 是"为了满足约束必须施加的前馈输入"引起的动力学漂移。
若约束在标称轨迹上已满足（$e=0$），则 $h=0$。

无约束情形（`:153-171`）走简化路径：只有 $u=P_u\tilde u$，$h$ 保持为 0。

---

## 4.6 从 Riccati 解构造控制器

**SLQ**（`ocs2_ddp/src/SLQ.cpp:127`）：

```cpp
projectedKm = -(ΔGm + P̃);           projectedKm -= B̃ᵀ·Sm;    // K̃ = -(P̃ + B̃ᵀSm + ΔGm)
projectedLv = -(ΔGv + r̃);           projectedLv -= B̃ᵀ·Sv;    // ℓ̃ = -(r̃ + B̃ᵀSv + ΔGv)

dstController.gainArray_[k]      = -CmProjected + Qu * projectedKm;     // K = -D†C + P_u K̃
dstController.biasArray_[k]      = ū - K·x̄;
dstController.deltaBiasArray_[k] = -EvProjected + Qu * projectedLv;     // Δb = -D†e + P_u ℓ̃
```

**推导**：投影坐标下的最优 $\delta\tilde u = \tilde K\delta x + \tilde\ell$，还原：

$$
\delta u = P_u(\tilde K\delta x+\tilde\ell) + P_x\delta x + u_0
= \underbrace{(P_uK̃ - D^\dagger C)}_{K}\delta x + \underbrace{(P_u\tilde\ell - D^\dagger e)}_{\Delta b}
$$

绝对坐标下：$u = \bar u + K(x-\bar x) + \Delta b = Kx + \underbrace{(\bar u - K\bar x)}_{b} + \Delta b$。

**iLQR**（`ILQR.cpp:162`）完全一样，只是 $\tilde K,\tilde\ell$ 在后向递推中已算好，
存在 `projectedKmTrajectoryStock_`、`projectedLvTrajectoryStock_` 里。

**两者唯一的算法差别**是 `computeHamiltonianHessian`：

```cpp
// SLQ.cpp:206
matrix_t SLQ::computeHamiltonianHessian(const ModelData& md, const matrix_t& Sm) const {
  return searchStrategyPtr_->augmentHamiltonianHessian(md, md.cost.dfduu);         // H = R
}
// ILQR.cpp:217
matrix_t ILQR::computeHamiltonianHessian(const ModelData& md, const matrix_t& Sm) const {
  matrix_t Hm = md.cost.dfduu;
  Hm.noalias() += (md.dynamics.dfdu.transpose() * Sm) * md.dynamics.dfdu;          // H = R + BᵀSm B
  return searchStrategyPtr_->augmentHamiltonianHessian(md, Hm);
}
```

与 §4.2.1、§4.3.2 的推导一致：连续时间的 $\mathcal{H}$ 对 $u$ 的二阶导只有 $R$
（因为 $\frac{\partial V}{\partial x}^\top f$ 对 $u$ 是线性的）；
离散时间因为 $V_{k+1}$ 的二次项被 $B\delta u$ 穿过，多了 $B^\top S_m'B$。

---

## 4.7 ⭐ 搜索策略

### 4.7.1 优点函数（merit function）

**文件**：`GaussNewtonDDP.cpp:500`

$$
\boxed{\;
M = \underbrace{J}_{\texttt{cost}}
 + \underbrace{\rho\sqrt{\textstyle\sum\|g\|^2}}_{\text{等式约束罚}}
 + \underbrace{\mathcal{L}_{eq}}_{\text{等式 AL}}
 + \underbrace{\mathcal{L}_{ineq}}_{\text{不等式 AL}}\;}
$$

```cpp
scalar_t merit = performanceIndex.cost;
merit += constraintPenaltyCoefficients_.penaltyCoeff * std::sqrt(performanceIndex.equalityConstraintsSSE);
merit += performanceIndex.equalityLagrangian;
merit += performanceIndex.inequalityLagrangian;
```

**为什么用 $\sqrt{\text{SSE}}$（即 $\ell_2$ 范数）而不是 SSE 本身？**
因为 $\ell_1$/$\ell_2$ 范数罚是**精确罚函数**（exact penalty）：
存在有限的 $\rho$ 使罚问题的极小点就是原约束问题的极小点。
而平方罚需要 $\rho\to\infty$ 才精确。见 Nocedal & Wright §17.2。

### 4.7.2 罚系数的自适应更新

**文件**：`GaussNewtonDDP.cpp:787`（初始化）与 `:807`（更新）

初始化：
$$
\rho_0 = \texttt{constraintPenaltyInitialValue\_},\qquad
\epsilon_0 = \rho_0^{-0.1}
$$

每次迭代后（`:807-817`）：

$$
\begin{cases}
\epsilon \leftarrow \epsilon\,/\,\rho^{0.9} & \text{若 } \|g\|^2_{SSE} < \epsilon\quad(\text{约束够好，只收紧容差})\\[4pt]
\rho \leftarrow \gamma\rho,\quad \epsilon\leftarrow\epsilon\,/\,\rho^{0.1} & \text{否则（约束不够好，加大罚）}
\end{cases}
$$

且 $\epsilon\ \ge$ `constraintTolerance_`（下界钳位）。

这是 **Conn–Gould–Toint 增广拉格朗日**的经典更新策略（LANCELOT 算法）：
指数 $0.9$ 与 $0.1$ 使得"罚系数增长时容差收紧得慢，罚系数不变时容差收紧得快"，
从而在可行性与最优性之间取得平衡。

### 4.7.3 并行 Armijo 回溯线搜索

**文件**：`ocs2_ddp/src/search_strategy/LineSearchStrategy.cpp`

**步长序列**：$\alpha_j = \alpha_{\max}\,c^{\,j}$，$j=0,1,2,\dots$，
其中 $c=$`contractionRate`（默认 0.5）。最大试探次数（`:71`）：

$$
n_{\max}=\left\lfloor\frac{\ln(\alpha_{\min}/\alpha_{\max})}{\ln c}\right\rfloor+1
$$

**接受准则**（`:235-236`）：

$$
\boxed{\;M(\alpha) < M(0) - c_1\,\alpha\,\Xi\;}
$$

其中 $c_1=$`armijoCoefficient`，$\Xi=$`unoptimizedControllerUpdateIS_`
= `computeControllerUpdateIS(unoptimizedController)`，即前馈更新量的积分平方：

$$
\Xi = \int_{t_0}^{t_f}\|\Delta b(t)\|^2\,\mathrm{d}t
$$

**⭐ 这是标准 Armijo 条件的一个 DDP 特化**。标准 Armijo 是
$f(x+\alpha p)\le f(x)+c_1\alpha\nabla f^\top p$，其中 $\nabla f^\top p<0$ 是方向导数。
DDP 中可以证明（Todorov & Li 2005; Tassa et al. 2012）：
沿 DDP 方向的预期代价下降为

$$
\Delta J(\alpha) = \alpha\sum_k \ell_k^{\!\top}G_v^{(k)} + \tfrac{\alpha^2}{2}\sum_k \ell_k^{\!\top}H_k\ell_k
$$

在投影坐标下 $H=I$、$\ell=-G_v$，故一阶项 $=-\alpha\sum\|\ell_k\|^2 = -\alpha\,\Xi$。
所以 $-\Xi$ 正是方向导数的代理量，代码用它是**理论正确**的。

**并行化**（`:188` 的 `lineSearchTask`）：

```cpp
const size_t alphaExp = alphaExpNext_++;              // 原子递增，无锁任务分发
const scalar_t stepLength = maxStepLength * pow(contractionRate, alphaExp);
if (stepLength < bestStepSize_) break;                 // 已有更大的可行步长，放弃
...
if (armijoCondition && stepLength > bestStepSize_) {   // 加锁更新最优
  bestStepSize_ = stepLength;
  swap(*bestSolutionRef_, workersSolution_[taskId]);
  terminateLinesearchTasks = 所有更大的 α 都已处理完;
}
if (terminateLinesearchTasks) {
  for (auto& rollout : rolloutRefStock_) rollout.abortRollout();   // 杀掉其他线程的积分
}
```

**巧妙之处**：由于 $\alpha$ 越大越好，一旦发现某个 $\alpha_j$ 满足条件
且所有 $\alpha_{j'}(j'<j)$ 都已验证失败，就可以立即中止其余线程。
`abortRollout()` 通过设置 `killIntegration_` 标志让积分器抛异常退出，
避免浪费算力在必然被丢弃的更小步长上。

**结果与单线程回溯线搜索完全等价**（代码注释 `:233` 明确指出），
只是把串行的"试 $\alpha_0$ 失败→试 $\alpha_1$..."变成并行。

**Riccati 修正**（`:294`）：线搜索策略的 $\Delta Q_m$ 用于 Hessian 修正：

```cpp
matrix_t Q_minus_PTRinvP = Q̃ - P̃ᵀP̃;      // 注意投影后 R=I，故 P̃ᵀR̃⁻¹P̃ = P̃ᵀP̃
deltaQm = Q_minus_PTRinvP;
hessian_correction::shiftHessian(strategy, deltaQm, multiple);
deltaQm -= Q_minus_PTRinvP;               // ΔQm = 修正量 - 原值 = 纯增量
deltaGv.setZero();  deltaGm.setZero();    // 线搜索不改 G
```

$Q-P^\top R^{-1}P$ 是**代价关于 $x$ 的 Schur 补**（消去 $u$ 后的有效状态曲率）。
只有它半正定，Riccati 递推才稳定。所以修正针对它而非 $Q$ 本身。

### 4.7.4 Levenberg–Marquardt 策略

**文件**：`ocs2_ddp/src/search_strategy/LevenbergMarquardtStrategy.cpp`

思想：不做线搜索（永远取 $\alpha=1$，`:75`），而是在**信赖域**意义下修正 Riccati。

**(a) 哈密顿量 Hessian 增广**（`:245`）：

$$
\boxed{\;H_{\text{aug}} = H + \mu\,B^{\!\top}B\;}
$$

```cpp
HmAug.noalias() += lmModule_.riccatiMultiple * dynamics.dfdu.transpose() * dynamics.dfdu;
```

**(b) Riccati 修正项**（`:230`）：

$$
\Delta Q_m = 0,\qquad
\Delta G_v = \mu\,B^{\!\top}h,\qquad
\Delta G_m = \mu\,B^{\!\top}A
$$

**⭐ 这三项的来源**：LM 相当于在代价里加了一项**动力学缺口的平方罚**

$$
\tilde l = l + \frac{\mu}{2}\big\|A\delta x+B\delta u+h\big\|^2
$$

展开：

$$
\frac{\mu}{2}\|\cdot\|^2 = \frac{\mu}{2}\Big(\delta u^{\!\top}B^{\!\top}B\delta u
+2\delta u^{\!\top}B^{\!\top}(A\delta x+h) + \|A\delta x+h\|^2\Big)
$$

对应地：
- $\partial^2/\partial u^2$ 增加 $\mu B^\top B$ → **(a)** ✓
- $\partial^2/\partial u\partial x$ 增加 $\mu B^\top A$ → $\Delta G_m$ ✓
- $\partial/\partial u$ 增加 $\mu B^\top h$ → $\Delta G_v$ ✓
- $\partial^2/\partial x^2$ 增加 $\mu A^\top A$ → 但代码里 $\Delta Q_m=0$

**为什么 $\Delta Q_m = 0$？** 因为 $\mu A^\top A$ 那一项只影响值函数的曲率估计，
不影响下降方向（$K$、$\ell$ 只由 $H$、$G_m$、$G_v$ 决定）。省掉它节省一次 $n_x^3$ 乘法，
且代码注释里 `deltaQm.setZero(...)` 明确如此。

**几何意义**：$\mu\|B\delta u\|^2$ 惩罚"输入偏离标称太远导致状态被推走太多"，
$\mu\to\infty$ 时 $\delta u\to0$（退化为不更新），$\mu\to0$ 时退化为纯 Gauss-Newton。
这正是 Levenberg–Marquardt 在最小二乘中的作用。

**(c) 信赖域半径自适应**（`:122-152`）：

$$
\rho = \frac{\text{实际下降}}{\text{预测下降}}
= \frac{M_{\text{prev}} - M_{\text{new}}}{M_{\text{prev}} - M_{\text{LQ 预测}}}
$$

```cpp
const auto actualReduction    = prevMerit - solution.performanceIndex.merit;
const auto expectedReduction  = solution.performanceIndex.merit - expectedCost;  // 由 Riccati 给出
const auto pho = reductionToPredictedReduction(actualReduction, expectedReduction);

if (pho < 0.25) {        // 模型不可信 → 收缩信赖域（增大 μ）
  ratio = max(1.0, ratio) * defaultRatio;
  μ = max(ratio * μ, defaultFactor);
} else if (pho > 0.75) { // 模型很准 → 扩大信赖域（减小 μ）
  ratio = min(1.0, ratio) / defaultRatio;
  μ = (ratio*μ > defaultFactor) ? ratio*μ : 0.0;    // 可以完全关掉
} else {
  ratio = 1.0;           // μ 不变
}
// 接受/拒绝
if (pho >= minAcceptedPho) { numSuccessiveRejections = 0; return true;  }
else                       { ++numSuccessiveRejections;   return false; }
```

阈值 0.25/0.75 是 Moré (1978) 提出的经典值，Nocedal & Wright §4.1 也采用。
`riccatiMultipleAdaptiveRatio` 是一个**二阶自适应**：
连续多次 $\rho<0.25$ 会让 ratio 本身也增大，实现"越不准就越快收缩"。

连续拒绝超过 `maxNumSuccessiveRejections` 则抛异常（`:177`）。

**LM vs 线搜索的取舍**：

| | 线搜索 | LM |
|---|---|---|
| 每次迭代 rollout 次数 | 1 ~ $n_{\max}$（并行） | 恰好 1 |
| 对不稳定系统 | 可能所有 $\alpha$ 都发散 | $\mu$ 大时天然保守 |
| 参数敏感性 | 低 | 中（需调 $\mu$ 的初值与比率） |
| 默认 | ✓ | |

---

## 4.8 ⭐ Hessian 修正的四种策略

**文件**：`ocs2_ddp/src/HessianCorrection.cpp:53`

目的：保证 $\hat Q := Q-P^\top R^{-1}P \succeq \lambda_{\min} I$。

### (1) `DIAGONAL_SHIFT`

$$
\hat Q \leftarrow \hat Q + \lambda_{\min} I
$$

最简单最快，但**不保证**结果 PSD（如果原来有小于 $-\lambda_{\min}$ 的特征值）。

### (2) `EIGENVALUE_MODIFICATION`

**文件**：`LinearAlgebra.cpp:52`

$$
\hat Q = V\Lambda V^{\!\top}\ \Longrightarrow\
\hat Q \leftarrow V\,\max(\Lambda,\lambda_{\min}I)\,V^{\!\top}
$$

```cpp
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(squareMatrix, Eigen::EigenvaluesOnly);
vector_t lambda = eig.eigenvalues();
for (j) if (lambda(j) < minEigenvalue) { hasNegative = true; lambda(j) = minEigenvalue; }
if (hasNegative) {
  eig.compute(squareMatrix, Eigen::ComputeEigenvectors);
  squareMatrix = eig.eigenvectors() * lambda.asDiagonal() * eig.eigenvectors().inverse();
} else {
  squareMatrix = 0.5 * (squareMatrix + squareMatrix.transpose());   // 只做对称化
}
```

**最精确**（改动量在 Frobenius 范数下最小），但需要完整特征分解，$O(n^3)$ 且常数大。
优化点：先只算特征值，只有确实需要修正时才算特征向量。

### (3) `GERSHGORIN_MODIFICATION`

**文件**：`LinearAlgebra.cpp:77`

用 **Gershgorin 圆盘定理**：矩阵的所有特征值落在圆盘
$\bigcup_i \{z: |z-a_{ii}|\le R_i\}$ 内，$R_i=\sum_{j\ne i}|a_{ij}|$。

若令 $a_{ii}\ge R_i+\lambda_{\min}$，则该圆盘的左端点 $\ge\lambda_{\min}$，
故所有特征值 $\ge\lambda_{\min}$。

```cpp
squareMatrix = 0.5 * (squareMatrix + squareMatrix.transpose());
for (i) {
  auto Ri = squareMatrix.col(i).cwiseAbs().sum() - std::abs(squareMatrix(i,i));  // 半径
  squareMatrix(i,i) = std::max(squareMatrix(i,i), Ri + minEigenvalue);
}
```

**$O(n^2)$，最快的保证 PSD 的方法**，但很保守（可能大幅改变矩阵）。
代码用列和而非行和——因为对称矩阵两者相同，而 Eigen 的列访问是连续内存，更快。

### (4) `CHOLESKY_MODIFICATION`

**文件**：`LinearAlgebra.cpp:90`

**修正 Cholesky 分解**（Gill–Murray–Wright）：对 $A-\lambda_{\min}I$ 做不完全 Cholesky
$P^\top(A+E)P = LL^\top$（$E$ 是算法自动选择的最小对角修正），
然后重构 $A\leftarrow LL^\top+\lambda_{\min}I$。

```cpp
A.diagonal().array() -= minEigenvalue;
sparse_matrix_t squareMatrix = 0.5*A.sparseView();  A.transposeInPlace();
squareMatrix += 0.5*A.sparseView();                       // 对称化
Eigen::IncompleteCholesky<scalar_t> incompleteCholesky(squareMatrix);
sparse_matrix_t M = incompleteCholesky.scalingS().asDiagonal().inverse()
                    * incompleteCholesky.matrixL();
LmTwisted.selfadjointView<Lower>() = M.selfadjointView<Lower>()
                    .twistedBy(incompleteCholesky.permutationP());   // 撤销置换
A = L * L.transpose();
A.diagonal().array() += minEigenvalue;
```

**精度与速度介于 (2) 和 (3) 之间**，且天然适合稀疏矩阵。
`Eigen::IncompleteCholesky` 内部会自动加对角偏移直到分解成功，
这正是 GMW 修正 Cholesky 的思想。

### 何时调用

1. **事件时刻与终端时刻**的 $\phi_{xx}$：`approximateOptimalControlProblem` 里
   （`GaussNewtonDDP.cpp:696, 725`），仅在 `LINE_SEARCH` 策略下
2. **中间时刻的 $\hat Q$**：`LineSearchStrategy::computeRiccatiModification`（`:294`）

LM 策略不用它们，因为 LM 通过 $\mu B^\top B$ 已经保证了 $H\succ0$。

---

## 4.9 主循环与数据管理

### 4.9.1 三组数据

```cpp
PrimalDataContainer nominalPrimalData_,  cachedPrimalData_;
DualDataContainer   nominalDualData_,    cachedDualData_;
PrimalSolution      optimizedPrimalSolution_;
DualSolution        optimizedDualSolution_;
ProblemMetrics      optimizedProblemMetrics_;
```

| 名称 | 含义 |
|---|---|
| `nominal*` | **当前线性化点**（LQ 近似围绕它展开） |
| `optimized*` | **本次迭代产出的新解**（线搜索/LM 的结果） |
| `cached*` | 上上次的数据，用 `swap` 循环复用内存，避免反复分配 |

迭代末尾（`GaussNewtonDDP.cpp:1080-1085`）：

```cpp
nominalDualData_.swap(cachedDualData_);
nominalPrimalData_.swap(cachedPrimalData_);
optimizedDualSolution_.swap(nominalDualData_.dualSolution);
optimizedPrimalSolution_.swap(nominalPrimalData_.primalSolution);
optimizedProblemMetrics_.swap(nominalPrimalData_.problemMetrics);
```

全部是 $O(1)$ 的指针交换。**MPC 场景下这个细节很重要**——
每个控制周期几十毫秒，内存分配会造成不可预测的延迟尖峰。

### 4.9.2 `PrimalDataContainer` / `DualDataContainer`

```cpp
struct PrimalDataContainer {
  PrimalSolution primalSolution;
  std::vector<ModelData> modelDataTrajectory;      // 每个时刻的 LQ
  std::vector<ModelData> modelDataEventTimes;      // 事件时刻的 LQ
  ModelData modelDataFinalTime;
  ProblemMetrics problemMetrics;
};
struct DualDataContainer {
  DualSolution dualSolution;
  std::vector<ScalarFunctionQuadraticApproximation> valueFunctionTrajectory;  // Sm,Sv,s
  std::vector<ModelData> projectedModelDataTrajectory;        // 投影后的 LQ
  std::vector<riccati_modification::Data> riccatiModificationTrajectory;
};
```

`riccati_modification::Data` 存投影矩阵与修正量：

```cpp
scalar_t time_;
matrix_t hamiltonianHessian_;        // H
matrix_t constraintRangeProjector_;  // D†
matrix_t constraintNullProjector_;   // P_u
matrix_t deltaQm_;                   // ΔQ
vector_t deltaGv_;                   // ΔGv
matrix_t deltaGm_;                   // ΔGm
```

### 4.9.3 初始化策略

`initializePrimalSolution()`（`:832`）的三级 fallback：

1. **用上次的控制器做 rollout**（`rolloutInitialController`）
   ——MPC 热启动的主路径
2. **用上次的状态-输入轨迹**（`extractInitialTrajectories`）
   ——控制器不可用（例如上次用了 feedforward）
3. **用 `Initializer`**（`rolloutInitializer`）
   ——前两者都失败或 horizon 延长出的新区间

返回值 `initialSolutionExists` 表示"是否不是纯 Initializer 生成"，
影响首次迭代的 `lqModelExpectedCost` 取值（`:1053`）：

```cpp
const auto lqModelExpectedCost = initialSolutionExists
    ? nominalDualData_.valueFunctionTrajectory.front().f   // Riccati 给的预测代价
    : performanceIndex_.merit;                              // 不可信，退回实际 merit
```

### 4.9.4 Riccati 求解的分区并行

`solveSequentialRiccatiEquationsImpl`（`:516`）：

- **首次迭代串行**（`totalNumIterations_ == 0`）：因为分区并行需要各分区的
  终端值函数作为初值，首次迭代时这些值未知
- **之后并行**：用上次迭代的 $S_m,S_v,s$ 作为各分区右端点的初值猜测

这是一个"用上次结果做预条件"的技巧，在 MPC 里非常有效——
相邻两个控制周期的值函数几乎相同。

### 4.9.5 SLQ 的 Riccati 积分细节

`ContinuousTimeRiccatiEquations` 继承 `OdeBase`，用**归一化时间** $z=-t$
（`:154` 的 `const scalar_t t = -z;`），把后向积分变成前向积分，
以便直接用 odeint 的正向步进器。

**状态向量打包**（`convert2Vector`/`convert2Matrix`，`:55`/`:86`）：
由于 $S_m$ 对称，只存上三角，节省近一半内存与带宽：

$$
\dim = \frac{n_x(n_x+1)}{2} + n_x + 1
$$

`s_vector_dim(nx)` 与 `riccati_matrix_dim(size)` 是这个映射的正逆。

**事件时刻的跳变**：`RiccatiTransversalityConditions.h` 给出值函数穿过事件时的条件。
对跳变映射 $x^+=j(x^-)$，链式法则给出

$$
S_m^- = A_j^{\!\top}S_m^+A_j + \phi^{(i)}_{xx},\quad
S_v^- = A_j^{\!\top}S_v^+ + \phi^{(i)}_x,\quad
s^- = s^+ + \phi^{(i)}
$$

其中 $A_j=\partial j/\partial x$，$\phi^{(i)}$ 是该事件的 pre-jump 代价。

---

## 4.10 风险敏感变体（iLEG）

`isRiskSensitive_` 为真时启用，对应 **Farshidian et al., "Risk sensitive iterative LEG"**。

### 连续时间（`ContinuousTimeRiccatiEquations.cpp:297`）

在标准 Riccati 上加噪声项。设过程噪声 $\mathrm{d}x = f\,\mathrm{d}t + \Sigma^{1/2}\mathrm{d}w$，
风险敏感代价 $\frac{1}{\sigma}\ln\mathbb{E}[e^{\sigma J}]$，则

$$
\begin{aligned}
-\dot S_m &\mathrel{+}= \sigma\,S_m\Sigma S_m\\
-\dot S_v &\mathrel{+}= \sigma\,S_m\Sigma S_v\\
-\dot s   &\mathrel{+}= \tfrac12\operatorname{tr}(\Sigma S_m) + \tfrac{\sigma}{2}S_v^{\!\top}\Sigma S_v
\end{aligned}
$$

```cpp
dSm.noalias() += riskSensitiveCoeff_ * Sm.transpose() * (dynamicsCovariance_ * Sm);
dSv.noalias() += riskSensitiveCoeff_ * (dynamicsCovariance_ * Sm).transpose() * Sv;
ds += 0.5 * (dynamicsCovariance_*Sm).trace() + 0.5*riskSensitiveCoeff_*Sv.dot(dynamicsCovariance_*Sv);
```

$\sigma>0$ 是**风险规避**（惩罚方差），$\sigma<0$ 是风险偏好，$\sigma\to0$ 退化为标准 LQG。

### 离散时间（`DiscreteTimeRiccatiEquations.cpp:159`）

先对下一时刻的值函数做**噪声预修正**，再调标准 iLQR 递推：

$$
\begin{aligned}
\tilde S_m' &= (I-S_m'\Sigma)^{-1}S_m'\\
\tilde S_v' &= (I-S_m'\Sigma)^{-1}S_v'\\
\tilde s'   &= s' + \sigma\,S_v^{\!\top}\Sigma S_v' - \tfrac{1}{2\sigma}\ln\det(I-S_m'\Sigma)
\end{aligned}
$$

```cpp
I_minus_Sm_Sigma_ = I - SmNext * dynamicsCovariance;
Eigen::LDLT<matrix_t> ldltSm(I_minus_Sm_Sigma_);
scalar_t det = ldltSm.vectorD().array().log().sum();     // ln det，用 LDLT 的 D 对角元
inv_I_minus_Sm_Sigma_.setIdentity();  ldltSm.solveInPlace(inv_I_minus_Sm_Sigma_);
SmNextStochastic_ = inv_I_minus_Sm_Sigma_ * SmNext;
...
computeMapILQR(..., SmNextStochastic_, SvNextStochastic_, sNextStochastic_, ...);
```

$(I-S_m\Sigma)^{-1}$ 来自**高斯积分的配方法**：
$\int e^{-\frac12\xi^\top(\Sigma^{-1}-S_m)\xi}\mathrm{d}\xi$ 的闭式解。
$\ln\det$ 项是高斯归一化常数。要求 $I-S_m\Sigma\succ0$，
否则风险敏感代价发散（"neurotic breakdown"）——
代码没有显式检查这一点，LDLT 会给出负的 $D$ 对角元导致 `log` 产生 NaN。

---

## 4.11 `ContinuousTimeLqr.h`：无限时域 LQR

**文件**：`ocs2_ddp/include/ocs2_ddp/ContinuousTimeLqr.h`

独立工具，解代数 Riccati 方程（ARE）：

$$
A^{\!\top}S + SA - (SB+P^{\!\top})R^{-1}(B^{\!\top}S+P) + Q = 0
$$

用于给 MPC 提供终端代价 $\phi(x)=\frac12 x^\top S x$（近似无穷时域尾部），
或作为局部稳定控制器的 fallback。

---

## 4.12 `DDP_Settings` 完全参考

**文件**：`ocs2_ddp/include/ocs2_ddp/DDP_Settings.h`

| 参数 | 默认 | 说明 |
|---|---|---|
| `algorithm_` | `SLQ` | `SLQ` 或 `ILQR` |
| `nThreads_` | 1 | 线程数 |
| `threadPriority_` | 99 | 实时优先级 |
| `maxNumIterations_` | 15 | 最大迭代次数（MPC 里常设 1~5） |
| `minRelCost_` | 1e-3 | 收敛判据：代价相对变化 |
| `constraintTolerance_` | 1e-3 | 收敛判据：约束 SSE |
| `displayInfo_` | false | 每次迭代打印详情 |
| `displayShortSummary_` | false | 结束时打印摘要 |
| `checkNumericalStability_` | true | 检查维度/PSD/NaN（**调试必开，生产可关**） |
| `debugPrintRollout_` | false | 打印 rollout 细节 |
| `absTolODE_` / `relTolODE_` | 1e-9 / 1e-6 | 积分器容差 |
| `maxNumStepsPerSecond_` | 10000 | 防止步长塌缩 |
| `timeStep_` | 1e-2 | iLQR 的固定步长 / 定步长积分器 |
| `backwardPassIntegratorType_` | `ODE45` | Riccati 积分器 |
| `constraintPenaltyInitialValue_` | 2.0 | $\rho_0$（必须 > 1） |
| `constraintPenaltyIncreaseRate_` | 2.0 | $\gamma$（必须 > 1） |
| `preComputeRiccatiTerms_` | true | 启用 reduced form Riccati |
| `useFeedbackPolicy_` | false | 输出 `LinearController`（true）还是 `FeedforwardController`（false） |
| `riskSensitiveCoeff_` | 0.0 | $\sigma$，非零启用 iLEG |
| `strategy_` | `LINE_SEARCH` | `LINE_SEARCH` 或 `LEVENBERG_MARQUARDT` |
| `lineSearch_` | — | 见 `StrategySettings.h` |
| `levenbergMarquardt_` | — | 同上 |

`ddp::loadSettings(filename, "ddp", verbose)` 从 Boost property_tree
格式的 `.info` 文件读取。

### `useFeedbackPolicy_` 的取舍

- `false`（默认）：MRT 只用前馈 $u(t)$，控制频率与 MPC 频率解耦，
  但对模型误差敏感
- `true`：MRT 用 $u=K(t)x+b(t)$，**高频反馈**（1kHz）跟踪低频 MPC（100Hz）的规划，
  抗扰动能力强得多

Grandia et al. (2019) "Feedback MPC for Torque-Controlled Legged Robots"
详细讨论了这个设计。四足机器人几乎必须开 `true`。

---

## 4.13 SLQ vs iLQR 对比小结

| | SLQ | iLQR |
|---|---|---|
| 时域 | 连续 | 离散 |
| Riccati | 微分方程，ODE 积分（自适应步长） | 差分方程，逐点递推 |
| $H$ | $R$ | $R+B^\top S_m'B$ |
| 时间网格 | 由积分器自适应决定 | 固定 `timeStep_` |
| 事件处理 | 积分中断 + 横截条件 | 索引跳转 |
| 精度 | 高（自适应） | 取决于 `timeStep_` |
| 速度 | 慢（ODE 求解开销） | 快 |
| 何时用 | 刚性系统、需要高精度 | 实时 MPC、动力学温和 |

代码复用度极高：`GaussNewtonDDP` 实现了 90% 的逻辑，
`SLQ`/`ILQR` 各自只实现 4 个虚函数
（`computeHamiltonianHessian`、`approximateIntermediateLQ`、
`calculateControllerWorker`、`riccatiEquationsWorker`）。

---

## 4.14 参考文献

见 [17 章](17-references.md) 的完整列表，本章直接相关的：

- **Farshidian et al. (2017)**, *Sequential Linear Quadratic Optimal Control for Nonlinear Switched Systems*, IFAC — **SLQ 原始论文**
- **Farshidian et al. (2017)**, *An Efficient Optimal Planning and Control Framework for Quadrupedal Locomotion*, ICRA
- **Sleiman, Farshidian, Hutter (2021)**, *Constraint Handling in Continuous-Time DDP-Based MPC*, ICRA — **本章的投影 + 增广拉格朗日方案**
- **Todorov & Li (2005)**, *A Generalized Iterative LQG Method*, ACC — iLQG/iLQR
- **Tassa, Erez, Todorov (2012)**, *Synthesis and Stabilization of Complex Behaviors through Online Trajectory Optimization*, IROS — DDP 正则化与线搜索
- **Jacobson & Mayne (1970)**, *Differential Dynamic Programming* — DDP 原典
- **Nocedal & Wright (2006)**, *Numerical Optimization*, 2nd ed. — §3.1（Armijo）、§4.1（信赖域）、§17.2（精确罚）、§19（内点法）
- **Gill, Murray, Wright (1981)**, *Practical Optimization* — 修正 Cholesky
- **Conn, Gould, Toint (1991)**, *A Globally Convergent Augmented Lagrangian Algorithm*, SIAM J. Numer. Anal. — 罚系数更新策略

**下一章**：多重打靶 SQP → [05 ocs2_sqp](05-ocs2-sqp.md)
