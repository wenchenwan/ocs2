# 08 · 约束处理：罚函数与增广拉格朗日

OCS2 提供**四种**约束处理机制，可以混用：

| 机制 | 适用 | 精确性 | 实现位置 |
|---|---|---|---|
| **零空间投影** | 状态-输入等式约束（$D$ 满行秩） | 精确 | [04](04-ocs2-ddp.md) §4.4、[03](03-ocs2-oc.md) §3.5.3 |
| **软约束（罚函数）** | 任意约束 | 近似（需大 $\mu$） | `ocs2_core/soft_constraint/` + `penalties/penalties/` |
| **增广拉格朗日** | 任意约束 | **渐近精确**（有限 $\rho$） | `ocs2_core/augmented_lagrangian/` + `penalties/augmented/` |
| **QP 内建 / 内点法** | 不等式约束 | 精确 | [05](05-ocs2-sqp.md)、[06](06-ocs2-ipm.md) |

本章推导第二、三类的全部数学。

---

## 8.1 为什么需要罚函数与增广拉格朗日

DDP（SLQ/iLQR）的 Riccati 递推**天然只能处理无约束问题**。
零空间投影解决了状态-输入等式约束，但：

- **仅状态的约束**（$g(x)=0$）无法投影——输入不出现在约束里
- **不等式约束**根本无法投影——激活集未知

所以 DDP 必须把这两类约束转成代价。SQP/IPM 虽然有 QP 求解器，
但对**软性偏好**（如"尽量远离关节限位"）用罚函数比硬约束更合适。

---

## 8.2 普通罚函数（软约束）

**接口**：`ocs2_core/penalties/penalties/PenaltyBase.h`

```cpp
class PenaltyBase {
  virtual scalar_t getValue(scalar_t t, scalar_t h) const = 0;
  virtual scalar_t getDerivative(scalar_t t, scalar_t h) const = 0;
  virtual scalar_t getSecondDerivative(scalar_t t, scalar_t h) const = 0;
};
```

对标量约束值 $h$ 返回罚值与导数。`MultidimensionalPenalty` 逐分量应用到向量约束。

### 8.2.1 二次罚（等式约束）

**文件**：`penalties/penalties/QuadraticPenalty.h`

$$
p(h)=\frac{\mu}{2}h^2,\qquad p'(h)=\mu h,\qquad p''(h)=\mu
$$

**特点**：光滑、Hessian 恒为常数 $\mu>0$（对 Riccati 友好）。
**缺点**：非精确罚——最优解满足 $\mu h^\star = -\lambda^\star$，
即 $h^\star = -\lambda^\star/\mu\ne0$。要精确需 $\mu\to\infty$，但那会让 Hessian 病态。

这正是增广拉格朗日要解决的问题（§8.3）。

### 8.2.2 ⭐ 松弛障碍罚（不等式约束）

**文件**：`penalties/penalties/RelaxedBarrierPenalty.h`

对 $h\ge0$：

$$
\boxed{\;p(h)=\begin{cases}
-\mu\ln h & h>\delta\\[6pt]
-\mu\ln\delta+\dfrac{\mu}{2}\left[\left(\dfrac{h-2\delta}{\delta}\right)^{2}-1\right] & h\le\delta
\end{cases}\;}
$$

**动机**：纯对数障碍 $-\mu\ln h$ 在 $h\le0$ 时无定义（$+\infty$）。
但 DDP 的线搜索可能产生临时违反约束的迭代点，一旦 $h<0$ 整个算法崩溃。
**松弛障碍**在 $h=\delta$ 处切换为二次延拓，使函数在整个 $\mathbb{R}$ 上有定义且 $C^2$。

#### ⭐ 验证 $C^2$ 连续性

记 $q(h)=-\mu\ln\delta+\frac{\mu}{2}\left[\left(\frac{h-2\delta}{\delta}\right)^2-1\right]$。

**值**：在 $h=\delta$：
$$
q(\delta)=-\mu\ln\delta+\frac{\mu}{2}\left[\left(\frac{\delta-2\delta}{\delta}\right)^2-1\right]
=-\mu\ln\delta+\frac{\mu}{2}\big[(-1)^2-1\big]=-\mu\ln\delta
$$
而 $-\mu\ln h\big|_{h=\delta}=-\mu\ln\delta$ ✓

**一阶导**：
$$
q'(h)=\frac{\mu}{2}\cdot 2\cdot\frac{h-2\delta}{\delta}\cdot\frac1\delta=\frac{\mu(h-2\delta)}{\delta^2}
$$
$$
q'(\delta)=\frac{\mu(\delta-2\delta)}{\delta^2}=-\frac{\mu}{\delta}
$$
而 $\dfrac{\mathrm{d}}{\mathrm{d}h}(-\mu\ln h)\Big|_{h=\delta}=-\dfrac{\mu}{\delta}$ ✓

**二阶导**：
$$
q''(h)=\frac{\mu}{\delta^2}\quad\text{（常数）},\qquad q''(\delta)=\frac{\mu}{\delta^2}
$$
而 $\dfrac{\mathrm{d}^2}{\mathrm{d}h^2}(-\mu\ln h)\Big|_{h=\delta}=\dfrac{\mu}{h^2}\Big|_{h=\delta}=\dfrac{\mu}{\delta^2}$ ✓

**三处全部匹配**——这个看似神秘的公式（尤其是 $h-2\delta$ 与 $-1$）
正是为了同时满足这三个条件而构造的。

#### 为什么是 $h-2\delta$ 而不是 $h-\delta$？

因为要匹配一阶导 $-\mu/\delta$。若用 $\frac{\mu}{2}\left(\frac{h-\delta}{\delta}\right)^2$，
则 $q'(\delta)=0\ne-\mu/\delta$。
偏移 $2\delta$ 使得二次函数的顶点在 $h=2\delta$，
从而在 $h=\delta$ 处有恰好正确的负斜率。

#### 性质

- $p''(h)=\mu/\delta^2>0$ 在延拓区，$p''(h)=\mu/h^2>0$ 在障碍区 → **全局凸** ✓
- $h\to0^+$ 时 $p$ 有限（$=q(0)=-\mu\ln\delta+\frac{\mu}{2}(4-1)=-\mu\ln\delta+\frac{3\mu}{2}$）
  → 不会数值溢出
- $h<0$ 时 $p$ 二次增长 → 强烈推回可行域

参考：Feller & Ebenbauer (2017), *Relaxed Logarithmic Barrier Function Based Model Predictive Control*, IEEE TAC。
OCS2 的 Grandia et al. (2019) *Feedback MPC for Torque-Controlled Legged Robots* 用它处理摩擦锥。

### 8.2.3 平方铰链罚（不等式约束）

**文件**：`penalties/penalties/SquaredHingePenalty.h`

$$
p(h)=\begin{cases}\dfrac{\mu}{2}(h-\delta)^2 & h<\delta\\[4pt] 0 & h\ge\delta\end{cases}
$$

$$
p'(h)=\begin{cases}\mu(h-\delta)&h<\delta\\0&\text{否则}\end{cases}\qquad
p''(h)=\begin{cases}\mu&h<\delta\\0&\text{否则}\end{cases}
$$

**$C^1$ 但非 $C^2$**（$p''$ 在 $h=\delta$ 跳变）。对 DDP 影响不大（Hessian 只需 PSD）。

**优点**：约束满足时罚为 0，不引入偏差；计算极快。
**缺点**：非精确罚（同二次罚）；$\delta>0$ 提供安全裕度但也引入保守性。

默认 $\mu=100$、$\delta=0.1$。

### 8.2.4 光滑绝对值罚（等式约束）

**文件**：`penalties/penalties/SmoothAbsolutePenalty.h`

$$
p(h)=\mu\sqrt{h^2+\delta^2}
$$
$$
p'(h)=\frac{\mu h}{\sqrt{h^2+\delta^2}},\qquad
p''(h)=\frac{\mu\delta^2}{(h^2+\delta^2)^{3/2}}
$$

**推导 $p''$**：
$$
p'=\mu h(h^2+\delta^2)^{-1/2}
$$
$$
p''=\mu(h^2+\delta^2)^{-1/2}+\mu h\cdot\left(-\tfrac12\right)(h^2+\delta^2)^{-3/2}\cdot 2h
=\mu\frac{(h^2+\delta^2)-h^2}{(h^2+\delta^2)^{3/2}}=\frac{\mu\delta^2}{(h^2+\delta^2)^{3/2}}\ \checkmark
$$

**这是 $\ell_1$ 罚 $\mu|h|$ 的光滑化**（也叫 pseudo-Huber）。误差界（头文件注释）：

$$
\big|\,|h|-\sqrt{h^2+\delta^2}\,\big|\le\delta\quad\forall h
$$

**证明**：设 $u=|h|$。$\sqrt{u^2+\delta^2}-u = \dfrac{\delta^2}{\sqrt{u^2+\delta^2}+u}$。
该式在 $u=0$ 取最大值 $\delta$，随 $u$ 单调递减到 0 ✓

**为什么重要**：$\ell_1$ 罚是**精确罚**（存在有限 $\mu$ 使罚问题的解 = 原问题的解，
条件是 $\mu>\|\lambda^\star\|_\infty$，见 Nocedal & Wright 定理 17.3）。
光滑化保留了这个性质（误差 $O(\delta)$），同时让 DDP 可用。

**注意 $p''\to0$ 当 $|h|\gg\delta$**：远离约束时 Hessian 趋于 0，
这与二次罚（$p''=\mu$ 恒定）行为很不同。好处是不会过度扭曲代价景观。

### 8.2.5 双边罚（箱约束）

**文件**：`penalties/penalties/DoubleSidedPenalty.h`

对 $l\le h\le u$：

$$
p_{box}(h)=p(h-l)+p(u-h)
$$
$$
p_{box}'(h)=p'(h-l)-p'(u-h),\qquad
p_{box}''(h)=p''(h-l)+p''(u-h)
$$

（链式法则：$\frac{\mathrm{d}}{\mathrm{d}h}(u-h)=-1$，故一阶导相减、二阶导相加）

配合 `StateInputSoftBoxConstraint` 可以一行代码给状态/输入加上下界，
无需写约束类。

---

## 8.3 ⭐ 增广拉格朗日

### 8.3.1 为什么二次罚不够

考虑 $\min f(x)$ s.t. $g(x)=0$。二次罚问题：

$$
\min_x\ f(x)+\frac{\rho}{2}g(x)^2
$$

一阶条件：$\nabla f+\rho\,g\,\nabla g=0$。而原问题的 KKT 是 $\nabla f+\lambda\nabla g=0$。
对比得 $\rho\,g^\star=\lambda^\star$，即

$$
g^\star=\frac{\lambda^\star}{\rho}\ne0\quad\text{除非}\ \rho\to\infty
$$

$\rho\to\infty$ 会让 Hessian 条件数 $\to\infty$，数值上不可行。

### 8.3.2 增广拉格朗日的思想

引入**乘子估计** $\lambda$，构造

$$
\boxed{\;L_A(x,\lambda)=f(x)-\lambda g(x)+\frac{\rho}{2}g(x)^2\;}
$$

（OCS2 用 $-\lambda g$ 的符号约定，见 `augmented/QuadraticPenalty.h:81`）

一阶条件：$\nabla f-\lambda\nabla g+\rho g\nabla g=0$，即

$$
\nabla f + (\rho g-\lambda)\nabla g=0
$$

与 KKT 对比得 $\lambda^\star_{\text{true}}=\lambda-\rho g^\star$。
于是**若 $\lambda$ 已经等于真实乘子，则 $g^\star=0$ 精确成立**，
无需 $\rho\to\infty$！

**对偶更新**（一阶乘子法 / Hestenes–Powell）：

$$
\lambda^{k+1}=\lambda^k-\rho\, g(x^{k+1})
$$

这是**对偶函数的梯度上升**：对偶函数 $d(\lambda)=\min_x L_A(x,\lambda)$
满足 $\nabla_\lambda d = -g(x^\star(\lambda))$，
步长取 $\rho$ 时收敛率最优（Bertsekas 1982）。

### 8.3.3 OCS2 的接口

**文件**：`ocs2_core/penalties/augmented/AugmentedPenaltyBase.h`

```cpp
class AugmentedPenaltyBase {
  virtual scalar_t getValue(scalar_t t, scalar_t l, scalar_t h) const = 0;
  virtual scalar_t getDerivative(scalar_t t, scalar_t l, scalar_t h) const = 0;
  virtual scalar_t getSecondDerivative(scalar_t t, scalar_t l, scalar_t h) const = 0;
  virtual scalar_t updateMultiplier(scalar_t t, scalar_t l, scalar_t h) const = 0;
  virtual scalar_t initializeMultiplier() const = 0;
};
```

比 `PenaltyBase` 多两个方法：`updateMultiplier`（对偶上升）与
`initializeMultiplier`（乘子初值）。

**数据流**（`ocs2_oc/src/oc_problem/OptimalControlProblemHelperFunction.cpp`）：

```
求解器迭代 i:
  ① approximateIntermediateLQ(..., multipliers, ...)   ← 用当前 λ 构造增广代价
  ② 求解 LQ 子问题 → 新的 primal solution
  ③ updateDualSolution(...) → 对每项调 updateMultiplier(t, λ, h_new) → 新的 λ
```

### 8.3.4 增广二次罚（等式约束）

**文件**：`augmented/QuadraticPenalty.h`

$$
p(h,\lambda)=-\lambda h+\frac{\rho}{2}h^2
$$
$$
\frac{\partial p}{\partial h}=-\lambda+\rho h,\qquad
\frac{\partial^2 p}{\partial h^2}=\rho
$$

**对偶更新**（`:85`）：

$$
\lambda^{+}=\lambda-\alpha\,\rho\,h
$$

```cpp
scalar_t updateMultiplier(scalar_t t, scalar_t l, scalar_t h) const override {
  return l - config_.stepSize * config_.scale * h;
}
scalar_t initializeMultiplier() const override { return 0.0; }
```

`stepSize` = $\alpha$（默认 0，即**不更新乘子**，退化为纯二次罚！）。
用户需显式设 $\alpha>0$ 才启用增广拉格朗日。$\alpha=1$ 对应标准的
$\lambda^+=\lambda-\rho h$。

### 8.3.5 ⭐ PHR 罚（不等式约束）

**文件**：`augmented/SlacknessSquaredHingePenalty.h`

**PHR** = Powell–Hestenes–Rockafellar。处理 $h\ge0$。

#### 从松弛变量推导

引入 $s\ge0$，把 $h\ge0$ 写成

$$
h-s=0,\qquad s\ge0
$$

对等式部分用增广拉格朗日：

$$
L_A(x,s,\lambda)=f(x)-\lambda(h-s)+\frac{\rho}{2}(h-s)^2
$$

**对 $s\ge0$ 显式最小化**（这是 PHR 的关键步骤）。
记 $\psi(s)=-\lambda(h-s)+\frac{\rho}{2}(h-s)^2$，

$$
\frac{\partial\psi}{\partial s}=\lambda-\rho(h-s)=0
\;\Longrightarrow\; s^\star=h-\frac{\lambda}{\rho}
$$

考虑 $s\ge0$ 约束：

$$
s^\star=\max\left(0,\ h-\frac{\lambda}{\rho}\right)
$$

**情形 A：$h\ge\lambda/\rho$**（约束"松"）→ $s^\star=h-\lambda/\rho$，此时 $h-s^\star=\lambda/\rho$：

$$
p=-\lambda\cdot\frac{\lambda}{\rho}+\frac{\rho}{2}\cdot\frac{\lambda^2}{\rho^2}
=-\frac{\lambda^2}{\rho}+\frac{\lambda^2}{2\rho}=-\frac{\lambda^2}{2\rho}
$$

**常数！** 与 $h$ 无关 → 导数为 0。

**情形 B：$h<\lambda/\rho$**（约束"紧"）→ $s^\star=0$，$h-s^\star=h$：

$$
p=-\lambda h+\frac{\rho}{2}h^2
$$

**合并**：

$$
\boxed{\;p(h,\lambda)=\begin{cases}
-\lambda h+\dfrac{\rho}{2}h^2 & h<\dfrac{\lambda}{\rho}\\[8pt]
-\dfrac{\lambda^2}{2\rho} & h\ge\dfrac{\lambda}{\rho}
\end{cases}\;}
$$

**这正是代码**（`SlacknessSquaredHingePenalty.h:85-91`）：

```cpp
scalar_t getValue(scalar_t t, scalar_t l, scalar_t h) const override {
  return (h < l / config_.scale) ? (-l*h + 0.5*config_.scale*h*h) : (-0.5*l*l/config_.scale);
}
scalar_t getDerivative(scalar_t t, scalar_t l, scalar_t h) const override {
  return (h < l / config_.scale) ? (-l + config_.scale*h) : 0.0;
}
scalar_t getSecondDerivative(scalar_t t, scalar_t l, scalar_t h) const override {
  return (h < l / config_.scale) ? config_.scale : 0.0;
}
```

#### 等价的紧凑形式

头文件注释给出

$$
p(h,\lambda)=\frac{1}{2\rho}\Big(\max\{0,\ \lambda-\rho h\}^2-\lambda^2\Big)
$$

**验证等价**：
- 若 $\lambda-\rho h>0$（即 $h<\lambda/\rho$）：
  $$
  \frac{1}{2\rho}\big[(\lambda-\rho h)^2-\lambda^2\big]
  =\frac{1}{2\rho}\big[\lambda^2-2\lambda\rho h+\rho^2h^2-\lambda^2\big]
  =-\lambda h+\frac{\rho}{2}h^2\ \checkmark
  $$
- 若 $\lambda-\rho h\le0$：$\frac{1}{2\rho}(0-\lambda^2)=-\frac{\lambda^2}{2\rho}$ ✓

#### $C^1$ 光滑性

在 $h=\lambda/\rho$ 处：$p'=-\lambda+\rho\cdot\frac{\lambda}{\rho}=0$，与右侧的 0 匹配 ✓
$p''$ 从 $\rho$ 跳到 0（$C^1$ 但非 $C^2$）。

#### 对偶更新

$$
\lambda^{+}=\max\Big(0,\ \max\big(\lambda-\alpha\rho h,\ (1-\alpha)\lambda\big)\Big)
$$

```cpp
scalar_t updateMultiplier(scalar_t t, scalar_t l, scalar_t h) const override {
  return std::max(0.0, std::max(l - config_.stepSize*config_.scale*h, (1.0 - config_.stepSize)*l));
}
```

**三层保护**：
1. **外层 $\max(0,\cdot)$**：保证 $\lambda\ge0$（不等式乘子的符号约束）
2. **$\lambda-\alpha\rho h$**：标准对偶上升（$h<0$ 时 $\lambda$ 增大 → 加强惩罚）
3. **$(1-\alpha)\lambda$**：**下界保护**——当 $h$ 很大（约束远离激活）时，
   标准更新会让 $\lambda$ 骤降到 0，导致下一轮突然违反约束时反应不及。
   $(1-\alpha)\lambda$ 让 $\lambda$ 至多按几何速率衰减，保留"记忆"。

默认 $\rho=10$、$\alpha=1$、$\lambda_0=0$。

### 8.3.6 ⭐ 光滑 PHR（修正松弛障碍）

**文件**：`augmented/ModifiedRelaxedBarrierPenalty.h`

PHR 只有 $C^1$，Hessian 跳变对 DDP 不友好。这个变体给出 $C^2$ 版本。

#### 构造

定义辅助量：

$$
v(h,\lambda)=\frac{\rho h}{\lambda},\qquad
w(\lambda)=\frac{\lambda^2}{\rho},\qquad
\frac{\partial v}{\partial h}=\frac{\rho}{\lambda}
$$

罚函数：

$$
\boxed{\;p(h,\lambda)=w(\lambda)\,\psi\big(v(h,\lambda)\big)
=\frac{\lambda^2}{\rho}\,\psi\!\left(\frac{\rho h}{\lambda}\right)\;}
$$

其中 $\psi$ 是**平移的松弛对数障碍**：

$$
\psi(v)=\begin{cases}
-\ln(1+v) & v>\bar\delta\\[4pt]
\tfrac12 c_2(v-\bar\delta)^2+c_1(v-\bar\delta)+c_0 & v\le\bar\delta
\end{cases}
$$

系数由 $C^2$ 匹配条件确定（`ModifiedRelaxedBarrierPenalty.h` 的 `QuadCoeff`）：

$$
c_0=-\ln(1+\bar\delta),\qquad
c_1=-\frac{1}{1+\bar\delta},\qquad
c_2=\frac{1}{(1+\bar\delta)^2}
$$

**验证**：
- $\psi(\bar\delta)$：二次式给 $c_0=-\ln(1+\bar\delta)$；对数式给 $-\ln(1+\bar\delta)$ ✓
- $\psi'(v)=-\frac{1}{1+v}$，在 $\bar\delta$ 处 $=-\frac{1}{1+\bar\delta}=c_1$；
  二次式导数在 $\bar\delta$ 处 $=c_1$ ✓
- $\psi''(v)=\frac{1}{(1+v)^2}$，在 $\bar\delta$ 处 $=\frac{1}{(1+\bar\delta)^2}=c_2$；
  二次式二阶导恒为 $c_2$ ✓

**定义域**：$\psi$ 在 $v>-1$ 上有定义（对数），
故要求 $\bar\delta>-1$（头文件注释明确指出这一点）。

#### 导数

链式法则：

$$
\frac{\partial p}{\partial h}=w(\lambda)\,\psi'(v)\,\frac{\rho}{\lambda},\qquad
\frac{\partial^2 p}{\partial h^2}=w(\lambda)\,\psi''(v)\left(\frac{\rho}{\lambda}\right)^{2}
$$

代码（`:86-113`）：

```cpp
scalar_t getValue(scalar_t t, scalar_t l, scalar_t h) const override {
  const scalar_t v = vFunc(l, h);                       // ρh/λ
  if (v > config_.relaxation) return -wFunc(l) * log(1.0 + v);
  else { const scalar_t vd = v - config_.relaxation;
         return wFunc(l) * (0.5*quadCoeff_.c2*vd*vd + quadCoeff_.c1*vd + quadCoeff_.c0); }
}
scalar_t getDerivative(...) const override {
  const scalar_t v = vFunc(l, h);
  if (v > config_.relaxation) return -wFunc(l) / (1.0 + v) * dvdhFunc(l);      // w·ψ'·(ρ/λ)
  else return wFunc(l) * (quadCoeff_.c2*(v-relaxation) + quadCoeff_.c1) * dvdhFunc(l);
}
scalar_t getSecondDerivative(...) const override {
  const scalar_t v = vFunc(l, h);  const scalar_t dvdh = dvdhFunc(l);
  if (v > config_.relaxation) return wFunc(l) / ((1.0+v)*(1.0+v)) * dvdh * dvdh;  // w·ψ''·(ρ/λ)²
  else return wFunc(l) * quadCoeff_.c2 * dvdh * dvdh;
}
```

#### 为什么这个缩放（$w=\lambda^2/\rho$、$v=\rho h/\lambda$）？

**关键性质**：$p$ 的行为与 $\lambda$ 自适应。考察 $\rho h\ll\lambda$（$v\approx0$）时：

$$
p\approx\frac{\lambda^2}{\rho}\left[-\ln(1+v)\right]
\approx\frac{\lambda^2}{\rho}\left[-v+\frac{v^2}{2}\right]
=\frac{\lambda^2}{\rho}\left[-\frac{\rho h}{\lambda}+\frac{\rho^2h^2}{2\lambda^2}\right]
=-\lambda h+\frac{\rho}{2}h^2
$$

**这正是 PHR 在紧约束区的形式！** 所以光滑 PHR 是 PHR 的 $C^2$ 光滑化，
在 $h$ 小的时候两者一致。

而 $h$ 大时（$v\gg1$），$p\approx-\frac{\lambda^2}{\rho}\ln v\to$ 缓慢增长的负值，
对应 PHR 的常数 $-\lambda^2/(2\rho)$ 的光滑版本。

#### 对偶更新

$$
\lambda^{+}=\max\left(\lambda_{\min},\ -\alpha\,\lambda\,\psi'(v)\cdot\frac{\rho}{\lambda}\cdot\frac{\lambda}{\rho}\right)
$$

代码（`:115-124`）：

```cpp
scalar_t updateMultiplier(scalar_t t, scalar_t l, scalar_t h) const override {
  const scalar_t v = vFunc(l, h);
  constexpr scalar_t lambdaMin = 1e-4;
  if (v > config_.relaxation)
    return std::max(lambdaMin, wFunc(l) * dvdhFunc(l) / (1 + v));
  else
    return std::max(lambdaMin, config_.stepSize * wFunc(l)
                    * (-quadCoeff_.c2*(v-relaxation) - quadCoeff_.c1) * dvdhFunc(l));
}
scalar_t initializeMultiplier() const override { return 1.0; }
```

在障碍区：$w\cdot\frac{\rho}{\lambda}\cdot\frac{1}{1+v}
=\frac{\lambda^2}{\rho}\cdot\frac{\rho}{\lambda}\cdot\frac{1}{1+v}
=\frac{\lambda}{1+v}=\frac{\lambda}{1+\rho h/\lambda}=\frac{\lambda^2}{\lambda+\rho h}$。

**这是 $-\partial p/\partial h$ 的值**，即"当前罚函数施加的等效约束力"。
用它作为新的乘子估计，正是对偶上升的自然选择
（KKT 条件要求 $\lambda^\star = \partial f/\partial h$）。

**$\lambda_{\min}=10^{-4}$ 硬编码**：因为 $v=\rho h/\lambda$ 与 $w=\lambda^2/\rho$
在 $\lambda\to0$ 时会除零。**这也意味着 $\lambda$ 永不为 0**，
即光滑 PHR 始终保持一个最小的障碍强度——这是它与 PHR 的一个实质差异。

**$\lambda_0=1.0$**（而非 0）也是因为除零。

### 8.3.7 增广光滑绝对值罚

**文件**：`augmented/SmoothAbsolutePenalty.h`

$$
p(h,\lambda)=-\lambda h+\mu\sqrt{h^2+\delta^2}
$$

即 §8.2.4 的罚 + 线性拉格朗日项。导数：

$$
\frac{\partial p}{\partial h}=-\lambda+\frac{\mu h}{\sqrt{h^2+\delta^2}},\qquad
\frac{\partial^2p}{\partial h^2}=\frac{\mu\delta^2}{(h^2+\delta^2)^{3/2}}
$$

对偶更新 $\lambda^+=\lambda-\alpha\mu h$，$\lambda_0=0$。

**特点**：结合了 $\ell_1$ 精确罚与增广拉格朗日。
理论上收敛最快（$\ell_1$ 已精确，乘子只需修正 $O(\delta)$ 的偏差），
但 $p''$ 在 $|h|$ 大时趋于 0，Riccati 可能需要额外的 Hessian 修正。

---

## 8.4 四种增广罚的对比

| 罚函数 | 约束类型 | 光滑性 | $p''$ 行为 | $\lambda_0$ | 适用 |
|---|---|---|---|---|---|
| `QuadraticPenalty` | $h=0$ | $C^\infty$ | 常数 $\rho$ | 0 | 等式约束的默认选择 |
| `SmoothAbsolutePenalty` | $h=0$ | $C^\infty$ | $\to0$ 当 $\|h\|$ 大 | 0 | 需要精确罚的等式约束 |
| `SlacknessSquaredHingePenalty` | $h\ge0$ | $C^1$ | $\rho$ 或 0（跳变） | 0 | 不等式约束的默认 |
| `ModifiedRelaxedBarrierPenalty` | $h\ge0$ | $C^2$ | 光滑过渡 | 1.0 | 需要 $C^2$ 的 DDP |

**选择建议**：
- **SQP/IPM**：只需 $C^1$，用 PHR 即可
- **SLQ/iLQR**：Riccati 需要 $C^2$ 的 Hessian，用 `ModifiedRelaxedBarrierPenalty`
- **等式约束优先用零空间投影**，投影不了（仅状态约束）才用增广拉格朗日

---

## 8.5 `MultidimensionalPenalty`

**文件**：`ocs2_core/penalties/MultidimensionalPenalty.h`

把标量罚函数应用到向量约束 $h\in\mathbb{R}^m$：

$$
P(h)=\sum_{i=1}^{m}p_i(h_i)
$$

支持**每个分量用不同的罚函数**（构造时传 `std::vector<std::unique_ptr<PenaltyBase>>`），
也支持所有分量共用一个。

二次近似的组装（配合约束的一二阶导）：

$$
\frac{\partial P}{\partial z}=\sum_i p_i'(h_i)\frac{\partial h_i}{\partial z}
= J_h^{\!\top}\,p'
$$
$$
\frac{\partial^2P}{\partial z^2}
=J_h^{\!\top}\operatorname{diag}(p'')J_h+\sum_i p_i'(h_i)\frac{\partial^2h_i}{\partial z^2}
$$

第一项是**高斯-牛顿部分**（PSD，因为 $p''\ge0$），
第二项需要约束的 Hessian（`ConstraintOrder::Quadratic`），**可能不定**。
这就是为什么很多约束只提供 `Linear` 阶——丢掉第二项换取 PSD 保证。

---

## 8.6 增广拉格朗日的完整数据流

以 DDP 为例：

```
迭代 i:
 ┌─ approximateOptimalControlProblem()
 │    └─ approximateIntermediateLQ(problem, t, x, u, multipliers[i], modelData)
 │         └─ 对每个 Lagrangian 容器:
 │              approx = collection->getQuadraticApproximation(t, x, u, multipliers.xxx, preComp)
 │                  内部: h = constraint->getQuadraticApproximation(...)
 │                        p = penalty->getValue/Derivative/SecondDerivative(t, λ, h)
 │                        组装 ∂P/∂z, ∂²P/∂z²
 │              modelData.cost += approx           ← 罚被吸收进代价
 │
 ├─ solveSequentialRiccatiEquations()               ← Riccati 完全不知道约束存在
 ├─ calculateController()
 ├─ takePrimalDualStep()
 │    └─ searchStrategy->run(...)                   ← 线搜索（merit 含 Lagrangian 项）
 │    └─ updateDualSolution(problem, primalSol, metrics, dualSol)
 │         └─ 对每项调 penalty->updateMultiplier(t, λ, h_new) → λ^{i+1}
 └─ updateConstraintPenalties()                     ← 只更新投影约束的 ρ（见 04 章 §4.7.2）
```

**注意两套罚系数**：
1. 增广拉格朗日各项的 $\rho$：**固定**（由用户在构造 penalty 时给定），
   靠乘子更新达到精确
2. `constraintPenaltyCoefficients_.penaltyCoeff`：**自适应**，
   用于 merit 函数里的 $\rho\sqrt{\text{SSE}}$ 项（[04 章](04-ocs2-ddp.md) §4.7.2）

---

## 8.7 实践建议

### 什么时候用哪种机制

| 约束 | 推荐机制 |
|---|---|
| 四足支撑腿零速度 $J_cv=0$ | **零空间投影**（$D$ 满行秩） |
| 摆动腿零接触力 $f=0$ | **零空间投影**（直接约束输入分量） |
| 摩擦锥 $\mu F_z-\|F_t\|\ge0$ | **软约束**（松弛障碍）或 **增广拉格朗日**（PHR） |
| 关节限位 $q_{\min}\le q\le q_{\max}$ | **软约束**（双边 + 松弛障碍） |
| 自碰撞距离 $d\ge d_{\min}$ | **软约束**（松弛障碍） |
| 末端位姿终端约束 | **增广拉格朗日**（等式，仅状态） |
| 力矩限 $\|\tau\|\le\tau_{\max}$ | **SQP/IPM 硬约束** 或 软约束 |

### 调参顺序

1. 先用**软约束**（罚函数）跑通，$\mu$ 从小到大试
2. 若约束违反不可接受，改用**增广拉格朗日**（$\rho$ 取软约束时的 $\mu$，$\alpha=1$）
3. 若仍不满足，检查是否问题本身不可行（约束冲突）
4. 硬约束是最后手段——它会显著增加求解时间且可能导致 QP 不可行

### 常见陷阱

- **`stepSize=0`**：`augmented::QuadraticPenalty` 的默认值，此时乘子不更新，
  增广拉格朗日退化为普通二次罚。这是最常见的"为什么增广拉格朗日没效果"的原因。
- **松弛参数 $\delta$ 太小**：松弛障碍在 $h\approx0$ 处 Hessian $\approx\mu/\delta^2$ 巨大，
  导致 Riccati 病态。建议 $\delta\ge0.01\times$（约束的典型量级）。
- **`ModifiedRelaxedBarrierPenalty` 的 $\bar\delta\le-1$**：定义域外，会产生 NaN。

---

## 8.8 参考文献

- **Hestenes (1969)**, *Multiplier and Gradient Methods*, JOTA
- **Powell (1969)**, *A Method for Nonlinear Constraints in Minimization Problems*
- **Rockafellar (1974)**, *Augmented Lagrange Multiplier Functions and Duality in Nonconvex Programming*, SIAM J. Control —— **PHR 罚的原始来源**
- **Bertsekas (1982)**, *Constrained Optimization and Lagrange Multiplier Methods* —— 增广拉格朗日的标准教材
- **Conn, Gould, Toint (1992)**, *LANCELOT: A Fortran Package for Large-Scale Nonlinear Optimization* —— 罚系数/容差更新策略
- **Feller & Ebenbauer (2017)**, *Relaxed Logarithmic Barrier Function Based Model Predictive Control of Linear Systems*, IEEE TAC 62(3) —— **松弛障碍**
- **Sleiman, Farshidian, Hutter (2021)**, *Constraint Handling in Continuous-Time DDP-Based Model Predictive Control*, ICRA —— **OCS2 的 PHR / 光滑 PHR 方案**（代码注释里的 "the corresponding paper" 指的就是它）
- **Grandia et al. (2019)**, *Feedback MPC for Torque-Controlled Legged Robots*, IROS —— 松弛障碍在四足上的应用
- **Nocedal & Wright (2006)**, *Numerical Optimization*, Ch. 17（罚函数与增广拉格朗日）

**下一章**：频域整形 → [09 Loopshaping](09-loopshaping.md)
