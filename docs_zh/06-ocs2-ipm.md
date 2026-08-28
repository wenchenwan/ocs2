# 06 · ocs2_ipm：多重打靶 + 原对偶内点法

`ocs2_ipm` 与 `ocs2_sqp` 共享全部多重打靶基础设施
（`multiple_shooting::Transcription`、`FilterLinesearch`、HPIPM），
区别只在**如何处理不等式约束**：

| | SQP | IPM |
|---|---|---|
| 不等式约束 | 交给 HPIPM 内部的内点法 | 在 OCS2 层面用障碍函数处理，喂给 HPIPM 的是**无不等式约束**的 QP |
| 松弛/对偶变量 | HPIPM 内部管理 | OCS2 显式管理（`slackStateIneq`、`dualStateIneq` 等） |
| 障碍参数 | HPIPM 内部 | OCS2 显式调度（`barrierParam`） |

**为什么要在外层再做一次内点法？**
因为这样可以在**非线性层面**做内点：每次迭代重新线性化约束，
松弛变量与对偶变量跨迭代保持连续，收敛更平滑。
而 SQP 把不等式交给 QP 后，每次 QP 都从头开始处理不等式，信息不复用。

设计上明确对标 **IPOPT**（代码注释多处写 "Conventions follows Ipopt"）。

---

## 6.1 ⭐ 原对偶内点法的数学

### 6.1.1 障碍问题

原问题：

$$
\min_w\ \Phi(w)\quad\text{s.t.}\quad g(w)=0,\ \ h(w)\ge0
$$

引入**松弛变量** $s>0$，把不等式变成等式 + 正性：

$$
\min_{w,s}\ \Phi(w)\quad\text{s.t.}\quad g(w)=0,\ \ h(w)-s=0,\ \ s>0
$$

再把 $s>0$ 用**对数障碍**吸收进代价：

$$
\boxed{\;\min_{w,s}\ \Phi(w)-\mu\sum_{i=1}^{m}\ln s_i
\quad\text{s.t.}\quad g(w)=0,\ \ h(w)-s=0\;}
$$

$\mu>0$ 是**障碍参数**。当 $\mu\to0^+$，障碍问题的解收敛到原问题的解
（中心路径 central path）。

### 6.1.2 KKT 条件

拉格朗日函数（$\nu$ 对应 $g$，$\lambda$ 对应 $h-s$）：

$$
\mathcal{L}=\Phi(w)-\mu\sum_i\ln s_i+\nu^{\!\top}g(w)-\lambda^{\!\top}\big(h(w)-s\big)
$$

（符号约定：OCS2 用 $-\lambda^\top(h-s)$，使得 $\lambda\ge0$ 对应 $h\ge0$ 的乘子）

一阶条件：

$$
\begin{aligned}
\nabla_w\mathcal{L}&=\nabla\Phi+\nabla g\,\nu-\nabla h\,\lambda=0 &&\text{(对偶可行性)}\\
\nabla_s\mathcal{L}&=-\mu S^{-1}\mathbf{1}+\lambda=0
 \;\Longleftrightarrow\; s_i\lambda_i=\mu\ \forall i &&\text{(松弛互补)}\\
g(w)&=0,\quad h(w)-s=0 &&\text{(原可行性)}\\
s&>0,\quad\lambda>0
\end{aligned}
$$

其中 $S=\operatorname{diag}(s)$。

**松弛互补条件 $s_i\lambda_i=\mu$** 是内点法的核心：
$\mu\to0$ 时退化为标准互补条件 $s_i\lambda_i=0$，
但保持 $\mu>0$ 让迭代点严格在可行域内部，避免了组合爆炸（哪些约束激活）。

### 6.1.3 ⭐ 凝聚（Condensing）：消去 $s$ 和 $\lambda$

**这是 OCS2 实现的关键一步**，让 QP 子问题不含额外变量。

对 KKT 系统做牛顿法。以 $s$ 与 $\lambda$ 的方程为例，线性化：

$$
\begin{aligned}
\Lambda\,\Delta s + S\,\Delta\lambda &= \mu\mathbf{1} - S\Lambda\mathbf{1} &&\text{(互补条件)}\\
\nabla h^{\!\top}\Delta w - \Delta s &= s - h &&\text{(}h-s=0\text{ 的线性化)}
\end{aligned}
$$

其中 $\Lambda=\operatorname{diag}(\lambda)$。从第二式解出

$$
\Delta s = \nabla h^{\!\top}\Delta w + h - s
$$

代入第一式：

$$
S\,\Delta\lambda = \mu\mathbf{1}-S\Lambda\mathbf{1}-\Lambda\big(\nabla h^{\!\top}\Delta w + h - s\big)
$$

$$
\Delta\lambda = S^{-1}\Big(\mu\mathbf{1}-\Lambda h\Big) - \lambda - S^{-1}\Lambda\,\nabla h^{\!\top}\Delta w
$$

（用了 $S\Lambda\mathbf{1}=\Lambda s$，故 $-S\Lambda\mathbf{1}+\Lambda s=0$）

把 $\Delta\lambda$ 代入对偶可行性方程 $\nabla^2\mathcal{L}\Delta w + \dots - \nabla h\,\Delta\lambda = \dots$：

$$
\boxed{
\begin{aligned}
\text{梯度增量：}&\quad \nabla h\,\Big(\frac{\lambda\odot h-\mu\mathbf{1}}{s}\Big)\;-\;\nabla h\,\lambda\\[4pt]
\text{Hessian 增量：}&\quad \nabla h\,\operatorname{diag}\!\Big(\frac{\lambda}{s}\Big)\,\nabla h^{\!\top}
\end{aligned}}
$$

**这正是代码 `ocs2_ipm/src/IpmHelpers.cpp:37` 的 `condenseIneqConstraints`**：

```cpp
void condenseIneqConstraints(scalar_t barrierParam, const vector_t& slack, const vector_t& dual,
                             const VectorFunctionLinearApproximation& ineqConstraint,
                             ScalarFunctionQuadraticApproximation& lagrangian) {
  // 对偶可行性：−∇h·λ
  lagrangian.dfdx.noalias() -= ineqConstraint.dfdx.transpose() * dual;
  lagrangian.dfdu.noalias() -= ineqConstraint.dfdu.transpose() * dual;

  // 凝聚系数
  const vector_t condensingLinearCoeff    = (dual.array()*ineqConstraint.f.array() - barrierParam) / slack.array();
  //                                        ↑ (λ⊙h − μ)/s
  const vector_t condensingQuadraticCoeff = dual.cwiseQuotient(slack);   // λ/s

  // 凝聚：加到梯度与 Hessian 上
  lagrangian.dfdx.noalias() += ineqConstraint.dfdx.transpose() * condensingLinearCoeff;
  const matrix_t Cq_dfdx = condensingQuadraticCoeff.asDiagonal() * ineqConstraint.dfdx;
  lagrangian.dfdxx.noalias() += ineqConstraint.dfdx.transpose() * Cq_dfdx;

  lagrangian.dfdu.noalias()  += ineqConstraint.dfdu.transpose() * condensingLinearCoeff;
  const matrix_t Cq_dfdu = condensingQuadraticCoeff.asDiagonal() * ineqConstraint.dfdu;
  lagrangian.dfduu.noalias() += ineqConstraint.dfdu.transpose() * Cq_dfdu;
  lagrangian.dfdux.noalias() += ineqConstraint.dfdu.transpose() * Cq_dfdx;
}
```

完全对应 ✓

**几何意义**：$\operatorname{diag}(\lambda/s)$ 是一个**自适应权重**：
- 约束接近激活（$s\to0$）→ 权重 $\to\infty$ → Hessian 在该方向变陡 → 步长自动变小
- 约束远离激活（$s$ 大）→ 权重 $\to0$ → 该约束几乎不影响子问题

这就是内点法"自动识别激活集"的机制。

**凝聚后 Hessian 保持 PSD**：因为 $\lambda>0$、$s>0$，
$\nabla h\operatorname{diag}(\lambda/s)\nabla h^\top\succeq0$ 是半正定的秩-$m$ 更新。
所以 HPIPM 收到的仍是一个凸 QP ✓

### 6.1.4 恢复 $\Delta s$ 与 $\Delta\lambda$

QP 解出 $\Delta w=(\delta x,\delta u)$ 后，回代：

`retrieveSlackDirection`（`IpmHelpers.cpp:70,83`）：

$$
\Delta s = h + \nabla h^{\!\top}\Delta w - s
$$

```cpp
vector_t slackDirection = stateInputIneqConstraints.f - slackStateInputIneq;   // h − s
slackDirection.noalias() += stateInputIneqConstraints.dfdx * dx;
slackDirection.noalias() += stateInputIneqConstraints.dfdu * du;
```

`retrieveDualDirection`（`IpmHelpers.cpp:95`）：

$$
\Delta\lambda = -\frac{\lambda\odot(s+\Delta s)-\mu\mathbf{1}}{s}
$$

```cpp
vector_t dualDirection = dual.cwiseProduct(slack + slackDirection);
dualDirection.array() -= barrierParam;
dualDirection.array() /= -slack.array();
```

**验证**：把 $\Delta s$ 代入 §6.1.3 的 $\Delta\lambda$ 表达式：

$$
\Delta\lambda = S^{-1}(\mu-\Lambda h)-\lambda-S^{-1}\Lambda\nabla h^{\!\top}\Delta w
= \frac{\mu - \lambda\odot h - \lambda\odot s - \lambda\odot\nabla h^\top\Delta w}{s}
$$

而 $s+\Delta s = h + \nabla h^\top\Delta w$，所以

$$
-\frac{\lambda\odot(s+\Delta s)-\mu}{s} = \frac{\mu-\lambda\odot h-\lambda\odot\nabla h^\top\Delta w}{s}
$$

两式相差 $-\lambda\odot s/s = -\lambda$。代码里这 $-\lambda$ 已经在
`condenseIneqConstraints` 的第一步（`dfdx -= ∇h·λ`）里体现，
所以是**一致的**——代码用的是 $\lambda^{new}=\lambda+\Delta\lambda$ 的直接表达式 ✓

### 6.1.5 ⭐ 分数到边界规则（Fraction-to-Boundary）

必须保证 $s+\alpha\Delta s>0$ 与 $\lambda+\alpha\Delta\lambda>0$。
标准 IPM 用**分数到边界规则**：

$$
\alpha^{\max}=\max\{\alpha\in(0,1]\ :\ v+\alpha\,\Delta v\ \ge\ (1-\tau)\,v\}
$$

$\tau\in(0,1)$ 是"允许消耗掉多少距离边界的裕度"（IPOPT 的 `tau_min`，OCS2 默认 0.995）。

**推导**：条件 $v_i+\alpha\Delta v_i\ge(1-\tau)v_i$ 即 $\alpha\Delta v_i\ge-\tau v_i$。
- 若 $\Delta v_i\ge0$：恒成立
- 若 $\Delta v_i<0$：$\alpha\le\dfrac{\tau v_i}{-\Delta v_i}=\dfrac{-\tau v_i}{\Delta v_i}$

故

$$
\alpha^{\max}=\min\left(1,\ \min_{i:\Delta v_i<0}\frac{-\tau v_i}{\Delta v_i}\right)
$$

**代码用了一个等价但更简洁的写法**（`IpmHelpers.cpp:103`）：

```cpp
scalar_t fractionToBoundaryStepSize(const vector_t& v, const vector_t& dv, scalar_t marginRate) {
  const vector_t invFractionToBoundary = (-1.0 / marginRate) * dv.cwiseQuotient(v);
  const auto alpha = invFractionToBoundary.maxCoeff();
  return alpha > 0.0 ? std::min(1.0/alpha, 1.0) : 1.0;
}
```

即先算 $\beta_i = \dfrac{-\Delta v_i}{\tau v_i}$，取 $\beta_{\max}=\max_i\beta_i$。
若 $\beta_{\max}>0$（存在 $\Delta v_i<0$，因 $v_i>0$），
则 $\alpha=\min(1/\beta_{\max},1)=\min\Big(\dfrac{\tau v_i}{-\Delta v_i},1\Big)$ ✓
若 $\beta_{\max}\le0$（所有 $\Delta v_i\ge0$），返回 1 ✓

**原对偶步长分离**：IPM 通常对 $w,s$ 用 $\alpha^{primal}$、对 $\lambda$ 用 $\alpha^{dual}$。
OCS2 的 `usePrimalStepSizeForDual`（默认 true）让两者取较小值，更保守但更稳。

---

## 6.2 初始化（严格照搬 IPOPT）

**文件**：`ocs2_ipm/include/ocs2_ipm/IpmInitialization.h`

### 松弛变量

$$
s_0 = (1+\kappa_s)\cdot\max\big(h(w_0),\ s_{\min}\big)
$$

```cpp
inline vector_t initializeSlackVariable(const vector_t& ineqConstraint,
                                        scalar_t initialSlackLowerBound, scalar_t initialSlackMarginRate) {
  return (1.0 + initialSlackMarginRate) * ineqConstraint.cwiseMax(initialSlackLowerBound);
}
```

- `initialSlackLowerBound` = IPOPT 的 `slack_bound_push`（默认 1e-4）
  ——保证即使 $h(w_0)<0$（初值不可行）也有 $s>0$
- `initialSlackMarginRate` = IPOPT 的 `slack_bound_frac`（默认 0.01）
  ——额外裕度，让初始点不贴着边界

### 对偶变量

$$
\lambda_0 = (1+\kappa_\lambda)\cdot\max\Big(\frac{\mu_0}{s_0},\ \lambda_{\min}\Big)
$$

```cpp
inline vector_t initializeDualVariable(const vector_t& slack, scalar_t barrierParam,
                                       scalar_t initialDualLowerBound, scalar_t initialDualMarginRate) {
  return (1.0 + initialDualMarginRate) * (barrierParam * slack.cwiseInverse()).cwiseMax(initialDualLowerBound);
}
```

对应 IPOPT 的 `bound_mult_init_method: mu-based`：
从互补条件 $s\lambda=\mu$ 反推 $\lambda=\mu/s$，
使初始点**恰好在中心路径上**。

---

## 6.3 障碍参数调度

**文件**：`ocs2_ipm/src/IpmSolver.cpp:898`

```cpp
scalar_t IpmSolver::updateBarrierParameter(scalar_t mu, const PerformanceIndex& baseline,
                                           const ipm::StepInfo& stepInfo) const {
  if (mu <= settings_.targetBarrierParameter) {
    return mu;                                    // 已到目标，不再降
  } else if (|Δmerit| < barrierReductionCostTol &&
             totalConstraintViolation < barrierReductionConstraintTol) {
    return std::min(mu * barrierLinearDecreaseFactor,               // 线性：μ ← κ·μ
                    std::pow(mu, barrierSuperlinearDecreasePower)); // 超线性：μ ← μ^θ
  } else {
    return mu;                                    // 当前 μ 下还没解好，保持
  }
}
```

**这是 IPOPT 的 "monotone" 更新策略**（Wächter & Biegler 2006 §3.2）：

$$
\mu_{k+1}=\min\Big(\kappa_\mu\,\mu_k,\ \mu_k^{\theta_\mu}\Big),
\qquad \kappa_\mu\in(0,1),\ \theta_\mu\in(1,2)
$$

**为什么取两者较小值？**

- **线性项** $\kappa_\mu\mu$：$\mu$ 大时主导（例如 $\mu=10^{-2}$，$\kappa_\mu=0.2$ → $2\times10^{-3}$；
  而 $\mu^{1.5}=10^{-3}$，更小，故取超线性）
- **超线性项** $\mu^{\theta_\mu}$：$\mu$ 小时主导（$\mu=10^{-6}$ → $\mu^{1.5}=10^{-9}$，
  而 $\kappa_\mu\mu=2\times10^{-7}$，更小，故取线性）

实际上：$\mu$ 大时超线性下降更快，$\mu$ 小时线性下降更快。
取 min 保证**始终按更快的那个走**，同时避免 $\mu$ 骤降导致的数值困难。

**关键约束**：只在"当前 $\mu$ 下的子问题已解得够好"时才降 $\mu$。
否则会跳过中心路径导致震荡。

默认值：$\mu_0=10^{-2}$、$\mu_{target}=10^{-4}$、$\kappa_\mu=0.2$、$\theta_\mu=1.5$。

> 注意 $\mu_{target}=10^{-4}$ **不算小**（IPOPT 默认降到 $10^{-11}$）。
> 这是 MPC 场景的有意折中：不追求高精度收敛，只要约束"够好"即可，换取实时性。

---

## 6.4 性能指标中的障碍项

**文件**：`ocs2_ipm/src/IpmPerformanceIndexComputation.cpp:41`

```cpp
PerformanceIndex computePerformanceIndex(const multiple_shooting::Transcription& transcription, scalar_t dt,
                                         scalar_t barrierParam,
                                         const vector_t& slackStateIneq, const vector_t& slackStateInputIneq) {
  auto performance = multiple_shooting::computePerformanceIndex(transcription, dt);

  // 障碍项计入 cost
  if (slackStateIneq.size() > 0)
    performance.cost -= dt * barrierParam * slackStateIneq.array().log().sum();
  if (slackStateInputIneq.size() > 0)
    performance.cost -= dt * barrierParam * slackStateInputIneq.array().log().sum();

  // h − s = 0 计入等式约束违反
  if (transcription.stateIneqConstraints.f.size() > 0)
    performance.equalityConstraintsSSE += dt * (transcription.stateIneqConstraints.f - slackStateIneq).squaredNorm();
  if (transcription.stateInputIneqConstraints.f.size() > 0)
    performance.equalityConstraintsSSE += dt * (transcription.stateInputIneqConstraints.f - slackStateInputIneq).squaredNorm();

  return performance;
}
```

**两个设计要点**：

1. **障碍项 $-\mu\sum\ln s_i$ 计入 `cost`**
   → 滤子线搜索的"代价"目标自动变成障碍目标 ✓
2. **$\|h-s\|^2$ 计入 `equalityConstraintsSSE`**
   → 滤子线搜索的"可行性"目标自动包含松弛一致性 ✓

于是**同一套 `FilterLinesearch` 代码可以直接复用**，不需要为 IPM 单独写线搜索。
这是很干净的复用设计。

`dualFeasibilitiesSSE` 则累加互补条件残差
（`evaluateComplementarySlackness`，即 $\|s\odot\lambda-\mu\mathbf{1}\|^2$）。

---

## 6.5 主循环

**文件**：`ocs2_ipm/src/IpmSolver.cpp:237` 起

```cpp
scalar_t barrierParam = settings_.initialBarrierParameter;
initializeSlackDualTrajectory(timeDiscretization, x, u, barrierParam,
                              slackStateIneq, dualStateIneq, slackStateInputIneq, dualStateInputIneq);

while (convergence == ipm::Convergence::FALSE) {
  // ① 构造 QP（含凝聚）
  const auto baselinePerformance = setupQuadraticSubproblem(
      timeDiscretization, initState, x, u, lmd, nu, barrierParam,
      slackStateIneq, dualStateIneq, slackStateInputIneq, dualStateInputIneq, metrics);

  // ② 解 QP（无不等式约束！），并恢复 Δs、Δλ、算最大步长
  const vector_t delta_x0 = initState - x[0];
  const auto deltaSolution = getOCPSolution(delta_x0, barrierParam,
      slackStateIneq, dualStateIneq, slackStateInputIneq, dualStateInputIneq);

  // ③ 滤子线搜索（原变量）
  const auto stepInfo = takePrimalStep(baselinePerformance, timeDiscretization, initState,
                                       deltaSolution, x, u, barrierParam, slackStateIneq, slackStateInputIneq, metrics);
  // ④ 更新对偶变量
  takeDualStep(deltaSolution, stepInfo, lmd, nu, dualStateIneq, dualStateInputIneq);

  // ⑤ 收敛检查
  convergence = checkConvergence(iter, barrierParam, baselinePerformance, stepInfo);

  // ⑥ 更新障碍参数
  barrierParam = updateBarrierParameter(barrierParam, baselinePerformance, stepInfo);
  ++iter;
}
```

### 步长计算（`getOCPSolution`，`:470`）

```cpp
scalar_array_t primalStepSizes(settings_.nThreads, 1.0);
scalar_array_t dualStepSizes(settings_.nThreads, 1.0);
// 并行地对每个节点：
deltaSlackStateIneq[i] = ipm::retrieveSlackDirection(stateIneqConstraints_[i], deltaXSol[i], barrierParam, slackStateIneq[i]);
deltaDualStateIneq[i]  = ipm::retrieveDualDirection(barrierParam, slackStateIneq[i], dualStateIneq[i], deltaSlackStateIneq[i]);
primalStepSizes[workerId] = std::min({primalStepSizes[workerId],
    ipm::fractionToBoundaryStepSize(slackStateIneq[i], deltaSlackStateIneq[i], settings_.fractionToBoundaryMargin), ...});
dualStepSizes[workerId] = std::min({dualStepSizes[workerId],
    ipm::fractionToBoundaryStepSize(dualStateIneq[i],  deltaDualStateIneq[i],  settings_.fractionToBoundaryMargin), ...});
// 全局归约
solution.maxPrimalStepSize = *std::min_element(primalStepSizes.begin(), primalStepSizes.end());
solution.maxDualStepSize   = *std::min_element(dualStepSizes.begin(),  dualStepSizes.end());
```

`maxPrimalStepSize` 是滤子线搜索的**起点**（而非 SQP 里的 1.0）——
保证第一次试探就已满足内点条件。

### 对偶步（`takeDualStep`，`:886`）

```cpp
if (settings_.computeLagrangeMultipliers) {
  multiple_shooting::incrementTrajectory(lmd, subproblemSolution.deltaLmdSol, stepInfo.primalStepSize, lmd);
  multiple_shooting::incrementTrajectory(nu,  subproblemSolution.deltaNuSol,  stepInfo.primalStepSize, nu);
}
const scalar_t dualStepSize = settings_.usePrimalStepSizeForDual
    ? std::min(stepInfo.primalStepSize, subproblemSolution.maxDualStepSize)
    : subproblemSolution.maxDualStepSize;
multiple_shooting::incrementTrajectory(dualStateIneq, subproblemSolution.deltaDualStateIneq, dualStepSize, dualStateIneq);
multiple_shooting::incrementTrajectory(dualStateInputIneq, subproblemSolution.deltaDualStateInputIneq, dualStepSize, dualStateInputIneq);
```

`lmd` 是**动力学约束的协态**（costate），`nu` 是状态-输入等式约束的乘子。
`computeLagrangeMultipliers` 为 false 时不更新它们——
此时 `dualFeasibilitiesSSE` 的记录不准，但**不影响算法正确性**
（设置注释里明确说了这一点），因为它们不参与原变量的更新。

---

## 6.6 收敛判据

**文件**：`IpmSolver.cpp:911`

```cpp
if ((iteration+1) >= settings_.ipmIteration)              return Convergence::ITERATIONS;
if (stepInfo.primalStepSize < settings_.alpha_min)        return Convergence::STEPSIZE;
if (|Δmerit| < costTol && θ < g_min)                      return Convergence::METRICS;
if (dx_norm < deltaTol && du_norm < deltaTol
    && barrierParam <= settings_.targetBarrierParameter)  return Convergence::PRIMAL;
return Convergence::FALSE;
```

**注意 `PRIMAL` 分支多了 `barrierParam <= targetBarrierParameter`**：
即使原变量不再变化，只要 $\mu$ 还没降到目标值就不算收敛——
因为降低 $\mu$ 会重新激活变化。这是 IPM 与 SQP 收敛判据的唯一区别。

---

## 6.7 `IpmSettings` 完全参考

| 参数 | 默认 | 说明 |
|---|---|---|
| `ipmIteration` | 10 | 最大迭代 |
| `deltaTol` / `costTol` | 1e-6 / 1e-4 | 收敛阈值 |
| `alpha_decay` / `alpha_min` | 0.5 / 1e-4 | 线搜索 |
| `g_max` / `g_min` / `armijoFactor` / `gamma_c` | 同 SQP | 滤子参数 |
| **`initialBarrierParameter`** | 1e-2 | $\mu_0$ |
| **`targetBarrierParameter`** | 1e-4 | $\mu_{target}$ |
| `barrierReductionCostTol` | 1e-2 | 降 $\mu$ 的代价变化阈值 |
| `barrierReductionConstraintTol` | 1e-2 | 降 $\mu$ 的约束违反阈值 |
| **`barrierLinearDecreaseFactor`** | 0.2 | $\kappa_\mu$ |
| **`barrierSuperlinearDecreasePower`** | 1.5 | $\theta_\mu$ |
| `initialSlackLowerBound` | 1e-4 | IPOPT `slack_bound_push` |
| `initialDualLowerBound` | 1e-4 | $\lambda$ 下界 |
| `initialSlackMarginRate` | 0.01 | IPOPT `slack_bound_frac` |
| `initialDualMarginRate` | 0.01 | $\lambda$ 裕度 |
| **`fractionToBoundaryMargin`** | 0.995 | IPOPT `tau_min` |
| `usePrimalStepSizeForDual` | true | 原对偶步长取较小者 |
| `computeLagrangeMultipliers` | false | 是否更新 $\nu,\lambda$（只影响日志） |
| `useFeedbackPolicy` / `createValueFunction` | true / false | 同 SQP |
| `dt` / `integratorType` | 0.01 / RK2 | 离散化 |
| `nThreads` / `threadPriority` | 4 / 50 | 并发 |
| `hpipmSettings` | — | QP 求解器设置 |

---

## 6.8 IPM vs SQP：如何选

| 场景 | 推荐 |
|---|---|
| 不等式约束少（<10 个/节点） | **SQP**——HPIPM 内部处理效率高 |
| 不等式约束多且大部分不激活 | **IPM**——凝聚后 QP 无不等式，且激活集信息跨迭代复用 |
| 需要严格内点（约束绝不可违反） | **IPM**——分数到边界规则保证 $h>0$ 始终成立 |
| 允许中间迭代不可行 | **SQP** |
| 参数调试成本敏感 | **SQP**——IPM 多 10 个障碍相关参数 |

**IPM 的一个实用优势**：由于始终 $h(w)>0$，
即使 MPC 中途被打断（超时），取当前迭代点也是**严格可行**的。
SQP 的中间迭代可能违反约束。这在安全关键场景很重要。

---

## 6.9 参考文献

- **Wächter & Biegler (2006)**, *On the Implementation of an Interior-Point Filter Line-Search Algorithm for Large-Scale Nonlinear Programming*, Mathematical Programming 106(1):25–57 — **IPOPT 论文，本模块的直接蓝本**
- **IPOPT 文档**：https://coin-or.github.io/Ipopt/OPTIONS.html （代码注释多处直接引用）
- **Nocedal & Wright (2006)**, *Numerical Optimization*, Ch. 19（Interior-Point Methods for Nonlinear Programming）
- **Boyd & Vandenberghe (2004)**, *Convex Optimization*, Ch. 11（Interior-point methods）——障碍法与中心路径的基础
- **Katayama & Ohtsuka**, *Efficient solution method based on inverse dynamics for interior point methods of nonlinear MPC*, ICRA 2022 —— OCS2 的 IPM 实现由 Sotaro Katayama 贡献，与其 robotoc 工作同源

**下一章**：一阶方法 PIPG → [07 ocs2_slp](07-ocs2-slp.md)
