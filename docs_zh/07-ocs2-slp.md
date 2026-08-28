# 07 · ocs2_slp：PIPG 一阶求解器

`ocs2_slp` 是 OCS2 里最"轻"的求解器：**不依赖任何外部 QP 库**，
用一阶原对偶方法（PIPG）求解多重打靶子问题。

**动机**：HPIPM/BLASFEO 是高度优化的 C 库，但

- 编译依赖重，难以移植到某些嵌入式平台
- 内点法每步需要因式分解，$O(Nn^3)$
- 无法在 GPU 上高效并行

PIPG 每步只做**矩阵-向量乘法**，$O(Nn^2)$，且完全并行友好。
代价是收敛速度慢（一阶方法），需要更多迭代。

**SLP = Sequential Linear Programming？**
名字有点误导——它其实是 SQP 的同构（同样的二次子问题），
只是子问题用 PIPG 而非 HPIPM 解。可以理解为 "Sequential LP-like Programming"，
或直接看成 "SQP with a first-order QP solver"。

---

## 7.1 ⭐ PIPG 算法

### 7.1.1 问题形式

**文件**：`ocs2_slp/include/ocs2_slp/pipg/PipgBounds.h` 的类注释

决策向量 $z=[u_0; x_1; u_1; \dots; u_n; x_{n+1}]$，问题：

$$
\min_z\ \tfrac12 z^{\!\top}Hz + h^{\!\top}z \qquad\text{s.t.}\qquad Gz=g
$$

（结构见 [03 章](03-ocs2-oc.md) §3.6）

### 7.1.2 鞍点形式

引入对偶变量 $w$，拉格朗日函数：

$$
\mathcal{L}(z,w)=\tfrac12 z^{\!\top}Hz+h^{\!\top}z+w^{\!\top}(Gz-g)
$$

最优解是鞍点：$\min_z\max_w\mathcal{L}$。

### 7.1.3 PIPG 迭代

**PIPG**（Proportional-Integral Projected Gradient，Yu, Elango & Açıkmeşe 2020）
是一个带"比例-积分"结构的原对偶方法：

$$
\boxed{
\begin{aligned}
v^{k} &= w^{k} + \beta\big(Gz^{k}-g\big) &&\text{(比例项：即时约束残差)}\\
z^{k+1} &= z^{k}-\alpha\big(Hz^{k}+h+G^{\!\top}v^{k}\big) &&\text{(原变量梯度下降)}\\
w^{k+1} &= w^{k}+\beta\big(Gz^{k+1}-g\big) &&\text{(积分项：对偶上升)}
\end{aligned}}
$$

**"比例-积分"的含义**：
- $w$ 累积历史约束残差 → **积分（I）项**
- $v = w+\beta(Gz-g)$ 在当前残差上加了即时项 → **比例（P）项**

这个结构来自控制论的 PI 控制器类比：$w$ 像积分器消除稳态误差，
$\beta(Gz-g)$ 像比例项提供快速响应。相比纯对偶上升（只有 I），
PI 结构显著改善了收敛速度。

### 7.1.4 步长

**文件**：`PipgBounds.h:66-68`

设 $\mu I\preceq H\preceq\lambda I$、$G^{\!\top}G\preceq\sigma I$，则

$$
\boxed{\;
\alpha_k=\frac{2}{(k+1)\mu+2\lambda},\qquad
\beta_k=\frac{(k+1)\mu}{2\sigma}\;}
$$

```cpp
scalar_t primalStepSize(size_t iteration) const {
  return 2.0 / ((static_cast<scalar_t>(iteration) + 1.0) * mu + 2.0 * lambda);
}
scalar_t dualStepSize(size_t iteration) const {
  return (static_cast<scalar_t>(iteration) + 1.0) * mu / (2.0 * sigma);
}
```

**性质**：
- $\alpha_k\to0$ 如 $O(1/k)$，$\beta_k\to\infty$ 如 $O(k)$
- $\alpha_k\beta_k\to \dfrac{1}{\sigma}$（常数）——保证原对偶更新的"平衡"
- $k=0$ 时 $\alpha_0=\dfrac{2}{\mu+2\lambda}\approx\dfrac{1}{\lambda}$，即标准梯度下降步长

这组步长来自 PIPG 论文的收敛性证明，保证 $O(1/k)$ 的遍历收敛率
（对强凸问题可达线性收敛）。

### 7.1.5 逐阶段展开的实现

**文件**：`ocs2_slp/src/pipg/PipgSolver.cpp:49`

代码**不组装全局 $H$、$G$**，而是把上面的迭代按时间阶段展开。
对阶段 $t$（`:152-187`）：

```cpp
// v_{t-1} = w_{t-1} + (β + β_last)(C·x_t − A·x_{t-1} − B·u_{t-1} − b)
V_[t-1] = W_[t-1] - (beta + betaLast) * b;
V_[t-1].array() += (beta + betaLast) * C.array() * X_[t].array();
V_[t-1].noalias() -= (beta + betaLast) * (A * X_[t-1]);
V_[t-1].noalias() -= (beta + betaLast) * (B * U_[t-1]);

// u_{t-1} ← u_{t-1} − α(R·u_{t-1} + P·x_{t-1} + r − Bᵀ·v_{t-1})
UNew_[t-1] = U_[t-1] - alpha * r;
UNew_[t-1].noalias() -= alpha * (R * U_[t-1]);
UNew_[t-1].noalias() -= alpha * (P * X_[t-1]);
UNew_[t-1].noalias() += alpha * (B.transpose() * V_[t-1]);

// x_t ← x_t − α(Q·x_t + q + C·v_{t-1} − Aᵀ_{next}·v_t + Pᵀ_{next}·u_t)
XNew_[t] = X_[t] - alpha * q;
XNew_[t].array() -= alpha * C.array() * V_[t-1].array();
XNew_[t].noalias() -= alpha * (Q * X_[t]);
if (t != N) {
  // 计算 v_t（下一阶段的对偶）并加入 x_t 的梯度
  vector_t VNext = W_[t] - (beta+betaLast)*bNext;
  VNext.array()  += (beta+betaLast) * CNext.array() * X_[t+1].array();
  VNext.noalias() -= (beta+betaLast) * (ANext * X_[t]);
  VNext.noalias() -= (beta+betaLast) * (BNext * U_[t]);
  XNew_[t].noalias() += alpha * (ANext.transpose() * VNext);
  XNew_[t].noalias() -= alpha * (PNext.transpose() * U_[t]);
}
```

**为什么 $x_t$ 的梯度里有 $A_{next}^\top v_t$？**
因为 $x_t$ 同时出现在两个动力学约束里：
- 作为 $t-1\to t$ 的**目标**：$C x_t - Ax_{t-1}-Bu_{t-1}-b=0$ → 贡献 $+C^\top v_{t-1}$
- 作为 $t\to t+1$ 的**源**：$Cx_{t+1}-A_{next}x_t-B_{next}u_t-b_{next}=0$ → 贡献 $-A_{next}^\top v_t$

$G^\top v$ 的这两块正是 $G$ 的双对角结构的转置。

**`scalingVectors` $C$ 是什么？**
原本动力学约束里 $x_{t}$ 的系数是 $-I$（见 §3.6 的 $G$ 矩阵）。
Ruiz 预条件（§7.4）会把它缩放成任意对角矩阵，
代码用向量 $C$ 存这个对角线，故 `C.array() * X_[t].array()` 是对角矩阵乘向量。

**`beta + betaLast`**：代码把 $v^k$ 与 $w^{k}$ 的更新合并，
用上一次的 $\beta_{k-1}$ 与本次 $\beta_k$ 之和，是一个等价的重排以减少一次遍历。

### 7.1.6 并行与同步

`updateVariablesTask`（`:99`）用 `std::atomic_int timeIndex` 分发阶段任务，
配合 `std::condition_variable` 做迭代间同步：

```cpp
while (keepRunning) {
  while ((t = timeIndex++) <= N) {  /* 更新阶段 t */ }
  // 等所有线程完成本轮，主线程检查收敛后放行下一轮
}
```

**这里有一个精妙的数据竞争规避**（`:128` 注释）：
> "Move the update of W to the front of the calculation of V to prevent data race."

因为 $v_{t}$ 需要 $w_{t}$ 的**旧值**，而 $w_t$ 的更新又需要 $x_t$ 的新值，
所以必须先用旧 $w$ 算完 $v$，再更新 $w$。

**内存复用技巧**（`:142-149`）：
```cpp
// UNew/XNew 此时存的是"上上次"的解，先用它们算差分再覆盖
UNew_[t-1] -= U_[t-1];
XNew_[t]   -= X_[t];
solutionSEArray[t-1] = UNew_[t-1].squaredNorm() + XNew_[t].squaredNorm();
```
避免额外分配"上次解"的缓冲区。

### 7.1.7 收敛判据

两个：
1. **原可行性**：$\|E^{-1}(Cx_t-Ax_{t-1}-Bu_{t-1}-b)\|_\infty <$ `absoluteTolerance`
   （$E^{-1}$ 是预条件的逆，把残差还原到原始尺度）
2. **解的变化**：$\sqrt{\text{SSE}(\Delta z)}/\|z\| <$ `relativeTolerance`

---

## 7.2 步长界的估计

PIPG 需要 $\mu$、$\lambda$、$\sigma$ 三个界。精确算它们需要特征值分解（太贵），
OCS2 用**Gershgorin 圆盘定理**给出上界。

### 7.2.1 $\lambda$：$H$ 的最大特征值上界

**文件**：`ocs2_slp/src/Helpers.cpp:53`

$$
\lambda_{\max}(H)\ \le\ \max_i\sum_j|H_{ij}|
$$

这是 Gershgorin 定理的直接推论（也是 $\|H\|_\infty$）。

```cpp
scalar_t hessianEigenvaluesUpperBound(const OcpSize& ocpSize,
                                      const std::vector<ScalarFunctionQuadraticApproximation>& cost) {
  const vector_t rowwiseAbsSumH = hessianAbsRowSum(ocpSize, cost);
  return rowwiseAbsSumH.maxCoeff();
}
```

`hessianAbsRowSum`（`:66`）逐阶段累加绝对值行和。由于 $H$ 是块对角的
（每个阶段一个 $\begin{bmatrix}Q&P^\top\\P&R\end{bmatrix}$ 块），
只需在块内求和：

```cpp
res.head(nu_0) = cost[0].dfduu.cwiseAbs().rowwise().sum();
for (k = 1; k < N; k++) {
  res.segment(curRow, nx_k)         = cost[k].dfdxx.cwiseAbs().rowwise().sum();
  res.segment(curRow, nx_k)        += cost[k].dfdux.transpose().cwiseAbs().rowwise().sum();  // Pᵀ 块
  res.segment(curRow+nx_k, nu_k)    = cost[k].dfdux.cwiseAbs().rowwise().sum();               // P 块
  res.segment(curRow+nx_k, nu_k)   += cost[k].dfduu.cwiseAbs().rowwise().sum();
  curRow += nx_k + nu_k;
}
res.tail(nx_N) = cost[N].dfdxx.cwiseAbs().rowwise().sum();
```

### 7.2.2 $\sigma$：$G^\top G$ 的最大特征值上界

**文件**：`Helpers.cpp:58`（入口）与 `:95`（核心）

关键观察（代码注释）：
> "since the G'G and GG' have exactly the same set of eigenvalues value: G G' < sigma I"

$G^\top G$ 与 $GG^\top$ 的**非零特征值完全相同**（这是 SVD 的直接推论）。
而 $GG^\top$ 的维度是约束数（$\approx Nn_x$），$G^\top G$ 的维度是变量数（$\approx N(n_x+n_u)$），
**前者更小**，且 $GG^\top$ 是**块三对角**的，可以逐块计算。

对第 $k$ 块行 $[\,\cdots\ -A_k\ -B_k\ C_k\ \cdots]$，$GG^\top$ 的对角块是

$$
(GG^{\!\top})_{kk}= C_k^2 + B_kB_k^{\!\top} + A_kA_k^{\!\top}
$$

（$C_k$ 是对角缩放，$C_k^2$ 即逐元素平方）

非对角块 $(GG^\top)_{k,k-1} = -A_k C_{k-1}$、$(GG^\top)_{k,k+1}=-C_kA_{k+1}^\top$。

```cpp
tempMatrixArray[k] = (scalingVectorsPtr == nullptr
                      ? matrix_t::Identity(nx_next, nx_next)
                      : (*scalingVectorsPtr)[k].cwiseProduct((*scalingVectorsPtr)[k]).asDiagonal().toDenseMatrix());
tempMatrixArray[k] += B * B.transpose();
if (k != 0) tempMatrixArray[k] += A * A.transpose();
absRowSumArray[k] = tempMatrixArray[k].cwiseAbs().rowwise().sum();      // 对角块的行和
if (k != 0)     absRowSumArray[k] += (A * C_{k-1}).cwiseAbs().rowwise().sum();      // 下对角块
if (k != N-1)   absRowSumArray[k] += (A_{next} * C_k).cwiseAbs().rowwise().sum();   // 上对角块
```

最后取全局最大值。**并行**：每个 $k$ 独立，用 `std::atomic_int` 分发。

### 7.2.3 $\mu$：$H$ 的最小特征值下界

**这个最难**——Gershgorin 给出的下界通常是负的（无用）。
OCS2 采取了**实用主义方案**（`SlpSolver.cpp:258`）：

```cpp
const auto muEstimated = [&]() {
  scalar_t maxScalingFactor = -1;
  for (auto& v : D) if (v.size() != 0) maxScalingFactor = std::max(maxScalingFactor, v.maxCoeff());
  return c * pipgSolver_.settings().lowerBoundH * maxScalingFactor * maxScalingFactor;
}();
```

即：$\mu = c\cdot\underline{h}\cdot(\max D)^2$，其中

- $\underline{h}$ = `lowerBoundH`，**用户指定的常数**（`PipgSettings`）
- $c$、$D$ 来自 Ruiz 预条件（缩放后 $H\to cDHD$，故下界也按 $c\,d^2$ 缩放）

**这是一个用户必须调的参数**。$\mu$ 取太大 → 步长过大 → 发散；
取太小 → 步长过小 → 收敛慢。

> **实践建议**：$\underline{h}$ 应设为代价 Hessian 中最小的正定块的特征值量级。
> 如果代价里有 $\frac{\epsilon}{2}\|u\|^2$ 这样的正则项，$\underline{h}\approx\epsilon$ 是合理起点。
> 若 $H$ 真的奇异（如某些状态不出现在代价里），PIPG 的收敛保证失效，
> 必须先加正则项。

---

## 7.3 `SlpSolver` 主循环

**文件**：`ocs2_slp/src/SlpSolver.cpp:159`

与 SQP 的结构**几乎逐行相同**（对照 [05 章](05-ocs2-sqp.md) §5.2）：

```cpp
while (convergence == slp::Convergence::FALSE) {
  const auto baselinePerformance = setupQuadraticSubproblem(timeDiscretization, initState, x, u, metrics);
  const vector_t delta_x0 = initState - x[0];
  const auto deltaSolution = getOCPSolution(delta_x0);      // ← 唯一的区别在这里
  const auto stepInfo = takeStep(baselinePerformance, timeDiscretization, initState, deltaSolution, x, u, metrics);
  convergence = checkConvergence(iter, baselinePerformance, stepInfo);
  ++iter;
}
```

`getOCPSolution`（`:239`）的五步：

```cpp
// ① 预条件
precondition::ocpDataInPlaceInParallel(threadPool_, delta_x0, pipgSolver_.size(),
                                       settings_.scalingIteration, dynamics_, cost_,
                                       D, E, scalingVectors, c);
// ② 估计 μ
const auto muEstimated = c * lowerBoundH * maxD * maxD;
// ③ 估计 λ
const auto lambdaScaled = slp::hessianEigenvaluesUpperBound(pipgSolver_.size(), cost_);
// ④ 估计 σ
const auto sigmaScaled = slp::GGTEigenvaluesUpperBound(threadPool_, pipgSolver_.size(),
                                                       dynamics_, nullptr, &scalingVectors);
// ⑤ 跑 PIPG
const pipg::PipgBounds pipgBounds{muEstimated, lambdaScaled, sigmaScaled};
pipgSolver_.solve(threadPool_, delta_x0, dynamics_, cost_, nullptr, scalingVectors, &EInv,
                  pipgBounds, deltaXSol, deltaUSol);
// ⑥ 反缩放 + 还原投影
precondition::descaleSolution(D, deltaXSol, deltaUSol);
multiple_shooting::remapProjectedInput(constraintsProjection_, deltaXSol, deltaUSol);
```

**注意 SLP 强制使用投影**（不像 SQP 可选）——
因为 PIPG 只处理等式约束 $Gz=g$（动力学），
状态-输入等式约束必须先被消去。

**Benchmark 计时器**：`preConditioning_`、`lambdaEstimation_`、`sigmaEstimation_`、`pipgSolverTimer_`
四个独立计时器，方便定位瓶颈。实测中预条件与界估计的开销可占总时间的 20~30%。

---

## 7.4 ⭐ Ruiz 均衡预条件

**文件**：`ocs2_oc/include/ocs2_oc/precondition/Ruzi.h`（文件名拼写为 "Ruzi"）
与 `ocs2_oc/src/precondition/Ruzi.cpp`

### 7.4.1 为什么需要预条件

一阶方法的收敛速度由**条件数** $\kappa(H)=\lambda_{\max}/\lambda_{\min}$ 决定。
机器人问题里代价矩阵常常量纲混杂（位置 $\sim$m、力 $\sim$100N、角度 $\sim$rad），
$\kappa$ 轻易达到 $10^6$，PIPG 会慢到不可用。

**Ruiz 均衡**（Ruiz 2001）通过对角缩放让矩阵各行各列的范数接近 1，
把 $\kappa$ 降到 $10^1\sim10^2$。

### 7.4.2 变换

对 $y := D^{-1}z$（即 $z=Dy$），问题变为

$$
\min_y\ \frac{c}{2}\,y^{\!\top}(DHD)\,y + c\,y^{\!\top}(Dh)
\qquad\text{s.t.}\qquad (EGD)\,y = Eg
$$

其中 $D$、$E$ 是对角正定矩阵，$c>0$ 是代价缩放标量。

**为什么这样变换保持等价**：
$z^\top Hz = y^\top DHDy$、$Gz-g = 0 \Leftrightarrow E(GDy-g)=0$（$E$ 可逆），
所以最优解满足 $z^\star = Dy^\star$ ✓

### 7.4.3 迭代

Ruiz 算法每轮：

$$
\begin{aligned}
\delta_i &\leftarrow \frac{1}{\sqrt{\|(\text{矩阵第 }i\text{ 列})\|_\infty}},&\quad D&\leftarrow D\operatorname{diag}(\delta)\\
\epsilon_j &\leftarrow \frac{1}{\sqrt{\|(\text{矩阵第 }j\text{ 行})\|_\infty}},&\quad E&\leftarrow E\operatorname{diag}(\epsilon)
\end{aligned}
$$

重复 `scalingIteration` 次（默认几次即可，Ruiz 证明了线性收敛到均衡态）。

代价缩放 $c$ 让 $\|DHD\|$ 的量级与 $\|EGD\|$ 相当。

### 7.4.4 实现要点

`ocpDataInPlaceInParallel` **原地**修改 `dynamics_` 与 `cost_`，
并输出 `D`、`E`、`scalingVectors`、`c`。

**`scalingVectors` 的由来**：预条件前，动力学约束里 $x_{k+1}$ 的系数是 $-I$。
缩放后变成 $-E_kD_{x_{k+1}}$，是一个**一般对角矩阵**。
代码用 `scalingVectors[k]` 存它的对角线，
PIPG 里所有 `C.array() * X_[t].array()` 都是在用它。

`descaleSolution(D, deltaXSol, deltaUSol)` 做 $z=Dy$ 的还原。
`EInv` 传给 PIPG 用于**在原始尺度下判断收敛**——
否则缩放后的残差小不代表原问题残差小。

---

## 7.5 `SlpSettings` 与 `PipgSettings`

### `SlpSettings`

与 `SqpSettings` 几乎相同（滤子参数、`dt`、`integratorType`、`nThreads`…），
额外的：

| 参数 | 说明 |
|---|---|
| `scalingIteration` | Ruiz 均衡迭代次数（默认 3~5） |
| `pipgSettings` | 见下 |

注意**没有** `projectStateInputEqualityConstraints`——SLP 强制投影。

### `PipgSettings`

| 参数 | 说明 |
|---|---|
| `nThreads` | PIPG 内部线程数 |
| `maxNumIterations` | PIPG 最大迭代，默认 **3000** |
| `absoluteTolerance` | 原残差 $\infty$-范数阈值，默认 1e-3 |
| `relativeTolerance` | 解变化的相对阈值，默认 1e-2 |
| `checkTerminationInterval` | 每隔几次迭代检查一次收敛（检查本身有开销） |
| **`lowerBoundH`** | $\underline h$，估计 $\mu$ 用，默认 **5e-6**（**最关键的调参**） |
| `displayShortSummary` | 打印迭代次数与残差 |

---

## 7.6 `SingleThreadPipg.h`

单线程参考实现，用于：
- 调试多线程版本的正确性
- 小问题（多线程开销大于收益）
- 理解算法（代码短得多，建议先读它）

---

## 7.7 SLP 的定位

| 维度 | SLP/PIPG | SQP/HPIPM |
|---|---|---|
| 外部依赖 | **无** | HPIPM + BLASFEO |
| 每次迭代复杂度 | $O(Nn^2)$ | $O(Nn^3)$ |
| 迭代次数 | 数百~数千 | 10~30 |
| 精度 | 中（一阶方法） | 高 |
| 并行度 | **极高**（纯 BLAS-2） | 中（Riccati 有串行依赖） |
| 调参难度 | **高**（$\underline h$、容差、缩放次数） | 低 |
| GPU 友好 | ✅ | ❌ |

**适用场景**：
- 嵌入式平台，无法编译 HPIPM
- 问题规模很大（长 horizon），$O(n^3)$ 不可接受
- 有 GPU/众核加速器
- 只需要中等精度的解（MPC 的 RTI 场景本来就不追求高精度）

**不适用**：
- 病态问题（即使 Ruiz 也救不回来）
- $H$ 接近奇异（$\mu$ 无法估计）
- 需要高精度的离线轨迹优化

---

## 7.8 参考文献

- **Yu, Elango, Açıkmeşe (2020)**, *Proportional-Integral Projected Gradient Method for Model Predictive Control*, arXiv:2009.06980 —— **PIPG 原始论文**，代码 `PipgSolver.h:50` 与 `PipgBounds.h:40` 直接给出了这个链接
- **Yu, Elango, Açıkmeşe, Topcu (2022)**, *Extrapolated Proportional-Integral Projected Gradient Method for Conic Optimization*, IEEE L-CSS
- **Ruiz, D. (2001)**, *A Scaling Algorithm to Equilibrate Both Rows and Columns Norms in Matrices*, Technical Report RAL-TR-2001-034 —— **Ruiz 均衡**，代码注释直接引用
- **Stellato et al. (2020)**, *OSQP: An Operator Splitting Solver for Quadratic Programs*, Math. Prog. Comp. —— 同为一阶方法，也用 Ruiz 预条件，可对照阅读
- **Chambolle & Pock (2011)**, *A First-Order Primal-Dual Algorithm for Convex Problems with Applications to Imaging*, JMIV —— PDHG，PIPG 的思想来源之一

**下一章**：约束处理的罚函数体系 → [08 约束与罚函数](08-constraints-penalties.md)
