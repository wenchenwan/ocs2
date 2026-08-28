# 05 · ocs2_sqp：多重打靶 SQP 与 HPIPM

`ocs2_sqp` 是 OCS2 里**最接近工业界主流 MPC 做法**的求解器：
多重打靶 + 稀疏 QP + 滤子线搜索。它也是四足机器人示例的默认求解器。

目录结构：

```
ocs2_sqp/
├── blasfeo_catkin/   BLASFEO（基础线性代数子程序，为小稠密矩阵优化）的 catkin 封装
├── hpipm_catkin/     HPIPM（高性能内点法）的 catkin 封装 + OCS2 接口
└── ocs2_sqp/         SqpSolver、SqpSettings、SqpLogging
```

---

## 5.1 SQP 的数学框架

### 5.1.1 NLP 形式

多重打靶把连续 OCP 变成有限维非线性规划。决策变量

$$
w = \big(x_0,u_0,x_1,u_1,\dots,u_{N-1},x_N\big)
$$

问题：

$$
\begin{aligned}
\min_w\quad & \Phi(w) := \sum_{k=0}^{N-1}\Delta t_k\, l(t_k,x_k,u_k) + \phi(x_N)\\
\text{s.t.}\quad
& x_0 = \hat x_0 &&\text{(初始状态)}\\
& x_{k+1} = F(t_k,x_k,u_k,\Delta t_k) &&k=0..N-1\quad\text{(动力学)}\\
& g_k(x_k,u_k)=0,\ \ h_k(x_k,u_k)\ge0 &&\text{(路径约束)}\\
& g_N(x_N)=0,\ \ h_N(x_N)\ge0 &&\text{(终端约束)}
\end{aligned}
$$

### 5.1.2 SQP 迭代

在当前迭代点 $w^{(i)}$ 处，SQP 求解**二次子问题**：

$$
\begin{aligned}
\min_{\delta w}\quad & \tfrac12\,\delta w^{\!\top}\mathcal{B}\,\delta w + \nabla\Phi^{\!\top}\delta w\\
\text{s.t.}\quad & \nabla g^{\!\top}\delta w + g = 0,\qquad \nabla h^{\!\top}\delta w + h \ge 0
\end{aligned}
$$

其中 $\mathcal{B}$ 理想上是拉格朗日函数的 Hessian $\nabla^2_{ww}\mathcal{L}$。

**⭐ OCS2 的选择**：$\mathcal{B}=\nabla^2\Phi$（**只用代价的 Hessian，忽略约束的曲率贡献**），
这称为 **Gauss-Newton SQP** 或 **constrained Gauss-Newton**。理由与 §2.3.3 相同：

- 约束的二阶导需要额外的 AD 开销
- $\nabla^2\Phi$ 若由 `StateInputGaussNewtonCostAd` 提供则天然 PSD，QP 一定可解
- 在 MPC 里只跑 1~3 次迭代，二阶收敛的优势体现不出来

代价：在约束高度弯曲处收敛慢。对机器人问题实测影响很小。

### 5.1.3 稀疏结构

关键在于：**多重打靶的 QP 具有块三对角 KKT 结构**。

$$
\mathcal{B}=\operatorname{blkdiag}\big(H_0,H_1,\dots,H_N\big),\qquad
H_k = \begin{bmatrix}Q_k & P_k^{\!\top}\\ P_k & R_k\end{bmatrix}
$$

动力学约束的 Jacobian 是**双对角块**：

$$
\begin{bmatrix}
A_0 & B_0 & -I\\
& & A_1 & B_1 & -I\\
& & & &\ddots
\end{bmatrix}
$$

通用稀疏 QP 求解器（OSQP、qpOASES）不利用这个结构，复杂度 $O(N^3 n^3)$。
**HPIPM 利用它做 Riccati 递推，复杂度降到 $O(N n^3)$**——这是它的核心优势。

---

## 5.2 `SqpSolver` 主循环

**文件**：`ocs2_sqp/ocs2_sqp/src/SqpSolver.cpp:183`

```cpp
void SqpSolver::runImpl(scalar_t initTime, const vector_t& initState, scalar_t finalTime) {
  // ① 时间离散化（含事件对齐）
  const auto& eventTimes = getReferenceManager().getModeSchedule().eventTimes;
  const auto timeDiscretization = timeDiscretizationWithEvents(initTime, finalTime, settings_.dt, eventTimes);

  // ② 注入参考轨迹
  for (auto& ocpDefinition : ocpDefinitions_)
      ocpDefinition.targetTrajectoriesPtr = &getReferenceManager().getTargetTrajectories();

  // ③ mode schedule 变了 → 调整上次的解
  if (!primalSolution_.timeTrajectory_.empty())
    trajectorySpread(primalSolution_.modeSchedule_, getReferenceManager().getModeSchedule(), primalSolution_);

  // ④ 初始化 x, u
  vector_array_t x, u;
  multiple_shooting::initializeStateInputTrajectories(initState, timeDiscretization,
                                                      primalSolution_, *initializerPtr_, x, u);

  // ⑤ SQP 主循环
  int iter = 0;
  sqp::Convergence convergence = sqp::Convergence::FALSE;
  while (convergence == sqp::Convergence::FALSE) {
    const auto baselinePerformance = setupQuadraticSubproblem(timeDiscretization, initState, x, u, metrics);
    const vector_t delta_x0 = initState - x[0];
    const auto deltaSolution = getOCPSolution(delta_x0);
    extractValueFunction(timeDiscretization, x);
    const auto stepInfo = takeStep(baselinePerformance, timeDiscretization, initState, deltaSolution, x, u, metrics);
    convergence = checkConvergence(iter, baselinePerformance, stepInfo);
    ++iter;  ++totalNumIterations_;
  }

  // ⑥ 组装解
  primalSolution_ = toPrimalSolution(timeDiscretization, std::move(x), std::move(u));
  problemMetrics_ = multiple_shooting::toProblemMetrics(timeDiscretization, std::move(metrics));
}
```

**注意 `delta_x0 = initState - x[0]`**：多重打靶允许 $x_0$ 与真实初始状态不同
（infeasible start），QP 的第一个约束是 $\delta x_0 = \hat x_0 - x_0$，
把这个"初始缺口"也交给 QP 一起消除。这是多重打靶相对单次打靶的一大优势。

---

## 5.3 QP 子问题的构造

**文件**：`SqpSolver.cpp:336` 的 `setupQuadraticSubproblem`

```cpp
std::atomic_int timeIndex{0};
auto parallelTask = [&](int workerId) {
  OptimalControlProblem& ocpDefinition = ocpDefinitions_[workerId];   // 每线程一份
  PerformanceIndex workerPerformance;                                  // 线程局部累加

  int i = timeIndex++;                       // 原子任务分发
  while (i < N) {
    if (time[i].event == AnnotatedTime::Event::PreEvent) {
      auto result = multiple_shooting::setupEventNode(ocpDefinition, time[i].time, x[i], x[i+1]);
      // 事件节点：无输入，用跳变映射
      ...
    } else {
      const scalar_t ti = getIntervalStart(time[i]);
      const scalar_t dt = getIntervalDuration(time[i], time[i+1]);
      auto result = multiple_shooting::setupIntermediateNode(ocpDefinition, sensitivityDiscretizer_,
                                                             ti, dt, x[i], x[i+1], u[i]);
      if (settings_.projectStateInputEqualityConstraints)
        multiple_shooting::projectTranscription(result, settings_.extractProjectionMultiplier);
      cost_[i] = std::move(result.cost);
      dynamics_[i] = std::move(result.dynamics);
      ...
    }
    i = timeIndex++;
  }
  if (i == N) {   // 恰好一个 worker 处理终端节点
    auto result = multiple_shooting::setupTerminalNode(ocpDefinition, getIntervalStart(time[N]), x[N]);
    ...
  }
  performance[workerId] += workerPerformance;
};
runParallel(std::move(parallelTask));

// 初始状态缺口计入性能指标
const vector_t initDynamicsViolation = initState - x.front();
metrics.front().dynamicsViolation += initDynamicsViolation;
performance.front().dynamicsViolationSSE += initDynamicsViolation.squaredNorm();
```

三点值得注意：

1. **无锁并行**：`std::atomic_int timeIndex` 做工作窃取式分发，负载自动均衡
2. **线程局部累加**：`workerPerformance` 避免 false sharing，最后一次性归约
3. **投影可选**：`projectStateInputEqualityConstraints` 决定是消去状态-输入等式约束（§3.5.3），
   还是把它们交给 HPIPM 处理

### 投影 vs 交给 QP

`getOCPSolution`（`:280`）里的分支：

```cpp
const bool hasStateInputConstraints = !ocpDefinitions_.front().equalityConstraintPtr->empty();
if (hasStateInputConstraints && !settings_.projectStateInputEqualityConstraints) {
  hpipmInterface_.resize(extractSizesFromProblem(dynamics_, cost_, &stateInputEqConstraints_));
  status = hpipmInterface_.solve(delta_x0, dynamics_, cost_, &stateInputEqConstraints_, ...);
} else {   // 无约束，或已投影 → 无约束 QP
  hpipmInterface_.resize(extractSizesFromProblem(dynamics_, cost_, nullptr));
  status = hpipmInterface_.solve(delta_x0, dynamics_, cost_, nullptr, ...);
}
```

**投影的好处**：QP 变量少 $n_c$ 个、无等式约束、HPIPM 走更快的代码路径。
**投影的代价**：需要 $D$ 满行秩；乘子要额外恢复（§3.5.4）。

对四足：每条支撑腿贡献 3 个零速度约束，4 条腿全支撑时消去 12 个变量，
QP 规模减小约 30%。默认开启投影。

---

## 5.4 HPIPM 接口

**文件**：`ocs2_sqp/hpipm_catkin/src/HpipmInterface.cpp`

### 5.4.1 HPIPM 是什么

**HPIPM**（High-Performance Interior-Point Method，Frison & Diehl）
是专门为**结构化最优控制 QP** 设计的求解器：

- 内部用 **Riccati 递推**（而非通用稀疏 LU）求解 KKT 系统
- 建立在 **BLASFEO** 之上——一个为"小到中等尺寸稠密矩阵"优化的 BLAS 实现，
  用手工 SIMD 内核与 panel-major 数据布局，在 $n<100$ 时比 OpenBLAS 快数倍
- 支持三种模式：`SPEED_ABS`、`SPEED`、`BALANCE`、`ROBUST`（精度/速度权衡）

### 5.4.2 接口做的事

```cpp
class HpipmInterface {
  void resize(OcpSize ocpSize);       // 重新分配 HPIPM 内部内存（维度变化时）
  hpipm_status solve(const vector_t& x0,
                     std::vector<VectorFunctionLinearApproximation>& dynamics,
                     std::vector<ScalarFunctionQuadraticApproximation>& cost,
                     std::vector<VectorFunctionLinearApproximation>* constraints,
                     vector_array_t& stateTrajectory, vector_array_t& inputTrajectory,
                     bool verbose);
  std::vector<ScalarFunctionQuadraticApproximation> getRiccatiCostToGo(dynamics0, cost0);
  matrix_array_t getRiccatiFeedback(dynamics0, cost0);
};
```

`resize()` 是**昂贵操作**（重新分配对齐内存），所以只在 `OcpSize` 变化时调用。
MPC 稳态运行时 horizon 与维度不变，`resize` 不会触发。

**所有约束映射为不等式**：HPIPM 的接口用 $\underline b\le Cx+Du\le\bar b$，
等式约束通过 $\underline b=\bar b$ 表达。

### 5.4.3 ⭐ 从 HPIPM 取回 Riccati 解

这是一个巧妙的设计。HPIPM 内部本来就做 Riccati 递推，
`getRiccatiCostToGo` 与 `getRiccatiFeedback` 把它的中间结果**直接暴露出来**，
于是 SQP 免费获得了：

- **值函数** $V_k(x)=\frac12 x^\top S_k x + s_k^\top x$（用于 MPC-Net、Deep Value MPC）
- **反馈增益** $K_k$（用于高频反馈控制）

`extractValueFunction`（`SqpSolver.cpp:311`）：

```cpp
if (settings_.createValueFunction) {
  valueFunction_ = hpipmInterface_.getRiccatiCostToGo(dynamics_[0], cost_[0]);
  for (int i = 0; i < time.size(); ++i)
    valueFunction_[i].dfdx.noalias() -= valueFunction_[i].dfdxx * x[i];   // 相对→绝对坐标
}
```

坐标修正的推导见 [04 章](04-ocs2-ddp.md) §4.2.5。

**为什么需要 `dynamics0`、`cost0`？** 因为 HPIPM 的 Riccati 从 $k=N$ 递推到 $k=1$，
$k=0$ 的值函数需要用第 0 阶段的数据额外补一步。

`toPrimalSolution`（`:321`）：

```cpp
if (settings_.useFeedbackPolicy) {
  matrix_array_t KMatrices = hpipmInterface_.getRiccatiFeedback(dynamics_[0], cost_[0]);
  if (settings_.projectStateInputEqualityConstraints)
    multiple_shooting::remapProjectedGain(constraintsProjection_, KMatrices);   // K = Pu·K̃ + Px
  return multiple_shooting::toPrimalSolution(time, modeSchedule, x, u, KMatrices);
}
```

---

## 5.5 ⭐ 滤子线搜索（Filter Line Search）

**文件**：`ocs2_oc/src/search_strategy/FilterLinesearch.cpp:34`

这是 SQP/IPM/SLP **三者共用**的步长接受准则，来自 **IPOPT**
（Wächter & Biegler 2006）。

### 5.5.1 动机：两个相互冲突的目标

约束优化里我们同时想要：
1. **降低代价** $\Phi$
2. **降低约束违反** $\theta := \|\text{constraint violation}\|$

传统做法是合成一个优点函数 $M=\Phi+\rho\theta$，但 $\rho$ 难调：
太小则不可行，太大则收敛慢，且可能出现 **Maratos 效应**
（好的步长因优点函数不降而被拒绝，导致收敛退化）。

**滤子方法**：不合成，而是把 $(\theta,\Phi)$ 当作**二维目标**，
只要在 Pareto 意义下"不被支配"就接受。

### 5.5.2 OCS2 的实现

约束违反的定义（`FilterLinesearch.h`）：

$$
\theta = \sqrt{\text{dynamicsViolationSSE}} + \sqrt{\text{equalityConstraintsSSE}}
$$

```cpp
static scalar_t totalConstraintViolation(const PerformanceIndex& p) {
  return std::sqrt(p.dynamicsViolationSSE + p.equalityConstraintsSSE);
}
```

`acceptStep()` 的三分支逻辑：

```cpp
if (θ_new > g_max) {
  // ① 违反太大：只看可行性
  accepted = θ_new < (1 - γ_c) * θ_old;
  stepType = CONSTRAINT;
}
else if (θ_new < g_min && θ_old < g_min && armijoDescentMetric < 0.0) {
  // ② 违反很小且有下降方向：走标准 Armijo
  accepted = Φ_new < Φ_old + armijoFactor * armijoDescentMetric;
  stepType = COST;
}
else {
  // ③ 中间地带：代价降 或 违反降，二者居一
  accepted = Φ_new < (Φ_old - γ_c*θ_old) || θ_new < (1-γ_c)*θ_old;
  stepType = DUAL;
}
```

（代码中 `Φ` 对应 `performance.merit`，`θ` 对应 `totalConstraintViolation`）

| 参数 | 默认 | 含义 |
|---|---|---|
| `g_max` | 1e6 | 违反上界，超过则进入"恢复模式"（分支 ①） |
| `g_min` | 1e-6 | 违反下界，低于则视为可行（分支 ②） |
| `gamma_c` | 1e-6 | 滤子的相对下降要求 |
| `armijoFactor` | 1e-4 | Armijo 系数 $c_1$ |

**三个分支的直观解释**：

- **①（远离可行域）**：先别管代价了，把约束降下来再说。
  要求违反至少下降 $\gamma_c$ 的相对量。
- **②（已经可行）**：这是无约束优化，用标准 Armijo。
- **③（接近可行）**：任一目标改善即可。
  `- γ_c*θ_old` 是一个小的"松弛裕度"，防止在 Pareto 前沿上原地打转。

### 5.5.3 Armijo 下降度量

```cpp
scalar_t armijoDescentMetric(const std::vector<ScalarFunctionQuadraticApproximation>& cost,
                             const vector_array_t& dx, const vector_array_t& du) {
  scalar_t metric = 0.0;
  for (int i = 0; i < cost.size(); i++) {
    if (cost[i].dfdx.size() > 0) metric += cost[i].dfdx.dot(dx[i]);
    if (cost[i].dfdu.size() > 0) metric += cost[i].dfdu.dot(du[i]);
  }
  return metric;
}
```

即 $\nabla\Phi^{\!\top}\delta w = \sum_k\big(q_k^\top\delta x_k + r_k^\top\delta u_k\big)$，
**代价梯度沿 QP 步方向的投影**。

注意：这里用的是 $\delta u$（remap 之后的真实输入增量），
而不是 $\delta\tilde u$，因为 `armijoDescentMetric` 在
`remapProjectedInput` **之前**调用（`SqpSolver.cpp:301` vs `:304`）。
—— 等等，仔细看代码：

```cpp
solution.armijoDescentMetric = armijoDescentMetric(cost_, deltaXSol, deltaUSol);   // :301
if (settings_.projectStateInputEqualityConstraints)
  multiple_shooting::remapProjectedInput(constraintsProjection_, deltaXSol, deltaUSol);  // :304
```

`cost_` 存的是**已投影的代价**（`projectTranscription` 改写了它），
`deltaUSol` 此时也还是 $\delta\tilde u$。所以两者坐标系一致 ✓
这是正确的：$\tilde\nabla\Phi^\top\delta\tilde w$ 与 $\nabla\Phi^\top\delta w$
在投影约束流形上相等。

### 5.5.4 回溯循环

**文件**：`SqpSolver.cpp:468` 的 `takeStep`

```cpp
scalar_t alpha = 1.0;
do {
  multiple_shooting::incrementTrajectory(u, du, alpha, uNew);   // u + α·du
  multiple_shooting::incrementTrajectory(x, dx, alpha, xNew);
  const PerformanceIndex performanceNew = computePerformance(time, initState, xNew, uNew, metricsNew);

  std::tie(stepAccepted, stepType) =
      filterLinesearch_.acceptStep(baseline, performanceNew, alpha * subproblemSolution.armijoDescentMetric);

  if (stepAccepted) {  x = move(xNew); u = move(uNew); return stepInfo; }
  else {
    alpha *= settings_.alpha_decay;                       // 默认 0.5
    if (alpha*deltaXnorm < deltaTol && alpha*deltaUnorm < deltaTol) break;  // 提前退出
  }
} while (alpha >= settings_.alpha_min);                   // 默认 1e-4
```

**注意 $\alpha$ 同时缩放 $\delta x$ 和 $\delta u$**——这与 DDP 不同！
DDP 只缩放前馈项（因为状态由积分产生），多重打靶里 $x$ 是独立变量，必须一起缩放。

**提前退出**：若 $\alpha\|\delta x\|$ 与 $\alpha\|\delta u\|$ 都小于 `deltaTol`，
说明步长已无意义，直接跳出而不必一路缩到 `alpha_min`——省下若干次昂贵的 `computePerformance`。

---

## 5.6 收敛判据

**文件**：`SqpSolver.cpp:563`

```cpp
enum class Convergence { FALSE, ITERATIONS, STEPSIZE, METRICS, PRIMAL };
```

| 返回值 | 条件 |
|---|---|
| `ITERATIONS` | `iteration + 1 >= sqpIteration` |
| `STEPSIZE` | `stepSize < alpha_min`（线搜索失败） |
| `METRICS` | $\|\Delta M\| <$ `costTol` **且** $\theta <$ `g_min` |
| `PRIMAL` | $\|\delta x\|<$ `deltaTol` **且** $\|\delta u\|<$ `deltaTol` |
| `FALSE` | 以上都不满足 |

**MPC 场景下通常 `sqpIteration = 1~3`**，即基本总是以 `ITERATIONS` 退出——
这是**实时迭代（Real-Time Iteration, RTI）**策略：
不追求单个周期内收敛，而是让 MPC 的滚动本身充当外层迭代。
Diehl et al. (2005) 证明了 RTI 在一定条件下的收敛性与稳定性。

---

## 5.7 `SqpSettings` 完全参考

**文件**：`ocs2_sqp/ocs2_sqp/include/ocs2_sqp/SqpSettings.h`

| 参数 | 默认 | 说明 |
|---|---|---|
| `sqpIteration` | 10 | 最大 SQP 迭代（MPC 常设 1~3） |
| `deltaTol` | 1e-6 | 原变量步长收敛阈值 |
| `costTol` | 1e-4 | 代价收敛阈值 |
| `alpha_decay` | 0.5 | 线搜索收缩率 |
| `alpha_min` | 1e-4 | 最小步长 |
| `g_max` / `g_min` | 1e6 / 1e-6 | 滤子阈值 |
| `armijoFactor` | 1e-4 | Armijo 系数 |
| `gamma_c` | 1e-6 | 滤子相对下降 |
| `projectStateInputEqualityConstraints` | true | 是否零空间投影 |
| `extractProjectionMultiplier` | false | 是否恢复投影乘子（需要 QR 而非 LU） |
| `useFeedbackPolicy` | true | 输出反馈控制器 |
| `createValueFunction` | false | 提取值函数 |
| `dt` | 0.01 | 时间离散步长 |
| `integratorType` | `RK2` | 灵敏度离散化方法 |
| `nThreads` | 4 | 线程数 |
| `threadPriority` | 50 | 实时优先级 |
| `hpipmSettings` | — | 见 `HpipmInterfaceSettings.h` |
| `enableLogging` | false | 记录每次迭代的详细日志 |
| `printSolverStatus` / `printSolverStatistics` / `printLinesearch` | false | 打印控制 |

---

## 5.8 `SqpLogging`：迭代日志

**文件**：`ocs2_sqp/ocs2_sqp/src/SqpLogging.cpp`

```cpp
struct LogEntry {
  int problemNumber, iteration;
  scalar_t time;
  scalar_t linearQuadraticApproximationTime, solveQpTime, linesearchTime;
  PerformanceIndex baselinePerformanceIndex;
  scalar_t totalConstraintViolationBaseline;
  StepInfo stepInfo;
  Convergence convergence;
};
```

可导出 CSV 供离线分析。调参时非常有用：
能看出时间花在 LQ 近似、QP 求解还是线搜索上，从而针对性优化。

典型分布（四足，12 状态 24 输入，horizon 1s，dt=15ms）：
LQ 近似 ~40%、QP ~35%、线搜索 ~25%。

---

## 5.9 测试用例（很好的学习材料）

| 文件 | 内容 |
|---|---|
| `test/testUnconstrained.cpp` | 无约束 LQ 问题，与解析解对比 |
| `test/testCircularKinematics.cpp` | 圆周运动跟踪，带状态-输入约束 |
| `test/testSwitchedProblem.cpp` | 带事件的切换系统 |
| `test/testValuefunction.cpp` | 验证 `getValueFunction` 与有限差分一致 |
| `hpipm_catkin/test/testHpipmInterface.cpp` | HPIPM 接口的独立测试 |

`ocs2_test_tools/ocs2_qp_solver` 提供了一个**稠密 QP 参考实现**
（把整个 KKT 系统组装成大矩阵直接求逆），用于交叉验证 HPIPM 的结果。
它慢但绝对可靠，是排查"HPIPM 是不是算错了"的黄金标准。

其 KKT 组装见 `ocs2_test_tools/ocs2_qp_solver/src/QpSolver.cpp:103`（约束矩阵）
与 `:162`（代价矩阵），结构与 [03 章](03-ocs2-oc.md) §3.6 的 `OcpToKkt` 相同。

---

## 5.10 SQP vs DDP：什么时候用哪个

| 维度 | SQP（多重打靶） | SLQ/iLQR（单次打靶） |
|---|---|---|
| 决策变量 | $x$ 和 $u$ 都是 | 只有 $u$ |
| 动力学 | 等式约束（允许缺口） | 恒满足 |
| 不可行初值 | ✅ 天然支持 | ❌ 必须先 rollout |
| 不稳定系统 | ✅ 数值稳定 | ❌ 误差指数放大 |
| 不等式约束 | ✅ QP 直接处理 | 需罚函数/增广拉格朗日 |
| 内存 | 大（存所有 $x_k$） | 小 |
| 依赖 | 需要 HPIPM+BLASFEO | 无外部 QP 依赖 |
| 时域 | 离散 | SLQ 连续 / iLQR 离散 |
| 事件处理 | 网格对齐 | 积分中断 |

**经验法则**：
- 有硬不等式约束（力矩限、摩擦锥硬约束）→ **SQP** 或 **IPM**
- 系统开环不稳定（人形、双足）→ **SQP**
- 追求最少依赖、连续时间精度 → **SLQ**
- 极端算力受限的嵌入式 → **iLQR** 或 **SLP**

---

## 5.11 参考文献

- **Bock & Plitt (1984)**, *A Multiple Shooting Algorithm for Direct Solution of Optimal Control Problems*, IFAC — 多重打靶原典
- **Frison & Diehl (2020)**, *HPIPM: a high-performance quadratic programming framework for model predictive control*, IFAC
- **Frison et al. (2018)**, *BLASFEO: Basic Linear Algebra Subroutines For Embedded Optimization*, ACM TOMS
- **Wächter & Biegler (2006)**, *On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming*, Math. Program. 106(1) — **滤子线搜索**（代码 `SqpSolver.cpp:476` 直接给出了这个链接）
- **Diehl, Bock, Schlöder (2005)**, *A Real-Time Iteration Scheme for Nonlinear Optimization in Optimal Feedback Control*, SIAM J. Control Optim.
- **Grandia et al. (2019)**, *Feedback MPC for Torque-Controlled Legged Robots*, IROS
- **Nocedal & Wright (2006)**, *Numerical Optimization*, Ch. 18（SQP）

**下一章**：内点法 → [06 ocs2_ipm](06-ocs2-ipm.md)
