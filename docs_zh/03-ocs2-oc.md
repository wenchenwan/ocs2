# 03 · ocs2_oc：最优控制问题与多重打靶转写

`ocs2_oc` 处在 `ocs2_core`（原语）与各求解器之间，负责：

1. 定义**问题的完整载体** `OptimalControlProblem`
2. 把问题在标称轨迹处**线性二次化**（`LinearQuadraticApproximator`）
3. **前向仿真**（`rollout/`）
4. **多重打靶转写**（`multiple_shooting/`）——SQP/IPM/SLP 共用
5. **参考信号管理**与求解器同步（`synchronized_module/`）
6. mode schedule 变化时的**轨迹时间轴调整**（`trajectory_adjustment/`）

91 个文件、约 1.26 万行。

---

## 3.1 `OptimalControlProblem`：问题的完整载体

**文件**：`ocs2_oc/include/ocs2_oc/oc_problem/OptimalControlProblem.h:48`

它是一个纯数据结构体，由 21 个 `unique_ptr` 与 1 个裸指针组成：

```cpp
struct OptimalControlProblem {
  /* 代价 */
  std::unique_ptr<StateInputCostCollection> costPtr;              // l(t,x,u)
  std::unique_ptr<StateCostCollection>      stateCostPtr;         // l(t,x)
  std::unique_ptr<StateCostCollection>      preJumpCostPtr;       // φ_i(x(t_i^-))
  std::unique_ptr<StateCostCollection>      finalCostPtr;         // φ(x(t_f))

  /* 软约束（罚函数形式） */
  std::unique_ptr<StateInputCostCollection> softConstraintPtr;
  std::unique_ptr<StateCostCollection>      stateSoftConstraintPtr;
  std::unique_ptr<StateCostCollection>      preJumpSoftConstraintPtr;
  std::unique_ptr<StateCostCollection>      finalSoftConstraintPtr;

  /* 硬等式约束 */
  std::unique_ptr<StateInputConstraintCollection> equalityConstraintPtr;      // g1(t,x,u)=0
  std::unique_ptr<StateConstraintCollection>      stateEqualityConstraintPtr; // g2(t,x)=0
  std::unique_ptr<StateConstraintCollection>      preJumpEqualityConstraintPtr;
  std::unique_ptr<StateConstraintCollection>      finalEqualityConstraintPtr;

  /* 硬不等式约束 */
  std::unique_ptr<StateInputConstraintCollection> inequalityConstraintPtr;     // h1(t,x,u)≥0
  std::unique_ptr<StateConstraintCollection>      stateInequalityConstraintPtr;// h2(t,x)≥0
  std::unique_ptr<StateConstraintCollection>      preJumpInequalityConstraintPtr;
  std::unique_ptr<StateConstraintCollection>      finalInequalityConstraintPtr;

  /* 增广拉格朗日（8 个，对应上面 8 类约束） */
  std::unique_ptr<StateInputAugmentedLagrangianCollection> equalityLagrangianPtr;
  std::unique_ptr<StateAugmentedLagrangianCollection>      stateEqualityLagrangianPtr;
  std::unique_ptr<StateInputAugmentedLagrangianCollection> inequalityLagrangianPtr;
  std::unique_ptr<StateAugmentedLagrangianCollection>      stateInequalityLagrangianPtr;
  std::unique_ptr<StateAugmentedLagrangianCollection>      preJumpEqualityLagrangianPtr;
  std::unique_ptr<StateAugmentedLagrangianCollection>      preJumpInequalityLagrangianPtr;
  std::unique_ptr<StateAugmentedLagrangianCollection>      finalEqualityLagrangianPtr;
  std::unique_ptr<StateAugmentedLagrangianCollection>      finalInequalityLagrangianPtr;

  /* 动力学与辅助 */
  std::unique_ptr<SystemDynamicsBase> dynamicsPtr;
  std::unique_ptr<PreComputation>     preComputationPtr;
  const TargetTrajectories*           targetTrajectoriesPtr;   // 由 ReferenceManager 注入
};
```

### 3.1.1 三条正交维度

理解这 24 个字段的关键是三条维度的笛卡尔积：

| 维度 | 取值 |
|---|---|
| **时刻类型** | intermediate（中间）· preJump（事件前）· final（终端） |
| **变量依赖** | state-only（$x$）· state-input（$x,u$） |
| **处理方式** | cost（代价）· softConstraint（罚）· constraint（硬约束）· Lagrangian（增广拉格朗日） |

注意几个"缺席"的组合是有物理意义的：
- **preJump/final 没有 state-input 项**：事件时刻与终端时刻没有输入变量
- **软约束没有增广拉格朗日版本**：软约束本身就是罚，不需要乘子

### 3.1.2 为什么用 `unique_ptr` 而不是值

因为这些是**多态容器**，且需要"空"语义（`empty()` 判断快速跳过）。
代价是深拷贝构造函数必须手写（`ocs2_oc/src/oc_problem/OptimalControlProblem.cpp`），
每个 `unique_ptr` 都要 `->clone()`。

**为什么需要深拷贝？** 因为多线程求解器要给每个线程一份独立副本
（`optimalControlProblemStock_`），而 `PreComputation` 有状态。

### 3.1.3 组装示例

以四足为例（`ocs2_robotic_examples/ocs2_legged_robot/src/LeggedRobotInterface.cpp`）：

```cpp
problemPtr_->costPtr->add("baseTrackingCost", getBaseTrackingCost(...));
problemPtr_->softConstraintPtr->add("frictionCone_LF", getFrictionConeSoftConstraint(0));
problemPtr_->equalityConstraintPtr->add("zeroForce_LF",  getZeroForceConstraint(0));
problemPtr_->equalityConstraintPtr->add("zeroVelocity_LF", getZeroVelocityConstraint(0));
problemPtr_->dynamicsPtr = std::make_unique<LeggedRobotDynamicsAD>(...);
problemPtr_->preComputationPtr = std::make_unique<LeggedRobotPreComputation>(...);
```

注意**摩擦锥用软约束、零力/零速度用硬约束**：
前者是不等式且允许轻微违反，后者是必须精确满足的等式（用零空间投影消去）。

---

## 3.2 `oc_data/`：解与指标的数据结构

### `PrimalSolution`

```cpp
struct PrimalSolution {
  scalar_array_t timeTrajectory_;
  vector_array_t stateTrajectory_;
  vector_array_t inputTrajectory_;
  size_array_t   postEventIndices_;    // 事件后第一个采样点的索引
  ModeSchedule   modeSchedule_;
  std::unique_ptr<ControllerBase> controllerPtr_;
};
```

**`postEventIndices_` 的语义**：在事件时刻，时间轴上会有**两个采样点共享同一时间值**
（$t_i^-$ 与 $t_i^+$），因为状态在此不连续。`postEventIndices_[j]` 记录第 $j$ 个事件的
$t_i^+$ 在数组中的下标。于是 `postEventIndices_[j] - 1` 就是 $t_i^-$ 的下标——
这个模式在 DDP 里到处出现（如 `GaussNewtonDDP.cpp:671`）。

### `DualSolution`

```cpp
struct DualSolution {
  scalar_array_t timeTrajectory;
  size_array_t   postEventIndices;
  std::vector<MultiplierCollection> intermediates;
  std::vector<MultiplierCollection> preJumps;
  MultiplierCollection              final;
};
```

对应上面 `OptimalControlProblem` 的 8 个 Lagrangian 容器。

### `ProblemMetrics` 与 `PerformanceIndex`

- `ProblemMetrics`：逐点的 `Metrics`（原始数据）
- `PerformanceIndex`：标量汇总（8 个字段，见 `PerformanceIndex.h:42`）

```cpp
scalar_t merit;                    // 优点函数（线搜索的目标）
scalar_t cost;                     // 纯代价
scalar_t dualFeasibilitiesSSE;     // 对偶可行性残差平方和
scalar_t dynamicsViolationSSE;     // 动力学缺口平方和（多重打靶特有）
scalar_t equalityConstraintsSSE;
scalar_t inequalityConstraintsSSE;
scalar_t equalityLagrangian;       // 等式约束的拉格朗日/罚贡献
scalar_t inequalityLagrangian;
```

`dynamicsViolationSSE` 是多重打靶的核心指标：
单次打靶（DDP）中状态由前向积分产生，动力学**恒等满足**，该项恒为 0；
多重打靶中 $x_k$ 是独立决策变量，$\|x_{k+1}-F(x_k,u_k)\|^2$ 就是"缺口"。

### `TimeDiscretization`

```cpp
struct AnnotatedTime {
  scalar_t time;
  enum class Event { None, PreEvent, PostEvent } event;
};
std::vector<AnnotatedTime> timeDiscretizationWithEvents(
    scalar_t initTime, scalar_t finalTime, scalar_t dt, const scalar_array_t& eventTimes);
```

把 $[t_0,t_f]$ 按步长 `dt` 均匀切分，**并在事件时刻插入成对的 PreEvent/PostEvent 节点**。
辅助函数：

```cpp
scalar_t getIntervalStart(const AnnotatedTime&);
scalar_t getIntervalDuration(const AnnotatedTime& start, const AnnotatedTime& end);
```

这保证了多重打靶的网格与事件对齐——事件不会落在区间内部，
否则该区间的动力学是不连续的，RK 积分会失效。

---

## 3.3 `LinearQuadraticApproximator`：把问题线性二次化

**文件**：`ocs2_oc/src/approximate_model/LinearQuadraticApproximator.cpp`

三个入口函数，分别对应三类时刻：

| 函数 | 行号 | 填充的 `ModelData` 字段 |
|---|---|---|
| `approximateIntermediateLQ` | `:41` | dynamics, dynamicsCovariance, cost, stateEqConstraint, stateInputEqConstraint |
| `approximatePreJumpLQ` | `:89` | dynamics（跳变映射！）, cost, stateEqConstraint |
| `approximateFinalLQ` | `:127` | cost, stateEqConstraint（dynamics 为空） |

### 3.3.1 中间时刻的完整流程

```cpp
// 1) 触发预计算
constexpr auto request = Request::Cost + Request::SoftConstraint
                       + Request::Constraint + Request::Dynamics + Request::Approximation;
preComputation.request(request, time, state, input);

// 2) 动力学线性化
modelData.dynamics = problem.dynamicsPtr->linearApproximation(time, state, input, preComp);
modelData.dynamicsBias.setZero(...);

// 3) 代价二次化（含软约束）
modelData.cost = ocs2::approximateCost(problem, time, state, input);

// 4) 硬约束线性化
modelData.stateEqConstraint      = problem.stateEqualityConstraintPtr->getLinearApproximation(...);
modelData.stateInputEqConstraint = problem.equalityConstraintPtr->getLinearApproximation(...);

// 5) ⭐ 增广拉格朗日项加到代价上
if (!problem.stateEqualityLagrangianPtr->empty()) {
  auto approx = problem.stateEqualityLagrangianPtr->getQuadraticApproximation(
                    time, state, multipliers.stateEq, preComp);
  modelData.cost.f    += approx.f;
  modelData.cost.dfdx += approx.dfdx;
  modelData.cost.dfdxx+= approx.dfdxx;
}
// ... 另外 3 类
```

**第 5 步是关键**：增广拉格朗日的罚项被**直接吸收进代价**，
于是 Riccati 递推完全不需要知道约束的存在。这是 OCS2 处理不等式约束的主要机制。

对于 state-input 类的拉格朗日项（`:76-83`），
因为它们同时依赖 $x,u$，用的是完整的 `modelData.cost += approx`（含 `dfdu`、`dfduu`、`dfdux`）。

### 3.3.2 `approximateCost` 内部

```cpp
auto cost = problem.costPtr->getQuadraticApproximation(t, x, u, targetTraj, preComp);
if (!problem.stateCostPtr->empty())      cost += (state-only 项，只有 x 相关字段);
if (!problem.softConstraintPtr->empty()) cost += ...;
if (!problem.stateSoftConstraintPtr->empty()) cost += ...;
```

即：**代价 = 显式代价 + 状态代价 + 软约束罚 + 状态软约束罚**。
再加上 §3.3.1 第 5 步的增广拉格朗日，构成 Riccati 看到的"有效代价"。

### 3.3.3 `ChangeOfInputVariables.h`：输入变量替换

**文件**：`ocs2_oc/include/ocs2_oc/approximate_model/ChangeOfInputVariables.h`

实现输入替换 $u = P_u\tilde u + P_x x + u_0$ 下的近似量变换。
这是零空间投影的数学核心，两个重载：

**向量函数**（动力学、约束）：
$$
f(x,u)\approx f + A\,\delta x + B\,\delta u
\;\Longrightarrow\;
\tilde f = f+B u_0,\quad \tilde A = A + B P_x,\quad \tilde B = B P_u
$$

**标量函数**（代价）：
$$
\begin{aligned}
\tilde f &= f + f_u^{\!\top} u_0 + \tfrac12 u_0^{\!\top} f_{uu} u_0\\
\tilde f_x &= f_x + P_x^{\!\top} f_u + P_x^{\!\top} f_{uu} u_0 + f_{ux}^{\!\top}u_0\\
\tilde f_u &= P_u^{\!\top}\big(f_u + f_{uu}u_0\big)\\
\tilde f_{xx} &= f_{xx} + P_x^{\!\top}f_{ux} + f_{ux}^{\!\top}P_x + P_x^{\!\top}f_{uu}P_x\\
\tilde f_{ux} &= P_u^{\!\top}\big(f_{ux} + f_{uu}P_x\big)\\
\tilde f_{uu} &= P_u^{\!\top} f_{uu} P_u
\end{aligned}
$$

**推导**：把 $u=P_u\tilde u+P_x\delta x+u_0$（其中 $\delta u = P_u\tilde u + P_x\delta x + u_0$）
代入二次型

$$
f + f_x^{\!\top}\delta x + f_u^{\!\top}\delta u + \tfrac12\delta x^{\!\top}f_{xx}\delta x
+ \delta u^{\!\top}f_{ux}\delta x + \tfrac12 \delta u^{\!\top}f_{uu}\delta u
$$

按 $\tilde u$ 与 $\delta x$ 的幂次收集同类项即得。例如 $\tilde u^\top\tilde u$ 项：
只有 $\frac12\delta u^\top f_{uu}\delta u$ 贡献 $\frac12 \tilde u^\top P_u^\top f_{uu}P_u\tilde u$，
故 $\tilde f_{uu}=P_u^\top f_{uu}P_u$。

调用点：`DDP_HelperFunctions.cpp:143` 的 `projectLQ()`，
以及 `multiple_shooting/Transcription.cpp:113-117` 的 `projectTranscription()`。

---

## 3.4 `rollout/`：前向仿真

```
RolloutBase
 ├─ TimeTriggeredRollout    事件时刻预先已知（由 ModeSchedule 给定）
 ├─ StateTriggeredRollout   事件由保护面零穿越触发
 └─ InitializerRollout      用 Initializer 生成轨迹（无控制器时）
```

### `RolloutBase::run` 的签名

```cpp
vector_t run(scalar_t initTime, const vector_t& initState, scalar_t finalTime,
             ControllerBase* controller, ModeSchedule& modeSchedule,
             scalar_array_t& timeTrajectory, size_array_t& postEventIndices,
             vector_array_t& stateTrajectory, vector_array_t& inputTrajectory);
```

注意 `modeSchedule` 是**非 const 引用**——`StateTriggeredRollout` 会**写回**
它实际检测到的事件时刻，因为事件时刻在状态触发模式下是仿真的结果而非输入。

### `RolloutSettings`

```cpp
IntegratorType integratorType = IntegratorType::ODE45;
scalar_t absTolODE, relTolODE;       // 自适应步长容差
size_t maxNumStepsPerSecond;         // 防止步长塌缩导致死循环
scalar_t timeStep;                   // 定步长积分器用
RootFinderType rootFindingAlgorithm; // 状态触发用
```

### ⭐ `StateTriggeredRollout` 的事件定位

问题：给定保护面 $g(t,x(t))\in\mathbb{R}^m$，找到最早的 $t^*$ 使某个分量变号。

流程（`ocs2_oc/src/rollout/StateTriggeredRollout.cpp`）：

1. 正常积分，`StateTriggeredEventHandler` 在每步后检查 $g$ 的符号
2. 检测到某分量在 $[t_a,t_b]$ 内变号 → 抛异常中断积分
3. 用 `RootFinder` 在 $[t_a,t_b]$ 内做**括号法求根**
4. 找到 $t^*$ 后积分到 $t^*$，施加跳变映射，记录事件，从 $x^+$ 继续

### ⭐ `RootFinder` 的四种算法

**文件**：`ocs2_oc/include/ocs2_oc/rollout/RootFinder.h`

都属于 **Regula Falsi（试位法）家族**。设括号 $[t_0,t_1]$，$f_0=g(t_0)$、$f_1=g(t_1)$ 异号。
基本试位法给出

$$
t_2 = t_1 - f_1\,\frac{t_1-t_0}{f_1-f_0}
$$

朴素试位法的问题是**一端会"卡住"**（stagnation）：若函数凸，某一端点永远不更新，
收敛退化为线性。四种改进都是**缩放被保留端点的函数值** $f_0\leftarrow \gamma f_0$：

| 算法 | 缩放因子 $\gamma$ | 收敛阶 |
|---|---|---|
| Regula Falsi | $1$（不缩放） | 线性 |
| **Illinois** | $\tfrac12$ | ≈1.442 |
| **Pegasus** | $\dfrac{f_1}{f_1+f_2}$ | ≈1.642 |
| **Anderson–Björck**（默认） | $\begin{cases}1-\dfrac{f_2}{f_1} & \text{若}>0\\ \tfrac12&\text{否则}\end{cases}$ | ≈1.7 |

其中 $f_2=g(t_2)$ 是新点的函数值。Anderson–Björck 的因子来自对函数二阶行为的估计，
是这一族里收敛最快的，故设为默认（`RootFinder.h:74`）。

参考文献见 [17 章](17-references.md) 的 Anderson & Björck (1973)、Dowell & Jarratt (1971, 1972)、Ford (1995)。

### `PerformanceIndicesRollout`

沿给定轨迹用梯形积分计算代价与约束违反量，得到 `PerformanceIndex`。

---

## 3.5 ⭐ `multiple_shooting/`：多重打靶转写

这是 SQP、IPM、SLP **三个求解器共享**的核心，理解它就理解了三者的公共骨架。

### 3.5.1 多重打靶 vs 单次打靶

**单次打靶**（DDP 用）：决策变量只有 $u_{0:N-1}$，状态由前向积分唯一确定。
- 优点：动力学恒满足，变量少
- 缺点：对初值敏感（"射出去"后期误差放大），不适合不稳定系统

**多重打靶**：决策变量是 $(x_{0:N}, u_{0:N-1})$，动力学作为**等式约束**

$$
x_{k+1} - F(t_k, x_k, u_k, \Delta t_k) = 0
$$

- 优点：数值条件好、可并行、天然支持不可行初值（infeasible start）
- 缺点：变量多，需要能处理稀疏 KKT 的 QP 求解器（HPIPM 正是为此设计）

参考：Bock & Plitt (1984) 的经典多重打靶论文。

### 3.5.2 `Transcription.h`：单个节点的转写

```cpp
struct Transcription {
  ScalarFunctionQuadraticApproximation cost;
  VectorFunctionLinearApproximation dynamics;
  ConstraintsSize constraintsSize;
  VectorFunctionLinearApproximation stateEqConstraints, stateInputEqConstraints;
  VectorFunctionLinearApproximation stateIneqConstraints, stateInputIneqConstraints;
  VectorFunctionLinearApproximation constraintsProjection;
  ProjectionMultiplierCoefficients projectionMultiplierCoefficients;
};
```

`setupIntermediateNode()`（`Transcription.cpp:40`）做四件事：

**(1) 动力学离散化 + 缺口**

```cpp
dynamics = sensitivityDiscretizer(*problem.dynamicsPtr, t, x, u, dt);
dynamics.f -= x_next;   // 变成 δx_{k+1} = A δx_k + B δu_k + b
```

即 `dynamics.f` 存的是**缺口（gap/defect）**

$$
b_k \;=\; F(t_k,x_k,u_k,\Delta t) - x_{k+1}
$$

线性化后的约束为 $\delta x_{k+1} = A_k\,\delta x_k + B_k\,\delta u_k + b_k$。

**(2) 代价的前向欧拉积分**

```cpp
cost = approximateCost(problem, t, x, u);
cost *= dt;
```

即 $\int_{t_k}^{t_{k+1}} l\,\mathrm{d}t\approx \Delta t\cdot l(t_k,x_k,u_k)$。

> **注意**：这是**一阶**积分，而动力学可以用 RK4。代价与动力学的精度不匹配是有意为之：
> 代价的一阶误差只影响收敛速度，不影响可行性；而动力学误差直接导致轨迹不可行。

**(3) 约束线性化**：四类约束各自调 `getLinearApproximation()`。

**(4) 事件节点与终端节点**：
- `setupEventNode()`（`:156`）用 `jumpMapLinearApproximation()`，且 `dynamics.dfdu` 置零（事件时刻无输入）
- `setupTerminalNode()`（`:125`）只有代价与约束，无动力学

### 3.5.3 ⭐ `projectTranscription`：消去状态-输入等式约束

**文件**：`Transcription.cpp:96`

给定线性化的状态-输入等式约束

$$
C\,\delta x + D\,\delta u + e = 0,\qquad D\in\mathbb{R}^{n_c\times n_u},\ \operatorname{rank}(D)=n_c
$$

我们要找出**参数化全部解**的形式。这是一个欠定线性系统，通解 = 特解 + 零空间。

#### QR 分解构造

对 $D^{\!\top}\in\mathbb{R}^{n_u\times n_c}$ 做 Householder QR
（`LinearAlgebra.cpp:160` 的 `qrConstraintProjection`）：

$$
D^{\!\top} = Q R = \begin{bmatrix}Q_1 & Q_2\end{bmatrix}\begin{bmatrix}R_1\\0\end{bmatrix}
= Q_1 R_1
$$

其中 $Q_1\in\mathbb{R}^{n_u\times n_c}$、$Q_2\in\mathbb{R}^{n_u\times(n_u-n_c)}$、
$R_1\in\mathbb{R}^{n_c\times n_c}$ 上三角非奇异，且 $Q^\top Q=I$。

**关键性质**：$D Q_2 = (Q_1R_1)^{\!\top}Q_2 = R_1^{\!\top}Q_1^{\!\top}Q_2 = 0$，
所以 $Q_2$ 的列张成 $\mathcal{N}(D)$（$D$ 的零空间）。

**伪逆**：定义 $D^{\dagger\top} := Q_1 R_1^{-\top}$，代码中
```cpp
const matrix_t pseudoInverse = R.solve(Q1.transpose());   // = R_1^{-1} Q_1^T
```
存的是 $R_1^{-1}Q_1^{\!\top}\in\mathbb{R}^{n_c\times n_u}$，即 $D^\top$ 的**左伪逆**：
$(R_1^{-1}Q_1^{\!\top})D^{\!\top}=R_1^{-1}Q_1^{\!\top}Q_1R_1=I$。

于是 $D\,(R_1^{-1}Q_1^{\!\top})^{\!\top} = D\,Q_1R_1^{-\top}$。注意
$D=R_1^\top Q_1^\top$，故 $D Q_1 R_1^{-\top} = R_1^\top Q_1^\top Q_1 R_1^{-\top}=R_1^\top R_1^{-\top}=I$。
即 $(R_1^{-1}Q_1^\top)^\top$ 是 $D$ 的右逆。

#### 通解

$$
\boxed{\;\delta u = \underbrace{Q_2}_{P_u}\,\delta\tilde u
 \;\underbrace{-\,(R_1^{-1}Q_1^{\!\top})^{\!\top} C}_{P_x}\,\delta x
 \;\underbrace{-\,(R_1^{-1}Q_1^{\!\top})^{\!\top} e}_{u_0}\;}
$$

验证：代入约束
$$
C\delta x + D\big(Q_2\delta\tilde u - D^{+}C\delta x - D^{+}e\big) + e
= C\delta x + 0 - C\delta x - e + e = 0\;\checkmark
$$
（其中 $D^{+}=(R_1^{-1}Q_1^\top)^\top$ 满足 $DD^+=I$、$DQ_2=0$）

代码对应（`LinearAlgebra.cpp:172-176`）：

```cpp
projectionTerms.dfdu = Q.rightCols(numInputs - numConstraints);   // P_u = Q2
projectionTerms.dfdx = -pseudoInverse.transpose() * constraint.dfdx;  // P_x
projectionTerms.f    = -pseudoInverse.transpose() * constraint.f;     // u_0
```

#### 应用替换

`projectTranscription()` 随后对动力学、代价、状态-输入不等式约束
统一施加 §3.3.3 的变量替换：

```cpp
changeOfInputVariables(dynamics, projection.dfdu, projection.dfdx, projection.f);
changeOfInputVariables(cost,     projection.dfdu, projection.dfdx, projection.f);
if (stateInputIneqConstraints.f.size() > 0)
  changeOfInputVariables(stateInputIneqConstraints, ...);
stateInputEqConstraints = VectorFunctionLinearApproximation();   // 约束已消失
```

**结果**：子问题变成关于 $\tilde u\in\mathbb{R}^{n_u-n_c}$ 的**无状态-输入等式约束** QP。
这既减小了 QP 规模，又避免了 QP 求解器处理等式约束的数值困难。

#### LU 版本

`luConstraintProjection`（`LinearAlgebra.cpp:183`）用 `Eigen::FullPivLU` 的
`kernel()` 与 `solve()` 得到同样的 $P_u,P_x,u_0$，速度略快但不提供伪逆。
`SqpSettings` 里 `extractProjectionMultiplier` 决定用哪个：
需要恢复乘子时用 QR（因为需要伪逆），否则用 LU
（`Transcription.cpp:104-111`）。

### 3.5.4 ⭐ `ProjectionMultiplierCoefficients`：恢复被消去的乘子

投影消去了约束，但用户可能想知道**约束力**（例如四足的接触力对应的乘子）。
`ProjectionMultiplierCoefficients::compute()`（`ProjectionMultiplierCoefficients.cpp:35`）
从 KKT 条件反解。

#### 推导

节点的拉格朗日函数（只保留与 $u$ 相关项）：

$$
\mathcal{L} = \tfrac12\delta u^{\!\top}R\,\delta u + \delta u^{\!\top}P\,\delta x + r^{\!\top}\delta u
 + \lambda^{\!\top}\big(C\delta x + D\delta u + e\big) + \nu^{\!\top}\big(\dots + B\delta u\big)
$$

其中 $\lambda$ 是状态-输入等式约束乘子，$\nu$ 是动力学（协态）乘子。
稳定性条件 $\partial\mathcal{L}/\partial\delta u = 0$：

$$
R\,\delta u + P\,\delta x + r + D^{\!\top}\lambda + B^{\!\top}\nu = 0
$$

左乘 $D$ 的左伪逆的转置（即 $R_1^{-1}Q_1^{\!\top}$，记作 $M$，满足 $MD^{\!\top}=I$）：

$$
\lambda = -M\big(R\,\delta u + P\,\delta x + r + B^{\!\top}\nu\big)
$$

现在把 $\delta u = P_u\delta\tilde u + P_x\delta x + u_0$ 代入：

$$
\lambda = -M\Big(R(P_u\delta\tilde u + P_x\delta x + u_0) + P\delta x + r + B^{\!\top}\nu\Big)
$$

整理成关于 $(\delta x, \delta\tilde u, \nu)$ 的仿射函数：

$$
\boxed{
\begin{aligned}
\lambda &= \Lambda_x\,\delta x + \Lambda_u\,\delta\tilde u + \Lambda_\nu\,\nu + \lambda_0\\[2pt]
\Lambda_x &= -M\,(P + R\,P_x)\\
\Lambda_u &= -M\,R\,P_u\\
\Lambda_\nu &= -M\,B^{\!\top}\\
\lambda_0 &= -M\,(r + R\,u_0)
\end{aligned}}
$$

对照代码（变量名 `semiprojectedCost_*` 表示"只做了部分替换的代价"）：

```cpp
vector_t semiprojectedCost_dfdu  = cost.dfdu  + cost.dfduu * projection.f;      // r + R·u0
matrix_t semiprojectedCost_dfdux = cost.dfdux + cost.dfduu * projection.dfdx;   // P + R·Px
matrix_t semiprojectedCost_dfduu = cost.dfduu * projection.dfdu;               // R·Pu

this->dfdx      = -pseudoInverse * semiprojectedCost_dfdux;   // Λ_x
this->dfdu      = -pseudoInverse * semiprojectedCost_dfduu;   // Λ_u
this->dfdcostate= -pseudoInverse * dynamics.dfdu.transpose(); // Λ_ν
this->f         = -pseudoInverse * semiprojectedCost_dfdu;    // λ_0
```

完全一致。字段名 `dfdcostate` 里的 "costate" 就是协态 $\nu$。

### 3.5.5 其余辅助文件

| 文件 | 作用 |
|---|---|
| `Initialization.h/.cpp` | `initializeStateInputTrajectories()`：用上次解 + Initializer 填初值 |
| `Helpers.h/.cpp` | `incrementTrajectory()`、`trajectoryNorm()`、`remapProjectedInput()`、`remapProjectedGain()`、`toPrimalSolution()` |
| `MetricsComputation.h/.cpp` | 从 `Transcription` 或直接从问题计算 `Metrics` |
| `PerformanceIndexComputation.h/.cpp` | `Transcription` → `PerformanceIndex` |
| `LagrangianEvaluation.h/.cpp` | 计算拉格朗日函数及其导数（IPM 用） |

**`remapProjectedInput`** 是投影的逆操作：QP 解出 $\delta\tilde u$ 后，
用 $\delta u = P_u\delta\tilde u + P_x\delta x + u_0$ 还原真实输入增量。
**`remapProjectedGain`** 对反馈增益做同样的还原：
若 QP 给出 $\delta\tilde u = \tilde K\delta x$，则

$$
\delta u = P_u\tilde K\,\delta x + P_x\,\delta x + u_0
\;\Longrightarrow\;
K = P_u\tilde K + P_x
$$

---

## 3.6 `oc_problem/` 的其余文件

### `OcpSize`

```cpp
struct OcpSize {
  int numStages;
  std::vector<int> numStates, numInputs;
  std::vector<int> numIneqConstraints;
};
```

HPIPM 与 PIPG 都需要预先知道每个阶段的维度以分配内存。
`extractSizesFromProblem()` 从 dynamics/cost 数组反推出 `OcpSize`。

### ⭐ `OcpToKkt`：稀疏 KKT 矩阵组装

**文件**：`ocs2_oc/include/ocs2_oc/oc_problem/OcpToKkt.h`

把逐阶段的数据拼成全局稀疏矩阵。决策向量定义为

$$
Z = \begin{bmatrix}u_0\\ x_1\\ u_1\\ \vdots\\ u_n\\ x_{n+1}\end{bmatrix}
$$

（注意 $x_0$ 不在里面，因为初始状态固定）

**约束矩阵** $G Z = g$：

$$
G=\begin{bmatrix}
-B_0 & I &&&&\\
& -A_1 & -B_1 & I &&\\
&&&\ddots&&\\
&&&& -A_n & -B_n & I\\
D_0 & 0 &&&&\\
& C_1 & D_1 & 0 &&\\
&&&\ddots&&\\
&&&& C_n & D_n & 0
\end{bmatrix},\qquad
g=\begin{bmatrix}A_0x_0+b_0\\ b_1\\ \vdots\\ b_n\\ -(C_0x_0+e_0)\\ -e_1\\ \vdots\\ -e_n\end{bmatrix}
$$

上半部分是动力学，下半部分是一般等式约束。

**Hessian** $H$ 与梯度 $h$：

$$
H=\begin{bmatrix}
R_0 &&&&\\
& Q_1 & P_1^{\!\top} &&\\
& P_1 & R_1 &&\\
&&&\ddots&\\
&&&& Q_{n+1}
\end{bmatrix},\qquad
h=\begin{bmatrix}P_0x_0+r_0\\ q_1\\ r_1\\ \vdots\\ q_{n+1}\end{bmatrix}
$$

于是原问题变成标准的等式约束 QP：

$$
\min_Z\ \tfrac12 Z^{\!\top}HZ + h^{\!\top}Z\quad\text{s.t.}\quad GZ=g
$$

PIPG（`ocs2_slp`）直接在这个形式上工作。
HPIPM 则**不组装全局矩阵**——它利用块三对角结构做 Riccati 递推，
这是它比通用 QP 求解器快一个数量级的原因。

### `precondition/Ruzi.h`

Ruiz 均衡预条件，见 [07 章](07-ocs2-slp.md) §7.4。

### `LoopshapingOptimalControlProblem.h`

把普通 OCP 包装成 loopshaping 增广 OCP，见 [09 章](09-loopshaping.md)。

---

## 3.7 `oc_solver/SolverBase.h`：求解器统一接口

```cpp
class SolverBase {
 public:
  void run(scalar_t initTime, const vector_t& initState, scalar_t finalTime);
  void run(scalar_t, const vector_t&, scalar_t, const ControllerBase* externalController);

  virtual void reset() = 0;
  virtual size_t getNumIterations() const = 0;
  virtual const PerformanceIndex& getPerformanceIndeces() const = 0;
  virtual void getPrimalSolution(scalar_t finalTime, PrimalSolution*) const = 0;
  virtual const DualSolution* getDualSolution() const;
  virtual ScalarFunctionQuadraticApproximation getValueFunction(scalar_t, const vector_t&) const;
  virtual ScalarFunctionQuadraticApproximation getHamiltonian(scalar_t, const vector_t&, const vector_t&);
  virtual vector_t getStateInputEqualityConstraintLagrangian(scalar_t, const vector_t&) const;

  void addSynchronizedModule(std::shared_ptr<SolverSynchronizedModule>);
  void setReferenceManager(std::shared_ptr<ReferenceManagerInterface>);
  void addSolverObserver(std::shared_ptr<SolverObserver>);
};
```

`run()` 是非虚的模板方法，内部：

```
preSolverRun  →  referenceManager->preSolverRun()
              →  每个 synchronizedModule->preSolverRun()
runImpl(...)  ←  纯虚，由具体求解器实现
postSolverRun →  每个 synchronizedModule->postSolverRun()
              →  每个 solverObserver->extract...()
```

**`getValueFunction` 与 `getHamiltonian` 的用途**：MPC-Net 需要哈密顿量的
二次近似作为训练损失（见 [14 章](14-mpcnet.md)），
`getValueFunction` 则用于 Deep Value MPC 一类的工作。

---

## 3.8 `synchronized_module/`：参考信号与观测

### `ReferenceManager`

管理 `ModeSchedule` 与 `TargetTrajectories`，**线程安全**：

```
外部线程 setTargetTrajectories() → 写入 buffer（加锁）
求解器线程 preSolverRun()        → buffer swap 到 active（无锁读）
```

这保证了求解过程中参考不会被中途改变（否则 LQ 近似与 rollout 会不一致）。

**装饰器模式**：`ReferenceManagerDecorator` 允许叠加行为，
如 `LoopshapingReferenceManager`（把参考扩展到增广状态空间）、
`RosReferenceManager`（从 ROS 话题接收参考）。

四足的 `SwitchedModelReferenceManager`
（`ocs2_legged_robot/reference_manager/SwitchedModelReferenceManager.h`）
在 `preSolverRun()` 里：
1. 从 `GaitSchedule` 取出当前 horizon 的接触序列 → 更新 `ModeSchedule`
2. 调用 `SwingTrajectoryPlanner` 生成摆动腿轨迹

这是"步态规划器 → MPC"的接口。

### `SolverObserver`

在每次求解后**提取指定约束项的乘子/指标**并回调，用于监控与可视化。
`ocs2_ros_interfaces/synchronized_module/SolverObserverRosCallbacks.h`
提供把它们发布到 ROS 话题的回调工厂。

---

## 3.9 ⭐ `trajectory_adjustment/TrajectorySpreading`

**问题**：MPC 每个周期 mode schedule 可能变化（步态相位推进、用户改变步态）。
上一次的解（primal + dual）是在**旧的 mode schedule** 下算的，
直接拿来热启动会导致"事件时刻对不上"。

**解法**：`TrajectorySpreading` 把旧轨迹的时间轴"拉伸/压缩"到新的事件时刻上。

```cpp
Status set(const ModeSchedule& oldModeSchedule,
           const ModeSchedule& newModeSchedule,
           const scalar_array_t& oldTimeTrajectory);

struct Status { bool willTruncate; bool willPerformTrajectorySpreading; };
```

算法思路：
1. 匹配新旧 mode sequence 的**公共子序列**
2. 对每对匹配的模态区间 $[t_i^{old},t_{i+1}^{old}]\to[t_i^{new},t_{i+1}^{new}]$，
   做仿射时间重映射
3. 无法匹配的尾部被截断（`willTruncate`）

之后 `adjustTrajectory()` 按新索引重排数据，
`adjustTimeTrajectory()` 改写时间戳，`getPostEventIndices()` 给出新的事件索引。

调用点：
- `LineSearchStrategy::computeSolution`（`LineSearchStrategy.cpp:103`）——调整 dual solution
- `SqpSolver::runImpl`（`SqpSolver.cpp:202`）——调整上次的 primal solution

`TrajectorySpreadingHelperFunctions.h` 提供对 `PrimalSolution`/`DualSolution` 的便捷封装。

---

## 3.10 `search_strategy/FilterLinesearch`

⭐ 见 [05 章](05-ocs2-sqp.md) §5.5——这是 SQP/IPM/SLP 共用的滤子线搜索。

---

## 3.11 小结

`ocs2_oc` 的价值在于**把"问题"与"算法"解耦得足够干净**：

- `OptimalControlProblem` 定义了什么是问题
- `LinearQuadraticApproximator` 定义了怎么把问题变成 LQ 数据
- `multiple_shooting::Transcription` 定义了怎么把问题变成 QP 数据
- `RolloutBase` 定义了怎么前向仿真

四个求解器分别选用其中一部分：
DDP 用 `LinearQuadraticApproximator` + `Rollout`；
SQP/IPM/SLP 用 `Transcription`（内部也调 `LinearQuadraticApproximator` 的 `approximateCost`）。

**下一章**：SLQ 与 iLQR 的完整数学 → [04 ocs2_ddp](04-ocs2-ddp.md)
