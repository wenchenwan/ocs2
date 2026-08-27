# 01 · 总览与整体架构

## 1.1 OCS2 是什么

OCS2 = **O**ptimal **C**ontrol for **S**witched **S**ystems，由 ETH Zürich 机器人系统实验室
（Robotic Systems Lab, RSL）开发的 C++ 最优控制工具箱。它的定位不是"又一个 MPC 库"，
而是围绕一个具体痛点设计的：

> **腿足机器人的动力学是分段的。** 四足小跑（trot）时，接触状态在 4 个支撑相之间周期性切换；
> 每次落足都会产生一个速度不连续（跳变映射）。经典 MPC 假设动力学光滑连续，遇到这类问题
> 要么把接触力当作平滑变量（丢失物理意义），要么退化成混合整数规划（实时性崩溃）。

OCS2 的做法是：**接触序列（mode schedule）由外部规划器给定**，在给定序列下动力学分段光滑，
于是可以用连续时间的 DDP/SQP 高效求解，并在跳变点上正确传播值函数与灵敏度。

### 支持的算法

| 算法 | Package | 时域 | 子问题求解 | 约束处理 |
|---|---|---|---|---|
| **SLQ** | `ocs2_ddp` | 连续时间 | 微分 Riccati（ODE 积分） | 投影 + 增广拉格朗日 / 罚函数 |
| **iLQR** | `ocs2_ddp` | 离散时间 | 差分 Riccati（递推） | 同上 |
| **SQP** | `ocs2_sqp` | 多重打靶 | HPIPM（结构化内点 QP） | QP 内建 + 投影 |
| **IPM** | `ocs2_ipm` | 多重打靶 | HPIPM（无约束 QP）+ 障碍项 | 原对偶内点法 |
| **SLP** | `ocs2_slp` | 多重打靶 | PIPG（一阶原对偶） | 投影 |
| PISOC | — | — | 路径积分（README 提及，本仓库未含独立实现） | — |

---

## 1.2 顶层架构：三层 + 两条正交轴

```
┌─────────────────────────────────────────────────────────────────────┐
│  应用层                                                              │
│  ocs2_robotic_examples/{ballbot,cartpole,quadrotor,legged_robot,…}  │
│  ocs2_*_ros  ← ROS 节点、可视化、键盘/交互标记指令                    │
│  ocs2_mpcnet ← 从 MPC 蒸馏策略                                       │
└───────────────────────────────┬─────────────────────────────────────┘
                                │ 继承 RobotInterface，装配 OptimalControlProblem
┌───────────────────────────────▼─────────────────────────────────────┐
│  建模层                                                              │
│  ocs2_pinocchio        URDF → 运动学/动力学 → 质心动力学/自碰撞       │
│  ocs2_robotic_tools    旋转表示、角速度映射、末端运动学接口           │
│  ocs2_perceptive       ESDF 距离场、障碍距离约束                     │
└───────────────────────────────┬─────────────────────────────────────┘
                                │
┌───────────────────────────────▼─────────────────────────────────────┐
│  求解层                                                              │
│                                                                      │
│    ┌──────────── 问题描述（与算法无关） ────────────┐                 │
│    │  ocs2_core : Types / Cost / Constraint /      │                 │
│    │              Dynamics / Penalty / Integrator  │                 │
│    │  ocs2_oc   : OptimalControlProblem            │                 │
│    │              LinearQuadraticApproximator      │                 │
│    │              Rollout / TimeDiscretization     │                 │
│    │              multiple_shooting::Transcription │                 │
│    └───────────────────┬───────────────────────────┘                 │
│                        │                                             │
│    ┌───────────────────▼───────────────────────────┐                 │
│    │  求解器（都实现 SolverBase 接口）              │                 │
│    │  ocs2_ddp  ocs2_sqp  ocs2_ipm  ocs2_slp       │                 │
│    └───────────────────┬───────────────────────────┘                 │
│                        │                                             │
│    ┌───────────────────▼───────────────────────────┐                 │
│    │  ocs2_mpc : MPC_BASE (滚动优化) + MRT (执行)   │                 │
│    └───────────────────────────────────────────────┘                 │
└─────────────────────────────────────────────────────────────────────┘
```

**两条正交轴**是理解 OCS2 的关键：

- **轴 A（问题）**：用户只需填 `OptimalControlProblem` 的各个字段。
  同一个问题可以喂给任意求解器。
- **轴 B（算法）**：求解器只依赖 `OptimalControlProblem` 的抽象接口，
  不知道背后是四足还是四旋翼。

这条分离在代码上体现为：`SolverBase`（`ocs2_oc/oc_solver/SolverBase.h`）
是唯一的求解器基类，`GaussNewtonDDP`、`SqpSolver`、`IpmSolver`、`SlpSolver` 全部继承它。

---

## 1.3 核心数据结构总览

| 结构体 | 定义位置 | 作用 |
|---|---|---|
| `OptimalControlProblem` | `ocs2_oc/oc_problem/OptimalControlProblem.h:48` | 问题的全部内容（代价、约束、动力学、乘子） |
| `PrimalSolution` | `ocs2_oc/oc_data/PrimalSolution.h` | 原变量解：时间/状态/输入轨迹 + 控制器 + mode schedule |
| `DualSolution` | `ocs2_oc/oc_data/DualSolution.h` | 对偶变量解：各类约束的拉格朗日乘子轨迹 |
| `ProblemMetrics` | `ocs2_oc/oc_data/ProblemMetrics.h` | 沿轨迹逐点的代价/约束违反量 |
| `PerformanceIndex` | `ocs2_oc/oc_data/PerformanceIndex.h:42` | 标量汇总：merit、cost、各类 SSE |
| `ModelData` | `ocs2_core/model_data/ModelData.h` | 单个时间点的 LQ 近似（动力学 + 代价 + 约束的一二阶导） |
| `ModeSchedule` | `ocs2_core/reference/ModeSchedule.h` | 事件时刻 + 模态序列 |
| `TargetTrajectories` | `ocs2_core/reference/TargetTrajectories.h` | 期望状态/输入轨迹（跟踪目标） |
| `Metrics` / `Multiplier` | `ocs2_core/model_data/{Metrics,Multiplier}.h` | 单点指标与乘子 |

### 线性/二次近似的统一表示

整个库只用两个结构体来承载所有导数信息
（`ocs2_core/include/ocs2_core/Types.h:78, 145, 234`）：

```cpp
struct ScalarFunctionLinearApproximation  { vector_t dfdx, dfdu; scalar_t f; };
struct ScalarFunctionQuadraticApproximation {
    matrix_t dfdxx, dfdux, dfduu; vector_t dfdx, dfdu; scalar_t f;
};
struct VectorFunctionLinearApproximation  { matrix_t dfdx, dfdu; vector_t f; };
struct VectorFunctionQuadraticApproximation { /* dfdxx 等为 matrix_array_t */ };
```

语义（见 `Types.h:74` 与 `Types.h:141` 的注释）：

$$
f(x,u)\;\approx\;f + \frac{\partial f}{\partial x}^{\!\top}\delta x
   + \frac{\partial f}{\partial u}^{\!\top}\delta u
   + \tfrac12 \delta x^{\!\top}\frac{\partial^2 f}{\partial x^2}\delta x
   + \delta u^{\!\top}\frac{\partial^2 f}{\partial u\partial x}\delta x
   + \tfrac12 \delta u^{\!\top}\frac{\partial^2 f}{\partial u^2}\delta u
$$

注意 `dfdux` 的形状是 $n_u\times n_x$（输入在左、状态在右），这个约定贯穿全库。

这两个结构体重载了 `operator+=`、`operator*=`，所以"多个代价项求和"
就是简单的 `cost += term.getQuadraticApproximation(...)`。

---

## 1.4 一次求解的完整数据流

以 SLQ 为例（`ocs2_ddp/src/GaussNewtonDDP.cpp:980` 的 `runImpl`）：

```
 输入: initTime, initState, finalTime
   │
   ├─ ① 初始化 primal：用上一次的控制器做前向 rollout
   │    （若无控制器则用 Initializer 生成操作点）
   │    GaussNewtonDDP::initializePrimalSolution()          :832
   │
   ├─ ② 初始化 dual + metrics
   │    initializeDualSolutionAndMetrics()
   │
   └─ while (未收敛):
        ├─ ③ 在标称轨迹处做 LQ 近似（多线程）
        │     approximateOptimalControlProblem()             :647
        │       → approximateIntermediateLQ()  各时刻
        │       → approximatePreJumpLQ()       事件时刻
        │       → approximateFinalLQ()         终端
        │
        ├─ ④ 后向传播解 Riccati
        │     solveSequentialRiccatiEquations()
        │       每步先 computeProjectionAndRiccatiModification()  :734
        │         ├─ computeHamiltonianHessian()   (SLQ/iLQR 不同)
        │         ├─ computeProjections()          约束零空间/值域投影
        │         ├─ projectLQ()                   把 LQ 投到零空间
        │         └─ searchStrategy->computeRiccatiModification()
        │       再积分/递推 Riccati 方程
        │
        ├─ ⑤ 由 Riccati 解生成新控制器
        │     calculateController() → calculateControllerWorker()
        │
        ├─ ⑥ 搜索策略：线搜索或 LM，得到新的 primal+dual
        │     takePrimalDualStep()                           :903
        │
        ├─ ⑦ 收敛判据
        │     searchStrategy->checkConvergence()
        │
        └─ ⑧ 未收敛则更新罚系数、交换 nominal/optimized 缓冲
              updateConstraintPenalties()                    :807
```

多重打靶类求解器（SQP/IPM/SLP）的循环骨架几乎一样，
只是把 ③④ 换成"构造 QP 子问题 + 调 QP 求解器"，⑥ 换成滤子线搜索。
对比 `ocs2_sqp/ocs2_sqp/src/SqpSolver.cpp:183` 与上面的流程，你会发现结构高度对称。

---

## 1.5 命名与代码风格约定

读代码前先记住这几条，能省很多时间：

| 约定 | 含义 |
|---|---|
| `xxxPtr_` 后缀下划线 | 成员变量 |
| `dfdx` / `dfdu` / `dfdxx` / `dfdux` / `dfduu` | 一阶/二阶导，`ux` 表示 $\partial^2/\partial u\partial x$ |
| `f` 字段 | 函数值本身（常数项） |
| `Sm` / `Sv` / `s` | 值函数的二次项矩阵 / 一次项向量 / 常数，即 $V=\frac12 x^\top S_m x + S_v^\top x + s$ |
| `Am` / `Bm` | 动力学的 $\partial f/\partial x$、$\partial f/\partial u$ |
| `Qm` / `Qv` / `q` | 代价的 $\partial^2 l/\partial x^2$、$\partial l/\partial x$、$l$ |
| `Rm` / `Rv` / `Pm` | 代价的 $\partial^2 l/\partial u^2$、$\partial l/\partial u$、$\partial^2 l/\partial u\partial x$ |
| `Hm` | 哈密顿量的 Hessian（对输入），SLQ 中 $=R_m$，iLQR 中 $=R_m+B^\top S_m B$ |
| `Km` / `Lv` | 反馈增益 / 前馈项 |
| `projected*` | 已经过约束零空间投影的量 |
| `nominal*` vs `optimized*` | 当前迭代的线性化点 vs 本次迭代产出的新解 |
| `*Stock_` | 每线程一份的资源数组（`nThreads_` 个副本） |
| `ad_scalar_t` / `ad_vector_t` | CppAD 自动微分标量/向量 |

**分隔注释**：源文件里大量出现的

```cpp
/******************************************************************************************************/
```

只是函数之间的视觉分隔，没有语义。

---

## 1.6 构建体系

OCS2 用 **catkin**（ROS1）组织，每个 package 有 `CMakeLists.txt` + `package.xml`。
`ocs2/package.xml` 是元包，列出所有子包依赖。

关键依赖：

| 外部依赖 | 用途 |
|---|---|
| Eigen 3.3+ | 全部线性代数 |
| Boost（system, filesystem, log, numeric/odeint） | 积分器、日志、文件 |
| **CppAD + CppADCodeGen** | 自动微分与 C 代码生成（`ocs2_thirdparty` 内联部分头文件） |
| **HPIPM + BLASFEO** | SQP/IPM 的结构化 QP 求解（`ocs2_sqp/{hpipm,blasfeo}_catkin` 是 catkin 封装） |
| **Pinocchio** | 刚体动力学（`ocs2_pinocchio`） |
| **HPP-FCL** | 碰撞检测（自碰撞约束） |
| ONNX Runtime | MPC-Net 策略推理 |
| PyTorch | MPC-Net 训练（Python 侧） |

CI 配置在 `.github/`，另有 `jenkins-pipeline`。

### 一个重要事实：`ocs2_ocs2` 与 `ocs2_frank_wolfe` 不参与实际构建

`ocs2_ocs2/CMakeLists.txt:47` 只编译 `src/lintTarget.cpp`，而这个文件
（`ocs2_ocs2/src/lintTarget.cpp`）里**所有 include 都被注释掉了**，只剩一个空 `main()`。
原因是这些头文件仍依赖已被删除的 `ocs2_core/Dimensions.h`（模板化的固定维度体系，
在 OCS2 转向动态维度 `vector_t` 后废弃）。

因此：**`GDDP`、`OCS2`、`GSLQPSolver`、`NumGDDP` 这些类现在编译不过**，
属于历史遗留代码。详见 [16-legacy-modules.md](16-legacy-modules.md)。

---

## 1.7 线程模型

OCS2 的实时性很大程度来自并行化。三处并行：

1. **LQ 近似并行**：每个时间点的动力学/代价线性化互相独立。
   用 `std::atomic_int` 做无锁任务分发
   （`ocs2_ddp/src/SLQ.cpp:90` 的 `nextTimeIndex_++` 模式）。
2. **线搜索并行**：不同步长 $\alpha$ 的 rollout 并行执行，
   先找到满足 Armijo 条件的最大 $\alpha$ 就中断其余线程
   （`ocs2_ddp/src/search_strategy/LineSearchStrategy.cpp:188`）。
3. **Riccati 分区并行**：把时间轴切成若干区间并行积分（首次迭代仍串行以保证正确初值，
   见 `GaussNewtonDDP.cpp:528`）。

基础设施是 `ocs2_core/thread_support/ThreadPool.h`，
以及 `Synchronized.h`（带锁的值封装）、`BufferedValue.h`（双缓冲）。

**每线程一份 `OptimalControlProblem` 副本**（`optimalControlProblemStock_`），
因为 `PreComputation` 是有状态的，不能跨线程共享。这解释了为什么
`OptimalControlProblem` 需要深拷贝构造函数
（`ocs2_oc/src/oc_problem/OptimalControlProblem.cpp`）。

---

## 1.8 与其他工具箱的关系

| 工具箱 | 与 OCS2 的关系 |
|---|---|
| **Control Toolbox (CT)** | 同为 ETH 出品，OCS2 的 SLQ 与 CT 的 SLQ 同源；CT 更通用，OCS2 更专注切换系统与机器人 |
| **Crocoddyl** | LAAS-CNRS，DDP 系；离散时间、与 Pinocchio 深度耦合；不处理切换/事件 |
| **acados** | 多重打靶 + HPIPM，与 OCS2 的 SQP 模块思路一致；acados 更偏代码生成与嵌入式 |
| **altro** | 增广拉格朗日 + iLQR；OCS2 的 AL 罚函数设计与其可对照 |
| **OpEn / ProxDDP** | 一阶/近端方法，可与 `ocs2_slp` 的 PIPG 对照 |

OCS2 的独特点：**连续时间 DDP + 事件驱动的切换系统 + 增广拉格朗日约束**这一组合。

---

## 1.9 下一步

- 想理解"问题怎么描述" → [02 ocs2_core](02-ocs2-core.md)
- 想理解"问题怎么组装" → [03 ocs2_oc](03-ocs2-oc.md)
- 想直接看算法数学 → [04 SLQ/iLQR](04-ocs2-ddp.md)
