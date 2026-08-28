# 18 · 代码索引与阅读路线图

本章是**双向索引**：按功能查文件，按文件查功能。
配合 `grep`/IDE 跳转使用。

> **路径约定**：本章的头文件路径是 **include 相对路径**（即 `#include <...>` 里写的那个），
> 例如 `ocs2_ddp/DDP_Settings.h` 的磁盘位置是
> `ocs2_ddp/include/ocs2_ddp/DDP_Settings.h`。
> `.cpp` 路径则是磁盘上的真实相对路径。

---

## 18.1 按功能查文件

### 我想改代价函数

| 需求 | 文件 |
|---|---|
| 写一个新的状态代价 | 继承 `ocs2_core/cost/StateCost.h` |
| 写一个新的状态-输入代价 | 继承 `ocs2_core/cost/StateInputCost.h` |
| 二次跟踪代价 | `ocs2_core/cost/QuadraticStateInputCost.h`（重载 `getStateInputDeviation`） |
| 用自动微分写代价 | `ocs2_core/cost/StateInputCostCppAd.h` |
| 高斯-牛顿（残差）形式 | `ocs2_core/cost/StateInputGaussNewtonCostAd.h` |
| 把代价加进问题 | `problem.costPtr->add("name", std::move(term))` |
| 看代价怎么被累加 | `ocs2_oc/src/approximate_model/LinearQuadraticApproximator.cpp` 的 `approximateCost` |

### 我想加约束

| 约束类型 | 文件 / 做法 |
|---|---|
| 状态-输入等式（会被投影消去） | 继承 `ocs2_core/constraint/StateInputConstraint.h`，加到 `equalityConstraintPtr` |
| 仅状态等式 | `StateConstraint.h` → `stateEqualityConstraintPtr` |
| 不等式（软约束） | 用 `StateInputSoftConstraint` 包一层罚函数 → `softConstraintPtr` |
| 不等式（增广拉格朗日） | `StateInputAugmentedLagrangian` → `inequalityLagrangianPtr` |
| 不等式（硬约束，SQP/IPM） | 加到 `inequalityConstraintPtr` |
| 箱式界（快捷方式） | `ocs2_core/soft_constraint/StateInputSoftBoxConstraint.h` |
| 线性约束 | `ocs2_core/constraint/LinearStateInputConstraint.h` |
| 用自动微分 | `StateInputConstraintCppAd.h` |

### 我想改罚函数

| 罚函数 | 文件 |
|---|---|
| 二次（等式） | `penalties/penalties/QuadraticPenalty.h` |
| 松弛障碍（不等式） | `penalties/penalties/RelaxedBarrierPenalty.h` |
| 平方铰链（不等式） | `penalties/penalties/SquaredHingePenalty.h` |
| 光滑绝对值（等式） | `penalties/penalties/SmoothAbsolutePenalty.h` |
| 箱式 | `penalties/penalties/DoubleSidedPenalty.h` |
| 增广二次 | `penalties/augmented/QuadraticPenalty.h` |
| PHR（增广，不等式） | `penalties/augmented/SlacknessSquaredHingePenalty.h` |
| 光滑 PHR（$C^2$） | `penalties/augmented/ModifiedRelaxedBarrierPenalty.h` |
| 向量化 | `penalties/MultidimensionalPenalty.h` |

### 我想改动力学

| 需求 | 文件 |
|---|---|
| 自动微分写动力学 | 继承 `ocs2_core/dynamics/SystemDynamicsBaseAD.h`，实现 `systemFlowMap` |
| 手写解析导数 | 继承 `SystemDynamicsBase.h`，实现 `computeFlowMap` + `linearApproximation` |
| 有限差分 | `dynamics/SystemDynamicsLinearizer.h` |
| 线性系统 | `dynamics/LinearSystemDynamics.h` |
| 跳变映射（事件） | 重载 `computeJumpMap` / `jumpMapLinearApproximation` |
| 保护面（状态触发） | 重载 `computeGuardSurfaces` |
| 从 URDF | `ocs2_pinocchio_interface/urdf.h` + `ocs2_centroidal_model/` |

### 我想调求解器

| 求解器 | 设置文件 | MPC 包装 |
|---|---|---|
| SLQ / iLQR | `ocs2_ddp/DDP_Settings.h` | `ocs2_ddp/GaussNewtonDDP_MPC.h` |
| SQP | `ocs2_sqp/SqpSettings.h` | `ocs2_sqp/SqpMpc.h` |
| IPM | `ocs2_ipm/IpmSettings.h` | `ocs2_ipm/IpmMpc.h` |
| SLP | `ocs2_slp/SlpSettings.h` + `pipg/PipgSettings.h` | `ocs2_slp/SlpMpc.h` |
| MPC 通用 | `ocs2_mpc/MPC_Settings.h` | — |
| Rollout | `ocs2_oc/rollout/RolloutSettings.h` | — |
| 线搜索 / LM | `ocs2_ddp/search_strategy/StrategySettings.h` | — |
| HPIPM | `hpipm_catkin/HpipmInterfaceSettings.h` | — |

### 我想理解某个算法步骤

| 步骤 | 文件:行 |
|---|---|
| DDP 主循环 | `ocs2_ddp/src/GaussNewtonDDP.cpp:980` |
| LQ 近似（DDP） | `ocs2_ddp/src/GaussNewtonDDP.cpp:647` |
| 约束投影 | `ocs2_ddp/src/GaussNewtonDDP.cpp:755` + `ocs2_core/src/misc/LinearAlgebra.cpp:129` |
| 投影 LQ | `ocs2_ddp/src/DDP_HelperFunctions.cpp:143` |
| 连续 Riccati | `ocs2_ddp/src/riccati_equations/ContinuousTimeRiccatiEquations.cpp:174` |
| 离散 Riccati | `ocs2_ddp/src/riccati_equations/DiscreteTimeRiccatiEquations.cpp:65` |
| 控制器生成（SLQ） | `ocs2_ddp/src/SLQ.cpp:127` |
| 控制器生成（iLQR） | `ocs2_ddp/src/ILQR.cpp:162` |
| 并行线搜索 | `ocs2_ddp/src/search_strategy/LineSearchStrategy.cpp:188` |
| LM 信赖域 | `ocs2_ddp/src/search_strategy/LevenbergMarquardtStrategy.cpp:122` |
| Hessian 修正 | `ocs2_ddp/src/HessianCorrection.cpp:53` |
| merit 函数 | `ocs2_ddp/src/GaussNewtonDDP.cpp:500` |
| 罚系数更新 | `ocs2_ddp/src/GaussNewtonDDP.cpp:807` |
| SQP 主循环 | `ocs2_sqp/ocs2_sqp/src/SqpSolver.cpp:183` |
| QP 构造 | `ocs2_sqp/ocs2_sqp/src/SqpSolver.cpp:336` |
| 滤子线搜索 | `ocs2_oc/src/search_strategy/FilterLinesearch.cpp:34` |
| 多重打靶转写 | `ocs2_oc/src/multiple_shooting/Transcription.cpp:40` |
| 转写投影 | `ocs2_oc/src/multiple_shooting/Transcription.cpp:96` |
| 投影乘子恢复 | `ocs2_oc/src/multiple_shooting/ProjectionMultiplierCoefficients.cpp:35` |
| IPM 凝聚 | `ocs2_ipm/src/IpmHelpers.cpp:37` |
| 分数到边界 | `ocs2_ipm/src/IpmHelpers.cpp:103` |
| 障碍参数更新 | `ocs2_ipm/src/IpmSolver.cpp:898` |
| PIPG 迭代 | `ocs2_slp/src/pipg/PipgSolver.cpp:49` |
| PIPG 步长 | `ocs2_slp/include/ocs2_slp/pipg/PipgBounds.h:66` |
| Ruiz 预条件 | `ocs2_oc/src/precondition/Ruzi.cpp` |
| 灵敏度离散化 | `ocs2_core/src/integration/SensitivityIntegratorImpl.cpp:130` |
| 状态触发 rollout | `ocs2_oc/src/rollout/StateTriggeredRollout.cpp` |
| 轨迹展开 | `ocs2_oc/src/trajectory_adjustment/TrajectorySpreading.cpp` |
| MRT 策略交换 | `ocs2_mpc/src/MRT_BASE.cpp:156` |

### 我想加一个新机器人

按这个顺序：

1. **准备 URDF**，确定浮动基类型
2. **定义状态/输入**：新建 `definitions.h`
3. **写动力学**：继承 `SystemDynamicsBaseAD`，或用 `ocs2_centroidal_model`
4. **写 `RobotInterface` 子类**：装配 `OptimalControlProblem`
   - 参考 `ocs2_cartpole/src/CartPoleInterface.cpp`（最简）
   - 或 `ocs2_legged_robot/src/LeggedRobotInterface.cpp`（最全）
5. **写配置** `config/mpc/task.info`
6. **写 ROS 节点**：参考 `ocs2_cartpole_ros/src/`
7. **写可视化**：继承 `DummyObserver`

模板：直接复制 `ocs2_cartpole` + `ocs2_cartpole_ros` 改名。

---

## 18.2 按文件查功能（核心文件速查）

### `ocs2_core`

| 文件 | 一句话 |
|---|---|
| `Types.h` | 全库类型别名 + 四个近似结构体 |
| `PreComputation.h` | 跨项共享计算的钩子 |
| `ComputationRequest.h` | 位掩码请求集合 |
| `NumericTraits.h` | 容差常数 |
| `dynamics/SystemDynamicsBase.h` | 动力学接口（流形 + 跳变 + 保护面） |
| `dynamics/SystemDynamicsBaseAD.h` | CppAD 版动力学 |
| `cost/StateInputCost.h` | 代价接口 |
| `cost/StateInputGaussNewtonCostAd.h` | 残差形式代价（PSD Hessian） |
| `constraint/StateInputConstraint.h` | 约束接口 |
| `soft_constraint/StateInputSoftConstraint.h` | 约束 + 罚 → 代价 |
| `augmented_lagrangian/*` | 增广拉格朗日接口 |
| `penalties/penalties/*` | 普通罚函数 |
| `penalties/augmented/*` | 增广罚函数 |
| `integration/Integrator.h` | odeint 封装（rollout 用） |
| `integration/SensitivityIntegrator.h` | 离散化 + 灵敏度（多重打靶用） |
| `integration/TrapezoidalIntegration.h` | 梯形积分（代价评估用） |
| `control/LinearController.h` | $u=Kx+b+\alpha\Delta b$ |
| `initialization/Initializer.h` | 无控制器时的轨迹生成 |
| `reference/ModeSchedule.h` | 事件时刻 + 模态序列 |
| `reference/TargetTrajectories.h` | 跟踪目标 |
| `model_data/ModelData.h` | 单时刻的完整 LQ 数据 |
| `model_data/Metrics.h` | 单时刻的评估值 |
| `misc/LinearAlgebra.h` | PSD 修正、约束投影、UUT 分解 |
| `misc/LinearInterpolation.h` | 共享查找的插值 |
| `misc/Collection.h` | 具名容器基类 |
| `misc/Benchmark.h` | 计时器 |
| `misc/LoadData.h` | `.info` 配置加载 |
| `automatic_differentiation/CppAdInterface.h` | AD → C 代码 → 编译 → 加载 |
| `thread_support/ThreadPool.h` | 线程池 |
| `loopshaping/LoopshapingDefinition.h` | 滤波器 + 状态/输入映射 |

### `ocs2_oc`

| 文件 | 一句话 |
|---|---|
| `oc_problem/OptimalControlProblem.h` | **问题的完整载体**（24 个字段） |
| `oc_problem/OcpSize.h` | 各阶段维度 |
| `oc_problem/OcpToKkt.h` | 稀疏 KKT 矩阵组装 |
| `oc_solver/SolverBase.h` | **求解器统一接口** |
| `oc_data/PrimalSolution.h` | 原变量解 |
| `oc_data/DualSolution.h` | 对偶变量解 |
| `oc_data/PerformanceIndex.h` | 标量性能汇总 |
| `oc_data/TimeDiscretization.h` | 事件对齐的时间网格 |
| `approximate_model/LinearQuadraticApproximator.h` | 问题 → LQ 数据 |
| `approximate_model/ChangeOfInputVariables.h` | 输入替换下的近似变换 |
| `multiple_shooting/Transcription.h` | 问题 → QP 数据 |
| `multiple_shooting/ProjectionMultiplierCoefficients.h` | 投影乘子恢复 |
| `multiple_shooting/Helpers.h` | 步长应用、remap、组装解 |
| `rollout/RolloutBase.h` | 前向仿真接口 |
| `rollout/StateTriggeredRollout.h` | 保护面触发的 rollout |
| `rollout/RootFinder.h` | 试位法家族 |
| `search_strategy/FilterLinesearch.h` | 滤子线搜索（三求解器共用） |
| `synchronized_module/ReferenceManager.h` | 线程安全的参考管理 |
| `synchronized_module/SolverObserver.h` | 求解后的指标提取 |
| `trajectory_adjustment/TrajectorySpreading.h` | mode schedule 变化时的轨迹适配 |
| `precondition/Ruzi.h` | Ruiz 均衡 |

### `ocs2_ddp`

| 文件 | 一句话 |
|---|---|
| `GaussNewtonDDP.h/.cpp` | **DDP 的 90% 逻辑**（主循环、投影、数据管理） |
| `SLQ.h/.cpp` | 连续时间特化（4 个虚函数） |
| `ILQR.h/.cpp` | 离散时间特化 |
| `riccati_equations/ContinuousTimeRiccatiEquations.h` | 微分 Riccati |
| `riccati_equations/DiscreteTimeRiccatiEquations.h` | 差分 Riccati |
| `riccati_equations/RiccatiModification.h` | 投影矩阵 + 修正项的容器 |
| `riccati_equations/RiccatiTransversalityConditions.h` | 事件处的值函数跳变 |
| `search_strategy/LineSearchStrategy.h` | 并行 Armijo 回溯 |
| `search_strategy/LevenbergMarquardtStrategy.h` | 信赖域 |
| `HessianCorrection.h` | 四种 PSD 修正 |
| `DDP_Data.h` | `PrimalDataContainer` / `DualDataContainer` |
| `DDP_HelperFunctions.h` | `projectLQ`、`incrementController` 等 |
| `ContinuousTimeLqr.h` | 代数 Riccati（无限时域） |
| `GaussNewtonDDP_MPC.h` | MPC 包装 |
| `unsupported/` | 遗留：哈密顿法 SLQ、BVP 求解器 |

---

## 18.3 术语对照表

| 英文 | 中文 | 说明 |
|---|---|---|
| flow map | 流形 / 连续动力学 | $\dot x=f(t,x,u)$ |
| jump map / reset map | 跳变映射 | $x^+=j(x^-)$ |
| guard surface | 保护面 | 零穿越触发事件 |
| mode schedule | 模态时序 | 事件时刻 + 模态序列 |
| pre-jump / post-jump | 事件前 / 事件后 | $t_i^-$ / $t_i^+$ |
| rollout | 前向仿真 | 用控制器积分动力学 |
| primal / dual solution | 原 / 对偶解 | 状态输入 / 乘子 |
| merit function | 优点函数 | 代价 + 约束罚 |
| line search | 线搜索 | 沿方向找步长 |
| filter line search | 滤子线搜索 | 代价与可行性的双目标接受准则 |
| trust region | 信赖域 | LM 的 $\mu$ 控制 |
| null-space projection | 零空间投影 | 消去等式约束 |
| range-space projector | 值域投影 | 加权伪逆 $D^\dagger$ |
| condensing | 凝聚 | 消去松弛/对偶变量 |
| fraction-to-boundary | 分数到边界 | 保持内点的步长规则 |
| barrier parameter | 障碍参数 | $\mu$ |
| central path | 中心路径 | $s\lambda=\mu$ 的解轨迹 |
| warm start | 热启动 | 用上次解初始化 |
| multiple shooting | 多重打靶 | 状态也是决策变量 |
| single shooting | 单次打靶 | 只有输入是决策变量 |
| defect / gap | 缺口 | $x_{k+1}-F(x_k,u_k)$ |
| costate | 协态 | 动力学约束的乘子 |
| transcription | 转写 | 连续问题 → 有限维 NLP |
| receding horizon | 滚动时域 | MPC 的核心机制 |
| real-time iteration (RTI) | 实时迭代 | 每周期只迭代 1~3 次 |
| centroidal dynamics | 质心动力学 | 牛顿-欧拉的 6 个方程 |
| centroidal momentum matrix | 质心动量矩阵 | $A(q)$，$\mathbf{h}=A\dot q$ |
| SRBD | 单刚体动力学 | 忽略腿质量的近似 |
| loopshaping | 频域整形 | 用滤波器塑造惩罚的频率特性 |
| trajectory spreading | 轨迹展开 | 时间轴重映射 |
| swing / stance | 摆动 / 支撑 | 腿的两个相位 |
| gait | 步态 | 周期性接触序列 |

---

## 18.4 常用 grep 模式

```bash
# 找某个类的定义
grep -rn "class GaussNewtonDDP" --include=*.h .

# 找某个函数的所有调用点
grep -rn "computeConstraintProjection" --include=*.cpp --include=*.h .

# 找所有 TODO
grep -rn "TODO\|FIXME\|HACK" --include=*.h --include=*.cpp .

# 找某个设置项在哪被读取
grep -rn "constraintPenaltyIncreaseRate_" .

# 列出某个包的所有公开接口
grep -n "virtual\|public:" ocs2_oc/include/ocs2_oc/oc_solver/SolverBase.h

# 找某个数学量的所有出现（如 Riccati 的 Sm）
grep -rn "\bSm\b" --include=*.cpp ocs2_ddp/

# 查一个 .info 配置项对应的 C++ 变量
grep -rn "loadPtreeValue.*minRelCost" .
```

---

## 18.5 从零到跑通的检查清单

```bash
# 1. 依赖
sudo apt install ros-noetic-desktop-full liburdfdom-dev liboctomap-dev libassimp-dev
# Pinocchio, HPP-FCL, ONNX Runtime 按官方文档装

# 2. 建 catkin workspace
mkdir -p ~/ocs2_ws/src && cd ~/ocs2_ws/src
git clone <this repo>
git clone https://github.com/leggedrobotics/ocs2_robotic_assets.git
git clone --recurse-submodules https://github.com/leggedrobotics/hpp-fcl.git
git clone https://github.com/leggedrobotics/pinocchio.git

# 3. 编译（先编简单的验证环境）
cd ~/ocs2_ws
catkin build ocs2_cartpole_ros

# 4. 跑
source devel/setup.bash
roslaunch ocs2_cartpole_ros cartpole.launch

# 5. 再编四足
catkin build ocs2_legged_robot_ros
roslaunch ocs2_legged_robot_ros legged_robot_sqp.launch
```

**首次运行会卡 1~5 分钟**（CppAD 编译），属正常。
详见 `ocs2_doc/docs/installation.rst`。

---

## 18.6 本文档章节 → 源码目录对照

| 章节 | 主要覆盖的目录 |
|---|---|
| [01 总览](01-overview-architecture.md) | 全局 |
| [02 core](02-ocs2-core.md) | `ocs2_core/` |
| [03 oc](03-ocs2-oc.md) | `ocs2_oc/` |
| [04 DDP](04-ocs2-ddp.md) | `ocs2_ddp/` |
| [05 SQP](05-ocs2-sqp.md) | `ocs2_sqp/` |
| [06 IPM](06-ocs2-ipm.md) | `ocs2_ipm/` |
| [07 SLP](07-ocs2-slp.md) | `ocs2_slp/`, `ocs2_oc/precondition/` |
| [08 约束](08-constraints-penalties.md) | `ocs2_core/{penalties,soft_constraint,augmented_lagrangian}/` |
| [09 Loopshaping](09-loopshaping.md) | `ocs2_core/loopshaping/`, `ocs2_oc/oc_problem/Loopshaping*` |
| [10 积分/Rollout](10-integration-rollout.md) | `ocs2_core/integration/`, `ocs2_oc/{rollout,trajectory_adjustment}/` |
| [11 MPC/ROS](11-mpc-ros.md) | `ocs2_mpc/`, `ocs2_ros_interfaces/`, `ocs2_python_interface/` |
| [12 机器人建模](12-robot-models.md) | `ocs2_pinocchio/`, `ocs2_robotic_tools/`, `ocs2_perceptive/` |
| [13 示例](13-robot-examples.md) | `ocs2_robotic_examples/`, `ocs2_raisim/` |
| [14 MPC-Net](14-mpcnet.md) | `ocs2_mpcnet/` |
| [15 工具](15-numerics-tools.md) | `ocs2_core/{misc,thread_support,automatic_differentiation}/`, `ocs2_test_tools/` |
| [16 遗留](16-legacy-modules.md) | `ocs2_ocs2/`, `ocs2_frank_wolfe/`, `ocs2_ddp/unsupported/` |
| [17 文献](17-references.md) | `ocs2_doc/docs/refs.bib` |

---

## 18.7 未被本文档单独展开的部分

为完整起见，列出本文档只做了概述而未逐文件展开的部分：

| 目录 | 原因 |
|---|---|
| `ocs2_msgs/` | 纯 ROS 消息定义，无逻辑（`mpc_state.msg`、`mpc_flattened_controller.msg` 等） |
| `ocs2_thirdparty/` | 内联的第三方头文件，不是 OCS2 的代码 |
| `ocs2_robotic_examples/ocs2_ballbot/generated/` | RobCoGen 自动生成的 20+ 文件（[13 章](13-robot-examples.md) §13.4 已说明其来源与作用） |
| 各包的 `test/` | 单元测试；[05 章](05-ocs2-sqp.md) §5.9、[09 章](09-loopshaping.md) §9.8 指出了其中最有学习价值的几个 |
| 各包的 `config/*.info` | 配置文件；参数含义在各章的 Settings 表格中给出 |
| `.github/`、`jenkins-pipeline` | CI 配置 |
| `ocs2_doc/tools/` | Sphinx/Doxygen 构建脚本 |

---

## 18.8 文档使用建议

- **查东西**：从本章的索引表出发
- **理解算法**：从对应章节的推导出发，再对照代码
- **改代码**：先看 §18.1 的"我想…"表格定位文件
- **调参**：各章末尾的 Settings 表 + [15 章](15-numerics-tools.md) §15.12 的调试清单
- **写论文**：[17 章](17-references.md) 的文献 + 引用建议

有疑问时，**代码是唯一的事实来源**。
本文档给出的行号基于撰写时的仓库状态，随代码演进可能偏移——
用函数名/类名 grep 比用行号更可靠。
