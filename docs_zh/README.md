# OCS2 源码完整导读（中文）

> 本目录是对 **OCS2（Optimal Control for Switched Systems）** 工具箱的一次系统性源码梳理。
> 目标读者是希望**逐行读懂这个仓库**的人：既讲清每个 package 的职责、类的继承关系与数据流，
> 也把算法背后的**数学推导完整写出来**（而不是只给结论），并给出**可追溯的参考文献**。
>
> 所有公式都与仓库中的具体代码位置一一对应，形如 `ocs2_ddp/src/SLQ.cpp:206`，可直接跳转阅读。

---

## 0. 这套文档怎么读

OCS2 的代码量约 **31 万行 C++**（不含 thirdparty），共 1600+ 源文件。
直接从 `ocs2_core` 一个文件一个文件读会迷路，因为它的核心抽象是"问题定义"和"求解器"的**正交分解**：

```
        用户描述"问题是什么"                  求解器决定"怎么解"
   ┌──────────────────────────────┐   ┌────────────────────────────────┐
   │  OptimalControlProblem       │   │  SLQ / iLQR   (ocs2_ddp)       │
   │   ├─ dynamicsPtr             │   │  SQP          (ocs2_sqp)       │
   │   ├─ costPtr / softConstraint│ ─▶│  IPM          (ocs2_ipm)       │
   │   ├─ equalityConstraintPtr   │   │  SLP/PIPG     (ocs2_slp)       │
   │   ├─ inequalityConstraintPtr │   │                                │
   │   └─ equalityLagrangianPtr   │   │  它们共享 ocs2_oc 里的          │
   │  (ocs2_oc/oc_problem)        │   │  LQ 近似 / rollout / 线搜索     │
   └──────────────────────────────┘   └────────────────────────────────┘
```

**建议阅读顺序**（也是本目录的编号顺序）：

| 序号 | 文档 | 内容 |
|---|---|---|
| 01 | [总览与整体架构](01-overview-architecture.md) | package 依赖图、数据流、命名约定、构建体系 |
| 02 | [ocs2_core：类型系统与问题建模原语](02-ocs2-core.md) | `Types.h`、代价/约束/动力学基类、PreComputation、集合类 |
| 03 | [ocs2_oc：最优控制问题与多重打靶转写](03-ocs2-oc.md) | `OptimalControlProblem`、LQ 近似、时间离散化、Transcription、投影 |
| 04 | [ocs2_ddp：SLQ 与 iLQR 的完整推导](04-ocs2-ddp.md) | ⭐ 连续/离散 Riccati 全推导、约束投影、Hessian 修正、线搜索/LM |
| 05 | [ocs2_sqp：多重打靶 SQP 与 HPIPM](05-ocs2-sqp.md) | ⭐ SQP 的 KKT 系统、滤子线搜索、值函数提取 |
| 06 | [ocs2_ipm：原对偶内点法](06-ocs2-ipm.md) | ⭐ 障碍函数、凝聚（condensing）、分数到边界规则 |
| 07 | [ocs2_slp：PIPG 一阶求解器](07-ocs2-slp.md) | ⭐ 比例-积分投影梯度法、Ruiz 均衡预条件、步长界估计 |
| 08 | [约束处理：罚函数与增广拉格朗日](08-constraints-penalties.md) | ⭐ 松弛障碍、平方铰链、PHR、光滑 PHR 的推导与对偶更新 |
| 09 | [Loopshaping：频域整形与系统增广](09-loopshaping.md) | ⭐ output/eliminate 两种模式的状态增广与代价变换 |
| 10 | [积分、Rollout 与切换系统](10-integration-rollout.md) | ⭐ 灵敏度离散化、状态触发 rollout、根搜索、轨迹展开 |
| 11 | [ocs2_mpc 与 ROS 接口](11-mpc-ros.md) | MPC/MRT 双线程架构、话题协议、可视化 |
| 12 | [机器人建模：Pinocchio 与质心动力学](12-robot-models.md) | ⭐ 质心动力学推导、SRBD 近似、自碰撞、末端运动学 |
| 13 | [机器人示例逐个拆解](13-robot-examples.md) | ⭐ ballbot / cartpole / quadrotor / 移动机械臂 / 四足 |
| 14 | [MPC-Net：从 MPC 蒸馏策略](14-mpcnet.md) | ⭐ 哈密顿量损失、混合专家网络、数据生成 |
| 15 | [数值与工程工具](15-numerics-tools.md) | CppAD 自动微分、线性代数、线程池、插值、日志 |
| 16 | [遗留与实验性模块](16-legacy-modules.md) | `ocs2_ocs2`(GDDP/切换时刻优化)、`ocs2_frank_wolfe`、RaiSim |
| 17 | [参考文献与延伸阅读](17-references.md) | 论文、教材、外部库的完整索引 |
| 18 | [代码索引与阅读路线图](18-code-index.md) | 按功能查文件、按文件查功能的双向索引 |

标 ⭐ 的章节含有完整的数学推导。

---

## 1. 一句话概括 OCS2

OCS2 求解如下**带切换的、约束的、连续时间最优控制问题**，并把它包成实时 MPC：

$$
\begin{aligned}
\min_{u(\cdot)}\quad & \phi\big(x(t_f)\big)+\int_{t_0}^{t_f} l\big(t,x(t),u(t)\big)\,\mathrm{d}t
   +\sum_{i\in\mathcal{E}} \phi_i\big(x(t_i^-)\big)\\
\text{s.t.}\quad
& \dot x(t)=f\big(t,x(t),u(t)\big), && t\notin\mathcal{E}\\
& x(t_i^+)=j\big(t_i,x(t_i^-)\big), && t_i\in\mathcal{E} \quad(\text{跳变/事件})\\
& x(t_0)=x_0\\
& g_1(t,x,u)=0,\qquad g_2(t,x)=0 && (\text{等式约束})\\
& h_1(t,x,u)\ge 0,\qquad h_2(t,x)\ge 0 && (\text{不等式约束})
\end{aligned}
$$

其中 $\mathcal{E}$ 是**事件时刻集合**（例如四足机器人落足/抬足的瞬间），
由 `ModeSchedule` 描述。"Switched Systems"正是指这类**分段动力学 + 跳变映射**的系统。

代码里这个问题的载体是 `ocs2::OptimalControlProblem`
（`ocs2_oc/include/ocs2_oc/oc_problem/OptimalControlProblem.h:48`）——
一个由 20 多个 `unique_ptr` 组成的结构体，每个指针对应上式的一项。

---

## 2. 仓库地图（按代码量与重要性排序）

| Package | 头文件+源文件 | 行数 | 职责 |
|---|---:|---:|---|
| `ocs2_core` | 273 | 29.6k | 类型、代价/约束/动力学抽象、罚函数、积分器、loopshaping、线程支持 |
| `ocs2_oc` | 91 | 12.6k | 最优控制问题定义、LQ 近似、rollout、多重打靶转写、参考管理 |
| `ocs2_ddp` | 55 | 11.6k | SLQ（连续时间 DDP）、iLQR（离散时间 DDP）、Riccati、搜索策略 |
| `ocs2_ocs2` | 22 | 4.4k | **（遗留）** GDDP：对切换时刻求梯度的双层优化 |
| `ocs2_ipm` | 18 | 3.7k | 多重打靶 + 原对偶内点法 |
| `ocs2_sqp` | 18 | 3.2k | 多重打靶 SQP + HPIPM QP 求解器绑定 |
| `ocs2_slp` | 20 | 2.7k | 顺序线性规划 + PIPG 一阶求解器 |
| `ocs2_robotic_tools` | 12 | 1.5k | 旋转表示、角速度映射、机器人接口基类 |
| `ocs2_test_tools` | 13 | 1.4k | `ocs2_qp_solver`：稠密 QP 参考实现（用于交叉验证） |
| `ocs2_mpc` | 15 | 1.4k | MPC 主循环、MRT（Model Reference Tracking）接口 |
| `ocs2_perceptive` | 14 | 1.4k | 距离场、三线性插值、末端到障碍距离约束 |
| `ocs2_frank_wolfe` | 10 | 1.3k | **（遗留）** Frank-Wolfe 线性化约束下的梯度下降 |
| `ocs2_python_interface` | 4 | 0.7k | pybind11 绑定 |
| `ocs2_pinocchio` | — | — | Pinocchio 接口、质心动力学、自碰撞、球近似 |
| `ocs2_ros_interfaces` | — | — | ROS 话题/服务、可视化、交互标记 |
| `ocs2_robotic_examples` | — | — | ballbot / cartpole / 双积分器 / 四旋翼 / 移动机械臂 / 四足 |
| `ocs2_mpcnet` | — | — | MPC-Net：从 MPC 蒸馏神经网络策略（C++ + PyTorch） |
| `ocs2_raisim` | — | — | RaiSim 物理引擎桥接 |
| `ocs2_thirdparty` | — | — | 内联的第三方头文件（Boost.Numeric.Odeint 等） |
| `ocs2_msgs` | — | — | ROS 消息定义 |
| `ocs2_doc` | — | — | Sphinx + Doxygen 官方文档源 |

---

## 3. 关于本文档的写作原则

1. **公式先于代码**：先写出连续/离散数学形式，再指出它落在哪一行代码上。
2. **推导不跳步**：Riccati、投影、内点法凝聚这些地方，中间步骤全部写出。
3. **符号与代码变量对齐**：例如 $S_m$ ↔ `Sm`、$\tilde u$ ↔ `projected input`，见每章开头的符号表。
4. **区分"在用"与"遗留"**：`ocs2_ocs2`、`ocs2_frank_wolfe` 已不参与实际构建，单独成章说明。
5. **引用可查**：每个算法都给出原始论文，见 [17-references.md](17-references.md)。

---

## 4. 快速上手路径

如果你只有一个下午：

1. 读 [01 总览](01-overview-architecture.md) 的架构图（10 分钟）
2. 读 [03 ocs2_oc](03-ocs2-oc.md) 的 `OptimalControlProblem` 一节（20 分钟）
3. 挑一个求解器精读：想理解 MPC 主流做法就读 [05 SQP](05-ocs2-sqp.md)；
   想理解 OCS2 的招牌算法就读 [04 SLQ](04-ocs2-ddp.md)（60–90 分钟）
4. 跑一遍 `ocs2_cartpole`，对照 [13 机器人示例](13-robot-examples.md)（30 分钟）

如果你要做机器人落地，[12 机器人建模](12-robot-models.md) 和
[13 示例](13-robot-examples.md) 是必读，因为 90% 的工程量在"怎么把机器人写成 OCP"。
