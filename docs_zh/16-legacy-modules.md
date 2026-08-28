# 16 · 遗留与实验性模块

本章记录仓库里**当前不参与实际构建**或处于实验状态的部分。
读源码时遇到它们不必困惑——它们是 OCS2 演化过程的化石层。

---

## 16.1 ⚠️ `ocs2_ocs2`：GDDP 与切换时刻优化

### 16.1.1 构建状态

**这个包目前只编译一个空的 `main()`。**

`ocs2_ocs2/CMakeLists.txt:47`：

```cmake
add_executable(${PROJECT_NAME}_lintTarget
  src/lintTarget.cpp
)
```

而 `src/lintTarget.cpp` 的全部内容是：

```cpp
// // sensitivity_equations
// #include <ocs2_ocs2/sensitivity_equations/BvpSensitivityEquations.h>
// ... （所有 include 都被注释掉）
// dummy target for clang toolchain
int main() { return 0; }
```

**原因**：这些头文件仍在 include 已被删除的 `ocs2_core/Dimensions.h`。

```cpp
// ocs2_ocs2/include/ocs2_ocs2/sensitivity_equations/BvpSensitivityEquations.h:37
#include <ocs2_core/Dimensions.h>     // ← 这个文件在仓库里不存在
```

`Dimensions<STATE_DIM, INPUT_DIM>` 是早期 OCS2 的**编译期固定维度**体系
（`state_vector_t = Eigen::Matrix<scalar_t, STATE_DIM, 1>`）。
后来 OCS2 全面改用**运行时动态维度**（`vector_t = Eigen::VectorXd`），
`Dimensions.h` 被删除，但 `ocs2_ocs2` 与 `ocs2_frank_wolfe` 未跟进迁移。

**验证**：`find . -name "Dimensions.h"` 返回空。

### 16.1.2 ⭐ 它本来做什么：双层优化

尽管代码不可用，**它的数学思想是 OCS2 的名字来源**，值得理解。

**问题**：前面所有章节都假设 `ModeSchedule`（事件时刻 $\tau=[\tau_1,\dots,\tau_K]$）
**由外部给定**。但事件时刻本身也应该被优化——
四足机器人的最优落足时刻是什么？

这是一个**双层优化（bilevel optimization）**：

$$
\boxed{
\begin{aligned}
\text{上层：}\quad &\min_{\boldsymbol{\tau}}\ \ J^\star(\boldsymbol{\tau})\qquad
\text{s.t.}\ \ \tau_1\le\tau_2\le\cdots\le\tau_K,\ \ \tau\in[t_0,t_f]\\[6pt]
\text{下层：}\quad &J^\star(\boldsymbol{\tau})=\min_{u(\cdot)}\ J\big(u;\boldsymbol{\tau}\big)
\quad\text{（这就是前面各章的 OCP）}
\end{aligned}}
$$

**难点**：上层需要 $\dfrac{\partial J^\star}{\partial\tau_i}$，
但 $J^\star$ 是一个优化问题的**最优值**，不能直接求导。

### 16.1.3 GDDP：梯度的计算

**GDDP** = **G**radient of **DDP**。核心是**灵敏度分析**。

#### 三类灵敏度方程

| 文件 | 计算什么 |
|---|---|
| `RolloutSensitivityEquations.h` | $\dfrac{\partial x(t)}{\partial\tau_i}$（状态轨迹对事件时刻的灵敏度） |
| `SensitivitySequentialRiccatiEquations.h` | $\dfrac{\partial S_m}{\partial\tau_i},\dfrac{\partial S_v}{\partial\tau_i},\dfrac{\partial s}{\partial\tau_i}$（值函数的灵敏度） |
| `BvpSensitivityEquations.h` | 两点边值问题形式的灵敏度（约束情形） |
| `BvpSensitivityErrorEquations.h` | 上面的误差补偿项 |

#### 灵敏度方程的推导思路

**关键工具：时间尺度变换**。把事件时刻的变化转化为动力学的"时间缩放"。

设第 $i$ 个子系统在 $[\tau_{i-1},\tau_i]$ 上生效。引入归一化时间 $\sigma\in[0,1]$：

$$
t = \tau_{i-1}+\sigma(\tau_i-\tau_{i-1})
$$

于是

$$
\frac{\mathrm{d}x}{\mathrm{d}\sigma}=(\tau_i-\tau_{i-1})\,f_i(t,x,u)
$$

现在事件时刻变成了**动力学的乘性系数**，对它求导就是普通的参数灵敏度：

$$
\frac{\partial}{\partial\tau_i}\left(\frac{\mathrm{d}x}{\mathrm{d}\sigma}\right)
=f_i + (\tau_i-\tau_{i-1})\,\frac{\partial f_i}{\partial x}\frac{\partial x}{\partial\tau_i}
$$

`computeEquivalentSystemMultiplier`（`GDDP.h:244`）算的就是这个"乘子"
$(\tau_i-\tau_{i-1})$ 的导数——对不同的 $(i,j)$ 组合有不同的值：

- 若 $j = i$（事件本身）：乘子为 $+1$
- 若 $j = i-1$（前一个事件）：乘子为 $-1$
- 其他：0

#### 代价的梯度

有了 $\dfrac{\partial x}{\partial\tau_i}$、$\dfrac{\partial S}{\partial\tau_i}$，
用链式法则组装：

$$
\frac{\partial J^\star}{\partial\tau_i}
=\underbrace{\frac{\partial J}{\partial x}\frac{\partial x}{\partial\tau_i}}_{\text{通过轨迹}}
+\underbrace{\frac{\partial J}{\partial \tau_i}\Big|_{x}}_{\text{显式依赖}}
$$

**包络定理的应用**：由于下层已最优（$\partial J/\partial u=0$），
$u$ 随 $\tau$ 的变化对 $J^\star$ 的一阶影响为零，
所以**不需要计算 $\partial u/\partial\tau$** ✓
这大幅简化了计算。

`GDDP.h` 里的 `nablaq`、`nablaQv`、`nablaRv`、`nablasHeuristics`、`nablaSvHeuristics`
（`:269-286`）就是各阶变分对 $\tau$ 的灵敏度。

`getCostFuntionDerivative`（`:172`，注意源码里的拼写笔误 "Funtion"）
返回最终的 $\nabla_\tau J^\star$。

### 16.1.4 上层求解：Frank-Wolfe

有了梯度，上层用**投影梯度**或 **Frank-Wolfe** 求解。
`OCS2.h` 类把两层串起来：

```cpp
template <size_t STATE_DIM, size_t INPUT_DIM>
class OCS2 {
  OCS2(rolloutPtr, systemDerivativesPtr, systemConstraintsPtr, costFunctionPtr,
       operatingTrajectoriesPtr, const SLQ_Settings&, referenceManagerPtr,
       heuristicsFunctionPtr, const GDDP_Settings&, const NLP_Settings&);
};
```

内部：
- `UpperLevelCost`（`upper_level_op/UpperLevelCost.h:47`）：
  实现 `NLP_Cost` 接口，`getCost(τ)` 内部跑一次 SLQ，`getCostDerivative(τ)` 跑一次 GDDP
- `UpperLevelConstraints`（`UpperLevelConstraints.h:39`）：
  实现 `NLP_Constraints`，编码 $\tau_1\le\cdots\le\tau_K$ 与 $\tau\in[t_0,t_f]$
- `GradientDescent`（来自 `ocs2_frank_wolfe`）：驱动上层迭代

其余相关类：
- `GSLQPSolver.h`：更早的 GSLQ（Gradient of SLQ）实现
- `NumGDDP.h`：**数值**版 GDDP（有限差分求 $\nabla_\tau J^\star$），用于验证解析 GDDP
- `FrankWolfeGDDP.h`：GDDP + Frank-Wolfe 的组合

### 16.1.5 为什么被弃用

1. **计算成本高**：每次上层迭代要跑一次完整的 SLQ + 一次 GDDP，
   实时性完全无法满足
2. **实践中不必要**：腿足机器人的接触时序由**步态规划器**给出更实用
   （周期步态、基于地形的落足点规划），效果好且快得多
3. **与动态维度重构冲突**：迁移成本高，收益低
4. **有替代方案**：现代做法是把切换时刻作为**决策变量嵌入多重打靶**
   （如 Winkler et al. 的 TOWR），或用相位参数化（Winkler et al. 2018）

**理论价值仍在**：如果你要研究切换时刻优化，
`sensitivity_equations/` 的推导仍然是很好的参考。

### 16.1.6 参考文献

- **Farshidian, Kamgarpour, Pardo, Buchli (2017)**,
  *Sequential Linear Quadratic Optimal Control for Nonlinear Switched Systems*, IFAC
  —— 包含 SLQ 与切换时刻灵敏度
- **Xu & Antsaklis (2004)**, *Optimal Control of Switched Systems Based on
  Parameterization of the Switching Instants*, IEEE TAC 49(1)
  —— **时间尺度变换技巧的来源**
- **Egerstedt, Wardi, Axelsson (2006)**, *Transition-Time Optimization for
  Switched-Mode Dynamical Systems*, IEEE TAC 51(1)
- **Johnson & Murphey (2011)**, *Second-Order Switching Time Optimization for
  Nonlinear Time-Varying Dynamic Systems*, IEEE TAC

---

## 16.2 ⚠️ `ocs2_frank_wolfe`

### 16.2.1 构建状态

与 `ocs2_ocs2` 同样的问题——依赖已删除的 `Dimensions.h`，且需要 **GLPK**
（GNU Linear Programming Kit）。

`ocs2_frank_wolfe/include/ocs2_frank_wolfe/FrankWolfeDescentDirection.h:35`：

```cpp
#include <glpk.h>
```

### 16.2.2 ⭐ Frank-Wolfe 算法

**问题**：在**凸紧集** $\mathcal{D}$ 上最小化可微函数 $f$：

$$
\min_{\theta\in\mathcal{D}}\ f(\theta)
$$

**Frank-Wolfe（条件梯度法）迭代**：

$$
\boxed{
\begin{aligned}
&\text{① 线性化：}\quad
   s^{k}=\arg\min_{s\in\mathcal{D}}\ \nabla f(\theta^{k})^{\!\top}s\qquad\text{（LP！）}\\
&\text{② 方向：}\quad d^{k}=s^{k}-\theta^{k}\\
&\text{③ 步长：}\quad \theta^{k+1}=\theta^{k}+\gamma_k\,d^{k},\qquad \gamma_k\in[0,1]
\end{aligned}}
$$

**核心优势**：
- 步骤 ① 是**线性规划**（在 $\mathcal{D}$ 上最小化线性函数），
  对多面体约束可以用单纯形法高效求解 → 这就是为什么需要 GLPK
- 迭代点 $\theta^{k+1}=(1-\gamma_k)\theta^k+\gamma_k s^k$ 是两个可行点的凸组合
  → **自动保持可行**，无需投影
- 对比投影梯度法：投影到复杂多面体可能比原问题还难

**收敛率**：$f(\theta^k)-f^\star = O(1/k)$（次线性）。
比投影梯度慢，但每步成本低得多。

**应用到切换时刻优化**：
$\mathcal{D}=\{\tau:\ t_0\le\tau_1\le\cdots\le\tau_K\le t_f\}$
是一个**单纯形**（多面体），LP 子问题的解一定在顶点上，
GLPK 秒解。

### 16.2.3 代码内容

| 文件 | 内容 |
|---|---|
| `FrankWolfeDescentDirection.h` | 上述算法的核心，调 GLPK 解 LP |
| `GradientDescent.h` | 完整的梯度下降框架（含线搜索） |
| `NLP_Cost.h` | 上层代价的接口：`getCost(θ)`、`getCostDerivative(θ)` |
| `NLP_Constraints.h` | 上层约束的接口：线性等式/不等式 |
| `NLP_Settings.h` | 设置 |

`FrankWolfeDescentDirection::run` 的额外功能：
**逐元素最大步长限制** $|\operatorname{diag}(e_v)\,d_v|\le1$，
防止某个方向变化过大（对应头文件注释的公式）。

### 16.2.4 参考文献

- **Frank & Wolfe (1956)**, *An Algorithm for Quadratic Programming*,
  Naval Research Logistics Quarterly 3(1-2):95–110 —— **原始论文**
- **Jaggi, M. (2013)**, *Revisiting Frank-Wolfe: Projection-Free Sparse Convex Optimization*,
  ICML —— 代码注释里的 `\cite jaggi13`
- **Lacoste-Julien & Jaggi (2015)**, *On the Global Linear Convergence of
  Frank-Wolfe Optimization Variants*, NeurIPS —— away-step 变体的线性收敛

---

## 16.3 `ocs2_ddp/unsupported/`

`ocs2_ddp` 里有一个 `unsupported/` 子目录：

| 文件 | 内容 |
|---|---|
| `SLQ_Hamiltonian.h` + `implementation/` | 基于**哈密顿系统**（而非 Riccati）的 SLQ 变体 |
| `MPC_OCS2.h` | 把 OCS2（双层优化）包成 MPC |
| `bvp_solver/BVPEquations.h` | 两点边值问题的方程 |
| `bvp_solver/SolveBVP.h` | BVP 求解器 |

### `SLQ_Hamiltonian` 的思路

不解 Riccati 方程，而直接积分**哈密顿正则方程**（庞特里亚金极小值原理）：

$$
\begin{aligned}
\dot x &= \frac{\partial\mathcal{H}}{\partial\lambda}=f(t,x,u)\\
\dot\lambda &= -\frac{\partial\mathcal{H}}{\partial x}
 = -\frac{\partial l}{\partial x}-\left(\frac{\partial f}{\partial x}\right)^{\!\top}\lambda\\
0 &= \frac{\partial\mathcal{H}}{\partial u}
 = \frac{\partial l}{\partial u}+\left(\frac{\partial f}{\partial u}\right)^{\!\top}\lambda
\end{aligned}
$$

边界条件：$x(t_0)=x_0$、$\lambda(t_f)=\dfrac{\partial\phi}{\partial x}(x(t_f))$
—— 这是一个**两点边值问题（TPBVP）**。

**与 Riccati 方法的关系**：
若假设 $\lambda(t)=S_m(t)\,\delta x(t)+S_v(t)$，
代入上面的方程组即可**导出 Riccati 方程**——
两者数学等价，Riccati 是"把 TPBVP 转化为初值问题"的技巧。

**为什么 unsupported**：
- TPBVP 需要打靶法或配点法求解，比 Riccati 的单次后向积分慢
- 数值稳定性差（协态方程是不稳定的后向系统）
- Riccati 直接给出**反馈增益**，哈密顿法只给开环解

保留它主要是历史与教学价值。

`bvp_solver/` 是配套的 BVP 求解器实现。

---

## 16.4 `ocs2_raisim`：实验性但可用

| 包 | 状态 |
|---|---|
| `ocs2_raisim_core` | 可用，但需要 RaiSim 授权（学术免费，商用付费） |
| `ocs2_raisim_ros` | 可用 |
| `ocs2_legged_robot_raisim` | 可用 |

**为什么单列**：RaiSim 不是开源的（需要单独申请许可证），
所以默认不在 CI 中构建。但它是评估 MPC 鲁棒性的重要工具
（见 [13 章](13-robot-examples.md) §13.8）。

`RaisimRollout` 实现 `RolloutBase` 接口，
把 rollout 从"用 MPC 自己的模型积分"换成"用真实物理引擎仿真"，
从而暴露模型误差的影响。

**参考**：Hwangbo, Lee, Hutter (2018),
*Per-Contact Iteration Method for Solving Contact Dynamics*, RA-L。

---

## 16.5 如何判断一个模块是否在用

三个检查：

**① 看 CMakeLists**

```bash
grep -n "add_library\|add_executable" <package>/CMakeLists.txt
```

若只有 `lintTarget`，说明没有实际编译产物。

**② 看是否在元包依赖里**

```bash
grep -n "exec_depend" ocs2/package.xml
```

注意：`ocs2_ocs2` **不在**元包依赖列表里，而 `ocs2_frank_wolfe` **在**——
但后者也编译不过。元包依赖列表本身并不完全可靠。

**③ 看头文件的 include 是否存在**

```bash
grep -rn "#include <ocs2_core/" <package>/include/ | \
  awk -F'[<>]' '{print $2}' | sort -u | \
  while read f; do [ -f "$(echo $f | sed 's|ocs2_core/|ocs2_core/include/ocs2_core/|')" ] || echo "MISSING: $f"; done
```

这能立刻发现 `ocs2_core/Dimensions.h` 这类问题。

---

## 16.6 与官方文档的对应

`ocs2_doc/docs/` 下的 Sphinx 文档（对应 https://leggedrobotics.github.io/ocs2/）：

| 文件 | 内容 |
|---|---|
| `intro.rst` | 项目简介 |
| `overview.rst` | 算法总览 |
| `installation.rst` | 安装步骤 |
| `getting-started.rst` | 快速上手 |
| `optimal_control_modules.rst` | 最优控制模块（对应本文档 02–08 章） |
| `from_urdf_to_ocp.rst` | URDF → OCP（对应本文档 12 章） |
| `robotic_examples.rst` | 机器人示例（对应本文档 13 章） |
| `mpcnet.rst` | MPC-Net（对应本文档 14 章） |
| `profiling.rst` | 性能分析 |
| `faq.rst` | 常见问题 |
| `refs.bib` | 参考文献（本文档 17 章包含它的全部条目） |

**官方文档不涉及**：数学推导、遗留模块、实现细节。
本文档的定位是补充这些。

---

## 16.7 迁移建议

如果你确实需要切换时刻优化：

### 方案 A：现代重写（推荐）

把 $\tau$ 作为**决策变量**加入多重打靶：
- 每个区间的时长 $\Delta t_k$ 变成决策变量
- 动力学变成 $x_{k+1}=F(x_k,u_k,\Delta t_k)$，
  灵敏度 $\partial F/\partial\Delta t_k = f(x,u)$（RK 的高阶版本更复杂但可 AD）
- 用 SQP/IPM 一起解

参考：Winkler et al. (2018) 的 TOWR、
以及 `CubicSpline::startTimeDerivative`（[13 章](13-robot-examples.md) §13.6.2）——
**OCS2 里已经有对时长求导的代码**，说明这个方向被考虑过。

### 方案 B：移植 GDDP

把 `sensitivity_equations/` 从 `Dimensions<>` 迁移到动态维度：
- `state_vector_t` → `vector_t`
- `state_matrix_t` → `matrix_t`
- 模板参数 `<STATE_DIM, INPUT_DIM>` → 构造函数参数
- `DDP_DataCollector` 需要适配新的 `PrimalDataContainer`/`DualDataContainer`

工作量中等（约 2000 行），但需要理解 GDDP 的全部数学。

### 方案 C：外部规划器

绝大多数情况下这才是对的做法：
用启发式或专门的接触规划器（如 MIP、采样法）给出 `ModeSchedule`，
MPC 只优化给定序列下的连续量。

四足的 `GaitSchedule` 就是这个思路的最简版本。

---

**下一章**：完整参考文献 → [17 参考文献](17-references.md)
