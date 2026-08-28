# 09 · Loopshaping：频域整形与系统增广

`ocs2_core/loopshaping/` 是 OCS2 里一个独立且优雅的子系统：
把**频域设计思想**引入时域最优控制。

**动机**（Grandia et al. 2019, *Frequency-Aware Model Predictive Control*）：

> 标准 MPC 用 $\frac12 u^\top Ru$ 惩罚输入幅值，这是一个**全频段均匀**的惩罚。
> 但实际系统往往只希望抑制**特定频段**：
> - 抑制高频输入（避免激励未建模的柔性模态、减少执行器磨损）
> - 允许低频大幅值输入（正常运动需要）
>
> 标准 MPC 做不到这一点——除非把 $R$ 调得很大，但那会同时抑制低频，让系统变迟钝。

**解法**：把一个**滤波器**串到输入通道上，惩罚**滤波后**的信号。
滤波器的频率响应决定了哪些频段被惩罚。

---

## 9.1 数学框架

### 9.1.1 滤波器

**文件**：`ocs2_core/loopshaping/LoopshapingFilter.h`

线性时不变滤波器，状态空间形式：

$$
\begin{aligned}
\dot x_f &= A_f\,x_f + B_f\,\eta\\
y_f &= C_f\,x_f + D_f\,\eta
\end{aligned}
$$

其中 $x_f\in\mathbb{R}^{n_f}$ 是滤波器状态，$\eta$ 是滤波器输入，$y_f$ 是输出。

传递函数 $W(s)=C_f(sI-A_f)^{-1}B_f+D_f$。

`Filter` 类存 $A_f,B_f,C_f,D_f$，并缓存它们的对角版本
（`a_`, `b_`, `c_`, `d_` 为 `Eigen::DiagonalMatrix`）。
若所有矩阵都是对角的，走**逐元素乘法**的快速路径
（`LoopshapingDefinition::isDiagonal()`），避免通用矩阵乘法。

`LoopshapingPropertyTree.h` 从 `.info` 配置文件构造滤波器，
支持直接给 $(A_f,B_f,C_f,D_f)$，也支持给传递函数系数
（通过 `ocs2_core/dynamics/TransferFunctionBase.h` 转成状态空间）。

### 9.1.2 增广系统

无论哪种模式，增广状态都是

$$
\boxed{\;\tilde x=\begin{bmatrix}x\\x_f\end{bmatrix}\in\mathbb{R}^{n_x+n_f}\;}
$$

（`LoopshapingDefinition::concatenateSystemAndFilterState`，
`getSystemState` 取前 $n_x$ 个，`getFilterState` 取后 $n_f$ 个）

增广输入 $\tilde u$ 的含义**取决于模式**。

---

## 9.2 ⭐ 两种模式

**文件**：`ocs2_core/loopshaping/LoopshapingDefinition.h:44`

```cpp
enum class LoopshapingType { outputpattern, eliminatepattern };
```

### 9.2.1 `outputpattern`（输出模式）

**结构**：

```
        ┌──────────────┐
 ũ ──┬─▶│  原系统 f    │──▶ ẋ
     │  └──────────────┘
     │  ┌──────────────┐
     └─▶│  滤波器 W    │──▶ y_f  （被惩罚的量）
        └──────────────┘
```

- **增广输入 = 原系统输入**：$\tilde u = u$
- **滤波器与系统并联**，滤波器输入也是 $u$
- **被惩罚的滤波输出**：$y_f = C_f x_f + D_f u$

代码（`LoopshapingDefinition.cpp:64-84, 92-111`）：

```cpp
// getSystemInput
case LoopshapingType::outputpattern:  systemInput = input;  break;

// getFilteredInput
case LoopshapingType::outputpattern:
  filteredInput = C_f·x_f + D_f·input;   break;
```

**增广动力学**：

$$
\dot{\tilde x}=\begin{bmatrix}f(t,x,\tilde u)\\ A_fx_f+B_f\tilde u\end{bmatrix}
$$

**增广代价**：

$$
\tilde l = l(t,x,\tilde u) + \tfrac12\,y_f^{\!\top}R_f\,y_f
= l + \tfrac12\big(C_fx_f+D_f\tilde u\big)^{\!\top}R_f\big(C_fx_f+D_f\tilde u\big)
$$

（`loopshapingCost(filteredInput)` = $\frac12 y_f^\top R_f y_f$，
$R_f$ 是 `LoopshapingDefinition::costMatrix()`）

**用途**：$W(s)$ 设计为**高通滤波器**，则高频输入被重罚、低频几乎不罚。
典型：$W(s)=\frac{s}{s+\omega_c}$，$\omega_c$ 以上的分量被完整惩罚。

### 9.2.2 `eliminatepattern`（消去模式）

**结构**：

```
        ┌──────────────┐    ┌──────────────┐
 ũ ────▶│  滤波器 W    │───▶│  原系统 f    │──▶ ẋ
        └──────────────┘ u  └──────────────┘
```

- **增广输入 = 滤波器输入**：$\tilde u = \eta$
- **滤波器串联在输入通道上**
- **原系统输入由滤波器输出给出**：$u = C_fx_f + D_f\tilde u$

代码（`LoopshapingDefinition.cpp:70-80, 104-106`）：

```cpp
// getSystemInput
case LoopshapingType::eliminatepattern:
  systemInput = C_f·x_f + D_f·input;   break;      // u = C_f x_f + D_f ũ

// getFilteredInput
case LoopshapingType::eliminatepattern:
  filteredInput = input;  break;                    // y_f = ũ
```

**增广动力学**：

$$
\dot{\tilde x}=\begin{bmatrix}f\big(t,x,\ C_fx_f+D_f\tilde u\big)\\ A_fx_f+B_f\tilde u\end{bmatrix}
$$

**增广代价**：

$$
\tilde l = l\big(t,x,\ C_fx_f+D_f\tilde u\big) + \tfrac12\,\tilde u^{\!\top}R_f\,\tilde u
$$

**用途**：$W(s)$ 设计为**低通滤波器**，则真实输入 $u$ 天然是平滑的
（高频分量被滤波器衰减），同时对 $\tilde u$ 的常规惩罚保持有效。

**名字的由来**："eliminate" 指真实输入 $u$ 被"消去"（不再是决策变量），
决策变量变成滤波器的输入 $\tilde u$。

### 9.2.3 对比

| | outputpattern | eliminatepattern |
|---|---|---|
| $\tilde u$ | $=u$（原输入） | $=\eta$（滤波器输入） |
| $u$ | $=\tilde u$ | $=C_fx_f+D_f\tilde u$ |
| $y_f$（被罚） | $=C_fx_f+D_f\tilde u$ | $=\tilde u$ |
| 滤波器拓扑 | 并联（旁路观测） | 串联（在通道上） |
| $W$ 典型设计 | 高通（罚高频） | 低通（平滑输入） |
| 输入约束 $u\in\mathcal{U}$ | 直接施加在 $\tilde u$ 上 | 需要通过 $C_fx_f+D_f\tilde u$ 表达 |

**两者的对偶性**：outputpattern 通过"惩罚高频"间接抑制高频；
eliminatepattern 通过"结构上无法产生高频"直接保证。
后者更强但限制了可达输入集。

---

## 9.3 平衡点求解

**文件**：`LoopshapingDefinition.cpp:141` 的 `getFilterEquilibrium`

MPC 初始化与操作点计算需要知道：给定期望的稳态 $u$，
滤波器状态 $x_f$ 与输入 $\eta$ 应该是多少？

稳态条件 $\dot x_f=0$：

$$
A_fx_f+B_f\eta=0,\qquad y_f=C_fx_f+D_f\eta
$$

**`findEquilibriumForInput(u, x_f, y_f)`**（outputpattern 用）：
已知滤波器输入 $\eta=u$，求 $x_f$ 与输出：

$$
x_f=-A_f^{-1}B_f\,u,\qquad y_f=C_fx_f+D_fu
$$

代码用 `Aqr_`（$A_f$ 的 ColPivHouseholderQR 预分解）求解。

**`findEquilibriumForOutput(y, x_f, \eta)`**（eliminatepattern 用）：
已知期望输出 $y=u$，求 $(x_f,\eta)$：

$$
\begin{bmatrix}A_f & B_f\\ C_f & D_f\end{bmatrix}
\begin{bmatrix}x_f\\ \eta\end{bmatrix}
=\begin{bmatrix}0\\ y\end{bmatrix}
$$

代码用 `ABCDqr_`（增广矩阵的 QR 预分解）求解。
**这个线性系统可能欠定或超定**（取决于滤波器维度），
ColPivHouseholderQR 给出最小二乘解。

**`findEquilibriumForOutputGivenState(y, x_f, \eta)`**：
$x_f$ 已知时只解 $\eta$：$D_f\eta = y - C_fx_f$，用 `Dqr_`。

**预分解的意义**：这些 QR 在 `Filter` 构造时算一次，
之后每次求平衡点只需回代，$O(n^2)$ 而非 $O(n^3)$。

---

## 9.4 装饰器体系

Loopshaping 用**装饰器模式**包装原问题的每一个组件。
目录结构完全镜像 `ocs2_core` 本身：

```
loopshaping/
├── LoopshapingDefinition.h        核心：类型 + 滤波器 + 状态/输入映射
├── LoopshapingFilter.h            滤波器状态空间
├── LoopshapingPreComputation.h    预计算增广量
├── LoopshapingPropertyTree.h      从配置文件加载
├── Loopshaping.h                  工厂函数汇总入口
├── dynamics/
│   ├── LoopshapingDynamics.h                基类
│   ├── LoopshapingDynamicsOutputPattern.h
│   ├── LoopshapingDynamicsEliminatePattern.h
│   └── LoopshapingFilterDynamics.h          纯滤波器仿真（MRT 用）
├── cost/
│   ├── LoopshapingCost.h / StateCost / StateInputCost
│   ├── LoopshapingCostOutputPattern.h
│   └── LoopshapingCostEliminatePattern.h
├── constraint/     （同样的 Output/Eliminate 分裂）
├── soft_constraint/
├── augmented_lagrangian/
└── initialization/LoopshapingInitializer.h
```

**每个装饰器做的事**：接收增广状态 $\tilde x$ 与输入 $\tilde u$，
- 用 `getSystemState()`、`getSystemInput()` 还原出 $(x,u)$
- 调用被包装的原对象
- 把结果的导数按链式法则**映射回增广空间**

### 9.4.1 导数映射（以 outputpattern 的代价为例）

原代价 $l(x,u)$，增广代价 $\tilde l(\tilde x,\tilde u)=l(x,\tilde u)+\frac12 y_f^\top R_fy_f$。
由于 outputpattern 下 $x = \tilde x_{1:n_x}$、$u=\tilde u$：

$$
\frac{\partial\tilde l}{\partial\tilde x}=\begin{bmatrix}\dfrac{\partial l}{\partial x}\\[6pt]
C_f^{\!\top}R_fy_f\end{bmatrix},\qquad
\frac{\partial\tilde l}{\partial\tilde u}=\frac{\partial l}{\partial u}+D_f^{\!\top}R_fy_f
$$

$$
\frac{\partial^2\tilde l}{\partial\tilde x^2}=
\begin{bmatrix}\dfrac{\partial^2l}{\partial x^2} & 0\\ 0 & C_f^{\!\top}R_fC_f\end{bmatrix},\quad
\frac{\partial^2\tilde l}{\partial\tilde u\partial\tilde x}=
\begin{bmatrix}\dfrac{\partial^2l}{\partial u\partial x} & D_f^{\!\top}R_fC_f\end{bmatrix},\quad
\frac{\partial^2\tilde l}{\partial\tilde u^2}=\frac{\partial^2l}{\partial u^2}+D_f^{\!\top}R_fD_f
$$

### 9.4.2 eliminatepattern 的额外复杂度

因为 $u=C_fx_f+D_f\tilde u$ 同时依赖 $\tilde x$ 的滤波器部分和 $\tilde u$，
链式法则更复杂：

$$
\frac{\partial\tilde l}{\partial\tilde x}=
\begin{bmatrix}\dfrac{\partial l}{\partial x}\\[6pt]
C_f^{\!\top}\dfrac{\partial l}{\partial u}\end{bmatrix},\qquad
\frac{\partial\tilde l}{\partial\tilde u}=D_f^{\!\top}\frac{\partial l}{\partial u}+R_f\tilde u
$$

二阶导需要展开 $\begin{bmatrix}I&0\\0&C_f\end{bmatrix}$ 与 $D_f$ 的所有组合，
这就是为什么 `LoopshapingCostEliminatePattern.h` 比 `OutputPattern` 版本长得多。

**对角优化**：`Filter` 缓存了 `diagCC_`（$C_f^\top R_fC_f$ 的对角版）、
`diagDC_`、`diagDD_`，当滤波器对角时直接用它们做逐元素运算，
避免 $O(n^3)$ 的矩阵乘法（`LoopshapingFilter.h:62-64`）。

---

## 9.5 求解器无关性

关键设计：**loopshaping 完全在问题层面完成**。

`ocs2_oc/oc_problem/LoopshapingOptimalControlProblem.h` 提供工厂函数：

```cpp
OptimalControlProblem create(const OptimalControlProblem& problem,
                             std::shared_ptr<LoopshapingDefinition> definition);
```

它把原问题的每个组件用对应装饰器包起来，返回一个**新的 `OptimalControlProblem`**。
之后**任何求解器**（SLQ/SQP/IPM/SLP）都可以直接用它，完全不知道 loopshaping 的存在。

配套的还有：
- `ocs2_oc/synchronized_module/LoopshapingReferenceManager.h`：把参考轨迹扩展到增广空间
- `ocs2_oc/oc_data/LoopshapingPrimalSolution.h`：从增广解中提取原系统的解
- `ocs2_robotic_tools/common/LoopshapingRobotInterface.h`：机器人接口的包装
- `ocs2_mpc/LoopshapingSystemObservation.h`：观测的增广
- `ocs2_ros_interfaces/mrt/LoopshapingDummyObserver.h`：可视化

---

## 9.6 MRT 侧的滤波器仿真

**文件**：`loopshaping/dynamics/LoopshapingFilterDynamics.h`

MPC 输出的是增广空间的策略 $\tilde u(t,\tilde x)$，
但执行器需要真实输入 $u$。所以控制侧必须：

1. 维护滤波器状态 $x_f$ 的估计（**实时积分滤波器动力学**）
2. 用 $\tilde x = [\hat x; x_f]$ 查询策略得到 $\tilde u$
3. 用 `getSystemInput(\tilde x, \tilde u)` 还原 $u$ 发给执行器

`LoopshapingFilterDynamics` 提供 `advance(dt, input)` 做步进积分。

**注意**：$x_f$ 是**控制器内部状态**，不可观测。
如果 MPC 与控制器的 $x_f$ 不同步（例如控制器重启），会产生瞬态。
实践中需要在切换时用 `getFilterEquilibrium()` 重新初始化。

---

## 9.7 设计示例

### 抑制高频（outputpattern + 高通）

一阶高通 $W(s)=\dfrac{s}{s+\omega_c}$ 的状态空间实现：

$$
A_f=-\omega_c,\quad B_f=1,\quad C_f=-\omega_c,\quad D_f=1
$$

验证：$W(s)=C_f(s-A_f)^{-1}B_f+D_f=\dfrac{-\omega_c}{s+\omega_c}+1=\dfrac{s}{s+\omega_c}$ ✓

$\omega\ll\omega_c$ 时 $|W|\approx\omega/\omega_c\to0$（低频不罚）；
$\omega\gg\omega_c$ 时 $|W|\to1$（高频全罚）。

### 平滑输入（eliminatepattern + 低通）

一阶低通 $W(s)=\dfrac{\omega_c}{s+\omega_c}$：

$$
A_f=-\omega_c,\quad B_f=\omega_c,\quad C_f=1,\quad D_f=0
$$

真实输入 $u = x_f$，其动力学 $\dot x_f=-\omega_c x_f+\omega_c\tilde u$
保证 $u$ 的变化率有界 → **输入天然平滑**。

**注意 $D_f=0$** 意味着 $\tilde u$ 对 $u$ 无直接影响（严格真传递函数）。
这会让 $\partial u/\partial\tilde u=0$，从而增广系统的 $\tilde B$ 在第一时刻为 0——
可能导致可控性问题。实践中常用 $W(s)=\frac{\omega_c}{s+\omega_c}+\epsilon$
（即 $D_f=\epsilon$ 小量）来避免。

---

## 9.8 测试

`ocs2_core/test/loopshaping/` 下有完整的单元测试：
`testLoopshapingDynamics.cpp`、`testLoopshapingCost.cpp`、
`testLoopshapingConstraint.cpp` 等，
用有限差分验证所有装饰器的导数映射正确性。

**这些测试是理解装饰器数学的最好材料**——
它们显式构造增广量并与数值导数对比。

---

## 9.9 参考文献

- **Grandia, Farshidian, Dosovitskiy, Ranftl, Hutter (2019)**,
  *Frequency-Aware Model Predictive Control*, IEEE RA-L 4(2):1517–1524
  —— **本模块的原始论文**
- **Skogestad & Postlethwaite (2005)**, *Multivariable Feedback Control: Analysis and Design*, 2nd ed.
  —— loopshaping 的经典教材（$H_\infty$ 语境下的加权函数设计）
- **Zhou, Doyle, Glover (1996)**, *Robust and Optimal Control* —— 加权函数与鲁棒性

**下一章**：积分与切换系统 → [10 积分/Rollout](10-integration-rollout.md)
