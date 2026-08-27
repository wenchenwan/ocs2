# 02 · ocs2_core：类型系统与问题建模原语

`ocs2_core` 是整个库的地基，**不含任何求解算法**，只提供：
类型定义、代价/约束/动力学的抽象接口、罚函数、积分器、loopshaping、并发工具。

273 个文件、约 3 万行。本章按目录逐个拆解。

---

## 2.1 `Types.h`：全库的类型词汇表

**文件**：`ocs2_core/include/ocs2_core/Types.h`

### 2.1.1 基础别名

```cpp
using scalar_t  = double;
using vector_t  = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
using matrix_t  = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
using vector_array_t  = std::vector<vector_t>;     // 一条轨迹
using vector_array2_t = std::vector<vector_array_t>; // 多条轨迹
```

> **历史注记**：早期 OCS2 用 `Dimensions<STATE_DIM, INPUT_DIM>` 模板做编译期固定维度，
> 现已全面转为动态维度。这是为什么遗留模块 `ocs2_ocs2` 编译不过——它还在 include
> 已被删除的 `ocs2_core/Dimensions.h`。

### 2.1.2 四个近似结构体

这是**全库最重要的四个类型**，所有导数信息都通过它们传递：

| 类型 | 表示的对象 | 字段 |
|---|---|---|
| `ScalarFunctionLinearApproximation` | 标量函数一阶近似 | `dfdx`(vec), `dfdu`(vec), `f`(scalar) |
| `ScalarFunctionQuadraticApproximation` | 标量函数二阶近似 | `dfdxx`, `dfdux`, `dfduu`(mat), `dfdx`, `dfdu`(vec), `f` |
| `VectorFunctionLinearApproximation` | 向量函数一阶近似 | `dfdx`, `dfdu`(mat), `f`(vec) |
| `VectorFunctionQuadraticApproximation` | 向量函数二阶近似 | `dfdxx`, `dfdux`, `dfduu` 为 `matrix_array_t`（每个输出分量一个 Hessian） |

标量二次近似的语义（`Types.h:141` 注释）：

$$
f(x,u)\approx \tfrac12\,\delta x^{\!\top} f_{xx}\,\delta x
 +\delta u^{\!\top} f_{ux}\,\delta x
 +\tfrac12\,\delta u^{\!\top} f_{uu}\,\delta u
 + f_x^{\!\top}\delta x + f_u^{\!\top}\delta u + f
$$

**关键形状约定**：`dfdux` 是 $n_u\times n_x$，即"行是输入、列是状态"。
写自定义代价项时搞反会在 `checkSize()` 处抛异常。

**代数运算**：这些结构体重载了 `operator+=`（逐字段相加）与 `operator*=`（整体缩放）。
于是"多个代价项累加"写起来就是：

```cpp
ScalarFunctionQuadraticApproximation cost;
for (auto& term : terms) cost += term.getQuadraticApproximation(t, x, u, targetTraj, preComp);
cost *= dt;   // 前向欧拉积分权重（多重打靶里用）
```

**尺寸/数值检查**：`checkSize()`、`checkBeingPSD()`（Types.h 尾部与
`ocs2_core/src/Types.cpp`）在 `checkNumericalStability_` 打开时被求解器调用，
定位"Hessian 非半正定""维度对不上"这类问题非常有效。

### 2.1.3 `NumericTraits.h`

提供 `numeric_traits::limitEpsilon<T>()`、`weakEpsilon<T>()` 等容差常数，
避免代码里散落魔法数。配合 `misc/Numerics.h` 的 `almost_eq`、`almost_ge`。

---

## 2.2 `dynamics/`：系统动力学

**继承链**：

```
OdeBase                                (integration/OdeBase.h)
  └─ ControlledSystemBase              带控制器的 ODE，可被积分器直接 rollout
       └─ SystemDynamicsBase           增加线性化接口
            ├─ SystemDynamicsBaseAD    用 CppAD 自动微分实现线性化
            ├─ SystemDynamicsLinearizer 用有限差分实现线性化
            └─ LinearSystemDynamics    直接给定 A、B 矩阵
```

### 2.2.1 `SystemDynamicsBase` 的接口

```cpp
virtual vector_t computeFlowMap(scalar_t t, const vector_t& x, const vector_t& u,
                                const PreComputation&) = 0;
virtual VectorFunctionLinearApproximation linearApproximation(...) = 0;

virtual vector_t computeJumpMap(scalar_t t, const vector_t& x, const PreComputation&);
virtual VectorFunctionLinearApproximation jumpMapLinearApproximation(...);

virtual vector_t computeGuardSurfaces(scalar_t t, const vector_t& x);
virtual matrix_t dynamicsCovariance(scalar_t t, const vector_t& x, const vector_t& u);
```

对应数学对象：

| 方法 | 数学 | 说明 |
|---|---|---|
| `computeFlowMap` | $\dot x = f(t,x,u)$ | 连续流形 |
| `linearApproximation` | $\delta\dot x = A\,\delta x + B\,\delta u + f$ | $A=\partial f/\partial x$, $B=\partial f/\partial u$ |
| `computeJumpMap` | $x^+ = j(t,x^-)$ | **跳变映射**，事件时刻状态不连续 |
| `jumpMapLinearApproximation` | $\delta x^+ = A_j\,\delta x^- + b_j$ | 跳变的线性化 |
| `computeGuardSurfaces` | $g(t,x)$ | **保护面**，零穿越触发事件（状态触发模式） |
| `dynamicsCovariance` | $\Sigma$ | 过程噪声协方差，仅风险敏感（iLEG）用 |

跳变映射的默认实现是恒等 $x^+=x^-$。四足机器人用它做落足冲击：
落足瞬间接触点速度被投影到零（碰撞模型）。

### 2.2.2 `SystemDynamicsBaseAD`：自动微分路径

用户只需实现

```cpp
virtual ad_vector_t systemFlowMap(ad_scalar_t t, const ad_vector_t& x,
                                  const ad_vector_t& u, const ad_vector_t& p) const = 0;
```

用 `ad_scalar_t` 写出动力学，`initialize()` 会：
1. 用 CppAD 记录计算图
2. 用 CppADCodeGen 生成 C 源码
3. 调用系统编译器编译成 `.so`
4. 后续直接 `dlopen` 调用——**运行时开销接近手写代码**

`ocs2_cartpole` 就是一个 30 行的完整示例
（`ocs2_robotic_examples/ocs2_cartpole/include/ocs2_cartpole/dynamics/CartPoleSystemDynamics.h`）。

编译产物默认落在 `/tmp/ocs2`（可通过 `libraryFolder` 改），
**第一次运行会卡几秒到几十秒，之后从磁盘加载**。这是新手最常见的困惑点。

### 2.2.3 `TransferFunctionBase.h`

SISO 传递函数 $\frac{b_0+b_1 s+\cdots}{a_0+a_1 s+\cdots}$ 到状态空间的转换，
供 loopshaping 模块构造滤波器用（见 [09 章](09-loopshaping.md)）。

---

## 2.3 `cost/`：代价函数

### 2.3.1 两类代价

```
StateCost        l(t, x, targetTrajectories, preComp)         状态代价（终端/事件/中间）
StateInputCost   l(t, x, u, targetTrajectories, preComp)      状态-输入代价（中间）
```

各自有：
- `getValue()` 返回标量
- `getQuadraticApproximation()` 返回 `ScalarFunctionQuadraticApproximation`

### 2.3.2 具体实现

| 类 | 说明 |
|---|---|
| `QuadraticStateCost` | $\frac12 (x-x_d)^\top Q (x-x_d)$，子类通过 `getStateDeviation()` 定义"偏差" |
| `QuadraticStateInputCost` | $\frac12\big[(x-x_d)^\top Q(x-x_d) + (u-u_d)^\top R(u-u_d)\big]$ |
| `StateCostCppAd` / `StateInputCostCppAd` | 用 CppAD 自动求二阶导 |
| `StateInputGaussNewtonCostAd` | ⭐ 见下 |

### 2.3.3 `StateInputGaussNewtonCostAd`：高斯-牛顿代价

**文件**：`ocs2_core/include/ocs2_core/cost/StateInputGaussNewtonCostAd.h`

用户定义**残差向量** $r(t,x,u,p)\in\mathbb{R}^m$，代价为

$$
l = \tfrac12\, r^{\!\top} r
$$

其精确 Hessian 是

$$
\nabla^2 l = J_r^{\!\top} J_r + \sum_{i=1}^{m} r_i \nabla^2 r_i,
\qquad J_r=\frac{\partial r}{\partial (x,u)}
$$

**高斯-牛顿近似**丢掉第二项：

$$
\nabla^2 l \;\approx\; J_r^{\!\top} J_r \;\succeq\; 0
$$

好处有两个：
1. 只需 $r$ 的**一阶**导数，自动微分开销减半；
2. $J_r^\top J_r$ **天然半正定**，Riccati 递推不会因为 Hessian 不定而发散。

代价是收敛阶从二次退化为超线性——但在 MPC 里每次只跑几次迭代，几乎没有影响。
这正是 OCS2 被称为 **Gauss-Newton DDP** 的原因（`GaussNewtonDDP` 类名即来源于此）。

### 2.3.4 `Collection<T>` 模式

**文件**：`ocs2_core/misc/Collection.h`

`StateCostCollection`、`StateInputCostCollection`、`StateConstraintCollection` 等
都继承自 `Collection<T>`，它是一个"具名项的容器"：

```cpp
collection.add("footPlacement", std::make_unique<MyCost>(...));
auto& term = collection.get<MyCost>("footPlacement");
```

`getValue()` 会遍历所有项求和，`getQuadraticApproximation()` 会累加所有二次近似。
这让"往问题里加一项代价"变成一行代码，且支持按名字取回做在线调参。

---

## 2.4 `constraint/`：约束

### 2.4.1 分类

OCS2 把约束按**依赖变量**和**约束类型**做 2×2 分类：

|  | 等式 $g=0$ | 不等式 $h\ge 0$ |
|---|---|---|
| **仅状态** $x$ | `StateConstraint` | `StateConstraint` |
| **状态-输入** $(x,u)$ | `StateInputConstraint` | `StateInputConstraint` |

再按**作用时刻**分为：中间（intermediate）、事件前（pre-jump）、终端（final）。
所以 `OptimalControlProblem` 里有 8 个约束容器指针。

### 2.4.2 `ConstraintOrder`

```cpp
enum class ConstraintOrder { Linear, Quadratic };
```

约束可以只提供线性近似（SQP/IPM 只需这个），也可以提供二次近似
（DDP 的增广拉格朗日需要，因为罚函数的 Hessian 里含 $\nabla^2 h$ 项）。

### 2.4.3 状态-输入等式约束的特殊地位

⚠️ 注意 `OptimalControlProblem.h:70` 的注释：

> `/** Intermediate equality constraints, full row rank w.r.t. inputs */`

**状态-输入等式约束 $D\,\delta u + C\,\delta x + e = 0$ 必须对输入满行秩**。
原因是 DDP 与 SQP 都用**零空间投影**来消去这类约束（见 [04](04-ocs2-ddp.md) §4.4
与 [03](03-ocs2-oc.md) §3.6），如果 $D$ 行亏秩，QR 分解会得到奇异的 $R$，
虽然 `setTriangularMinimumEigenvalues()` 做了兜底，但解会失去意义。

典型用法：四足机器人的**支撑腿零速度约束** $J_c\,v = 0$，
每条支撑腿贡献 3 行，且这些行关于关节速度输入是满秩的。

---

## 2.5 `soft_constraint/`：软约束

把约束 $h(x,u)\ge0$ 转成代价项 $p(h(x,u))$：

| 类 | 说明 |
|---|---|
| `StateSoftConstraint` | 包装 `StateConstraint` + 罚函数 |
| `StateInputSoftConstraint` | 包装 `StateInputConstraint` + 罚函数 |
| `StateInputSoftBoxConstraint` | 直接对状态/输入的**箱式界**做罚，无需写约束类 |

链式法则给出软约束的一二阶导（`ocs2_core/src/soft_constraint/StateInputSoftConstraint.cpp`）：

$$
\frac{\partial p}{\partial z}=p'(h)\,\frac{\partial h}{\partial z},\qquad
\frac{\partial^2 p}{\partial z^2}
 = p''(h)\,\frac{\partial h}{\partial z}\frac{\partial h}{\partial z}^{\!\top}
 + p'(h)\,\frac{\partial^2 h}{\partial z^2}
$$

其中 $z=(x,u)$。第二项需要约束的二阶导——这就是 `ConstraintOrder::Quadratic` 的用途。
若约束只提供线性近似，则该项被丢弃（相当于对软约束也做高斯-牛顿近似）。

罚函数 $p(\cdot)$ 的具体形式见 [08 章](08-constraints-penalties.md)。

---

## 2.6 `augmented_lagrangian/`：增广拉格朗日

与软约束的区别：**罚函数额外依赖拉格朗日乘子 $\lambda$**。

```cpp
class StateInputAugmentedLagrangianInterface {
  virtual ScalarFunctionQuadraticApproximation
      getQuadraticApproximation(t, x, u, const Multiplier&, preComp);
  virtual void updateLagrangian(t, x, u, Multiplier&);   // 对偶上升
};
```

数据流：
1. 求解器在 LQ 近似时把当前乘子传进来，得到增广代价
   （`ocs2_oc/src/approximate_model/LinearQuadraticApproximator.cpp:64-83`）
2. 求解器求解后调用 `updateDualSolution()` 做对偶更新
   （`ocs2_oc/src/oc_problem/OptimalControlProblemHelperFunction.cpp`）

推导见 [08 章](08-constraints-penalties.md)。

---

## 2.7 `penalties/`：罚函数库

两套并列的体系：

```
penalties/penalties/           普通罚函数 p(h)         → 软约束用
  PenaltyBase
  ├─ QuadraticPenalty          μ/2 · h²                （等式）
  ├─ RelaxedBarrierPenalty     -μ ln(h)，h<δ 时二次延拓（不等式）
  ├─ SquaredHingePenalty       μ/2 · (h-δ)²  if h<δ    （不等式）
  ├─ SmoothAbsolutePenalty     μ√(h²+δ²)              （等式，L1 光滑化）
  └─ DoubleSidedPenalty        p(h-l) + p(u-h)         （箱约束）

penalties/augmented/           增广罚函数 p(h, λ)      → 增广拉格朗日用
  AugmentedPenaltyBase
  ├─ QuadraticPenalty              -λh + ρ/2·h²
  ├─ SlacknessSquaredHingePenalty  PHR 罚（不等式）
  ├─ ModifiedRelaxedBarrierPenalty 光滑 PHR（不等式）
  └─ SmoothAbsolutePenalty         -λh + μ√(h²+δ²)
```

`MultidimensionalPenalty` 把标量罚函数逐分量应用到向量约束上，
并组装出向量约束的二次近似。

完整数学在 [08 章](08-constraints-penalties.md)。

---

## 2.8 `integration/`：数值积分

### 2.8.1 两套并存的积分设施

**(a) `Integrator` / `IntegratorBase`** —— 基于 Boost.Numeric.Odeint，用于 **rollout**（前向仿真）

`steppers.h` 定义了可用的步进器：

| `IntegratorType` | 说明 |
|---|---|
| `EULER` | 显式欧拉，1 阶 |
| `RK4` | 经典 4 阶 Runge-Kutta，定步长 |
| `ODE45` | Dormand-Prince 5(4)，**自适应步长**（默认） |
| `RK5_VARIABLE` | 5 阶变步长 |
| `ADAMS_BASHFORTH` / `ADAMS_BASHFORTH_MOULTON` | 多步法 |
| `BULIRSCH_STOER` | 外推法，高精度 |

`eigenIntegration.h` 是让 Eigen 向量能被 odeint 使用的适配层
（定义 `norm_inf`、`operator/` 等 odeint 要求的代数运算）。

`RungeKuttaDormandPrince5.h` 是自实现的 DP5，用于对 odeint 的行为做精确控制。

**(b) `SensitivityIntegrator`** —— 用于 **多重打靶的离散化 + 灵敏度传播**

这是 SQP/IPM/SLP 的核心工具，见 §2.8.2。

**(c) `TrapezoidalIntegration.h`** —— 梯形法则，用于沿轨迹积分代价：

$$
\int_{t_0}^{t_N} l\,\mathrm{d}t \approx \sum_{k=0}^{N-1}\frac{t_{k+1}-t_k}{2}\big(l_k+l_{k+1}\big)
$$

`ocs2_oc/rollout/PerformanceIndicesRollout.h` 用它算 `PerformanceIndex`。

### 2.8.2 ⭐ 灵敏度离散化的推导

**文件**：`ocs2_core/src/integration/SensitivityIntegratorImpl.cpp`

多重打靶需要把连续动力学 $\dot x = f(t,x,u)$ 变成离散映射

$$
x_{k+1}=F(t_k,x_k,u_k,\Delta t)
$$

并给出其 Jacobian $A_k=\partial F/\partial x_k$、$B_k=\partial F/\partial u_k$。
朴素做法是对整个 RK 步做自动微分，但那会重复计算。OCS2 的做法是
**把灵敏度方程与状态方程一起离散化**。

#### 连续灵敏度方程

定义状态转移灵敏度 $S_x(t)=\partial x(t)/\partial x_k$、输入灵敏度 $S_u(t)=\partial x(t)/\partial u_k$。
对 $\dot x = f(t,x,u)$ 两边关于 $x_k$、$u_k$ 求导（$u$ 在区间内为常值，零阶保持）：

$$
\dot S_x = \frac{\partial f}{\partial x}S_x,\qquad S_x(t_k)=I
$$
$$
\dot S_u = \frac{\partial f}{\partial x}S_u + \frac{\partial f}{\partial u},\qquad S_u(t_k)=0
$$

#### 用同一个 RK 格式积分

对增广系统 $(x, S_x, S_u)$ 施加 RK4，各级 $k_i$ 的**状态部分**照常：

$$
\begin{aligned}
\kappa_1 &= f(t,\;x,\;u)\\
\kappa_2 &= f(t+\tfrac{\Delta t}{2},\;x+\tfrac{\Delta t}{2}\kappa_1,\;u)\\
\kappa_3 &= f(t+\tfrac{\Delta t}{2},\;x+\tfrac{\Delta t}{2}\kappa_2,\;u)\\
\kappa_4 &= f(t+\Delta t,\;x+\Delta t\,\kappa_3,\;u)
\end{aligned}
$$

**灵敏度部分**：记 $\mathcal{A}_i=\partial f/\partial x$、$\mathcal{B}_i=\partial f/\partial u$
在第 $i$ 级的取值点处的值。由链式法则，第 $i$ 级 $\kappa_i$ 对 $x_k$ 的导数为

$$
\frac{\partial \kappa_1}{\partial x_k}=\mathcal{A}_1,\qquad
\frac{\partial \kappa_2}{\partial x_k}=\mathcal{A}_2\Big(I+\tfrac{\Delta t}{2}\frac{\partial\kappa_1}{\partial x_k}\Big),\ \dots
$$

代码里做了一个**变量复用技巧**：把 `k_i.dfdx` 原地累加成 $\partial\kappa_i/\partial x_k$。
对照 `SensitivityIntegratorImpl.cpp:155-160`：

```cpp
matrix_t tmp = dt_halve * k2.dfdx * k1.dfdx;  k2.dfdx += tmp;   // A2 + (dt/2) A2 A1
tmp.noalias() = dt_halve * k3.dfdx * k2.dfdx; k3.dfdx += tmp;
tmp.noalias() = dt * k4.dfdx * k3.dfdx;       k4.dfdx += tmp;
```

即递推

$$
\hat{\mathcal{A}}_1=\mathcal{A}_1,\quad
\hat{\mathcal{A}}_2=\mathcal{A}_2+\tfrac{\Delta t}{2}\mathcal{A}_2\hat{\mathcal{A}}_1,\quad
\hat{\mathcal{A}}_3=\mathcal{A}_3+\tfrac{\Delta t}{2}\mathcal{A}_3\hat{\mathcal{A}}_2,\quad
\hat{\mathcal{A}}_4=\mathcal{A}_4+\Delta t\,\mathcal{A}_4\hat{\mathcal{A}}_3
$$

输入灵敏度同理（`:148-150`）：

$$
\hat{\mathcal{B}}_1=\mathcal{B}_1,\quad
\hat{\mathcal{B}}_2=\mathcal{B}_2+\tfrac{\Delta t}{2}\mathcal{A}_2\hat{\mathcal{B}}_1,\quad
\hat{\mathcal{B}}_3=\mathcal{B}_3+\tfrac{\Delta t}{2}\mathcal{A}_3\hat{\mathcal{B}}_2,\quad
\hat{\mathcal{B}}_4=\mathcal{B}_4+\Delta t\,\mathcal{A}_4\hat{\mathcal{B}}_3
$$

最后按 RK4 权重 $(\frac16,\frac13,\frac13,\frac16)$ 组合（`:164-167`）：

$$
\boxed{
\begin{aligned}
A_k &= I+\Delta t\Big(\tfrac16\hat{\mathcal{A}}_1+\tfrac13\hat{\mathcal{A}}_2+\tfrac13\hat{\mathcal{A}}_3+\tfrac16\hat{\mathcal{A}}_4\Big)\\
B_k &= \Delta t\Big(\tfrac16\hat{\mathcal{B}}_1+\tfrac13\hat{\mathcal{B}}_2+\tfrac13\hat{\mathcal{B}}_3+\tfrac16\hat{\mathcal{B}}_4\Big)\\
b_k &= x_k+\Delta t\Big(\tfrac16\kappa_1+\tfrac13\kappa_2+\tfrac13\kappa_3+\tfrac16\kappa_4\Big)
\end{aligned}}
$$

代码里 `k1.dfdx.diagonal().array() += 1.0;` 就是加上那个 $I$。

**RK2（中点法）** 与 **欧拉** 是同样推导的降阶版本。
`selectDynamicsSensitivityDiscretization()`（`SensitivityIntegrator.cpp:57`）
按 `SensitivityIntegratorType` 返回对应函数指针。

**为什么不直接自动微分整个 RK 步？** 因为那样需要对 $f$ 求 $O(\text{stage}^2)$ 次导数，
而上面的做法每级只求一次 $\partial f/\partial x$、$\partial f/\partial u$，
再做矩阵乘法组合，成本低得多。

### 2.8.3 事件处理

- `SystemEventHandler`：时间触发事件（预先知道 $t_i$）
- `StateTriggeredEventHandler`：状态触发事件（保护面 $g(t,x)$ 零穿越）

积分器在检测到事件时抛出特定异常，rollout 捕获后做跳变、再从新初值继续积分。
`killIntegration_` 标志用于线搜索并行时提前中止无用的 rollout
（`ocs2_oc/rollout/StateTriggeredRollout.h:64`）。

---

## 2.9 `control/`：控制器表示

```
ControllerBase
 ├─ FeedforwardController        u(t) = u_ff(t)          （插值查表）
 ├─ LinearController            u(t) = K(t)·x + b(t)     （时变线性反馈）
 └─ StateBasedLinearController  按状态而非时间索引增益
```

`LinearController` 是 DDP 的原生输出，字段：

```cpp
matrix_array_t gainArray_;       // K_k
vector_array_t biasArray_;       // b_k
vector_array_t deltaBiasArray_;  // Δb_k  ← 线搜索用
```

**`deltaBiasArray_` 的作用**：线搜索时新控制器为

$$
u = K_k\,x + b_k + \alpha\,\Delta b_k
$$

即**只对前馈项做步长缩放，反馈增益全量接受**。
这是 DDP 的标准做法（Todorov & Li 2005），因为反馈增益的改变不影响标称轨迹上的一阶最优性，
只影响偏离时的修正行为，缩放它没有意义。

实现见 `incrementController()`（`ocs2_ddp/src/DDP_HelperFunctions.cpp`）。

`ControllerAdjustmentBase` 用于 mode schedule 变化时调整控制器时间轴。

---

## 2.10 `initialization/`：初始化器

当没有可用控制器时（第一次 MPC、或 horizon 延长后的新区间），
需要给出一个"操作点"来填充轨迹：

| 类 | 策略 |
|---|---|
| `DefaultInitializer` | $u=0$，$x_{k+1}=x_k$ |
| `OperatingPoints` | 用户给定的时间-状态-输入轨迹，按时间插值 |

四足机器人用自定义 `LeggedRobotInitializer`：
把机器人重量均分到支撑腿上作为初始接触力
（`ocs2_robotic_examples/ocs2_legged_robot/src/initialization/LeggedRobotInitializer.cpp`），
这比 $u=0$ 好得多（$u=0$ 意味着机器人自由落体）。

---

## 2.11 `reference/`：参考信号

### `ModeSchedule`

```cpp
struct ModeSchedule {
  scalar_array_t eventTimes;   // 大小 N-1
  size_array_t   modeSequence; // 大小 N
};
```

不变式：`eventTimes` 严格递增，`modeSequence.size() == eventTimes.size() + 1`。
`modeAtTime(t)` 用二分查找返回 $t$ 所在的模态。

四足的 mode 是接触状态的位掩码：`modeNumber2StanceLeg()` 把整数
解码成 4 个 bool（`ocs2_legged_robot/gait/MotionPhaseDefinition.h`）。

### `TargetTrajectories`

```cpp
struct TargetTrajectories {
  scalar_array_t timeTrajectory;
  vector_array_t stateTrajectory;
  vector_array_t inputTrajectory;
};
```

`getDesiredState(t)` / `getDesiredInput(t)` 做线性插值，超出范围时钳位到端点。

---

## 2.12 `PreComputation`：跨项共享计算

**文件**：`ocs2_core/include/ocs2_core/PreComputation.h`

**动机**：四足机器人的代价、约束、动力学**都需要正运动学**。
如果每一项各算一次 `pinocchio::forwardKinematics`，开销会翻好几倍。

**机制**：

```cpp
enum class Request { Dynamics=1, Cost=2, Constraint=4, SoftConstraint=8, Approximation=16 };
```

是**位掩码**，`RequestSet` 用 `operator+` 做并集。求解器在计算前调用：

```cpp
constexpr auto request = Request::Cost + Request::SoftConstraint
                       + Request::Constraint + Request::Dynamics + Request::Approximation;
preComputation.request(request, t, x, u);
```

用户的 `PreComputation` 子类据此决定算什么。之后每个代价/约束项的 getter
都收到**同一个** `PreComputation&`，用 `ocs2::cast<Derived>(preComp)` 取回结果。

三个回调对应三类时刻：`request()`（中间）、`requestPreJump()`（事件前）、`requestFinal()`（终端）。

实例：`LeggedRobotPreComputation`
（`ocs2_robotic_examples/ocs2_legged_robot/src/LeggedRobotPreComputation.cpp`）
预算正运动学 + 摆动腿参考轨迹。

**这也是为什么 `OptimalControlProblem` 必须每线程一份**——`PreComputation` 有状态。

---

## 2.13 `misc/`：杂项工具

| 文件 | 内容 |
|---|---|
| `LinearAlgebra.h/.cpp` | ⭐ PSD 修正、约束投影、UUT 分解，见 [04 章](04-ocs2-ddp.md) §4.4/§4.6 |
| `LinearInterpolation.h` | 轨迹插值，返回 `(index, alpha)` 对以便复用查找结果 |
| `Lookup.h` | 二分查找区间索引 |
| `Collection.h` | 具名容器基类（§2.3.4） |
| `Numerics.h` | `almost_eq`、`almost_ge` 等容差比较 |
| `LoadData.h`, `LoadStdVectorOfPair.h` | 从 `.info`（Boost property_tree）读配置 |
| `Benchmark.h` | `RepeatedTimer`：max/avg/last 计时，求解器用它输出性能报告 |
| `Log.h` | Boost.Log 封装 |
| `Display.h` | 容器的格式化打印 |
| `randomMatrices.h` | 生成随机 PSD 矩阵，单元测试用 |
| `LTI_Equations.h` | 线性时不变系统的解析解 |
| `LinearFunction.h` | 分段线性函数 |

### `LinearInterpolation` 的性能设计

```cpp
auto indexAlpha = LinearInterpolation::timeSegment(t, timeArray);      // 查一次
auto A = LinearInterpolation::interpolate(indexAlpha, dataArray, accessorA);
auto B = LinearInterpolation::interpolate(indexAlpha, dataArray, accessorB);
```

Riccati 积分里一个时刻要插值十几个量，共享 `indexAlpha` 避免重复二分查找。
`ModelDataLinearInterpolation.h` 定义了访问 `ModelData` 各字段的 accessor
（`model_data::dynamics_dfdx` 等），配合上面的模式使用。

---

## 2.14 `automatic_differentiation/`

| 文件 | 内容 |
|---|---|
| `Types.h` | `ad_base_t = CppAD::cg::CG<double>`、`ad_scalar_t = CppAD::AD<ad_base_t>` |
| `CppAdInterface.h` | 封装"记录→生成 C 代码→编译→加载"全流程 |
| `CppAdSparsity.h` | 稀疏模式工具，只生成非零导数条目 |
| `FiniteDifferenceMethods.h` | 有限差分，用于校验 AD 结果 |

`CppAdInterface` 的三个入口：

```cpp
void createModels(ApproximationOrder, bool verbose);   // 强制重新生成
void loadModels(bool verbose);                          // 只加载
void loadModelsIfAvailable(ApproximationOrder, bool);   // 有则加载、无则生成
```

默认编译选项 `{"-O3","-g","-march=native","-mtune=native","-ffast-math"}`
（`CppAdInterface.h:71,83`）。注意 `-march=native` 意味着**生成的库不可跨机器分发**。

---

## 2.15 `thread_support/`

| 文件 | 内容 |
|---|---|
| `ThreadPool.h` | 固定大小线程池，`runParallel(task, N)` 阻塞直到全部完成 |
| `Synchronized.h` | `Synchronized<T>`：值 + mutex，`lock()` 返回 RAII 句柄 |
| `BufferedValue.h` | 双缓冲，写线程写 buffer、读线程 swap |
| `SetThreadPriority.h` | 设置实时调度优先级（`SCHED_FIFO`） |
| `ExecuteAndSleep.h` | 固定周期执行（MPC 主循环限频用） |

`Synchronized` + `BufferedValue` 是 MPC/MRT 双线程通信的基础：
MPC 线程算完策略写入 buffer，控制线程在自己的节拍上 `updatePolicy()` 尝试 swap，
**用 `try_lock` 而非 `lock`**，保证控制线程永不阻塞
（`ocs2_mpc/src/MRT_BASE.cpp:157`）。

---

## 2.16 `loopshaping/`

单独成章：[09 · Loopshaping](09-loopshaping.md)。

---

## 2.17 `model_data/`

| 类型 | 内容 |
|---|---|
| `ModelData` | 单个时刻的完整 LQ 数据：`time`, `stateDim`, `inputDim`, `dynamics`, `dynamicsBias`, `dynamicsCovariance`, `cost`, `stateEqConstraint`, `stateInputEqConstraint` |
| `Metrics` | 单个时刻的评估值：`cost`, `dynamicsViolation`, 各类约束值与拉格朗日值 |
| `Multiplier` | `{scalar_t penalty; vector_t lagrangian;}` |
| `MultiplierCollection` | 按约束类别分组的乘子集合 |
| `ModelDataLinearInterpolation.h` | 对 `std::vector<ModelData>` 的字段级插值访问器 |

`ModelData` 是 DDP 的工作单元：LQ 近似阶段并行填充，
Riccati 阶段串行（分区并行）消费。

---

## 2.18 小结与阅读建议

`ocs2_core` 的设计可以概括为**三层抽象**：

1. **数据层**（`Types.h`、`model_data/`）：只有数值，无行为
2. **模型层**（`cost/`、`constraint/`、`dynamics/`）：纯虚接口 + 少量具体实现
3. **工具层**（`integration/`、`misc/`、`thread_support/`）：算法无关的数值/工程设施

读源码时的最短路径：
`Types.h` → `dynamics/SystemDynamicsBase.h` → `cost/StateInputCost.h` →
`constraint/StateInputConstraint.h` → `PreComputation.h` → `misc/Collection.h`。

理解了这 6 个文件，`ocs2_core` 的其余部分都是它们的具体化。

**下一章**：这些原语如何组装成一个完整的最优控制问题 → [03 ocs2_oc](03-ocs2-oc.md)
