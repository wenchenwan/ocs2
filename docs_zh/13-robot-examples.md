# 13 · 机器人示例逐个拆解

`ocs2_robotic_examples/` 下有 6 个完整示例，按复杂度递增。
每个都是"如何把一个机器人写成 `OptimalControlProblem`"的教科书。

**每个示例的目录结构**（以 cartpole 为例）：

```
ocs2_cartpole/          核心库：动力学、参数、Interface
  ├── include/ocs2_cartpole/
  │   ├── definitions.h              STATE_DIM / INPUT_DIM
  │   ├── CartPoleParameters.h       物理参数（从 .info 读）
  │   ├── CartPoleInterface.h        RobotInterface 实现
  │   └── dynamics/CartPoleSystemDynamics.h
  ├── src/
  ├── config/mpc/task.info           所有配置（代价矩阵、求解器设置…）
  └── test/
ocs2_cartpole_ros/      ROS 层：节点、可视化、launch
  ├── src/CartpoleDummyVisualization.cpp
  ├── src/CartPoleMpcNode.cpp
  ├── src/CartPoleDummyNode.cpp
  └── launch/cartpole.launch
```

**`RobotInterface` 基类**（`ocs2_robotic_tools/common/RobotInterface.h`）：

```cpp
class RobotInterface {
  virtual const OptimalControlProblem& getOptimalControlProblem() const = 0;
  virtual const Initializer& getInitializer() const = 0;
  virtual std::shared_ptr<ReferenceManagerInterface> getReferenceManagerPtr() const = 0;
};
```

---

## 13.1 双积分器（`ocs2_double_integrator`）

**最简单的示例，从这里开始读。**

$$
\dim x = 2,\qquad \dim u = 1
$$

$$
x=\begin{bmatrix}p\\ v\end{bmatrix},\qquad
\dot x = \begin{bmatrix}0&1\\0&0\end{bmatrix}x+\begin{bmatrix}0\\ 1/m\end{bmatrix}u
$$

即 $\ddot p = u/m$（牛顿第二定律）。

**实现**：直接用 `ocs2_core` 的 `LinearSystemDynamics`，
无需自己写动力学类。代价用 `QuadraticStateInputCost`。

**这个例子的价值**：
- 问题是**严格线性二次的**，最优解有解析形式（LQR）
- 可以用来**验证求解器正确性**：SLQ/SQP/IPM/SLP 应该都在**一次迭代**内收敛到同一个解
- `test/` 下有与解析解对比的单元测试

有 Python 绑定（`DoubleIntegratorPyBindings.h`），是学习 `ocs2_python_interface` 的最佳起点。

---

## 13.2 倒立摆小车（`ocs2_cartpole`）

$$
\dim x = 4,\qquad \dim u = 1
$$

$$
x=\begin{bmatrix}\theta\\ p\\ \dot\theta\\ \dot p\end{bmatrix}
$$

（$\theta$ 是摆角，$p$ 是小车位置；**注意顺序是角度在前**）

### ⭐ 动力学推导

**文件**：`ocs2_cartpole/include/ocs2_cartpole/dynamics/CartPoleSystemDynamics.h:52`

用拉格朗日方法。设小车质量 $m_c$、摆质量 $m_p$、摆半长 $l$（质心到铰点）、
摆关于铰点的转动惯量（含 Steiner 项）$I_\theta = I_{cm}+m_pl^2$。

**动能**：小车 $\frac12m_c\dot p^2$；摆的质心速度为
$(\dot p + l\dot\theta\cos\theta,\ -l\dot\theta\sin\theta)$，故

$$
T=\tfrac12m_c\dot p^2+\tfrac12m_p\big[(\dot p+l\dot\theta\cos\theta)^2+(l\dot\theta\sin\theta)^2\big]+\tfrac12I_{cm}\dot\theta^2
$$

展开并合并：

$$
T=\tfrac12(m_c+m_p)\dot p^2+m_pl\,\dot p\dot\theta\cos\theta+\tfrac12(I_{cm}+m_pl^2)\dot\theta^2
$$

**势能**：$V = m_pgl\cos\theta$（$\theta=0$ 为竖直向上）

拉格朗日方程 $\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial L}{\partial\dot q}-\frac{\partial L}{\partial q}=Q$
给出（略去 $\dot\theta^2$ 的科氏项按代码约定放到右端）：

$$
\boxed{\;
\underbrace{\begin{bmatrix}
I_\theta & m_pl\cos\theta\\
m_pl\cos\theta & m_c+m_p
\end{bmatrix}}_{M(\theta)}
\begin{bmatrix}\ddot\theta\\ \ddot p\end{bmatrix}
=
\underbrace{\begin{bmatrix}
m_pgl\sin\theta\\
u+m_pl\dot\theta^2\sin\theta
\end{bmatrix}}_{b(\theta,\dot\theta,u)}\;}
$$

**逐行对应代码**：

```cpp
const ad_scalar_t cosTheta = cos(state(0));
const ad_scalar_t sinTheta = sin(state(0));

Eigen::Matrix<ad_scalar_t, 2, 2> I;
I << param_.poleSteinerMoi_,                    param_.poleMass_ * param_.poleHalfLength_ * cosTheta,
     param_.poleMass_ * param_.poleHalfLength_ * cosTheta,  param_.cartMass_ + param_.poleMass_;

Eigen::Matrix<ad_scalar_t, 2, 1> rhs(
    param_.poleMass_ * param_.poleHalfLength_ * param_.gravity_ * sinTheta,
    input(0) + param_.poleMass_ * param_.poleHalfLength_ * pow(state(2), 2) * sinTheta);

ad_vector_t stateDerivative(STATE_DIM);
stateDerivative << state.tail<2>(), I.inverse() * rhs;
```

`poleSteinerMoi_` = $I_\theta = I_{cm}+m_pl^2$（Steiner = 平行轴定理）。

**注意 `I.inverse()`**：这是 $2\times2$ 矩阵，Eigen 用解析公式求逆，很快。
但因为是在 CppAD 表达式里，会被记录进计算图并生成对应的 C 代码。

### 为什么这是好例子

- **非线性但低维**：可以画出完整的相图，直观理解 DDP 的收敛过程
- **强非线性**：摆的旋转（swing-up）需要非平凡的最优策略
- **有输入约束**：`task.info` 里配了输入的软约束（力限制）
- **文档来源**：头文件注释给出了推导的参考链接

`ocs2_cartpole_ros/launch/cartpole.launch` 启动后可在 RViz 里看到 swing-up 过程。

---

## 13.3 四旋翼（`ocs2_quadrotor`）

$$
\dim x = 12,\qquad \dim u = 4
$$

$$
x=\begin{bmatrix}\mathbf{p}\ (3)\\ \boldsymbol{\theta}_{xyz}\ (3)\\ \dot{\mathbf{p}}\ (3)\\ \boldsymbol{\omega}_{body}\ (3)\end{bmatrix},
\qquad
u=\begin{bmatrix}F_z\\ M_x\\ M_y\\ M_z\end{bmatrix}
$$

**注意状态的混合表示**：位置/姿态用**欧拉角**（XYZ），
速度用**世界系线速度** + **机体系角速度**。这是航空领域的常规约定。

### 关键实现

**文件**：`ocs2_quadrotor/src/QuadrotorSystemDynamics.cpp:38`

```cpp
Eigen::Matrix<scalar_t,3,1> eulerAngle = state.segment<3>(3);
Eigen::Matrix<scalar_t,3,3> T = getMappingFromLocalAngularVelocityToEulerAnglesXyzDerivative<scalar_t>(eulerAngle);
Eigen::Matrix<scalar_t,3,1> eulerAngleDerivatives = T * state.segment<3>(9);   // θ̇ = T·ω_body
```

这里用的是 **XYZ 欧拉角 + 机体系角速度**的映射
（与 [12 章](12-robot-models.md) §12.1.2 的 ZYX + 世界系角速度不同，
但推导思路一致）。

**平移动力学**：

$$
m\ddot{\mathbf{p}} = R(\boldsymbol{\theta})\begin{bmatrix}0\\0\\F_z\end{bmatrix}
- \begin{bmatrix}0\\0\\mg\end{bmatrix}
$$

代码里 `stateDerivative(6) = Fz * t2 * t4;`（$t_2=1/m$，$t_4=\sin\theta$）等
是这个式子的展开——**由 MATLAB 符号工具生成**，
所以有 `t2, t3, ..., t13` 这样的公共子表达式变量。

**旋转动力学**（欧拉方程，机体系）：

$$
I\dot{\boldsymbol{\omega}} + \boldsymbol{\omega}\times(I\boldsymbol{\omega}) = \mathbf{M}
$$

四旋翼的惯量矩阵对角且 $I_{xx}=I_{yy}=$ `Thxxyy_`，$I_{zz}=$ `Thzz_`，
所以展开式里有 $I_{zz}^2$（`t12`）这样的项。

### 线性化

**文件**：`QuadrotorSystemDynamics.cpp:113`

这个例子**手写了解析导数**（而非用 CppAD），
用 `JacobianOfAngularVelocityMapping` 处理欧拉角映射的导数：

```cpp
jacobianOfAngularVelocityMapping_ = JacobianOfAngularVelocityMapping(eulerAngle, angularVelocity).transpose();
```

**这是学习"如何手写复杂系统解析导数"的好材料**——
但也说明了为什么大多数示例用 CppAD：手写容易出错且难维护。

---

## 13.4 Ballbot（`ocs2_ballbot`）

$$
\dim x = 10,\qquad \dim u = 3
$$

**Ballbot** = 站在一个球上的倒立摆机器人（ETH 的 Rezero）。
三个全向轮驱动球，机身靠倾斜保持平衡。

$$
x=\begin{bmatrix}\mathbf{q}\ (5)\\ \dot{\mathbf{q}}\ (5)\end{bmatrix},\qquad
\mathbf{q}=\begin{bmatrix}x_{ball}\\ y_{ball}\\ \psi\\ \theta_p\\ \phi\end{bmatrix}
$$

（球在地面的 2D 位置 + 机身的 ZYX 欧拉角）

$$
u = \begin{bmatrix}\tau_1\\ \tau_2\\ \tau_3\end{bmatrix}\quad\text{（三个全向轮的力矩）}
$$

### ⭐ 执行器映射矩阵

**文件**：`ocs2_ballbot/src/dynamics/BallbotSystemDynamics.cpp:18`

动力学形式：

$$
M(\mathbf{q})\ddot{\mathbf{q}} + h(\mathbf{q},\dot{\mathbf{q}}) = S^{\!\top}(\mathbf{q})\,\boldsymbol{\tau}
$$

其中 $S^\top\in\mathbb{R}^{5\times3}$ 是**执行器映射矩阵**，
把 3 个轮力矩映射到 5 个广义力。

代码里 `S_transposed` 的每个元素都是显式的三角函数表达式，例如：

```cpp
const ad_scalar_t c1 = (sqrt_2 * ballRadius_) / (4.0 * wheelRadius_);
S_transposed(3, 0) = 2.0 * c1 * sroll;
S_transposed(3, 1) = c1 * (2.0*sroll + sqrt_3*croll);
S_transposed(3, 2) = c1 * (2.0*sroll - sqrt_3*croll);
S_transposed(4, 0) =  2.0 * c1;
S_transposed(4, 1) = -c1;
S_transposed(4, 2) = -c1;
```

**$\sqrt{2}$、$\sqrt{3}$ 的来源**：三个全向轮**均匀分布在球面上**，
相邻夹角 120°（$\cos120°=-\frac12$、$\sin120°=\frac{\sqrt3}{2}$ → $\sqrt3$），
且轮轴与竖直方向成特定倾角（→ $\sqrt2$）。
这些常数编码了机械设计的几何。

**$\frac{R_{ball}}{R_{wheel}}$**：轮与球的传动比。

### 刚体动力学：RobCoGen 生成

`generated/` 目录下的 20+ 个文件是 **RobCoGen** 生成的：

```
inertia_properties.h/.impl.h    连杆惯性参数
transforms.h/.impl.h            坐标变换
jacobians.h/.impl.h             Jacobian
forward_dynamics.h/.impl.h      前向动力学（ABA 算法）
inverse_dynamics.h/.impl.h      逆动力学（RNEA 算法）
jsim.h/.impl.h                  关节空间惯量矩阵（M）
```

**RobCoGen**（Frigerio, Buchli, Caldwell 2016）从机器人描述文件
生成**针对特定机器人优化的 C++ 代码**：
展开所有循环、消除零元素运算、内联常数。
比通用库（Pinocchio/RBDL）快 2~5 倍，代价是不灵活。

调用（`:68-75`）：

```cpp
using trait_t = typename iit::rbd::tpl::TraitSelector<ad_scalar_t>::Trait;
iit::Ballbot::dyn::tpl::InertiaProperties<trait_t> inertias;
iit::Ballbot::tpl::MotionTransforms<trait_t> transforms;
iit::Ballbot::dyn::tpl::ForwardDynamics<trait_t> forward_dyn(inertias, transforms);
forward_dyn.fd(qdd, state.head<5>(), state.tail<5>(), new_input);
```

**注意 `trait_t`**：RobCoGen 的模板化设计让同一份代码可以用 `double` 或 `ad_scalar_t` 实例化，
从而**整个前向动力学都可以被 CppAD 自动微分** ✓

这是 OCS2 里唯一用 RobCoGen 的例子（历史遗留，其余都用 Pinocchio）。

**参考**：Minniti, Farshidian, Grandia, Hutter (2019),
*Whole-Body MPC for a Dynamically Stable Mobile Manipulator*, RA-L。

---

## 13.5 移动机械臂（`ocs2_mobile_manipulator`）

**最"工程化"的示例**：支持 4 种底盘类型、任意 URDF、自碰撞、末端约束。

### 13.5.1 四种模型

**文件**：`ocs2_mobile_manipulator/include/ocs2_mobile_manipulator/ManipulatorModelInfo.h:42`

```cpp
enum class ManipulatorModelType {
  DefaultManipulator = 0,                   // 固定基机械臂（直接用 URDF）
  WheelBasedMobileManipulator = 1,          // + 可驱动的 XY-Yaw（差速/阿克曼底盘）
  FloatingArmManipulator = 2,               // + 不可驱动的 XYZ-RPY（虚拟浮动基，用于规划）
  FullyActuatedFloatingArmManipulator = 3,  // + 可驱动的 XYZ-RPY（全向浮动，如无人机机械臂）
};
```

对应四个动力学类：

| 类 | 状态 | 输入 | 动力学 |
|---|---|---|---|
| `DefaultManipulatorDynamics` | $q_{arm}$ | $\dot q_{arm}$ | $\dot x = u$（纯积分器） |
| `WheelBasedMobileManipulatorDynamics` | $[x,y,\theta,q_{arm}]$ | $[v,\omega,\dot q_{arm}]$ | ⭐ 非完整约束 |
| `FloatingArmManipulatorDynamics` | $q_{arm}$（基座固定为参数） | $\dot q_{arm}$ | $\dot x=u$ |
| `FullyActuatedFloatingArmManipulatorDynamics` | $[\mathbf{p},\boldsymbol{\theta},q_{arm}]$ | 全部速度 | $\dot x=u$ |

### ⭐ 13.5.2 非完整底盘的动力学

**文件**：`src/dynamics/WheelBasedMobileManipulatorDynamics.cpp:48`

```cpp
ad_vector_t dxdt(info_.stateDim);
const auto theta = state(2);
const auto v = input(0);      // 机体系前进速度
dxdt << cos(theta) * v, sin(theta) * v, input(1), input.tail(info_.armDim);
```

即**独轮车模型**：

$$
\begin{bmatrix}\dot x\\ \dot y\\ \dot\theta\\ \dot q_{arm}\end{bmatrix}
=\begin{bmatrix}v\cos\theta\\ v\sin\theta\\ \omega\\ \dot q_{arm}\end{bmatrix}
$$

**非完整约束**：$\dot x\sin\theta-\dot y\cos\theta=0$（不能横向移动）。

注意这个约束**没有被显式写成约束**，而是**隐含在动力学参数化里**——
因为输入只有 $(v,\omega)$，横向速度物理上不可能产生。
这比"给全向模型加约束"高效得多（少一个状态-输入约束，少一次投影）。

**代价**：底盘的机动性受限，MPC 需要规划出"倒车-转向"这类机动。
Giftthaler et al. (2017), *Efficient Kinematic Planning for Mobile Manipulators
with Non-Holonomic Constraints Using Optimal Control* 讨论了这个问题。

### 13.5.3 约束与代价

**文件**：`src/MobileManipulatorInterface.cpp`

```cpp
problem.stateSoftConstraintPtr->add("endEffector", getEndEffectorConstraint(...));      // 末端位姿跟踪
problem.finalSoftConstraintPtr->add("finalEndEffector", getEndEffectorConstraint(...)); // 终端位姿
problem.stateSoftConstraintPtr->add("selfCollision", getSelfCollisionConstraint(...));  // 自碰撞
problem.softConstraintPtr->add("jointLimits", getJointLimitSoftConstraint(...));        // 关节限位
problem.costPtr->add("inputCost", getQuadraticInputCost(...));                          // 输入正则
```

**`EndEffectorConstraint`**（`constraint/EndEffectorConstraint.h`）：

$$
h(x)=\begin{bmatrix}\mathbf{p}_{ee}(x)-\mathbf{p}_{ref}\\
\text{quaternionDistance}\big(\mathbf{q}_{ee}(x),\mathbf{q}_{ref}\big)\end{bmatrix}\in\mathbb{R}^6
$$

位置误差用欧氏距离，姿态误差用四元数距离（见 [12 章](12-robot-models.md) §12.1.5）。
配合 `StateSoftConstraint` + 二次罚，等价于一个位姿跟踪代价。

**`MobileManipulatorPinocchioMapping`**：把 OCS2 状态映射到 Pinocchio 的 $(q,v)$。
四种模型的映射不同（有无浮动基、浮动基是否可驱动），
这个类把差异隔离在一处。

### 13.5.4 支持的机器人

`config/` 下有 7 个机器人的配置：
Franka Emika Panda、Kinova J2N6/J2N7、Mabi Mobile、PR2、Ridgeback-UR5、
以及一个通用的 dummy。

**只需换 URDF + 改 `task.info`** 就能支持新机器人——
这是 `ocs2_mobile_manipulator` 设计得好的地方。

---

## 13.6 ⭐ 四足机器人（`ocs2_legged_robot`）

**最完整的示例**，也是 OCS2 的旗舰应用。以 ANYmal 为例：

$$
\dim x = 24,\qquad \dim u = 24
$$

$$
x = \begin{bmatrix}\mathbf{h}/m\ (6)\\ \mathbf{p}_{base}\ (3)\\ \boldsymbol{\theta}_{base}^{zyx}\ (3)\\ \mathbf{q}_j\ (12)\end{bmatrix},
\qquad
u = \begin{bmatrix}\mathbf{f}_{LF},\mathbf{f}_{RF},\mathbf{f}_{LH},\mathbf{f}_{RH}\ (12)\\ \dot{\mathbf{q}}_j\ (12)\end{bmatrix}
$$

（质心动力学，见 [12 章](12-robot-models.md) §12.3）

### 13.6.1 步态系统

#### `Gait`

**文件**：`gait/Gait.h`

```cpp
struct Gait {
  scalar_t duration;                    // 一个步态周期的时长
  std::vector<scalar_t> eventPhases;    // (0,1) 内的相位切换点，大小 N-1
  std::vector<size_t> modeSequence;     // 模态序列，大小 N
};
```

**用"相位"而非"时间"参数化**，使得步态可以任意缩放周期。

辅助函数：

```cpp
scalar_t wrapPhase(scalar_t phase);                  // 归一化到 [0,1)
int      getModeIndexFromPhase(scalar_t, const Gait&);
size_t   getModeFromPhase(scalar_t, const Gait&);
scalar_t timeLeftInGait(scalar_t phase, const Gait&);
scalar_t timeLeftInMode(scalar_t phase, const Gait&);
bool     isValidGait(const Gait&);                   // 检查不变式
```

#### 典型步态

| 步态 | `modeSequence` | 说明 |
|---|---|---|
| **stance** | `[15]` | 四腿全支撑（15 = 0b1111） |
| **trot** | `[9, 6]` | 对角腿交替（9=0b1001=LF+RH，6=0b0110=RF+LH） |
| **trot with flight** | `[9, 0, 6, 0]` | 中间插入腾空相（0 = 无腿触地） |
| **pace** | `[12, 3]` | 同侧腿交替 |
| **static walk** | `[14, 13, 11, 7]` | 每次抬一条腿，始终三腿支撑 |

配置在 `config/command/gait.info`。

#### `GaitSchedule`

```cpp
class GaitSchedule {
  ModeSchedule getModeSchedule(scalar_t lowerBoundTime, scalar_t upperBoundTime);
  void insertModeSequenceTemplate(const ModeSequenceTemplate&, scalar_t startTime, scalar_t finalTime);
  void setGaitAfterTime(const ModeSequenceTemplate&, scalar_t startTime);
};
```

维护一个**时间轴上的模态序列**，MPC 每次向它请求 $[t, t+T]$ 内的 `ModeSchedule`。
用户切换步态时调 `insertModeSequenceTemplate`，
从指定时刻开始把新步态循环铺开。

#### `LegLogic`

从 `ModeSchedule` 提取每条腿的**接触序列**与**相位时间**：

```cpp
feet_array_t<std::vector<ContactTiming>> extractContactTimings(const ModeSchedule&);
scalar_t getTimeOfNextLiftOff(scalar_t time, const std::vector<ContactTiming>&);
scalar_t getTimeOfNextTouchDown(scalar_t time, const std::vector<ContactTiming>&);
```

摆动腿轨迹规划需要知道"这条腿什么时候抬起、什么时候落下"。

### 13.6.2 ⭐ 摆动腿轨迹规划

#### `CubicSpline`

**文件**：`foot_planner/CubicSpline.cpp:36`

给定起点 $(t_0,p_0,v_0)$ 与终点 $(t_1,p_1,v_1)$，构造三次多项式。

用**归一化时间** $\tau=\dfrac{t-t_0}{\Delta t}\in[0,1]$：

$$
p(\tau)=c_3\tau^3+c_2\tau^2+c_1\tau+c_0
$$

四个边界条件：

$$
p(0)=p_0,\quad p(1)=p_1,\quad
\frac{\mathrm{d}p}{\mathrm{d}t}\Big|_0=v_0,\quad \frac{\mathrm{d}p}{\mathrm{d}t}\Big|_1=v_1
$$

注意 $\frac{\mathrm{d}p}{\mathrm{d}t}=\frac{1}{\Delta t}\frac{\mathrm{d}p}{\mathrm{d}\tau}$，故

$$
p'(0)=v_0\Delta t,\qquad p'(1)=v_1\Delta t
$$

解得（$\Delta p = p_1-p_0$、$\Delta v = v_1-v_0$）：

$$
\boxed{
\begin{aligned}
c_0 &= p_0\\
c_1 &= v_0\,\Delta t\\
c_2 &= 3\Delta p - (2v_0+v_1)\Delta t\\
c_3 &= -2\Delta p + (v_0+v_1)\Delta t
\end{aligned}}
$$

**验证代码**（`CubicSpline.cpp:44-52`）：

```cpp
scalar_t dp = end.position - start.position;    // Δp
scalar_t dv = end.velocity - start.velocity;    // Δv
dc0_ = 0.0;
dc1_ = start.velocity;                          // v0
dc2_ = -(3.0 * start.velocity + dv);            // −(3v0 + Δv) = −(2v0 + v1)
dc3_ = (2.0 * start.velocity + dv);             // 2v0 + Δv = v0 + v1
c0_ = dc0_ * dt_ + start.position;              // p0                     ✓
c1_ = dc1_ * dt_;                               // v0·Δt                  ✓
c2_ = dc2_ * dt_ + 3.0 * dp;                    // 3Δp − (2v0+v1)Δt       ✓
c3_ = dc3_ * dt_ - 2.0 * dp;                    // −2Δp + (v0+v1)Δt       ✓
```

完全一致。

**`dc*_` 的用途**：它们是 $c_i$ 对 $\Delta t$ 的偏导数：

$$
\frac{\partial c_0}{\partial\Delta t}=0,\quad
\frac{\partial c_1}{\partial\Delta t}=v_0,\quad
\frac{\partial c_2}{\partial\Delta t}=-(2v_0+v_1),\quad
\frac{\partial c_3}{\partial\Delta t}=v_0+v_1
$$

用于 `startTimeDerivative` / `finalTimeDerivative`——
**对摆动时长求导**，这在优化步态时序（switching time optimization）时需要。

`startTimeDerivative(t)` 的推导（`:76-81`）：
$t_0$ 变化时，$\Delta t = t_1-t_0$ 减小，且归一化时间 $\tau=\frac{t-t_0}{\Delta t}$ 也变：

$$
\frac{\partial p}{\partial t_0}
=\underbrace{\frac{\partial p}{\partial\tau}\frac{\partial\tau}{\partial t_0}}_{\text{通过 }\tau}
+\underbrace{\sum_i\frac{\partial p}{\partial c_i}\frac{\partial c_i}{\partial\Delta t}\frac{\partial\Delta t}{\partial t_0}}_{\text{通过系数}}
$$

$\frac{\partial\Delta t}{\partial t_0}=-1$，$\frac{\partial\tau}{\partial t_0}=-\frac{t_1-t}{\Delta t^2}$，故

```cpp
scalar_t dCoff = -(dc3_*tn*tn*tn + dc2_*tn*tn + dc1_*tn + dc0_);   // 系数项
scalar_t dTn = -(t1_ - t) / (dt_ * dt_);                            // ∂τ/∂t0
return velocity(t) * dt_ * dTn + dCoff;                             // 注意 velocity 已含 1/dt
```

#### `SplineCpg`

**文件**：`foot_planner/SplineCpg.h`

**CPG** = Central Pattern Generator。把摆动分成**两段三次样条**：

```
        midHeight
           ╱╲
          ╱  ╲
   ──────╱    ╲──────
  liftOff  midTime  touchDown
```

```cpp
SplineCpg(CubicSpline::Node liftOff, scalar_t midHeight, CubicSpline::Node touchDown);
```

中点 $(t_{mid}, h_{mid}, 0)$ 的速度设为 0（最高点）。
两段各自是三次样条，整体 $C^1$ 连续。

**为什么不用一段五次样条？**
一段五次可以满足 6 个条件（起终点的位置/速度/加速度），
但**无法直接指定中间高度**。两段三次给了对轨迹形状的直接控制，
这在越障时很重要。

#### `SwingTrajectoryPlanner`

**文件**：`foot_planner/SwingTrajectoryPlanner.h:40`

```cpp
struct Config {
  scalar_t liftOffVelocity = 0.0;
  scalar_t touchDownVelocity = 0.0;
  scalar_t swingHeight = 0.1;
  scalar_t swingTimeScale = 0.15;   // 短于此时长的摆动相，高度和速度按比例缩小
};

void update(const ModeSchedule&, scalar_t terrainHeight);
void update(const ModeSchedule&, const feet_array_t<scalar_array_t>& liftOffHeightSequence,
            const feet_array_t<scalar_array_t>& touchDownHeightSequence);   // ← 感知版本

scalar_t getZvelocityConstraint(size_t leg, scalar_t time) const;
scalar_t getZpositionConstraint(size_t leg, scalar_t time) const;
```

流程：
1. 从 `ModeSchedule` 提取每条腿的接触时序（`LegLogic`）
2. 对每个摆动相构造一个 `SplineCpg`
3. MPC 每步查询 $z_{ref}(t)$ 与 $\dot z_{ref}(t)$ 作为约束目标

**`swingTimeScale` 的作用**：快速步态下摆动相可能只有 0.1s，
此时 0.1m 的抬腿高度需要极大的垂直加速度，物理上不可行。
该参数让高度随摆动时长线性缩放：

$$
h_{\text{actual}} = h_{\text{nominal}}\cdot\min\left(1,\ \frac{T_{swing}}{T_{scale}}\right)
$$

**感知版本**：接受每次抬脚/落脚的**地形高度序列**，
使摆动轨迹贴合不平地面（`ocs2_perceptive_anymal` 用它）。

### 13.6.3 约束体系

**文件**：`src/LeggedRobotInterface.cpp:181-190`

对每条腿 $i$：

```cpp
// 不等式约束（可选：硬约束或软约束）
if (useHardFrictionConeConstraint) {
  problemPtr_->inequalityConstraintPtr->add(footName + "_frictionCone", getFrictionConeConstraint(i, μ));
} else {
  problemPtr_->softConstraintPtr->add(footName + "_frictionCone",
                                      getFrictionConeSoftConstraint(i, μ, barrierPenaltyConfig));
}
// 硬等式约束（用零空间投影消去）
problemPtr_->equalityConstraintPtr->add(footName + "_zeroForce",      getZeroForceConstraint(i));
problemPtr_->equalityConstraintPtr->add(footName + "_zeroVelocity",   getZeroVelocityConstraint(...));
problemPtr_->equalityConstraintPtr->add(footName + "_normalVelocity", getNormalVelocityConstraint(...));
```

#### ⭐ 摩擦锥约束

**文件**：`constraint/FrictionConeConstraint.h` / `.cpp`

$$
\boxed{\;h(u)=\mu\big(F_z + F_{grip}\big)-\sqrt{F_x^2+F_y^2+\varepsilon}\ \ge\ 0\;}
$$

**三个设计要点**：

**(1) 正则项 $\varepsilon$**（`regularization`，默认 25.0）

标准摩擦锥 $\mu F_z\ge\sqrt{F_x^2+F_y^2}$ 在 $F_x=F_y=0$ 处**梯度不存在**
（$\sqrt{\cdot}$ 在 0 处不可导），Hessian 发散。

加 $\varepsilon$ 后：
- 处处可导 ✓
- 零穿越点从 $F_z=0$ 变成 $F_z=\frac{\sqrt\varepsilon}{\mu}$，
  即**提供了一个抛物线形的安全裕度**

代码注释明确指出了后一点：
> "when Fx = Fy = 0 the constraint zero-crossing will be at Fz = 1/frictionCoefficient * sqrt(regularization) instead of Fz = 0"

$\varepsilon=25$、$\mu=0.7$ 时零穿越在 $F_z=\frac{5}{0.7}\approx7.1$N——
即要求法向力至少 7N，防止"零力接触"这种病态解。

**(2) 抓握力 $F_{grip}$**

允许"拉"地面（磁吸、爪、真空吸盘）。$F_{grip}>0$ 把锥的顶点下移。

**(3) 地形法向**

`setSurfaceNormalInWorld(n)` 设置局部地形法向，
约束在**局部切平面坐标系**下表达：

```cpp
const vector3_t localForce = t_R_w * forcesInWorldFrame;   // 世界系 → 地形系
return coneConstraint(localForce);
```

`t_R_w` 是从世界系到地形系的旋转。

##### 导数推导

**局部力的一阶导**（`:129`）：

$$
\frac{\partial h}{\partial F_x}=-\frac{F_x}{\sqrt{F_x^2+F_y^2+\varepsilon}},\quad
\frac{\partial h}{\partial F_y}=-\frac{F_y}{\sqrt{\cdot}},\quad
\frac{\partial h}{\partial F_z}=\mu
$$

**二阶导**：记 $s=F_x^2+F_y^2+\varepsilon$，$\|F_t\|=\sqrt s$。

$$
\frac{\partial^2h}{\partial F_x^2}
=-\frac{\partial}{\partial F_x}\left(\frac{F_x}{\sqrt s}\right)
=-\frac{\sqrt s - F_x\cdot\frac{F_x}{\sqrt s}}{s}
=-\frac{s-F_x^2}{s^{3/2}}
=-\frac{F_y^2+\varepsilon}{s^{3/2}}
$$

$$
\frac{\partial^2h}{\partial F_x\partial F_y}
=-\frac{\partial}{\partial F_y}\left(\frac{F_x}{\sqrt s}\right)
= F_x\cdot\frac{F_y}{s^{3/2}}=\frac{F_xF_y}{s^{3/2}}
$$

**完全对应代码**（`:141-149`）：

```cpp
coneDerivatives.d2Cone_dF2(0,0) = -(F_y_square + regularization) / F_tangent_square_pow32;
coneDerivatives.d2Cone_dF2(0,1) = localForces.x() * localForces.y() / F_tangent_square_pow32;
coneDerivatives.d2Cone_dF2(1,1) = -(F_x_square + regularization) / F_tangent_square_pow32;
// 其余为 0（因为 h 对 F_z 是线性的）
```

**Hessian 是负定的**（$-\frac{F_y^2+\varepsilon}{s^{3/2}}<0$），
说明 $h$ 是**凹函数** → $h\ge0$ 定义了一个**凸集** ✓
（这与摩擦锥是二阶锥的事实一致）

##### 世界系导数

$$
\frac{\partial h}{\partial u}=\frac{\partial h}{\partial F^{local}}\cdot{}^tR_w,
\qquad
\frac{\partial^2h}{\partial u^2}={}^tR_w^{\!\top}\frac{\partial^2h}{\partial (F^{local})^2}\,{}^tR_w
$$

（`:167-178` 的 `computeConeConstraintDerivatives`）

`hessianDiagonalShift`（默认 1e-6）在需要严格凸二次近似时加到对角上。

#### 其他约束

| 约束 | 数学 | 作用于 |
|---|---|---|
| `ZeroForceConstraint` | $\mathbf{f}_i=0$ | **摆动腿**（3 行，直接约束输入分量，$D$ 是选择矩阵） |
| `ZeroVelocityConstraintCppAd` | $\mathbf{v}_{ee,i}(x,u)=0$ | **支撑腿**（3 行，末端速度为零） |
| `NormalVelocityConstraintCppAd` | $v_{ee,i}^z = \dot z_{ref}(t)$ | **摆动腿**（1 行，跟踪摆动轨迹） |
| `EndEffectorLinearConstraint` | $A\,\mathbf{v}_{ee}+b=0$ | 通用的末端线性约束基类 |

**这些都是状态-输入等式约束，会被零空间投影消去**（[04 章](04-ocs2-ddp.md) §4.4）。

四腿全支撑时：$4\times3=12$ 行零速度约束
→ 投影后输入维度从 24 降到 12。

**注意零力与零速度是互斥的**：`isActive(time)` 根据当前模态决定：

```cpp
bool ZeroForceConstraint::isActive(scalar_t time) const {
  return !referenceManagerPtr_->getContactFlags(time)[contactPointIndex_];   // 摆动腿才激活
}
bool ZeroVelocityConstraint::isActive(scalar_t time) const {
  return referenceManagerPtr_->getContactFlags(time)[contactPointIndex_];    // 支撑腿才激活
}
```

于是**约束的数量随模态变化**——这正是切换系统的特征。

### 13.6.4 代价

**文件**：`cost/LeggedRobotQuadraticTrackingCost.h`

$$
l = \tfrac12(x-x_{ref})^{\!\top}Q(x-x_{ref}) + \tfrac12(u-u_{ref})^{\!\top}R(u-u_{ref})
$$

$u_{ref}$ 的构造很有讲究：`LeggedRobotStateInputQuadraticCost::getStateInputDeviation`
（`cost/LeggedRobotQuadraticTrackingCost.h:61`）调用 `weightCompensatingInput`
（`common/utils.h:63`）：

$$
\mathbf{f}_{i,ref}^z=\frac{mg}{n_{\text{stance}}},\qquad \mathbf{f}_{i,ref}^{x,y}=0
$$

即**重力补偿输入**：把机器人重量均分到支撑腿。
这让"什么都不做"的代价接近 0，避免了 MPC 一开始就要"对抗"代价函数。
同一个函数也被 `LeggedRobotInitializer`（`initialization/LeggedRobotInitializer.cpp:58`）
用作初始输入猜测。

### 13.6.5 `SwitchedModelReferenceManager`

**文件**：`reference_manager/SwitchedModelReferenceManager.h`

```cpp
class SwitchedModelReferenceManager : public ReferenceManager {
  void setModeSchedule(const ModeSchedule&) override;
  contact_flag_t getContactFlags(scalar_t time) const;
  const std::shared_ptr<GaitSchedule>& getGaitSchedule() { return gaitSchedulePtr_; }
  const std::shared_ptr<SwingTrajectoryPlanner>& getSwingTrajectoryPlanner() { return swingTrajectoryPtr_; }
 protected:
  void modifyReferences(scalar_t initTime, scalar_t finalTime, const vector_t& initState,
                        TargetTrajectories&, ModeSchedule&) override;
};
```

`modifyReferences` 在每次求解前被调用：
1. 从 `GaitSchedule` 取 $[t_0,t_f]$ 的 `ModeSchedule`
2. 用它更新 `SwingTrajectoryPlanner`
3. 约束类通过 `getContactFlags(t)` 与 `getZpositionConstraint(leg,t)` 读取结果

**这是"步态规划器 → MPC"的完整接口**，也是理解四足 MPC 架构的关键。

### 13.6.6 `LeggedRobotPreComputation`

预计算并缓存：
- Pinocchio 的正运动学（`forwardKinematics` + `updateFramePlacements`）
- 各末端的 Jacobian
- 当前时刻的摆动腿参考（$z_{ref}$、$\dot z_{ref}$）

所有约束/代价项共享，避免重复计算。见 [02 章](02-ocs2-core.md) §2.12。

### 13.6.7 ROS 层

`ocs2_legged_robot_ros/` 提供：

| 节点 | 用途 |
|---|---|
| `legged_robot_ddp_mpc` / `legged_robot_sqp_mpc` / `legged_robot_ipm_mpc` | MPC 节点（三种求解器） |
| `legged_robot_dummy` | 无物理引擎的仿真 + RViz 可视化 |
| `legged_robot_target` | 键盘控制目标位姿 |
| `legged_robot_gait_command` | 键盘切换步态 |
| `legged_robot_pose_command` | 交互标记控制 |

`LeggedRobotVisualizer` 画出：机器人模型、接触力箭头、
支撑多边形、摆动腿轨迹、CoM 轨迹。

---

## 13.7 `ocs2_perceptive_anymal`

在 `ocs2_legged_robot` 基础上加入感知：

- 用 `ocs2_perceptive` 的距离场约束避障
- 用感知得到的地形高度更新 `SwingTrajectoryPlanner`
- `ocs2_anymal_*` 子包提供 ANYmal 特定的模型与配置

**参考**：Gaertner, Bjelonic, Farshidian, Hutter (2021),
*Collision-Free MPC for Legged Robots in Static and Dynamic Scenes*, ICRA。

---

## 13.8 `ocs2_raisim`：物理引擎桥接

| 包 | 内容 |
|---|---|
| `ocs2_raisim_core` | `RaisimRollout`（用 RaiSim 做 rollout）、状态转换接口 |
| `ocs2_raisim_ros` | ROS 集成与可视化 |
| `ocs2_legged_robot_raisim` | 四足在 RaiSim 中的仿真 |

**`RaisimRollout` 的意义**：用**真实物理引擎**代替 MPC 自己的模型做 rollout。
于是可以评估**模型误差**的影响——
MPC 用简化的质心模型规划，RaiSim 用全身刚体动力学 + 接触求解器仿真。

这是评估 MPC 鲁棒性的标准做法。`MRT_ROS_Dummy_Loop` 做不到这一点
（因为它用的就是 MPC 自己的模型）。

---

## 13.9 示例对比总表

| 示例 | $\dim x$ | $\dim u$ | 动力学 | 约束 | 切换 | 导数 |
|---|---:|---:|---|---|---|---|
| double integrator | 2 | 1 | 线性 | 无 | 无 | 解析 |
| cartpole | 4 | 1 | 非线性 | 输入软约束 | 无 | CppAD |
| quadrotor | 12 | 4 | 非线性 | 无 | 无 | **手写解析** |
| ballbot | 10 | 3 | RobCoGen RBD | 无 | 无 | CppAD |
| mobile manipulator | 可变 | 可变 | 运动学/非完整 | 自碰撞、末端、限位 | 无 | CppAD |
| **legged robot** | 24 | 24 | 质心动力学 | 摩擦锥、零力/零速、法向速度 | **有** | CppAD |

**学习路径建议**：
double integrator（30 分钟）→ cartpole（1 小时）→
mobile manipulator（3 小时）→ legged robot（1~2 天）。

---

## 13.10 参考文献

- **Minniti, Farshidian, Grandia, Hutter (2019)**, *Whole-Body MPC for a Dynamically Stable Mobile Manipulator*, RA-L 4(4) —— **Ballbot**
- **Giftthaler, Farshidian, Sandy, Stadelmann, Buchli (2017)**, *Efficient Kinematic Planning for Mobile Manipulators with Non-Holonomic Constraints Using Optimal Control*, ICRA —— **移动机械臂**
- **Farshidian, Neunert, Winkler, Rey, Buchli (2017)**, *An Efficient Optimal Planning and Control Framework for Quadrupedal Locomotion*, ICRA —— **四足**
- **Grandia, Farshidian, Ranftl, Hutter (2019)**, *Feedback MPC for Torque-Controlled Legged Robots*, IROS —— 摩擦锥与松弛障碍
- **Sleiman, Farshidian, Minniti, Hutter (2021)**, *A Unified MPC Framework for Whole-Body Dynamic Locomotion and Manipulation*, RA-L 6(3)
- **Gaertner, Bjelonic, Farshidian, Hutter (2021)**, *Collision-Free MPC for Legged Robots in Static and Dynamic Scenes*, ICRA
- **Frigerio, Buchli, Caldwell, Semini (2016)**, *RobCoGen: a code generator for efficient kinematics and dynamics of articulated robots*, JOSER —— **Ballbot 用的代码生成器**
- **Hwangbo, Lee, Hutter (2018)**, *Per-Contact Iteration Method for Solving Contact Dynamics*, RA-L —— **RaiSim** 的接触求解
- **Winkler, Bellicoso, Hutter, Buchli (2018)**, *Gait and Trajectory Optimization for Legged Systems through Phase-based End-Effector Parameterization*, RA-L —— 相位参数化的步态

**下一章**：MPC-Net → [14 MPC-Net](14-mpcnet.md)
