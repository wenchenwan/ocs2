# 12 · 机器人建模：Pinocchio 与质心动力学

从 URDF 到最优控制问题的桥梁。本章覆盖：

- `ocs2_robotic_tools`：旋转表示与角速度映射
- `ocs2_pinocchio_interface`：Pinocchio 封装与末端运动学
- `ocs2_centroidal_model`：⭐ 质心动力学（腿足机器人的核心）
- `ocs2_self_collision`：自碰撞约束
- `ocs2_sphere_approximation`：球近似
- `ocs2_perceptive`：距离场与感知约束

---

## 12.1 ⭐ 旋转表示与角速度映射

**目录**：`ocs2_robotic_tools/common/`

这一小节的数学是所有浮动基机器人建模的基础，**极易出错**，值得完整推导。

### 12.1.1 ZYX 欧拉角与旋转矩阵

**文件**：`RotationTransforms.h`

OCS2 用 **ZYX 内旋欧拉角**（yaw-pitch-roll，$\theta=[\psi,\theta_p,\phi]^\top$）：

$$
R(\psi,\theta_p,\phi)=R_z(\psi)\,R_y(\theta_p)\,R_x(\phi)
$$

$$
R_z=\begin{bmatrix}c_\psi&-s_\psi&0\\ s_\psi&c_\psi&0\\ 0&0&1\end{bmatrix},\quad
R_y=\begin{bmatrix}c_{\theta}&0&s_{\theta}\\ 0&1&0\\ -s_{\theta}&0&c_{\theta}\end{bmatrix},\quad
R_x=\begin{bmatrix}1&0&0\\ 0&c_\phi&-s_\phi\\ 0&s_\phi&c_\phi\end{bmatrix}
$$

主要函数：

```cpp
matrix3_t getRotationMatrixFromZyxEulerAngles(const vector3_t& eulerAngles);
vector3_t getEulerAnglesFromRotationMatrix(const matrix3_t& R);
void makeEulerAnglesUnique(vector3_t& eulerAngles);        // 归一化到主值区间
matrix3_t getMappingFromEulerAnglesZyxDerivativeToGlobalAngularVelocity(const vector3_t& eulerAngles);
```

### 12.1.2 ⭐ 欧拉角速率 → 角速度的映射

**这是整个建模里最容易搞错的地方**：$\dot\theta\ne\omega$。

角速度 $\omega$ 是三个连续旋转的角速度**在世界系下的叠加**：

$$
\omega = \dot\psi\,\mathbf{e}_z + \dot\theta_p\,R_z(\psi)\mathbf{e}_y
 + \dot\phi\,R_z(\psi)R_y(\theta_p)\mathbf{e}_x
$$

**推导**：
- 第一次旋转 $R_z(\psi)$ 绕**世界系** $z$ 轴 → 贡献 $\dot\psi\,\mathbf{e}_z$
- 第二次旋转 $R_y(\theta_p)$ 绕**已被 $R_z$ 转过的** $y$ 轴 → 该轴在世界系下是 $R_z\mathbf{e}_y$
- 第三次旋转 $R_x(\phi)$ 绕**已被 $R_zR_y$ 转过的** $x$ 轴 → 该轴是 $R_zR_y\mathbf{e}_x$

展开各项：

$$
R_z(\psi)\mathbf{e}_y=\begin{bmatrix}-s_\psi\\ c_\psi\\ 0\end{bmatrix},\qquad
R_z(\psi)R_y(\theta_p)\mathbf{e}_x=\begin{bmatrix}c_\psi c_\theta\\ s_\psi c_\theta\\ -s_\theta\end{bmatrix}
$$

于是

$$
\boxed{\;\omega = T(\theta)\,\dot\theta,\qquad
T(\theta)=\begin{bmatrix}
0 & -s_\psi & c_\psi c_{\theta}\\
0 & c_\psi & s_\psi c_{\theta}\\
1 & 0 & -s_{\theta}
\end{bmatrix}\;}
$$

（列的顺序对应 $\dot\theta=[\dot\psi,\dot\theta_p,\dot\phi]^\top$）

**逆映射**（`getEulerAnglesZyxDerivativesFromGlobalAngularVelocity`）：

$$
\dot\theta = T^{-1}(\theta)\,\omega
$$

$$
\det T = \cdots = -c_\theta
$$

**⚠️ 万向节死锁**：$\theta_p=\pm\pi/2$ 时 $c_\theta=0$，$T$ 奇异，$T^{-1}$ 不存在。
物理含义：俯仰 90° 时 yaw 与 roll 退化为同一个自由度。

这是 OCS2 用 ZYX 欧拉角的**根本限制**：
四足机器人俯仰不会到 90°，所以没问题；
但做空翻、体操动作的人形机器人必须换四元数或旋转矩阵表示。

### 12.1.3 ⭐ 映射的导数（Jacobian）

**文件**：`RotationDerivativesTransforms.h`、`AngularVelocityMapping.h`

线性化动力学需要 $\partial(T\dot\theta)/\partial\theta$。
`JacobianOfAngularVelocityMapping(eulerAngle, angularVelocity)` 返回

$$
\frac{\partial}{\partial\theta}\Big(T(\theta)\,\omega\Big)
\quad\text{（注意是对固定 }\omega\text{ 求 }\theta\text{ 的导数）}
$$

四旋翼的 `linearApproximation`（`QuadrotorSystemDynamics.cpp:118`）用它：

```cpp
jacobianOfAngularVelocityMapping_ = JacobianOfAngularVelocityMapping(eulerAngle, angularVelocity).transpose();
```

**其他函数**：
- `getGlobalAngularVelocityFromEulerAnglesZyxDerivatives`
- `getGlobalAngularAccelerationFromEulerAnglesZyxDerivatives`：
  $\dot\omega = T\ddot\theta + \dot T\dot\theta$，需要 $\dot T$
- `getMappingFromLocalAngularVelocityToEulerAnglesXyzDerivative`：
  XYZ 欧拉角版本（四旋翼用）

### 12.1.4 `SkewSymmetricMatrix.h`

$$
[\,a\,]_\times=\begin{bmatrix}0&-a_3&a_2\\ a_3&0&-a_1\\ -a_2&a_1&0\end{bmatrix},
\qquad [\,a\,]_\times b = a\times b
$$

在质心动力学的力矩项与 Jacobian 里大量出现。

### 12.1.5 四元数的处理

`RotationTransforms.h` 还提供：

```cpp
vector4_t quaternionDistance(const quaternion_t& q, const quaternion_t& qRef);
matrix_t  quaternionDistanceJacobian(const quaternion_t& q, const quaternion_t& qRef);
```

**四元数距离**用于末端姿态跟踪代价：由于 $q$ 与 $-q$ 表示同一姿态，
不能直接用 $\|q-q_{ref}\|$。OCS2 用

$$
d(q,q_{ref}) = \text{（}q_{ref}^{-1}\otimes q\text{ 的虚部）}\times\operatorname{sign}(\text{实部})
$$

这在 $\pm q$ 下一致，且在 $q=q_{ref}$ 处光滑。

`ocs2_mobile_manipulator` 的末端姿态约束用它。

---

## 12.2 `ocs2_pinocchio_interface`

### 12.2.1 `PinocchioInterface`

**文件**：`ocs2_pinocchio_interface/include/ocs2_pinocchio_interface/PinocchioInterface.h`

```cpp
template <typename SCALAR>
class PinocchioInterfaceTpl {
  const Model& getModel() const;
  Data& getData();
  PinocchioInterfaceCppAd toCppAd() const;      // ← 关键：标量类型转换
};
```

薄封装 Pinocchio 的 `Model`（静态：连杆、关节、惯性参数）
与 `Data`（动态：每次计算的缓存）。

**`toCppAd()`**：把 `double` 模型转成 `CppAD::AD<CG<double>>` 模型，
从而可以对整条 Pinocchio 计算链做自动微分。这是 OCS2 能"自动"得到
复杂机器人动力学导数的关键。

`urdf.h` 提供从 URDF 文件/字符串构造：

```cpp
PinocchioInterface getPinocchioInterfaceFromUrdfFile(const std::string& urdfFilePath);
PinocchioInterface getPinocchioInterfaceFromUrdfFile(const std::string& path, const JointModel& rootJoint);
```

`rootJoint` 决定浮动基类型：
- 不给 → 固定基（机械臂）
- `JointModelFreeFlyer` → 6DoF 浮动基（四足、人形）
- `JointModelComposite`（PX+PY+RZ）→ 平面移动基（轮式底盘）

### 12.2.2 `EndEffectorKinematics`

**文件**：`ocs2_robotic_tools/end_effector/EndEffectorKinematics.h`（接口）

```cpp
template <typename SCALAR>
class EndEffectorKinematics {
  virtual std::vector<vector3_t> getPosition(const vector_t& state) const = 0;
  virtual std::vector<vector3_t> getVelocity(const vector_t& state, const vector_t& input) const = 0;
  virtual std::vector<vector3_t> getOrientationError(const vector_t& state,
                                                     const std::vector<quaternion_t>& refs) const = 0;
  virtual std::vector<VectorFunctionLinearApproximation> getPositionLinearApproximation(...) const = 0;
  // ...
};
```

两个实现：

| 类 | 导数来源 | 特点 |
|---|---|---|
| `PinocchioEndEffectorKinematics` | Pinocchio 的解析 Jacobian | 快，需要 `PreComputation` 先调 `computeJointJacobians` |
| `PinocchioEndEffectorKinematicsCppAd` | CppAD 自动微分 | 首次编译慢，之后同样快；无需手写导数链 |

**`PinocchioStateInputMapping`**：把 OCS2 的 $(x,u)$ 映射到 Pinocchio 的 $(q,v)$。
不同机器人的状态定义不同（质心动力学 vs 关节空间），
这个接口把差异隔离开。

---

## 12.3 ⭐ 质心动力学（Centroidal Dynamics）

**目录**：`ocs2_pinocchio/ocs2_centroidal_model/`

这是四足/人形机器人 MPC 的**标准建模方法**。

### 12.3.1 为什么用质心动力学

全身刚体动力学：

$$
M(q)\ddot q + C(q,\dot q)\dot q + g(q) = S^{\!\top}\tau + \sum_i J_{c,i}^{\!\top}f_i
$$

$M$ 是 $n\times n$（$n=18$ 对四足），求逆很贵，且 $\tau$ 与 $f$ 都是决策变量。

**关键观察**：前 6 行（浮动基部分）**不含 $\tau$**（$S$ 的前 6 列为 0），
即浮动基的运动**完全由外力决定**。这 6 个方程就是牛顿-欧拉方程。

**质心动力学**只保留这 6 个方程，把关节运动降级为运动学量：

$$
\boxed{\;\dot{\mathbf{h}}=\begin{bmatrix}\dot{\mathbf{l}}\\ \dot{\mathbf{k}}\end{bmatrix}
=\begin{bmatrix}m\mathbf{g}+\sum_i \mathbf{f}_i\\
\sum_i(\mathbf{p}_i-\mathbf{p}_{com})\times \mathbf{f}_i+\sum_j\boldsymbol{\tau}_j\end{bmatrix}\;}
$$

其中：
- $\mathbf{l}$ = 线动量、$\mathbf{k}$ = 关于质心的角动量
- $\mathbf{f}_i$ = 第 $i$ 个接触点的力（世界系）
- $\mathbf{p}_i$ = 接触点位置、$\mathbf{p}_{com}$ = 质心位置
- $\boldsymbol{\tau}_j$ = 6DoF 接触的力矩

**推导**（牛顿-欧拉）：
- 线动量定理：$\dot{\mathbf{l}}=\sum\mathbf{F}_{ext}=m\mathbf{g}+\sum_i\mathbf{f}_i$
- 角动量定理（关于质心）：$\dot{\mathbf{k}}=\sum\mathbf{M}_{ext}
  =\sum_i(\mathbf{p}_i-\mathbf{p}_{com})\times\mathbf{f}_i+\sum_j\boldsymbol{\tau}_j$
  （重力对质心无力矩 ✓）

### 12.3.2 OCS2 的状态与输入定义

**文件**：`ocs2_centroidal_model/PinocchioCentroidalDynamics.h:39`

$$
x=\begin{bmatrix}\mathbf{l}/m\\ \mathbf{k}/m\\ \mathbf{p}_{base}\\ \boldsymbol{\theta}_{base}^{zyx}\\ \mathbf{q}_j\end{bmatrix}
\in\mathbb{R}^{6+6+n_j},\qquad
u=\begin{bmatrix}\mathbf{f}_{1..n_c}\\ \boldsymbol{\tau}_{1..n_w}\\ \dot{\mathbf{q}}_j\end{bmatrix}
$$

**三个设计选择**：

1. **归一化动量** $\mathbf{h}/m$：让状态各分量量级接近（否则质量 50kg 的机器人
   动量数值远大于位置），改善数值条件
2. **关节速度作为输入**：而非关节力矩。这让动力学对输入是**线性**的
   （见下），大幅简化 MPC。真实力矩由底层 WBC（whole-body controller）
   通过逆动力学算出
3. **动量在质心系表达**：质心系是"原点在 CoM、姿态与世界系对齐"的坐标系

**维度**（ANYmal 四足）：$n_j=12$，$n_c=4$，$n_w=0$
→ $\dim x = 6+6+12=24$，$\dim u = 12+12=24$。

### 12.3.3 流形

**文件**：`PinocchioCentroidalDynamics.cpp:54`

```cpp
vector_t PinocchioCentroidalDynamics::getValue(scalar_t time, const vector_t& state, const vector_t& input) {
  vector_t f(info.stateDim);
  f << getNormalizedCentroidalMomentumRate(interface, info, input),   // 前 6 维
       mapping_.getPinocchioJointVelocity(state, input);              // 后 6+n_j 维
  return f;
}
```

**前 6 维**（`ModelHelperFunctions.cpp:167`）：

```cpp
Eigen::Matrix<SCALAR_T, 6, 1> centroidalMomentumRate;
centroidalMomentumRate << info.robotMass * gravityVector, Zero3;      // 重力
for (i in 3DoF contacts) {
  centroidalMomentumRate.head<3>() += f_i;                            // ∑ f_i
  centroidalMomentumRate.tail<3>() += (p_i - p_com).cross(f_i);       // ∑ r_i × f_i
}
for (i in 6DoF contacts) {
  centroidalMomentumRate.head<3>() += f_i;
  centroidalMomentumRate.tail<3>() += (p_i - p_com).cross(f_i) + τ_i;
}
centroidalMomentumRate /= info.robotMass;                             // 归一化
```

**⭐ 注意：这对输入 $u$ 是线性的！**
$\dot{\mathbf{h}}/m$ 是 $\mathbf{f}_i$、$\boldsymbol{\tau}_i$ 的线性组合，
系数只依赖状态（通过 $\mathbf{p}_i-\mathbf{p}_{com}$）。
这是质心模型对 MPC 友好的关键性质。

**后 6+n_j 维**：广义坐标的导数 $\dot q$。

### 12.3.4 ⭐ 从动量反解基座速度

**文件**：`CentroidalModelPinocchioMapping.cpp:84`

状态里存的是动量 $\mathbf{h}$，但 $\dot q$ 需要基座速度 $v_b$。
两者的关系由**质心动量矩阵（Centroidal Momentum Matrix, CMM）** $A(q)$ 给出：

$$
\mathbf{h} = A(q)\,\dot q = \underbrace{A_b(q)}_{6\times6}\,v_b + \underbrace{A_j(q)}_{6\times n_j}\,\dot q_j
$$

（Orin & Goswami 2008 定义了 CMM；Pinocchio 的 `computeCentroidalMap` 算它，存在 `data.Ag`）

**反解**：

$$
\boxed{\;v_b = A_b^{-1}\big(\mathbf{h}-A_j\dot q_j\big)\;}
$$

代码：

```cpp
const auto& A = getCentroidalMomentumMatrix(*pinocchioInterfacePtr_);    // data.Ag
const Eigen::Matrix<SCALAR,6,6> Ab = A.template leftCols<6>();
const auto Ab_inv = computeFloatingBaseCentroidalMomentumMatrixInverse(Ab);

Eigen::Matrix<SCALAR,6,1> momentum = info.robotMass * getNormalizedMomentum(state, info);
if (info.centroidalModelType == CentroidalModelType::FullCentroidalDynamics)
  momentum.noalias() -= A.rightCols(info.actuatedDofNum) * jointVelocities;

vPinocchio.head<6>().noalias() = Ab_inv * momentum;
vPinocchio.tail(info.actuatedDofNum) = jointVelocities;
```

**注意 SRBD 模式跳过了 $-A_j\dot q_j$ 项**——见 §12.3.6。

#### $A_b^{-1}$ 的解析形式

$A_b$ 有特殊结构（Wensing & Orin 2016）。以基座速度 $(\dot{\mathbf{p}}_{base},\omega)$
为自变量，并记 $\mathbf{c}=\mathbf{p}_{base}-\mathbf{p}_{com}$（**CoM→base**，与代码约定一致）：

$$
A_b=\begin{bmatrix}m I_3 & m[\,\mathbf{c}\,]_\times\\ 0 & \bar I_{com}\end{bmatrix}
$$

$\bar I_{com}$ 是关于质心的惯性张量。
（左下块为 0，因为关于质心的角动量不依赖基座线速度——这正是"质心系"定义的结果）

分块求逆：

$$
A_b^{-1}=\begin{bmatrix}\frac{1}{m}I_3 & -[\,\mathbf{c}\,]_\times\,\bar I_{com}^{-1}\\
0 & \bar I_{com}^{-1}\end{bmatrix}
$$

**验证**（右上块）：
$$
mI\cdot\big(-[c]_\times\bar I^{-1}\big)+m[c]_\times\cdot\bar I^{-1}
=-m[c]_\times\bar I^{-1}+m[c]_\times\bar I^{-1}=0\ \checkmark
$$
其余块显然，故 $A_bA_b^{-1}=I_6$ ✓

`computeFloatingBaseCentroidalMomentumMatrixInverse` 用这个解析式，
只需对 $3\times3$ 的 $\bar I_{com}$ 求逆，比 $6\times6$ 通用求逆快得多。

### 12.3.5 线性化

**文件**：`PinocchioCentroidalDynamics.cpp:68`

```cpp
VectorFunctionLinearApproximation getLinearApproximation(scalar_t t, const vector_t& x, const vector_t& u) {
  auto dynamics = VectorFunctionLinearApproximation::Zero(info.stateDim, info.stateDim, info.inputDim);
  dynamics.f = getValue(t, x, u);

  computeNormalizedCentroidalMomentumRateGradients(x, u);     // ← 见下

  matrix_t dfdq = matrix_t::Zero(info.stateDim, info.generalizedCoordinatesNum);
  matrix_t dfdv = matrix_t::Zero(info.stateDim, info.generalizedCoordinatesNum);
  dfdq.topRows<6>() << normalizedLinearMomentumRateDerivativeQ_, normalizedAngularMomentumRateDerivativeQ_;
  dfdv.bottomRows(info.generalizedCoordinatesNum).setIdentity();

  std::tie(dynamics.dfdx, dynamics.dfdu) = mapping_.getOcs2Jacobian(x, dfdq, dfdv);   // 链式法则

  // 显式输入依赖（接触力/力矩）
  dynamics.dfdu.topRows<3>()      += normalizedLinearMomentumRateDerivativeInput_;
  dynamics.dfdu.middleRows(3, 3)  += normalizedAngularMomentumRateDerivativeInput_;
  return dynamics;
}
```

#### 动量变化率的梯度

**文件**：`PinocchioCentroidalDynamics.cpp:97`

对 3DoF 接触：

$$
\dot{\mathbf{k}}/m \ni \frac{1}{m}(\mathbf{p}_i-\mathbf{p}_{com})\times\mathbf{f}_i
$$

**对 $q$ 求导**（$\mathbf{f}_i$ 固定）：

$$
\frac{\partial}{\partial q}\Big[\frac{1}{m}(\mathbf{p}_i-\mathbf{p}_{com})\times\mathbf{f}_i\Big]
=-\frac{1}{m}[\,\mathbf{f}_i\,]_\times\,J_{i}^{com}
$$

其中 $J_i^{com}=\dfrac{\partial(\mathbf{p}_i-\mathbf{p}_{com})}{\partial q}$，
用了 $a\times b=-b\times a$ 与 $\frac{\partial(a\times b)}{\partial a}=-[b]_\times$。

代码（`:110-115`）：

```cpp
const Vector3 contactForceInWorldFrame = getContactForces(input, i, info);
f_hat = skewSymmetricMatrix(contactForceInWorldFrame) / info.robotMass;      // [f]_× / m
const auto J = getTranslationalJacobianComToContactPointInWorldFrame(interface, info, i);
normalizedAngularMomentumRateDerivativeQ_.noalias() -= f_hat * J;            // −[f]_× J / m ✓
normalizedLinearMomentumRateDerivativeInput_.block<3,3>(0, 3*i).diagonal().array() = 1.0/info.robotMass;
```

$J_i^{com}$ 的计算（`ModelHelperFunctions.cpp:151`）：

$$
J_i^{com}=J_i - J_{com},\qquad J_{com}=\frac{1}{m}A_{lin}(q)
$$

```cpp
pinocchio::getFrameJacobian(model, data, frameIdx, pinocchio::LOCAL_WORLD_ALIGNED, J_full);
Matrix3x J_com = getCentroidalMomentumMatrix(interface).topRows<3>() / info.robotMass;
return (J_full.topRows<3>() - J_com);
```

**为什么 $J_{com}=A_{lin}/m$？** 因为 $\mathbf{l}=m\dot{\mathbf{p}}_{com}=A_{lin}\dot q$，
故 $\dot{\mathbf{p}}_{com}=\frac{A_{lin}}{m}\dot q$，即 $J_{com}=A_{lin}/m$ ✓

**对 $u$ 求导**：

$$
\frac{\partial(\dot{\mathbf{l}}/m)}{\partial\mathbf{f}_i}=\frac{1}{m}I_3,\qquad
\frac{\partial(\dot{\mathbf{k}}/m)}{\partial\mathbf{f}_i}=\frac{1}{m}[\,\mathbf{p}_i-\mathbf{p}_{com}\,]_\times
$$

（$a\times b$ 对 $b$ 的导数是 $[a]_\times$）

#### `getOcs2Jacobian`：链式法则

**文件**：`CentroidalModelPinocchioMapping.cpp:113`

Pinocchio 给出的是 $\partial f/\partial q$、$\partial f/\partial v$，
但 OCS2 需要 $\partial f/\partial x$、$\partial f/\partial u$。
由于 $q=q(x)$、$v=v(x,u)$：

$$
\frac{\partial f}{\partial x}=\frac{\partial f}{\partial q}\frac{\partial q}{\partial x}
+\frac{\partial f}{\partial v}\frac{\partial v}{\partial x},\qquad
\frac{\partial f}{\partial u}=\frac{\partial f}{\partial v}\frac{\partial v}{\partial u}
$$

其中 $\partial v/\partial x$ 来自 §12.3.4 的 $v_b=A_b^{-1}(\mathbf{h}-A_j\dot q_j)$：

$$
\frac{\partial v_b}{\partial(\mathbf{h}/m)}=m\,A_b^{-1},\qquad
\frac{\partial v_b}{\partial \dot q_j}=-A_b^{-1}A_j,\qquad
\frac{\partial \dot q_j}{\partial u}=[\,0\ \ I\,]
$$

代码 `:131`：
```cpp
floatingBaseVelocitiesDerivativeState.leftCols(6) = info.robotMass * Ab_inv;   // m·Ab⁻¹ ✓
```

还有一项 $\partial v_b/\partial q$（因为 $A_b$、$A_j$ 依赖 $q$），
用 Pinocchio 的 `dHdq`、`dFda` 等导数（`:137-143`）。

### 12.3.6 ⭐ SRBD 近似

**文件**：`ModelHelperFunctions.cpp:61`

```cpp
enum class CentroidalModelType { FullCentroidalDynamics, SingleRigidBodyDynamics };
```

**SRBD（单刚体动力学）假设**：
> 腿的质量与惯量相对躯干可忽略，且**关节运动不产生显著动量**。

于是：

$$
A(q)\approx\begin{bmatrix}A_b^{nom}(\theta_{base}) & 0\end{bmatrix}
$$

即 $A_j\approx0$，且 $A_b$ 只依赖基座姿态（用名义构型 $q^{nom}$ 的惯量）。

代码构造 $A_b$（`:62-74`）：

```cpp
const vector3_t eulerAnglesZyx = q.segment<3>(3);
const matrix3_t mappingZyx = getMappingFromEulerAnglesZyxDerivativeToGlobalAngularVelocity(eulerAnglesZyx);  // T(θ)
const matrix3_t R = getRotationMatrixFromZyxEulerAngles(eulerAnglesZyx);      // R(θ)
const vector3_t c = R * info.comToBasePositionNominal;                        // CoM→base 向量（世界系）
const matrix3_t c_hat = skewSymmetricMatrix(c);

matrix6_t Ab = matrix6_t::Zero();
Ab.topLeftCorner<3,3>().diagonal().array() = info.robotMass;                   // m·I
Ab.topRightCorner<3,3>().noalias() = info.robotMass * c_hat * mappingZyx;      // m·[c]_× T
Ab.bottomRightCorner<3,3>().noalias() = (R * info.centroidalInertiaNominal) * (R.transpose() * mappingZyx);
                                                                              // R Ī Rᵀ T
data.Ag.leftCols<6>() = Ab;
data.com[0] = q.head<3>() - c;                                                // p_com = p_base − c
```

**数学形式**：

$$
A_b=\begin{bmatrix}
m I_3 & m\,[\,\mathbf{c}\,]_\times\,T(\theta)\\
0 & R\,\bar I^{nom}\,R^{\!\top}\,T(\theta)
\end{bmatrix}
$$

**逐块解释**：
- $(1,1)$ **$mI$**：$\mathbf{l}=m\dot{\mathbf{p}}_{com}$，而 $\dot{\mathbf{p}}_{com}=\dot{\mathbf{p}}_{base}+\omega\times\mathbf{c}$，
  线速度部分贡献 $mI$
- $(1,2)$ **$m[\mathbf{c}]_\times T$**：代码里 `comToBasePositionNominal` 是 **CoM→base** 向量，
  经 $R$ 旋转后 $\mathbf{c}=\mathbf{p}_{base}-\mathbf{p}_{com}$
  （与 `data.com[0] = q.head<3>() - c` 一致）。于是
  $$
  \dot{\mathbf{p}}_{com}=\dot{\mathbf{p}}_{base}+\omega\times(\mathbf{p}_{com}-\mathbf{p}_{base})
  =\dot{\mathbf{p}}_{base}-\omega\times\mathbf{c}
  =\dot{\mathbf{p}}_{base}+\mathbf{c}\times\omega
  =\dot{\mathbf{p}}_{base}+[\,\mathbf{c}\,]_\times\,T\dot\theta
  $$
  两边乘 $m$ 即得该块
- $(2,1)$ **$0$**：角动量（关于质心）不依赖基座线速度 ✓
- $(2,2)$ **$R\bar I^{nom}R^\top T$**：$\mathbf{k}=I_{world}\omega$，
  其中 $I_{world}=R\bar I^{nom}R^\top$（惯量张量的坐标变换），再乘 $T$

**SRBD 的代价**：
- 忽略了腿摆动产生的角动量（对高速运动有影响）
- 惯量固定为名义构型（蹲下/站立时惯量实际变化）

**SRBD 的收益**：
- 无需 `computeCentroidalMap`（Pinocchio 里较贵的调用）
- $A_j=0$ 使 §12.3.4 的反解简化为 $v_b = A_b^{-1}\mathbf{h}$
- 导数计算大幅简化

**实测**：ANYmal 的 trot/walk 用 SRBD 与 Full 几乎无差别，
但用 SRBD 快 30~40%。跳跃、快跑时差别开始显现。

### 12.3.7 `CentroidalModelRbdConversions`

**文件**：`ocs2_centroidal_model/CentroidalModelRbdConversions.h`

MPC 的质心状态 ↔ 全身 RBD 状态的双向转换。

**RBD → Centroidal**（`:99`）：

```cpp
vector_t computeCentroidalStateFromRbdModel(const vector_t& rbdState) {
  // rbdState = [θ_base(zyx), p_base, q_j, ω_base(world), v_base(world), q̇_j]
  qPinocchio << p_base, θ_base, q_j;
  vPinocchio << v_base,
                getEulerAnglesZyxDerivativesFromGlobalAngularVelocity(θ_base, ω_base),  // ω → θ̇
                q̇_j;
  updateCentroidalDynamics(pinocchioInterface_, info, qPinocchio);
  const auto& A = getCentroidalMomentumMatrix(pinocchioInterface_);
  getNormalizedMomentum(state, info).noalias() = A * vPinocchio / info.robotMass;    // h/m = A·v/m
  getGeneralizedCoordinates(state, info) = qPinocchio;
  return state;
}
```

**Centroidal → RBD**（`:128`）：反向，用 §12.3.4 的 $v_b=A_b^{-1}(\mathbf{h}-A_j\dot q_j)$。

**加速度计算**（`:85`）：给定期望的 $\dot{\mathbf{h}}$ 与 $\ddot q_j$，求 $\ddot q_b$：

$$
\dot{\mathbf{h}} = \dot A\dot q + A\ddot q
= \dot A\dot q + A_b\ddot q_b + A_j\ddot q_j
\;\Longrightarrow\;
\ddot q_b = A_b^{-1}\big(\dot{\mathbf{h}}-\dot A\dot q - A_j\ddot q_j\big)
$$

```cpp
Vector6 centroidalMomentumRate = info.robotMass * getNormalizedCentroidalMomentumRate(pinocchioInterface_, info, input);
centroidalMomentumRate.noalias() -= Adot * vPinocchio;
centroidalMomentumRate.noalias() -= Aj * jointAccelerations.head(info.actuatedDofNum);
const Vector6 qbaseDdot = Ab_inv * centroidalMomentumRate;
```

**这是 WBC（全身控制器）的入口**：MPC 给出 $\dot{\mathbf{h}}$ 与接触力，
WBC 用这个转换得到关节加速度，再用逆动力学算力矩。

### 12.3.8 `AccessHelperFunctions`

类型安全的状态/输入分量访问：

```cpp
auto getNormalizedMomentum(state, info);        // → state.head<6>()
auto getGeneralizedCoordinates(state, info);    // → state.tail(6 + nj)
auto getBasePose(state, info);
auto getJointAngles(state, info);
auto getContactForces(input, i, info);          // → input.segment<3>(3*i)
auto getContactTorques(input, i, info);
auto getJointVelocities(input, info);
```

**全部是 `Eigen::Block` 表达式，零拷贝**。
用它们而非手写 `segment<3>(3*i)` 可以避免索引错误——
质心模型的输入布局（先所有 3DoF 力，再所有 6DoF 力+力矩，最后关节速度）
很容易搞错。

---

## 12.4 `ocs2_self_collision`：自碰撞避免

**目录**：`ocs2_pinocchio/ocs2_self_collision/`

### 12.4.1 `PinocchioGeometryInterface`

封装 Pinocchio 的 `GeometryModel`（碰撞几何）与 **HPP-FCL**（碰撞检测）。

```cpp
class PinocchioGeometryInterface {
  PinocchioGeometryInterface(const PinocchioInterface&, const std::vector<std::string>& collisionLinkPairs,
                             const std::vector<std::pair<size_t,size_t>>& collisionObjectPairs);
  std::vector<hpp::fcl::DistanceResult> computeDistances(const PinocchioInterface&) const;
  size_t getNumCollisionPairs() const;
};
```

用户指定要检查的**连杆对**（如"左前小腿 vs 躯干"），
`computeDistances` 返回每对的**有符号距离**与最近点。

### 12.4.2 ⭐ 距离约束的梯度

**文件**：`ocs2_self_collision/src/SelfCollision.cpp`

约束 $d_k(q)\ge d_{\min}$，需要 $\partial d/\partial q$。

**关键几何事实**：设最近点对为 $\mathbf{p}_1$（物体 1 上）与 $\mathbf{p}_2$（物体 2 上），
单位法向 $\mathbf{n}=\dfrac{\mathbf{p}_2-\mathbf{p}_1}{\|\mathbf{p}_2-\mathbf{p}_1\|}$，则

$$
\boxed{\;\frac{\partial d}{\partial q}=\mathbf{n}^{\!\top}\big(J_{\mathbf{p}_2}-J_{\mathbf{p}_1}\big)\;}
$$

其中 $J_{\mathbf{p}_i}$ 是最近点的平移 Jacobian。

**为什么？**（包络定理 / Danskin 定理）
$d=\min_{\mathbf{x}\in\mathcal{O}_1,\mathbf{y}\in\mathcal{O}_2}\|\mathbf{y}-\mathbf{x}\|$
是一个最小值函数。由包络定理，对参数 $q$ 求导时**不需要考虑最优点 $\mathbf{p}_1,\mathbf{p}_2$ 随 $q$ 的移动**
（那部分的贡献为零，因为在最优点处目标对 $\mathbf{x},\mathbf{y}$ 的偏导为零，
或位于边界的法向），只需对 $q$ 的显式依赖求导：

$$
\frac{\partial d}{\partial q}=\frac{\partial}{\partial q}\|\mathbf{p}_2(q)-\mathbf{p}_1(q)\|
=\frac{(\mathbf{p}_2-\mathbf{p}_1)^{\!\top}}{\|\mathbf{p}_2-\mathbf{p}_1\|}\Big(J_{\mathbf{p}_2}-J_{\mathbf{p}_1}\Big)
=\mathbf{n}^{\!\top}(J_{\mathbf{p}_2}-J_{\mathbf{p}_1})\ \checkmark
$$

**注意**：这个梯度在 $d=0$（穿透）时法向定义会翻转，
且当最近点对**跳变**时（如物体绕过棱边）梯度不连续。
实践中通过保持足够的安全裕度（$d_{\min}>0$）规避。

### 12.4.3 两条实现路径

| 类 | 导数 | 特点 |
|---|---|---|
| `SelfCollisionConstraint` | 用上面的解析公式 + Pinocchio Jacobian | 快，但依赖 FCL 的最近点信息 |
| `SelfCollisionConstraintCppAd` | 全 CppAD | 更通用；但 FCL 无法 AD，所以**最近点对被"冻结"**为参数 |

**`SelfCollisionCppAd` 的技巧**：
把上一次 FCL 算出的最近点在**各自连杆局部系**的坐标作为 CppAD 的**参数**，
于是距离函数变成"两个固定局部点之间的距离"，这是纯运动学，可以 AD。
只要机器人构型变化不大，这个近似很准。

`ocs2_self_collision_visualization` 提供 RViz 可视化：画出碰撞对与最近点连线。

---

## 12.5 `ocs2_sphere_approximation`

**目录**：`ocs2_pinocchio/ocs2_sphere_approximation/`

**动机**：FCL 的精确碰撞检测（凸包-凸包）较慢且梯度不光滑。
用**球集合**近似连杆几何后：

- 球-球距离 = $\|\mathbf{c}_1-\mathbf{c}_2\|-r_1-r_2$，**解析且光滑**
- 球-点云/距离场距离 = 查表 $-r$，极快

```cpp
class SphereApproximation {
  SphereApproximation(const hpp::fcl::CollisionGeometry& geometry, size_t objectId,
                      scalar_t maxExcess, scalar_t shrinkRatio);
  const std::vector<vector3_t>& getSphereCentersToObjectCenter() const;
  const std::vector<scalar_t>& getSphereRadii() const;
};
```

`maxExcess` 控制近似质量：球集合的并集不超出原几何体 `maxExcess` 距离。
`shrinkRatio` 控制球半径的收缩（更保守 = 更安全但更保守）。

`PinocchioSphereInterface` 管理整个机器人的球近似，
`PinocchioSphereKinematics` / `...CppAd` 提供球心位置的运动学与导数。

---

## 12.6 `ocs2_perceptive`：感知约束

**目录**：`ocs2_perceptive/`

### 12.6.1 距离变换（Distance Transform）

```cpp
class DistanceTransformInterface {
  virtual scalar_t getValue(const vector3_t& position) const = 0;
  virtual std::pair<scalar_t, vector3_t> getLinearApproximation(const vector3_t& p) const = 0;
};
```

给定 3D 点，返回到最近障碍的**有符号距离**及其梯度。
典型实现是 **ESDF**（Euclidean Signed Distance Field），
由感知模块（如 Voxblox）从点云生成。

`ComputeDistanceTransform.h` 提供从占据栅格计算 EDT 的算法
（Felzenszwalb & Huttenlocher 的 $O(n)$ 平方距离变换）。

### 12.6.2 ⭐ 三线性插值

**文件**：`ocs2_perceptive/interpolation/TrilinearInterpolation.h`

ESDF 存储在离散网格上，查询任意点需要插值。
三线性插值同时给出**值与梯度**，且梯度是解析的（对 MPC 至关重要）。

设查询点在体素内的归一化坐标 $(\alpha,\beta,\gamma)\in[0,1]^3$，
八个顶点值 $v_{ijk}$（$i,j,k\in\{0,1\}$），则

$$
f(\alpha,\beta,\gamma)=\sum_{i,j,k\in\{0,1\}} v_{ijk}\,
w_i(\alpha)\,w_j(\beta)\,w_k(\gamma)
$$

其中 $w_0(t)=1-t$、$w_1(t)=t$。

**梯度**（对归一化坐标）：

$$
\frac{\partial f}{\partial\alpha}=\sum_{i,j,k}v_{ijk}\,w_i'(\alpha)w_j(\beta)w_k(\gamma),
\qquad w_0'=-1,\ w_1'=+1
$$

再除以体素尺寸得到世界系梯度：

$$
\nabla_{\mathbf{p}}f=\begin{bmatrix}
\frac{1}{h_x}\partial_\alpha f\\
\frac{1}{h_y}\partial_\beta f\\
\frac{1}{h_z}\partial_\gamma f\end{bmatrix}
$$

`BilinearInterpolation.h` 是 2D 版本（用于高程图 elevation map）。

**注意**：三线性插值的梯度在体素边界**不连续**（$C^0$ 但非 $C^1$）。
对 DDP 会产生 Hessian 的跳变。若需要 $C^1$，应使用三次样条插值，
但代价高得多。OCS2 的选择是接受这个不连续，靠罚函数的光滑性补偿。

### 12.6.3 `EndEffectorDistanceConstraint`

约束末端（或球近似的球心）到障碍的距离：

$$
h(x) = d\big(\mathbf{p}_{ee}(x)\big) - d_{\min}\ \ge\ 0
$$

梯度由链式法则：

$$
\frac{\partial h}{\partial x}=\nabla_{\mathbf{p}}d\cdot J_{ee}(x)
$$

其中 $\nabla_{\mathbf{p}}d$ 来自距离场插值，$J_{ee}$ 来自运动学。

`EndEffectorDistanceConstraintCppAd` 是 CppAD 版本，
把距离场的值与梯度作为**参数**传入（因为距离场本身无法 AD）。

`ocs2_perceptive_anymal` 是它在 ANYmal 上的应用示例。

---

## 12.7 小结：从 URDF 到 OCP 的完整链条

```
URDF 文件
   │ getPinocchioInterfaceFromUrdfFile(path, JointModelFreeFlyer())
   ▼
PinocchioInterface (Model + Data)
   │ createCentroidalModelInfo(interface, type, nominalJointAngles, threeDofContactNames, sixDofContactNames)
   ▼
CentroidalModelInfo  ←── 定义了 stateDim, inputDim, 接触点索引, 名义惯量
   │
   ├─▶ PinocchioCentroidalDynamicsAD  ──▶ problem.dynamicsPtr
   │
   ├─▶ PinocchioEndEffectorKinematicsCppAd
   │     ├─▶ EndEffectorLinearConstraint (零速度)  ──▶ problem.equalityConstraintPtr
   │     └─▶ SwingTrajectoryPlanner 的参考
   │
   ├─▶ FrictionConeConstraint  ──▶ problem.softConstraintPtr
   ├─▶ ZeroForceConstraint     ──▶ problem.equalityConstraintPtr
   │
   ├─▶ SelfCollisionConstraintCppAd  ──▶ problem.stateSoftConstraintPtr
   │
   └─▶ LeggedRobotPreComputation  ──▶ problem.preComputationPtr
```

全部在 `LeggedRobotInterface` 的构造函数里完成
（`ocs2_robotic_examples/ocs2_legged_robot/src/LeggedRobotInterface.cpp`）。

---

## 12.8 参考文献

- **Orin, D. & Goswami, A. (2008)**, *Centroidal Momentum Matrix of a Humanoid Robot: Structure and Properties*, IROS —— **CMM 的定义**
- **Orin, Goswami, Lee (2013)**, *Centroidal Dynamics of a Humanoid Robot*, Autonomous Robots 35(2) —— 质心动力学系统性阐述
- **Wensing, P. & Orin, D. (2016)**, *Improved Computation of the Humanoid Centroidal Dynamics and Application for Whole-Body Control*, Int. J. Humanoid Robotics —— $A_b$ 的结构与高效求逆
- **Dai, Valenzuela, Tedrake (2014)**, *Whole-body Motion Planning with Centroidal Dynamics and Full Kinematics*, Humanoids
- **Carpentier, J. et al. (2019)**, *The Pinocchio C++ Library*, SII —— Pinocchio
- **Featherstone, R. (2008)**, *Rigid Body Dynamics Algorithms* —— 刚体动力学的标准教材
- **Pardo, Möller, Neunert, Winkler, Buchli (2016)**, *Evaluating Direct Transcription and Nonlinear Optimization Methods for Robot Motion Planning*, RA-L
- **Sleiman, Farshidian, Minniti, Hutter (2021)**, *A Unified MPC Framework for Whole-Body Dynamic Locomotion and Manipulation*, RA-L
- **Felzenszwalb & Huttenlocher (2012)**, *Distance Transforms of Sampled Functions*, Theory of Computing —— EDT 算法
- **Pan, Chitta, Manocha (2012)**, *FCL: A General Purpose Library for Collision and Proximity Queries*, ICRA
- **Gaertner, Bjelonic, Farshidian, Hutter (2021)**, *Collision-Free MPC for Legged Robots in Static and Dynamic Scenes*, ICRA —— 球近似 + 距离场

**下一章**：机器人示例逐个拆解 → [13 机器人示例](13-robot-examples.md)
