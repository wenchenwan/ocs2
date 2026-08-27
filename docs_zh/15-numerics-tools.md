# 15 · 数值与工程工具

本章覆盖 `ocs2_core` 里那些"不属于任何算法但处处被用到"的基础设施。
它们决定了 OCS2 能否在实时系统上跑起来。

---

## 15.1 ⭐ 自动微分（CppAD + CppADCodeGen）

### 15.1.1 为什么需要

机器人动力学的导数：
- **手写**：几百行公式，极易出错，改模型就要重写（见四旋翼的例子）
- **有限差分**：$O(n)$ 次函数求值，精度差（截断误差 vs 舍入误差的权衡）
- **自动微分**：精确到机器精度，且计算量与函数本身同阶

OCS2 选择 **CppAD + CppADCodeGen**：不仅自动微分，还**生成 C 代码并编译**，
使运行时性能接近手写。

### 15.1.2 类型体系

**文件**：`ocs2_core/automatic_differentiation/Types.h`

```cpp
using ad_base_t   = CppAD::cg::CG<scalar_t>;      // CodeGen 标量（记录表达式树）
using ad_scalar_t = CppAD::AD<ad_base_t>;         // AD 包装
using ad_vector_t = Eigen::Matrix<ad_scalar_t, Eigen::Dynamic, 1>;
using ad_matrix_t = Eigen::Matrix<ad_scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
```

**两层嵌套的含义**：
- `CG<double>` 不做数值计算，而是**记录操作序列**（构建计算图）
- `AD<CG<double>>` 在这个记录之上再做自动微分的 taping

结果：一次"求值"实际上是构建了一棵**符号表达式树**，
CppADCodeGen 再把它翻译成 C 源码。

### 15.1.3 `CppAdInterface` 的四阶段流程

**文件**：`ocs2_core/automatic_differentiation/CppAdInterface.h`

```
① Taping（记录）
   ad_vector_t x(nx), p(np), y;
   CppAD::Independent(x);          // 声明自变量
   adFunction(x, p, y);            // 用户函数，用 ad_scalar_t 运算
   CppAD::ADFun<ad_base_t> fun(x, y);   // 构建计算图

② 代码生成
   CppAD::cg::ModelCSourceGen<scalar_t> cgen(fun, modelName);
   cgen.setCreateJacobian(true);
   cgen.setCreateHessian(true);
   cgen.setCreateSparseJacobian(true);   // ← 稀疏模式
   → 生成 .c 文件

③ 编译
   CppAD::cg::GccCompiler<scalar_t> compiler;
   compiler.setCompileFlags({"-O3","-g","-march=native","-mtune=native","-ffast-math"});
   → 生成 .so

④ 加载与调用
   dynamicLib_ = std::make_unique<CppAD::cg::LinuxDynamicLib<scalar_t>>(libraryPath);
   model_ = dynamicLib_->model(modelName);
   model_->ForwardZero(x, y);       // 求值
   model_->SparseJacobian(x, jac);  // 一阶导
   model_->SparseHessian(x, hess);  // 二阶导
```

三个入口：

```cpp
void createModels(ApproximationOrder order, bool verbose);      // 强制重新生成
void loadModels(bool verbose);                                   // 只加载已有
void loadModelsIfAvailable(ApproximationOrder order, bool verbose);  // 有则加载、无则生成
```

`ApproximationOrder`：`Zero`（只求值）、`First`（+Jacobian）、`Second`（+Hessian）。
**只生成需要的阶数**——Hessian 的代码量与编译时间远大于 Jacobian。

### 15.1.4 参数化函数

```cpp
CppAdInterface(ad_parameterized_function_t adFunction, size_t variableDim, size_t parameterDim,
               std::string modelName, std::string folderName, std::vector<std::string> compileFlags);
```

区分**变量** $x$（要求导）与**参数** $p$（不求导，运行时可变）。

典型用途：
- 参数 = 目标轨迹（每步变化，但不需要对它求导）
- 参数 = 上一次的碰撞最近点（见 [12 章](12-robot-models.md) §12.4.3）
- 参数 = 距离场的局部值与梯度（[12 章](12-robot-models.md) §12.6.3）

**这个设计让"不可微的外部数据"可以进入 AD 函数** ✓

### 15.1.5 稀疏性

**文件**：`automatic_differentiation/CppAdSparsity.h`

机器人的 Jacobian/Hessian 通常很稀疏
（例如末端位置只依赖该腿的关节角，与其他腿无关）。

```cpp
SparsityPattern getJacobianVariableSparsity(int rangeDim, int variableDim);
SparsityPattern getHessianVariableSparsity(int variableDim);
```

CppADCodeGen 只为**非零条目**生成代码，
可以把生成的代码量减少一个数量级，编译更快、运行更快。

### 15.1.6 实践要点

**⚠️ 首次运行会卡住**

第一次运行会 taping + 生成 C 代码 + 调 gcc 编译，
四足机器人可能要 **1~5 分钟**。之后从磁盘加载，秒级。

```cpp
// ModelSettings 里的开关
bool recompileLibrariesCppAd = true;   // 改成 false 跳过重新生成
std::string modelFolderCppAd = "/tmp/ocs2";
```

**开发时开 `true`（改了模型自动重编），部署时关 `false`。**

**⚠️ `-march=native` 不可跨机器**

生成的 `.so` 用了编译机的指令集。
交叉编译或分发到不同 CPU 时必须改编译选项。

**⚠️ 磁盘缓存的失效**

OCS2 用**模型名**做缓存键，**不检查函数内容是否变了**。
所以改了动力学但没改模型名，会加载到旧的库！
排查"改了代码没生效"时先删 `/tmp/ocs2`。

**验证工具**：`FiniteDifferenceMethods.h`

```cpp
bool checkSystemDynamicsDerivative(SystemDynamicsBase& system, scalar_t t,
                                   const vector_t& x, const vector_t& u,
                                   scalar_t tolerance, bool verbose);
```

用有限差分校验 AD 结果。**写新模型时必跑**。

---

## 15.2 ⭐ 线性代数工具

**文件**：`ocs2_core/misc/LinearAlgebra.h` / `.cpp`

已在 [04 章](04-ocs2-ddp.md) 详细推导的：
- `makePsdEigenvalue` / `makePsdGershgorin` / `makePsdCholesky`（§4.8）
- `computeInverseMatrixUUT` / `computeConstraintProjection`（§4.4）
- `qrConstraintProjection` / `luConstraintProjection`（[03 章](03-ocs2-oc.md) §3.5.3）

本节补充其余工具。

### `setTriangularMinimumEigenvalues`

```cpp
void setTriangularMinimumEigenvalues(matrix_t& Lr, scalar_t minEigenValue = ...);
```

把三角矩阵的对角元绝对值下限设为 `minEigenValue`，**保持符号**：

```cpp
if (eigenValue < 0.0) eigenValue = std::min(-minEigenValue, eigenValue);
else                  eigenValue = std::max( minEigenValue, eigenValue);
```

防止 QR 分解得到近奇异的 $R$ 时求逆爆炸。
**保号很重要**：直接取绝对值会破坏分解的符号结构，
导致后续的投影矩阵不正确。

### `symmetricEigenSolver`

```cpp
std::pair<vector_t, matrix_t> symmetricEigenSolver(const matrix_t& mat);
```

对称矩阵的特征值与特征向量，`Eigen::SelfAdjointEigenSolver` 的封装。

### `rank` / `eigenvalues` / `nullspaceProjection`

诊断与测试用的辅助函数。

---

## 15.3 插值

### 15.3.1 `LinearInterpolation`

**文件**：`ocs2_core/misc/LinearInterpolation.h`

**核心设计：分离"查找"与"插值"**

```cpp
// ① 查找一次
std::pair<int, scalar_t> indexAlpha = LinearInterpolation::timeSegment(t, timeArray);
//                       ↑ index      ↑ alpha ∈ [0,1]

// ② 复用查找结果，插值多个量
auto A = LinearInterpolation::interpolate(indexAlpha, dataArray, accessorA);
auto B = LinearInterpolation::interpolate(indexAlpha, dataArray, accessorB);
auto C = LinearInterpolation::interpolate(indexAlpha, dataArray, accessorC);
```

插值公式：

$$
v(t) = \alpha\,v_{i} + (1-\alpha)\,v_{i+1},
\qquad \alpha = \frac{t_{i+1}-t}{t_{i+1}-t_i}
$$

（注意 OCS2 的 $\alpha$ 是"距离右端点的归一化距离"，与常见约定相反）

**为什么这样设计？** Riccati 积分的每个时间点要插值十几个量
（$A$、$B$、$Q$、$q$、$R$、$r$、$P$、$\Delta Q$、$\Delta G_m$、$\Delta G_v$…）。
共享 `indexAlpha` 避免了十几次二分查找，
在 $N=1000$ 的轨迹上能省 20%~30% 的后向传播时间。

**边界处理**：$t$ 超出范围时钳位到端点（不外插）。

### 15.3.2 `ModelDataLinearInterpolation`

**文件**：`ocs2_core/model_data/ModelDataLinearInterpolation.h`

定义访问 `ModelData` 各字段的 accessor：

```cpp
namespace model_data {
  inline const matrix_t& dynamics_dfdx(const std::vector<ModelData>& v, size_t i) { return v[i].dynamics.dfdx; }
  inline const matrix_t& dynamics_dfdu(...);
  inline const vector_t& dynamicsBias(...);
  inline const matrix_t& cost_dfdxx(...);
  // ...
}
```

配合 `LinearInterpolation::interpolate(indexAlpha, modelDataArray, model_data::dynamics_dfdx)` 使用。

`RiccatiModificationInterpolation.h` 是 `riccati_modification::Data` 的对应版本。

### 15.3.3 `Lookup`

```cpp
namespace lookup {
  int findIndexInTimeArray(const scalar_array_t& timeArray, scalar_t time);
  int findActiveIntervalInTimeArray(const scalar_array_t& timeArray, scalar_t time);
  int findBoundedActiveIntervalInTimeArray(const scalar_array_t& timeArray, scalar_t time);
}
```

二分查找的不同变体，区别在边界与越界时的行为。
`ModeSchedule::modeAtTime` 用它。

---

## 15.4 并发工具

### 15.4.1 `ThreadPool`

**文件**：`ocs2_core/thread_support/ThreadPool.h`

```cpp
class ThreadPool {
  ThreadPool(size_t nThreads, int priority = 0);
  template <typename Functor> void runParallel(Functor taskFunction, size_t N);
  size_t numThreads() const;
};
```

`runParallel(task, N)` **阻塞直到 N 个任务全部完成**。
`task` 的签名是 `void(int workerId)`，`workerId` 用于索引每线程资源。

**为什么不用 `std::async` / OpenMP？**
- `std::async` 每次可能创建新线程（不可预测的延迟）
- OpenMP 引入额外依赖，且线程优先级难控制
- 自实现可以设置实时优先级（`SCHED_FIFO`）

### 15.4.2 无锁任务分发模式

OCS2 大量使用这个模式：

```cpp
std::atomic_int timeIndex{0};
auto task = [&](int workerId) {
  int i;
  while ((i = timeIndex++) < N) {      // 原子递增，天然负载均衡
    doWork(i, perThreadResource[workerId]);
  }
};
threadPool.runParallel(task, nThreads);
```

**优点**：
- 无锁（`atomic_int` 的 `fetch_add`）
- **自动负载均衡**：快的线程自动多领任务
- 无需预先划分工作范围

出现位置：`SLQ::approximateIntermediateLQ`、
`SqpSolver::setupQuadraticSubproblem`、
`slp::GGTAbsRowSumInParallel`、`PipgSolver::solve` 等。

### 15.4.3 `Synchronized<T>`

```cpp
template <typename T>
class Synchronized {
  LockedPtr lock();               // RAII 句柄，析构时解锁
  const T& operator*() const;
};
```

值 + mutex 的封装，避免"忘记加锁"。
`ReferenceManager` 用它保护 buffer。

### 15.4.4 `BufferedValue<T>`

双缓冲：写线程写 buffer，读线程 swap。
配合 `Synchronized` 实现 MPC/MRT 的策略传递。

### 15.4.5 `SetThreadPriority`

```cpp
void setThreadPriority(int priority, std::thread& thread);
```

设置 `SCHED_FIFO` 实时调度。需要：
- root 权限，或
- `/etc/security/limits.conf` 里配 `rtprio`

**失败时只打印警告，不抛异常**——
使得程序在没有实时权限的开发机上也能跑（只是不保证时序）。

### 15.4.6 `ExecuteAndSleep`

```cpp
template <typename Functor>
void executeAndSleep(Functor task, scalar_t frequency);
```

固定周期执行：测量 `task` 的耗时，睡眠剩余时间。
用于 `mpcDesiredFrequency_` 的限频。

---

## 15.5 配置文件加载

**文件**：`ocs2_core/misc/LoadData.h`

基于 **Boost.PropertyTree**，格式是 `.info`：

```
ddp
{
  algorithm                     SLQ
  nThreads                      3
  maxNumIterations              10
  minRelCost                    0.1
  constraintTolerance           1e-3
}

Q
{
  scaling 1e+0
  (0,0)  100.0    ; 状态 0
  (1,1)  100.0    ; 状态 1
}
```

API：

```cpp
template <typename T> void loadPtreeValue(const pt::ptree&, T& value, const std::string& key, bool verbose);
template <typename T> void loadCppDataType(const std::string& file, const std::string& key, T& value);
void loadEigenMatrix(const std::string& file, const std::string& key, Eigen::MatrixXd& matrix);
void loadStdVector(...);
void loadStdVectorOfPair(...);   // LoadStdVectorOfPair.h
```

**`loadEigenMatrix` 的语法**：
- `scaling` 是整体缩放因子
- `(i,j) value` 指定单个元素，未指定的为 0
- 对角矩阵只需写 `(i,i)`

**`verbose` 参数**：打印"加载了什么值"或"用了默认值"。
**强烈建议开着**——大多数"参数没生效"的问题都是键名拼错，
verbose 会显示它用了默认值。

**注意**：`.info` 的注释符号是 `;`，不是 `#`。

---

## 15.6 性能测量

**文件**：`ocs2_core/misc/Benchmark.h`

```cpp
namespace benchmark {
class RepeatedTimer {
  void startTimer();
  void endTimer();
  scalar_t getTotalInMilliseconds() const;
  scalar_t getMaxIntervalInMilliseconds() const;
  scalar_t getAverageInMilliseconds() const;
  scalar_t getLastIntervalInMilliseconds() const;
  int getNumTimedIntervals() const;
  void reset();
};
}
```

用 `std::chrono::steady_clock`。求解器里到处是：

```cpp
linearQuadraticApproximationTimer_.startTimer();
approximateOptimalControlProblem();
linearQuadraticApproximationTimer_.endTimer();
```

`getBenchmarkingInfo()` 汇总输出：

```
DDP Benchmarking
    Initialization    :   0.5 [ms] (2.1%)
    LQ Approximation  :   4.8 [ms] (40.3%)
    Backward Pass     :   3.1 [ms] (26.1%)
    Compute Controller:   0.3 [ms] (2.5%)
    Search Strategy   :   3.2 [ms] (26.9%)
```

**调优的第一步永远是看这个表**。
`ocs2_doc/docs/profiling.rst` 有更详细的性能分析指南
（包括 perf、valgrind 的用法）。

---

## 15.7 数值比较

**文件**：`ocs2_core/misc/Numerics.h`

```cpp
namespace numerics {
  template <typename T> bool almost_eq(T x, T y, T prec = ...);
  template <typename T> bool almost_le(T x, T y, T prec = ...);
  template <typename T> bool almost_ge(T x, T y, T prec = ...);
}
```

实现基于**相对误差 + 绝对误差**的混合判据：

$$
|x-y| \le \varepsilon\,\max\big(|x|,|y|,1\big)
$$

**为什么不能直接 `==`？** 浮点数的舍入误差。
例如线搜索里 `stepLength >= minStepLength` 用 `almost_ge`，
否则 $\alpha_{\min}$ 那一步可能因为 $10^{-16}$ 的误差被跳过。

`NumericTraits.h` 提供容差常数：

```cpp
namespace numeric_traits {
  template <typename T> constexpr T limitEpsilon();   // ~1e-9 for double
  template <typename T> constexpr T weakEpsilon();    // 更宽松
}
```

---

## 15.8 显示与日志

### `Display.h`

```cpp
template <typename T> std::string toDelimitedString(const std::vector<T>& v, const std::string& delim = ", ");
std::ostream& operator<<(std::ostream&, const std::vector<T>&);
```

容器的格式化输出，调试时打印轨迹用。

### `Log.h`

Boost.Log 的封装：

```cpp
namespace log {
  void init(const Settings& settings, std::ostream* console = &std::cout);
  Settings loadSettings(const std::string& fileName, const std::string& fieldName);
}
```

支持日志级别、文件输出、格式化。

### `CommandLine.h`

简单的命令行参数解析（`getCommandLineOption`）。

---

## 15.9 测试辅助

### `randomMatrices.h`

```cpp
template <typename MatrixType> MatrixType generateSPDmatrix(int n);          // 随机对称正定
template <typename MatrixType> MatrixType generateSPDmatrix(int n, scalar_t maxEigenvalue);
matrix_t generateFullRowRankmatrix(size_t m, size_t n);                       // 满行秩
```

生成 SPD 矩阵的方法：$M = A^\top A + \epsilon I$，
其中 $A$ 是随机矩阵。保证正定且条件数可控。

单元测试里大量使用——例如测试 Riccati 递推时需要随机但合法的 $Q$、$R$。

### `ocs2_qp_solver`（`ocs2_test_tools`）

**文件**：`ocs2_test_tools/ocs2_qp_solver/`

一个**稠密 QP 参考实现**：把整个 KKT 系统组装成大矩阵直接求解。

```cpp
std::pair<vector_array_t, vector_array_t> solveLinearQuadraticProblem(
    const std::vector<LinearQuadraticStage>& lqApproximation, const vector_t& dx0);
```

组装（`src/QpSolver.cpp:103, 162`）：

$$
\begin{bmatrix}H & G^{\!\top}\\ G & 0\end{bmatrix}
\begin{bmatrix}w\\ \lambda\end{bmatrix}
=\begin{bmatrix}-g\\ -b\end{bmatrix}
$$

用 `Eigen::LDLT` 或 `FullPivLU` 直接求解这个鞍点系统。

**用途**：
- **交叉验证**：SLQ/SQP/IPM/SLP 在同一个 LQ 问题上应给出相同结果
- **调试**：怀疑 HPIPM 算错时用它对照
- **教学**：代码短，KKT 结构一目了然

**局限**：$O((N(n_x+n_u))^3)$，只能用于小问题。

`QpDiscreteTranscription.h` 提供从 `OptimalControlProblem` 到 `LinearQuadraticStage` 的转写，
`testProblemsGeneration.h` 生成随机测试问题。

---

## 15.10 `Collection<T>`：具名容器

**文件**：`ocs2_core/misc/Collection.h`

```cpp
template <typename T>
class Collection {
  void add(std::string name, std::unique_ptr<T> term);
  template <typename Derived> Derived& get(const std::string& name);
  bool empty() const;
  size_t size() const;
  const std::vector<std::unique_ptr<T>>& getTerms() const;
};
```

**设计要点**：
- **按名字存取**：`get<FrictionConeConstraint>("LF_frictionCone").setSurfaceNormalInWorld(n)`
  —— 在线调整某一项的参数
- **`empty()` 快速跳过**：求解器里大量 `if (!problem.xxxPtr->empty())`，
  避免为空容器走完整路径
- **`clone()` 深拷贝**：多线程需要每线程一份

派生类（`StateCostCollection`、`StateConstraintCollection` 等）
在其上实现 `getValue()`（遍历求和）与 `getQuadraticApproximation()`（遍历累加）。

**`getTermsSize(t)`**：返回每一项的约束维度（考虑 `isActive`），
多重打靶的 `ConstraintsSize` 用它。

---

## 15.11 `ocs2_thirdparty`

内联的第三方头文件，避免版本冲突：

| 内容 | 说明 |
|---|---|
| `boost/numeric/odeint` | 部分 odeint 头文件（OCS2 需要的版本可能与系统装的不同） |
| `cppad/cg` 相关 | CppADCodeGen 的部分头文件 |
| 其他 | 小型 header-only 库 |

**这是 catkin 生态的常见做法**：
把不稳定或版本敏感的依赖内联进来，保证可复现构建。

---

## 15.12 调试清单

遇到问题时按这个顺序检查：

### 求解器不收敛

1. 打开 `displayInfo_` / `printSolverStatus`，看每次迭代的 `PerformanceIndex`
2. 打开 `checkNumericalStability_`，看是否有维度/PSD 错误
3. 看 `equalityConstraintsSSE` 是否在下降——不降说明约束不可行或投影失败
4. 看线搜索的步长——总是取到 `alpha_min` 说明下降方向不对
5. 用 `ocs2_qp_solver` 对照单次 LQ 子问题的解

### 导数不对

1. `FiniteDifferenceMethods::checkSystemDynamicsDerivative`
2. 检查 `dfdux` 的形状（$n_u\times n_x$，不是 $n_x\times n_u$）
3. 删 `/tmp/ocs2` 强制重新生成 CppAD 库
4. 检查 `PreComputation` 是否请求了 `Request::Approximation`

### 性能不达标

1. 看 `getBenchmarkingInfo()` 的分解
2. `nThreads_` 是否设对（不要超过物理核数）
3. `checkNumericalStability_` 在生产环境应关闭
4. `useFeedbackPolicy_ = false` 可以省掉增益矩阵的计算与传输
5. SQP：`sqpIteration` 降到 1~3（RTI）
6. SLQ：`preComputeRiccatiTerms_ = true`（reduced form）

### 事件/切换有问题

1. `debugPrintRollout_` 看事件检测
2. 检查 `postEventIndices_` 的长度与时间轴的重复值
3. `getNumFunctionCalls()` 异常大 = 步长塌缩

---

## 15.13 参考文献

- **Bell, B. (2023)**, *CppAD: A Package for Differentiation of C++ Algorithms*,
  https://coin-or.github.io/CppAD/
- **Leal, A. (2018)**, *CppADCodeGen*, https://github.com/joaoleal/CppADCodeGen
- **Griewank & Walther (2008)**, *Evaluating Derivatives: Principles and Techniques of
  Algorithmic Differentiation*, 2nd ed. —— **自动微分的标准教材**
- **Guennebaud, Jacob et al. (2010)**, *Eigen v3*, https://eigen.tuxfamily.org
- **Ahnert & Mulansky (2011)**, *Odeint – Solving Ordinary Differential Equations in C++*,
  AIP Conf. Proc.
- **Higham, N. (2002)**, *Accuracy and Stability of Numerical Algorithms*, 2nd ed.
  —— 浮点比较、条件数、稳定性

**下一章**：遗留模块 → [16 遗留模块](16-legacy-modules.md)
