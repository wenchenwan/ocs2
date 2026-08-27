# 10 · 积分、Rollout 与切换系统

本章聚焦 OCS2 处理**混合动力系统（hybrid system）**的机制：
事件、跳变、保护面、以及围绕它们的数值积分设施。

这是 OCS2 区别于其他 MPC 库的核心能力之一。

---

## 10.1 切换系统的数学模型

$$
\begin{cases}
\dot x = f_{\sigma(t)}(t,x,u), & t\notin\mathcal{E}\\[4pt]
x(t_i^+) = j_i\big(t_i,x(t_i^-)\big), & t_i\in\mathcal{E}
\end{cases}
$$

- $\sigma(t)\in\{1,\dots,M\}$ 是**模态信号**（哪个子系统在生效）
- $\mathcal{E}=\{t_1,t_2,\dots\}$ 是**事件时刻集合**
- $j_i$ 是**跳变映射**（reset map），状态在事件处不连续

### 事件的两种触发方式

| 方式 | 触发条件 | 代码 |
|---|---|---|
| **时间触发** | $t=t_i$，$t_i$ 由 `ModeSchedule` 预先给定 | `TimeTriggeredRollout` |
| **状态触发** | $g_k(t,x)=0$ 且符号从 + 变 −（零穿越） | `StateTriggeredRollout` |

**四足机器人用时间触发**：步态规划器（`GaitSchedule`）预先给出落足/抬足时刻。
**弹跳球、碰撞问题用状态触发**：接触时刻取决于轨迹本身。

---

## 10.2 `ModeSchedule`：模态序列的表示

**文件**：`ocs2_core/include/ocs2_core/reference/ModeSchedule.h`

```cpp
struct ModeSchedule {
  scalar_array_t eventTimes;   // 大小 N-1，严格递增
  size_array_t   modeSequence; // 大小 N
};
```

不变式：`modeSequence.size() == eventTimes.size() + 1`。

时间轴划分：

```
mode[0]      mode[1]      mode[2]           mode[N-1]
────────┬────────────┬────────────┬── … ──┬─────────
      t[0]         t[1]         t[2]    t[N-2]
```

`modeAtTime(t)` 用 `std::upper_bound` 二分查找，$O(\log N)$。

### 四足机器人的模态编码

`ocs2_legged_robot/gait/MotionPhaseDefinition.h`：

```cpp
// mode 是 4 位二进制掩码，bit i 表示第 i 条腿是否触地
// 例如 mode 15 = 0b1111 = STANCE（四腿全支撑）
// mode  9 = 0b1001 = LF 与 RH 触地 = trot 的一个相位
contact_flag_t modeNumber2StanceLeg(size_t modeNumber);
size_t stanceLeg2ModeNumber(const contact_flag_t& stanceLegs);
```

于是"模态"与"接触状态"一一对应，
每个模态对应一组不同的约束（哪些腿零速度、哪些腿零力）。

---

## 10.3 积分器体系

### 10.3.1 `OdeBase` 与 `ControlledSystemBase`

```cpp
class OdeBase {
  virtual vector_t computeFlowMap(scalar_t t, const vector_t& x) = 0;
  virtual vector_t computeJumpMap(scalar_t t, const vector_t& x) { return x; }
  virtual vector_t computeGuardSurfaces(scalar_t t, const vector_t& x) { return vector_t(0); }
  int getNumFunctionCalls() const;   // 性能统计
};

class ControlledSystemBase : public OdeBase {
  void setController(ControllerBase* controller);
  vector_t computeFlowMap(scalar_t t, const vector_t& x) final {
    const vector_t u = controllerPtr_->computeInput(t, x);
    return computeFlowMap(t, x, u);      // 转发到带 u 的版本
  }
};
```

`getNumFunctionCalls()` 用于诊断：如果某次 rollout 的调用次数异常多，
说明积分器在某处步长塌缩（通常是动力学不连续或数值刚性）。

### 10.3.2 步进器（`steppers.h`）

基于 **Boost.Numeric.Odeint**，`ocs2_thirdparty` 内联了部分 odeint 头文件
以避免版本冲突。

| `IntegratorType` | odeint 步进器 | 阶 | 步长 |
|---|---|---|---|
| `EULER` | `euler` | 1 | 定 |
| `RK4` | `runge_kutta4` | 4 | 定 |
| `ODE45` | `runge_kutta_dopri5` + 控制器 | 5(4) | **自适应** |
| `RK5_VARIABLE` | `runge_kutta_cash_karp54` | 5(4) | 自适应 |
| `ADAMS_BASHFORTH` | `adams_bashforth` | 可配 | 定 |
| `ADAMS_BASHFORTH_MOULTON` | 预测-校正 | 可配 | 定 |
| `BULIRSCH_STOER` | `bulirsch_stoer` | 自适应阶 | 自适应 |

`eigenIntegration.h` 提供 odeint 所需的 Eigen 适配：
```cpp
// odeint 需要向量类型支持这些运算
namespace boost::numeric::odeint {
  template<...> struct vector_space_norm_inf<Eigen::Matrix<...>> { ... };
}
```

`RungeKuttaDormandPrince5.h` 是自实现的 DP5，
提供比 odeint 更细粒度的步长控制与事件检测钩子。

### 10.3.3 `Observer`

积分过程中**记录轨迹**的回调：

```cpp
class Observer {
  void observe(const vector_t& state, scalar_t time);   // 每个接受的步调用
};
```

配合 odeint 的 `integrate_adaptive` 或 `integrate_times` 使用。
`integrate_times` 用于需要**指定输出时刻**的场景（MPC 的固定网格输出）。

### 10.3.4 `SystemEventHandler` 与 `StateTriggeredEventHandler`

```cpp
class SystemEventHandler {
  std::atomic_bool killIntegration_{false};
  virtual int checkEvent(const vector_t& state, scalar_t time);
  // 返回非零 → 抛异常中断积分
};
```

`killIntegration_` 是**线搜索并行的关键**：
`LineSearchStrategy` 找到最优步长后调 `abortRollout()` 置位这个标志，
其他线程的积分会在下一步检查时抛异常退出，不再浪费算力。

`StateTriggeredEventHandler` 额外检查保护面符号：

```cpp
int checkEvent(const vector_t& state, scalar_t time) override {
  const vector_t guard = systemPtr_->computeGuardSurfaces(time, state);
  for (i) if (guard[i] < 0 && lastGuard_[i] > 0) return i+1;   // 零穿越
  lastGuard_ = guard;
  return 0;
}
```

---

## 10.4 ⭐ `TimeTriggeredRollout`

**文件**：`ocs2_oc/src/rollout/TimeTriggeredRollout.cpp`

**算法**：

```
输入: [t0, tf], x0, controller, modeSchedule
输出: timeTrajectory, stateTrajectory, inputTrajectory, postEventIndices

1. 从 modeSchedule 取出 [t0, tf] 内的事件时刻 {τ_1, …, τ_K}
2. 切分区间: [t0, τ_1], [τ_1, τ_2], …, [τ_K, tf]
3. for 每个区间 [a, b]:
     3a. 积分 ẋ = f(t,x,u(t,x)) 从 a 到 b，记录轨迹
     3b. if b 是事件时刻:
           记录 postEventIndices.push_back(当前轨迹长度)
           x ← computeJumpMap(b, x)
           把 (b, x^+) 也压入轨迹      ← 时间值重复！
4. 计算 inputTrajectory = controller->computeInput(t, x) 逐点
```

**关键细节：事件时刻有两个采样点**

```
时间:  … t_{k-1}   τ    τ    t_{k+2} …
索引:  …   k-1    k    k+1    k+2   …
状态:  …   x     x^-   x^+     x    …
                  ↑     ↑
        postEventIndices[j] = k+1
```

于是 `postEventIndices[j] - 1` 是 $x^-$ 的索引，`postEventIndices[j]` 是 $x^+$ 的索引。
这个约定在整个代码库里一致使用（例如 `GaussNewtonDDP.cpp:671`
`const size_t preEventIndex = postEventIndices_[timeIndex] - 1;`）。

**为什么不合并成一个点？** 因为 LQ 近似需要在 $x^-$ 处线性化跳变映射
（`approximatePreJumpLQ`），而 Riccati 后向传播需要在 $x^+$ 处的值函数。
两个点承载不同的信息。

---

## 10.5 ⭐ `StateTriggeredRollout`

**文件**：`ocs2_oc/src/rollout/StateTriggeredRollout.cpp`

比时间触发复杂得多，因为事件时刻**未知**。

### 算法

```
1. 从 t0 开始积分，StateTriggeredEventHandler 每步检查保护面符号
2. 检测到 g_i 在 [t_a, t_b] 内变号 → 抛异常，积分中断
3. 用 RootFinder 在 [t_a, t_b] 内求根:
     3a. setInitBracket({t_a, t_b}, {g(t_a), g(t_b)})
     3b. 循环:
           t_new = display()              // 试位法给出新试探点
           积分到 t_new，算 g(t_new)
           updateBracket(t_new, g(t_new)) // 更新括号
         直到 |t_b - t_a| < eps 或 |g| < eps
4. 积分到 τ*，施加跳变映射，记录事件
5. 更新 modeSchedule（写回实际检测到的事件时刻！）
6. 从 x^+ 继续，回到步骤 1
```

**注意步骤 5**：`RolloutBase::run` 的 `ModeSchedule&` 是**非 const 引用**，
状态触发模式下事件时刻是仿真的**输出**而非输入。

### ⭐ `RootFinder`：试位法家族的完整推导

**文件**：`ocs2_oc/include/ocs2_oc/rollout/RootFinder.h`

#### 基本试位法（Regula Falsi）

给定括号 $[t_0,t_1]$，$f_0=g(t_0)$、$f_1=g(t_1)$ 异号。
过两点作直线，取其零点：

$$
\frac{f_1-f_0}{t_1-t_0}=\frac{0-f_1}{t_2-t_1}
\;\Longrightarrow\;
\boxed{\;t_2=t_1-f_1\frac{t_1-t_0}{f_1-f_0}\;}
$$

更新括号：若 $f_2$ 与 $f_1$ 异号，新括号 $[t_1,t_2]$；否则 $[t_0,t_2]$。

#### 停滞问题

若 $g$ 是凸函数且根在右侧，则每次迭代 $t_0$ 端点**永不更新**，
$f_0$ 保持不变。收敛退化为线性，且收敛因子接近 1（非常慢）。

#### 改进：缩放保留端点的函数值

所有改进算法都是同一个思路：当端点 $t_0$ 被保留时，
把 $f_0$ 乘以缩放因子 $\gamma\in(0,1)$，"人为拉近"该端点，
迫使下一次试探点向该侧移动。

| 算法 | $\gamma$ | 说明 |
|---|---|---|
| **Illinois** | $\tfrac12$ | 最简单，固定减半 |
| **Pegasus** | $\dfrac{f_1}{f_1+f_2}$ | 用新旧函数值的比例 |
| **Anderson–Björck** | $\begin{cases}m,&m>0\\ \tfrac12,&\text{否则}\end{cases}$，$m=1-\dfrac{f_2}{f_1}$ | 基于二阶信息 |

**Anderson–Björck 的推导思路**：
$m = 1-f_2/f_1$ 估计了函数在该区间的"弯曲程度"。
若 $f_2$ 与 $f_1$ 同号且 $|f_2|<|f_1|$（正常收敛），则 $m\in(0,1)$，
缩放量与收敛速度自适应匹配。若 $m\le0$（$f_2$ 比 $f_1$ 还大，异常），
退回 Illinois 的保守值 $\frac12$。

**收敛阶**（Ford 1995 的分析）：

| 算法 | 收敛阶 |
|---|---|
| Regula Falsi | 1（线性） |
| Illinois | ≈1.442 |
| Pegasus | ≈1.642 |
| Anderson–Björck | ≈1.7 |

（对比：二分法 = 1，牛顿法 = 2，割线法 ≈1.618）

**OCS2 默认用 Anderson–Björck**（`RootFinder.h:74`），
因为它在保持括号性（保证收敛）的前提下速度最快。

**为什么不用牛顿法？** 因为需要 $\partial g/\partial t$，
而保护面对时间的导数需要 $\frac{\partial g}{\partial x}\dot x$，
额外的雅可比计算不划算；且牛顿法不保证括号，可能跳出区间。

---

## 10.6 `InitializerRollout`

**文件**：`ocs2_oc/include/ocs2_oc/rollout/InitializerRollout.h`

不做积分，直接调 `Initializer::compute(t, x, tNext, u, xNext)` 逐步生成轨迹。
用于：
- MPC 第一次运行（无历史控制器）
- horizon 延长产生的新时间区间

四足的 `LeggedRobotInitializer` 把机器人重量按支撑腿数均分作为初始接触力：

$$
f_i^z = \frac{mg}{n_{\text{stance}}},\qquad f_i^{x,y}=0
$$

这比 $u=0$（自由落体）好得多——初始 rollout 至少不会立刻发散。

---

## 10.7 ⭐ `TrajectorySpreading`：mode schedule 变化时的轨迹适配

**文件**：`ocs2_oc/include/ocs2_oc/trajectory_adjustment/TrajectorySpreading.h`

### 问题

MPC 每个周期：
- 时间前进了 $\Delta t_{mpc}$
- 步态相位推进，可能有事件被"消费掉"，horizon 末尾出现新事件
- 用户可能改变步态（trot → walk）

上一次的解（primal + dual）是在**旧 mode schedule** 下的。
如果直接拿来热启动：
- 事件时刻对不上 → LQ 近似在错误的位置线性化跳变
- 乘子轨迹的事件索引错位 → 增广拉格朗日项张冠李戴

### 解法

```cpp
Status set(const ModeSchedule& oldModeSchedule,
           const ModeSchedule& newModeSchedule,
           const scalar_array_t& oldTimeTrajectory);
```

算法：
1. **匹配模态序列**：找出新旧 `modeSequence` 的最长公共前缀
   （从当前时刻开始，因为过去的模态已经确定）
2. **时间重映射**：对每对匹配的模态区间
   $[\tau_i^{old},\tau_{i+1}^{old}]\to[\tau_i^{new},\tau_{i+1}^{new}]$，
   做**仿射映射**

   $$
   t^{new}=\tau_i^{new}+\frac{\tau_{i+1}^{new}-\tau_i^{new}}{\tau_{i+1}^{old}-\tau_i^{old}}\big(t^{old}-\tau_i^{old}\big)
   $$

   实现上不做插值，而是**重新分配采样点的索引**（"spreading"）：
   把旧轨迹的采样点按比例映射到新的事件划分上。
3. **截断**：无法匹配的尾部（新序列出现了旧序列没有的模态）被丢弃，
   `Status::willTruncate = true`

### 接口

```cpp
template <typename T> void adjustTrajectory(std::vector<T>& trajectory) const;
template <typename T> std::vector<T> extractEventsArray(const std::vector<T>& array) const;
void adjustTimeTrajectory(scalar_array_t& timeTrajectory) const;
const size_array_t& getPostEventIndices() const;
```

模板设计使得同一个 spreading 策略可以应用到状态、输入、乘子、ModelData 等任意数据。

`TrajectorySpreadingHelperFunctions.h` 提供对 `PrimalSolution`/`DualSolution` 的
一站式封装：

```cpp
TrajectorySpreading::Status trajectorySpread(const ModeSchedule& oldMS, const ModeSchedule& newMS,
                                             PrimalSolution& primalSolution);
void trajectorySpread(const TrajectorySpreading&, const DualSolution& src, DualSolution& dst);
```

### 调用点

| 位置 | 用途 |
|---|---|
| `LineSearchStrategy::computeSolution`（`LineSearchStrategy.cpp:103`） | 线搜索时调整 dual solution |
| `LevenbergMarquardtStrategy::run`（`:88`） | 同上 |
| `SqpSolver::runImpl`（`SqpSolver.cpp:202`） | 调整上次的 primal solution 用于热启动 |
| `SlpSolver::runImpl`、`IpmSolver::runImpl` | 同上 |

---

## 10.8 `TimeDiscretization`：多重打靶的网格

**文件**：`ocs2_oc/include/ocs2_oc/oc_data/TimeDiscretization.h`

```cpp
struct AnnotatedTime {
  scalar_t time;
  enum class Event { None, PreEvent, PostEvent } event;
};

std::vector<AnnotatedTime> timeDiscretizationWithEvents(
    scalar_t initTime, scalar_t finalTime, scalar_t dt, const scalar_array_t& eventTimes);
```

**算法**：
1. 生成均匀网格 $t_0, t_0+\Delta t, t_0+2\Delta t,\dots$
2. 对每个事件时刻 $\tau_i$：
   - 找到它落在哪个区间
   - **插入一对节点** $(\tau_i, \texttt{PreEvent})$ 与 $(\tau_i, \texttt{PostEvent})$
   - 如果 $\tau_i$ 离某个网格点很近（小于某阈值），**合并**以避免超短区间

**为什么必须对齐？** 因为多重打靶的每个区间内要做 RK 积分。
如果事件落在区间内部，该区间的动力学是不连续的，RK 会给出完全错误的结果。

**辅助函数**：

```cpp
scalar_t getIntervalStart(const AnnotatedTime& t);
    // PostEvent 节点的"区间起点"是它自己（事件后立即开始）
scalar_t getIntervalDuration(const AnnotatedTime& start, const AnnotatedTime& end);
    // PreEvent 到 PostEvent 的"时长"为 0
```

在 `SqpSolver::setupQuadraticSubproblem` 里：

```cpp
if (time[i].event == AnnotatedTime::Event::PreEvent) {
  auto result = multiple_shooting::setupEventNode(ocpDefinition, time[i].time, x[i], x[i+1]);
  // 跳变映射，dt 不参与
} else {
  const scalar_t ti = getIntervalStart(time[i]);
  const scalar_t dt = getIntervalDuration(time[i], time[i+1]);
  auto result = multiple_shooting::setupIntermediateNode(..., ti, dt, x[i], x[i+1], u[i]);
}
```

---

## 10.9 `TrapezoidalIntegration`

**文件**：`ocs2_core/include/ocs2_core/integration/TrapezoidalIntegration.h`

沿离散轨迹积分标量或向量：

$$
\int_{t_0}^{t_N}\!\! v(t)\,\mathrm{d}t\ \approx\
\sum_{k=0}^{N-1}\frac{t_{k+1}-t_k}{2}\big(v_k+v_{k+1}\big)
$$

用于计算 `PerformanceIndex` 的积分项（代价、约束 SSE）。
比前向欧拉精度高一阶（$O(\Delta t^2)$ vs $O(\Delta t)$），
且对非均匀网格（含事件）自然适用。

**注意**：这与多重打靶转写里的代价离散化（前向欧拉，见 [03](03-ocs2-oc.md) §3.5.2）
**不一致**。这是有意的：
- **转写用欧拉**：与 QP 的线性化保持一致，梯度精确对应
- **评估用梯形**：更准确地反映真实代价

若两者不一致导致线搜索行为异常，可以检查这一点。

---

## 10.10 `PerformanceIndicesRollout`

**文件**：`ocs2_oc/include/ocs2_oc/rollout/PerformanceIndicesRollout.h`

沿给定的 $(t,x,u)$ 轨迹计算 `PerformanceIndex`，
用梯形积分处理中间项，逐点累加事件项与终端项。

DDP 的线搜索每试一个步长都要调它一次，是**热点函数**。

---

## 10.11 事件对 Riccati 的影响

### 横截条件

**文件**：`ocs2_ddp/include/ocs2_ddp/riccati_equations/RiccatiTransversalityConditions.h`

值函数在事件处的传播。设跳变映射线性化 $\delta x^+ = A_j\delta x^- + b_j$，
事件代价 $\phi^{(i)}(x^-)$，则

$$
V^-(\delta x^-)=\phi^{(i)}(\delta x^-)+V^+\big(A_j\delta x^-+b_j\big)
$$

展开：

$$
\boxed{
\begin{aligned}
S_m^- &= \phi^{(i)}_{xx} + A_j^{\!\top}S_m^+A_j\\
S_v^- &= \phi^{(i)}_x + A_j^{\!\top}\big(S_v^+ + S_m^+b_j\big)\\
s^-   &= \phi^{(i)} + s^+ + S_v^{+\top}b_j + \tfrac12 b_j^{\!\top}S_m^+b_j
\end{aligned}}
$$

**这与 [04 章](04-ocs2-ddp.md) §4.3.2 的离散 Riccati 形式完全一致**——
因为跳变映射就是一个"没有输入的离散动力学步"。

在 SLQ 里，Riccati 积分到事件时刻时暂停，施加上述跳变，再继续积分
（`ContinuousTimeRiccatiEquations.cpp:135` 的 `computeJumpMap`）。

在 iLQR 里，事件节点在递推中被特殊处理
（`ILQR.cpp:227` 的 `riccatiEquationsWorker`，用 `finalValueTemp` 保存 pre-jump 值）。

---

## 10.12 实用建议

### 调试事件相关的问题

1. **打开 `debugPrintRollout_`**：会打印每个事件的检测时刻与状态跳变
2. **检查 `postEventIndices_`**：长度应等于 horizon 内的事件数
3. **检查时间轴是否有重复值**：事件时刻必须有两个采样点
4. **`getNumFunctionCalls()` 异常大**：通常是保护面在零附近抖动
   （Zeno 现象），需要加滞回或最小驻留时间

### 状态触发的陷阱

- **保护面必须光滑**：`RootFinder` 假设区间内单调，不连续的保护面会失败
- **多个保护面同时穿越**：代码只处理第一个，可能漏事件
- **Zeno 行为**：事件间隔趋于 0（如弹跳球停止时），积分器会卡死。
  实践中需要在动力学层面加"最小驻留时间"或检测到 Zeno 后切换模态

### 时间触发的陷阱

- **事件时刻必须在 horizon 内**：`ModeSchedule` 应覆盖 $[t_0, t_f]$
- **`dt` 与事件间隔的关系**：若两个事件间隔 < `dt`，多重打靶会产生退化区间。
  `timeDiscretizationWithEvents` 有合并逻辑，但极端情况仍需注意

---

## 10.13 参考文献

- **Anderson, N. & Björck, Å. (1973)**, *A New High Order Method of Regula Falsi Type for Computing a Root of an Equation*, BIT 13:253–264
- **Dowell, M. & Jarratt, P. (1971)**, *A Modified Regula Falsi Method for Computing the Root of an Equation*, BIT 11:168–174 —— **Illinois 方法**
- **Dowell, M. & Jarratt, P. (1972)**, *The "Pegasus" Method for Computing the Root of an Equation*, BIT 12:503–508
- **Ford, J.A. (1995)**, *Improved Algorithms of Illinois-Type for the Numerical Solution of Nonlinear Equations*, Tech. Report, Univ. of Essex
  （以上四篇在 `RootFinder.h:51-63` 的注释中被直接引用）
- **Hairer, Nørsett, Wanner (1993)**, *Solving Ordinary Differential Equations I: Nonstiff Problems* —— Dormand-Prince 与自适应步长
- **Farshidian et al. (2017)**, *Sequential Linear Quadratic Optimal Control for Nonlinear Switched Systems*, IFAC —— 切换系统的 SLQ 与横截条件
- **Van der Schaft & Schumacher (2000)**, *An Introduction to Hybrid Dynamical Systems* —— 混合系统理论、Zeno 现象

**下一章**：MPC 与 ROS → [11 MPC/ROS](11-mpc-ros.md)
