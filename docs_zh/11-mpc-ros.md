# 11 · ocs2_mpc 与 ROS 接口

前面几章讲的是"给定 $[t_0,t_f]$ 与 $x_0$ 求一次最优轨迹"。
本章讲**如何把它变成实时闭环控制器**。

---

## 11.1 MPC 的双线程架构

OCS2 的实时架构基于一个关键观察：

> **MPC 求解慢（10~100ms），控制回路快（1~5ms）。**
> 两者必须解耦，否则控制回路会被求解阻塞。

```
┌────────────────────────────────────────────────────────────┐
│  MPC 线程（慢，~10-100 Hz）                                 │
│                                                             │
│   observation ──▶ MPC_BASE::run(t, x)                       │
│                     └─▶ solver->run(t, x, t+T)              │
│                            （SLQ/SQP/IPM/SLP）              │
│                     └─▶ PrimalSolution                      │
│                              │                              │
│                              ▼  写入 buffer（加锁）          │
│                     ┌──────────────────┐                    │
│                     │  policy buffer   │                    │
│                     └──────────────────┘                    │
└──────────────────────────────│─────────────────────────────┘
                               │ try_lock + swap（非阻塞）
┌──────────────────────────────▼─────────────────────────────┐
│  控制线程（快，~100-1000 Hz）                               │
│                                                             │
│   MRT_BASE::updatePolicy()      ← 尝试取新策略（失败就用旧的）│
│   MRT_BASE::evaluatePolicy(t, x, …)  → u = K(t)x + b(t)     │
│                                  ──▶ 发给执行器             │
│   MRT_BASE::setCurrentObservation(obs)  ──▶ 送回 MPC 线程   │
└────────────────────────────────────────────────────────────┘
```

**MRT** = **M**odel **R**eference **T**racking，
指"跟踪 MPC 给出的参考模型（策略）"的那一侧。

---

## 11.2 `MPC_BASE`

**文件**：`ocs2_mpc/include/ocs2_mpc/MPC_BASE.h:44`

```cpp
class MPC_BASE {
 public:
  explicit MPC_BASE(mpc::Settings mpcSettings);
  virtual bool run(scalar_t currentTime, const vector_t& currentState);
  virtual SolverBase* getSolverPtr() = 0;
  scalar_t getTimeHorizon() const { return mpcSettings_.timeHorizon_; }
 protected:
  virtual void calculateController(scalar_t initTime, const vector_t& initState, scalar_t finalTime) = 0;
  bool isFirstMpcRun() const { return initRun_; }
 private:
  bool initRun_ = true;
  const mpc::Settings mpcSettings_;
  benchmark::RepeatedTimer mpcTimer_;
};
```

`run()` 的实现（`MPC_BASE.cpp:53`）非常短：

```cpp
bool MPC_BASE::run(scalar_t currentTime, const vector_t& currentState) {
  if (!initRun_ && currentTime >= getSolverPtr()->getFinalTime()) {
    std::cerr << "WARNING: The MPC time-horizon is smaller than the MPC starting time.\n";
    return false;
  }
  const scalar_t finalTime = currentTime + mpcSettings_.timeHorizon_;
  if (mpcSettings_.debugPrint_) mpcTimer_.startTimer();
  calculateController(currentTime, currentState, finalTime);   // ← 纯虚，子类实现
  initRun_ = false;
  if (mpcSettings_.debugPrint_) { mpcTimer_.endTimer(); /* 打印 max/avg/last */ }
  return true;
}
```

**滚动时域**：每次调用都用 $[t, t+T]$ 作为新的优化区间，
$T=$ `timeHorizon_` 固定。这就是 receding horizon。

### 具体实现

| 类 | 文件 | 求解器 |
|---|---|---|
| `GaussNewtonDDP_MPC` | `ocs2_ddp/include/ocs2_ddp/GaussNewtonDDP_MPC.h` | SLQ 或 iLQR |
| `SqpMpc` | `ocs2_sqp/ocs2_sqp/include/ocs2_sqp/SqpMpc.h` | SQP |
| `IpmMpc` | `ocs2_ipm/include/ocs2_ipm/IpmMpc.h` | IPM |
| `SlpMpc` | `ocs2_slp/include/ocs2_slp/SlpMpc.h` | SLP |

它们的 `calculateController` 几乎相同（以 DDP 为例，`GaussNewtonDDP_MPC.h:69`）：

```cpp
void calculateController(scalar_t initTime, const vector_t& initState, scalar_t finalTime) override {
  if (settings().coldStart_) ddpPtr_->reset();   // 冷启动：丢弃上次的解
  ddpPtr_->run(initTime, initState, finalTime);
}
```

**`coldStart_` 的取舍**：
- `false`（默认）：**热启动**，用上次的解初始化 → 通常 1~3 次迭代即收敛
- `true`：每次从零开始 → 慢得多，但避免"陷在上次的局部极小里"

对周期性任务（步行）热启动效果极好；
对突变任务（急停、避障）偶尔冷启动能跳出局部解。

---

## 11.3 `mpc::Settings`

**文件**：`ocs2_mpc/include/ocs2_mpc/MPC_Settings.h`

| 参数 | 默认 | 说明 |
|---|---|---|
| `timeHorizon_` | 1.0 | 预测时域 $T$（秒） |
| `solutionTimeWindow_` | -1 | 发布给 MRT 的解的时间窗口；**负值表示发布整个 horizon** |
| `debugPrint_` | false | 打印每次求解的耗时统计 |
| `coldStart_` | false | 是否每次重置求解器 |
| `mpcDesiredFrequency_` | -1 | MPC 目标频率（Hz）；**负值表示"尽可能快"** |
| `mrtDesiredFrequency_` | 100.0 | MRT 目标频率（Hz），用于 dummy 仿真循环 |

### `solutionTimeWindow_` 的作用

MPC 算出的是 $[t, t+T]$ 上的完整轨迹，但通常只需要发布前面一小段：

```cpp
// MPC_ROS_Interface.cpp:220
scalar_t finalTime = mpcInitObservation.time + mpc_.settings().solutionTimeWindow_;
if (mpc_.settings().solutionTimeWindow_ < 0) {
  finalTime = mpc_.getSolverPtr()->getFinalTime();     // 发布全部
}
```

**为什么要截断？**
1. **减少通信量**：ROS 消息更小，序列化更快
2. **避免使用过期的远期规划**：MPC 会在下个周期重新规划，
   远期部分本来就会被丢弃

**为什么默认发布全部（-1）？**
因为如果 MPC 突然变慢（一次求解超时），MRT 需要有足够长的策略撑到下次更新。
截断太短会导致"策略用完了"。

### `mpcDesiredFrequency_` 的取舍

- **负值（默认）**：MPC 线程一算完立刻开始下一次
  → 最大化利用算力，但 MPC 频率不确定（随问题难度波动）
- **正值**：用 `ExecuteAndSleep` 限频到指定值
  → 频率稳定，可预测；剩余算力留给其他任务

四足机器人通常设为负值（能算多快算多快）。

---

## 11.4 `MRT_BASE`

**文件**：`ocs2_mpc/include/ocs2_mpc/MRT_BASE.h:58`

```cpp
class MRT_BASE {
 public:
  virtual void setCurrentObservation(const SystemObservation& observation) = 0;

  bool updatePolicy();                                    // 尝试从 buffer 取新策略
  const PrimalSolution& getPolicy() const;
  const CommandData& getCommand() const;
  const PerformanceIndex& getPerformanceIndices() const;

  void evaluatePolicy(scalar_t t, const vector_t& x, vector_t& mpcState, vector_t& mpcInput, size_t& mode);
  void rolloutPolicy(scalar_t t, const vector_t& x, const scalar_t& dt,
                     vector_t& mpcState, vector_t& mpcInput, size_t& mode);
  void initRollout(const RolloutBase* rolloutPtr);

  bool initialPolicyReceived() const;
  void reset();
 protected:
  virtual void modifyActiveSolution(const CommandData&, PrimalSolution&) {}
  // buffer / active 双份数据
  std::unique_ptr<CommandData>     activeCommandPtr_,           bufferCommandPtr_;
  std::unique_ptr<PrimalSolution>  activePrimalSolutionPtr_,    bufferPrimalSolutionPtr_;
  std::unique_ptr<PerformanceIndex> activePerformanceIndicesPtr_, bufferPerformanceIndicesPtr_;
  std::mutex bufferMutex_;
  std::atomic_bool newPolicyInBuffer_;
};
```

### ⭐ `updatePolicy()`：非阻塞的策略交换

**文件**：`ocs2_mpc/src/MRT_BASE.cpp:156`

```cpp
bool MRT_BASE::updatePolicy() {
  std::unique_lock<std::mutex> lock(bufferMutex_, std::try_to_lock);   // ← 关键：try_lock
  if (lock.owns_lock()) {
    mrtTrylockWarningCount_ = 0;
    if (newPolicyInBuffer_) {
      activeCommandPtr_.swap(bufferCommandPtr_);                       // O(1) 指针交换
      activePrimalSolutionPtr_.swap(bufferPrimalSolutionPtr_);
      activePerformanceIndicesPtr_.swap(bufferPerformanceIndicesPtr_);
      newPolicyInBuffer_ = false;
      modifyActiveSolution(*activeCommandPtr_, *activePrimalSolutionPtr_);
      return true;
    }
    return false;      // buffer 里没新东西
  } else {
    ++mrtTrylockWarningCount_;
    if (mrtTrylockWarningCount_ > mrtTrylockWarningThreshold_)
      std::cerr << "[MRT_BASE::updatePolicy] failed to lock the policyBufferMutex for "
                << mrtTrylockWarningCount_ << " consecutive times.\n";
    return false;      // 锁被 MPC 线程占着，继续用旧策略
  }
}
```

**`try_to_lock` 是整个设计的核心**：
控制线程**绝不阻塞**。拿不到锁就继续用上一个策略——
反正 MPC 策略在短时间内变化不大，用旧的一个周期完全可以接受。

**连续拿不到锁的告警**：如果连续 $N$ 次 try_lock 失败，
说明 MPC 线程持锁时间过长（可能是在锁内做了拷贝而非交换），需要排查。

**`modifyActiveSolution` 钩子**：子类可以在策略切换时做额外处理，
如 loopshaping 的滤波器状态同步。

### `evaluatePolicy` vs `rolloutPolicy`

```cpp
// evaluatePolicy: 直接查表 + 插值
void MRT_BASE::evaluatePolicy(scalar_t t, const vector_t& x, vector_t& mpcState, vector_t& mpcInput, size_t& mode) {
  mpcInput = activePrimalSolutionPtr_->controllerPtr_->computeInput(t, x);   // u = K(t)x + b(t)
  mpcState = LinearInterpolation::interpolate(t, timeTrajectory_, stateTrajectory_);
  mode = activePrimalSolutionPtr_->modeSchedule_.modeAtTime(t);
}

// rolloutPolicy: 用策略做一步前向仿真
void MRT_BASE::rolloutPolicy(scalar_t t, const vector_t& x, const scalar_t& dt, ...) {
  rolloutPtr_->run(t, x, t + dt, controllerPtr, modeSchedule, timeTraj, postEventIdx, stateTraj, inputTraj);
  mpcState = stateTraj.back();
  mpcInput = inputTraj.back();
  mode = modeSchedule.modeAtTime(t + dt);
}
```

| | `evaluatePolicy` | `rolloutPolicy` |
|---|---|---|
| 用途 | **真实机器人**：直接算控制量 | **仿真**：模拟机器人状态演化 |
| 开销 | 极小（一次插值 + 矩阵乘） | 一次积分 |
| 状态来源 | 传感器（$x$ 是测量值） | 仿真积分结果 |

`MRT_ROS_Dummy_Loop` 用 `rolloutPolicy` 做无物理引擎的"假机器人"。

### 具体实现

| 类 | 场景 |
|---|---|
| `MPC_MRT_Interface` | **单进程**：MPC 与控制在同一程序的两个线程 |
| `MRT_ROS_Interface` | **跨进程**：通过 ROS 话题通信 |
| `MpcnetDummyLoopRos` | MPC-Net 策略的执行循环 |

---

## 11.5 `MPC_MRT_Interface`：单进程双线程

**文件**：`ocs2_mpc/include/ocs2_mpc/MPC_MRT_Interface.h:50`

```cpp
class MPC_MRT_Interface final : public MRT_BASE {
 public:
  explicit MPC_MRT_Interface(MPC_BASE& mpc);
  void setCurrentObservation(const SystemObservation& currentObservation) override;
  void advanceMpc();                                    // ← 在 MPC 线程调用
  matrix_t getLinearFeedbackGain(scalar_t time);
  ScalarFunctionQuadraticApproximation getValueFunction(scalar_t time, const vector_t& state);
  vector_t getStateInputEqualityConstraintLagrangian(scalar_t time, const vector_t& state);
};
```

典型用法：

```cpp
MPC_MRT_Interface mpcInterface(mpc);
// ---- MPC 线程 ----
std::thread mpcThread([&]() {
  while (running) {
    mpcInterface.advanceMpc();      // 内部：run() + 更新 buffer
  }
});
// ---- 控制线程 ----
while (running) {
  mpcInterface.setCurrentObservation(observation);
  if (mpcInterface.initialPolicyReceived()) {
    mpcInterface.updatePolicy();
    mpcInterface.evaluatePolicy(t, x, optimalState, optimalInput, mode);
    sendToActuators(optimalInput);
  }
}
```

`advanceMpc()` 内部（`MPC_MRT_Interface.cpp`）：
1. 读当前 observation（加锁）
2. 调 `mpc_.run(t, x)`
3. 从 solver 取出 `PrimalSolution`、`CommandData`、`PerformanceIndex`
4. 写入 buffer（加锁），置 `newPolicyInBuffer_ = true`

**这是延迟最低的方案**（无序列化、无网络），
适合真实机器人的板载控制器。

---

## 11.6 `ocs2_ros_interfaces`：跨进程

### 11.6.1 `MPC_ROS_Interface`

**文件**：`ocs2_ros_interfaces/src/mpc/MPC_ROS_Interface.cpp`

MPC 侧的 ROS 节点封装：

```
订阅:  <robot>_mpc_observation      (ocs2_msgs/mpc_observation)
发布:  <robot>_mpc_policy           (ocs2_msgs/mpc_flattened_controller)
服务:  <robot>_mpc_reset            (ocs2_msgs/reset)
```

**"flattened controller"**：`LinearController` 的增益矩阵与偏置被展平成一维数组
传输，接收端重建。见 `ocs2_ros_interfaces/common/RosMsgConversions.h`。

工作流程：
1. `mpcObservationCallback` 收到观测 → 存入 buffer
2. `spin()` 循环里：取 buffer → `mpc_.run()` → 发布策略
3. `solutionTimeWindow_` 决定发布多长的轨迹（`:220`）

### 11.6.2 `MRT_ROS_Interface`

MRT 侧：订阅策略话题，发布观测话题。
收到策略消息后反序列化成 `PrimalSolution`，写入 buffer。

**支持 ROS 的独立 spinner 线程**，
使得 `updatePolicy()` 的语义与单进程版本完全一致。

### 11.6.3 `MRT_ROS_Dummy_Loop`

**文件**：`ocs2_ros_interfaces/include/ocs2_ros_interfaces/mrt/MRT_ROS_Dummy_Loop.h`

无物理引擎的"假机器人"：用 `rolloutPolicy` 积分 MPC 给出的动力学模型，
把结果当作"真实状态"发布出去。

```cpp
void MRT_ROS_Dummy_Loop::run(const SystemObservation& initObservation,
                             const TargetTrajectories& initTargetTrajectories);
```

**用途**：
- 快速验证 MPC 是否收敛、轨迹是否合理
- RViz 可视化
- 不依赖 Gazebo/RaiSim，启动快

**局限**：因为用的就是 MPC 自己的模型，**不存在模型误差**，
所以看不出鲁棒性问题。真实评估需要 RaiSim（`ocs2_raisim`）或 Gazebo。

`DummyObserver` 接口让用户在每步注入可视化逻辑
（各机器人的 `*Visualizer` 类都实现它）。

### 11.6.4 指令接口

| 类 | 用途 |
|---|---|
| `TargetTrajectoriesRosPublisher` | 通用的目标轨迹发布 |
| `TargetTrajectoriesKeyboardPublisher` | 键盘输入 → 目标（终端交互） |
| `TargetTrajectoriesInteractiveMarker` | RViz 交互标记 → 目标（鼠标拖拽） |

`RosReferenceManager` 是 `ReferenceManagerDecorator` 的实现：
订阅 ROS 话题，把收到的 `TargetTrajectories`/`ModeSchedule` 注入求解器。

### 11.6.5 可视化工具

`ocs2_ros_interfaces/visualization/VisualizationHelpers.h`：
生成 RViz marker 的辅助函数（箭头、球、线条、坐标系）。
`VisualizationColors.h` 定义配色。

`multiplot/` 目录下是 `rqt_multiplot` 的配置文件，
用于实时绘制求解器的性能指标。

### 11.6.6 `SolverObserverRosCallbacks`

**文件**：`ocs2_ros_interfaces/synchronized_module/SolverObserverRosCallbacks.h`

工厂函数，创建把**指定约束项的乘子/指标**发布到 ROS 话题的回调。
用法：

```cpp
auto observer = SolverObserver::ConstraintTermObserver(
    SolverObserver::Type::Intermediate, "frictionCone_LF",
    ros::createConstraintCallback(nodeHandle, "metric_topic", ...));
solver.addSolverObserver(observer);
```

于是可以在 rqt_plot 里实时看某条腿的摩擦锥违反量——
**调试机器人 MPC 的利器**。

---

## 11.7 `SystemObservation` 与 `CommandData`

```cpp
struct SystemObservation {
  size_t   mode = 0;      // 当前模态（接触状态）
  scalar_t time = 0.0;
  vector_t state;
  vector_t input;         // 上一时刻实际施加的输入（用于热启动）
};

struct CommandData {
  SystemObservation  mpcInitObservation;   // 这次 MPC 求解的初始观测
  TargetTrajectories mpcTargetTrajectories; // 这次求解用的目标
};
```

`CommandData` 与 `PrimalSolution` 一起传给 MRT，
使得控制侧知道"这个策略是基于什么时刻、什么目标算出来的"，
可以据此判断策略是否过期。

---

## 11.8 `ocs2_python_interface`

**文件**：`ocs2_python_interface/include/ocs2_python_interface/PythonInterface.h`

pybind11 绑定，暴露：

```python
interface.setObservation(t, x, u)
interface.advanceMpc()
interface.getMpcSolution()          # → (t_array, x_array, u_array)
interface.computeFlowMap(t, x, u)
interface.getLinearFeedbackGain(t)
interface.getValueFunction(t, x)
interface.getHamiltonian(t, x, u)   # ← MPC-Net 训练用
```

`PythonInterface.h` 是通用基类，各机器人有自己的绑定
（如 `ocs2_ballbot/BallbotPyBindings.h`）。

`ocs2_mpcnet_core/MpcnetPybindMacros.h` 提供宏，
一行代码生成一个机器人的完整 Python 模块。

**用途**：
- 快速原型验证（Jupyter 里调参）
- MPC-Net 的数据生成（Python 训练循环调 C++ MPC）
- 与 RL 框架（Gym）集成

---

## 11.9 实时性工程要点

### 线程优先级

```cpp
ocs2::setThreadPriority(settings.threadPriority, thread);   // SCHED_FIFO
```

MPC 线程通常设 50~99，控制线程设更高（99）。
需要 `CAP_SYS_NICE` 权限或 `rtprio` 的 ulimit 配置。

### 内存分配

**实时循环里绝不能有 malloc**。OCS2 的做法：
- 预分配所有轨迹缓冲（`resize` 只在维度变化时调用）
- 用 `swap` 而非拷贝交换数据
- Eigen 的 `noalias()` 避免临时对象

`GaussNewtonDDP` 的 `nominal/cached/optimized` 三缓冲轮转
（[04 章](04-ocs2-ddp.md) §4.9.1）就是这个原则的体现。

### 时序诊断

打开 `debugPrint_` 后 `MPC_BASE` 会打印：

```
### MPC Benchmarking
###   Maximum : 45.2[ms].
###   Average : 12.3[ms].
###   Latest  : 11.8[ms].
```

**关注 Maximum 而非 Average**——
实时系统的瓶颈是最坏情况。若 Maximum 远大于 Average，
通常是某次迭代触发了重新分配内存或线搜索退化。

各求解器的 `getBenchmarkingInfo()` 给出更细的分解
（LQ 近似 / 后向传播 / 线搜索各占多少）。

---

## 11.10 参考文献

- **Farshidian, Jelavic, Satapathy, Giftthaler, Buchli (2017)**,
  *Real-Time Motion Planning of Legged Robots: A Model Predictive Control Approach*, Humanoids
- **Grandia, Farshidian, Ranftl, Hutter (2019)**,
  *Feedback MPC for Torque-Controlled Legged Robots*, IROS —— MRT 侧高频反馈的必要性
- **Rawlings, Mayne, Diehl (2017)**, *Model Predictive Control: Theory, Computation, and Design*, 2nd ed. —— MPC 稳定性理论
- **Diehl, Bock, Schlöder (2005)**, *A Real-Time Iteration Scheme for Nonlinear Optimization in Optimal Feedback Control*, SIAM J. Control Optim.

**下一章**：机器人建模 → [12 机器人模型](12-robot-models.md)
