# 14 · MPC-Net：从 MPC 蒸馏神经网络策略

`ocs2_mpcnet/` 实现 **MPC-Net**：用 MPC 作为"专家"训练神经网络策略。

**动机**：
- MPC 每步要解一个优化问题，**算力需求高**（几十 ms，需要多核 CPU）
- 神经网络推理是**一次前向传播**（几十 μs，可在 MCU 上跑）
- 但纯 RL 训练慢、样本效率低、且难以保证约束

**MPC-Net 的思路**：不用 RL 的奖励信号，而用 **MPC 的哈密顿量**作为损失。
这是"first principles guided policy search"——由最优控制原理直接指导策略学习。

**目录结构**：

```
ocs2_mpcnet/
├── ocs2_mpcnet_core/          核心：C++ 数据生成 + Python 训练框架
│   ├── include/ocs2_mpcnet_core/
│   │   ├── MpcnetDefinitionBase.h        观测/动作变换的机器人特定定义
│   │   ├── MpcnetInterfaceBase.h         Python 绑定的入口
│   │   ├── MpcnetPybindMacros.h          一行生成 Python 模块的宏
│   │   ├── control/                      MpcnetOnnxController、BehavioralController
│   │   ├── rollout/                      数据生成与策略评估
│   │   └── dummy/                        ROS 执行循环
│   └── python/ocs2_mpcnet_core/
│       ├── mpcnet.py                     训练主循环
│       ├── loss/{hamiltonian,behavioral_cloning,cross_entropy}.py
│       ├── policy/{linear,nonlinear,mixture_of_*}.py
│       ├── memory/circular.py            回放缓冲
│       ├── config.py, helper.py
├── ocs2_ballbot_mpcnet/       Ballbot 的 MPC-Net
└── ocs2_legged_robot_mpcnet/  四足的 MPC-Net
```

---

## 14.1 ⭐ 核心思想：用哈密顿量作损失

### 14.1.1 为什么不用行为克隆

**行为克隆（BC）**：$\min_\pi\ \mathbb{E}\big[\|\pi(x)-u^{MPC}(x)\|_R^2\big]$

问题：
1. **不知道 MPC 的输入哪些维度重要**——BC 对所有维度一视同仁，
   但某些输入的偏差代价很大，某些几乎无所谓
2. **分布漂移**：策略偏离专家轨迹后没有纠正信号
3. **需要精确的专家动作**——但 MPC 的解本身只是数值近似

### 14.1.2 哈密顿量损失

**核心洞察**：最优控制的必要条件是**哈密顿量对输入取极小**（庞特里亚金极小值原理）：

$$
u^\star(t,x)=\arg\min_u\ \mathcal{H}(t,x,u)
= \arg\min_u\Big\{l(t,x,u)+\frac{\partial V}{\partial x}^{\!\top}f(t,x,u)\Big\}
$$

所以**直接让策略去最小化哈密顿量**，而不是去模仿 $u^\star$：

$$
\boxed{\;\min_\pi\ \mathbb{E}_{(t,x)\sim\mathcal{D}}\Big[\mathcal{H}\big(t,x,\pi(x)\big)\Big]\;}
$$

**优点**：
- 损失函数**自带各维度的相对重要性**（由 $\mathcal{H}$ 的曲率给出）
- 即使策略偏离专家，$\mathcal{H}$ 仍给出正确的下降方向
- 不需要专家的精确解，只需要哈密顿量的局部近似

### 14.1.3 二次近似

$\mathcal{H}$ 本身难以在线求值（需要值函数梯度）。
MPC-Net 用 MPC 提供的**二次近似**：

$$
\boxed{
\begin{aligned}
\mathcal{H}(x,u)\;\approx\;&
\tfrac12\,\delta x^{\!\top}\mathcal{H}_{xx}\,\delta x
+\delta u^{\!\top}\mathcal{H}_{ux}\,\delta x
+\tfrac12\,\delta u^{\!\top}\mathcal{H}_{uu}\,\delta u\\
&+\mathcal{H}_x^{\!\top}\delta x+\mathcal{H}_u^{\!\top}\delta u+\mathcal{H}
\end{aligned}}
$$

其中 $\delta x=x_{query}-x_{nominal}$、$\delta u=u_{query}-u_{nominal}$，
展开点 $(x_{nominal},u_{nominal})$ 是 MPC 求解时的标称轨迹点。

这正是 `HamiltonianLoss` 的文档字符串
（`loss/hamiltonian.py:42`）里写的公式。

**这些系数从哪来？** `SolverBase::getHamiltonian(t, x, u)`
（[03 章](03-ocs2-oc.md) §3.7），各求解器实现它：
DDP 从 Riccati 解的 $S_m,S_v$ 构造 $\frac{\partial V}{\partial x}$，
SQP 从 HPIPM 的 Riccati 中间结果构造。

**训练时**：策略输出 $u_{query}=\pi(x)$，
代入上式得到标量损失，反向传播即可。
**每个数据点的损失是解析的二次型，无需再调 MPC** ✓

### 14.1.4 与 BC 的对比

`loss/behavioral_cloning.py:42` 也提供了 BC 损失：

$$
BC(u)=\delta u^{\!\top}R\,\delta u
$$

**这可以看作哈密顿量损失的退化情形**：
若令 $\mathcal{H}_{uu}=2R$、其余项为 0，两者等价。
即 BC 假设"所有输入维度的重要性由固定的 $R$ 给出"，
而哈密顿量损失让**这个权重随状态自适应**。

代码保留 BC 是为了做消融实验对照。

---

## 14.2 混合专家网络（Mixture of Experts）

### 14.2.1 动机

四足机器人的最优策略是**分段的**：
trot 的 LF-RH 支撑相与 RF-LH 支撑相，最优控制律完全不同。
单个网络要拟合这种分段函数很困难（需要很大容量，且在切换处震荡）。

**混合专家**：每个模态一个"专家"网络，一个"门控"网络决定权重。

### 14.2.2 结构

**文件**：`policy/mixture_of_linear_experts.py:43`

```python
class MixtureOfLinearExpertsPolicy(BasePolicy):
    def __init__(self, config):
        self.gating_net = torch.nn.Sequential(
            torch.nn.Linear(self.observation_dimension, self.expert_number),
            torch.nn.Softmax(dim=1)
        )
        self.expert_nets = torch.nn.ModuleList(
            [_LinearExpert(i, self.observation_dimension, self.action_dimension)
             for i in range(self.expert_number)]
        )
```

输出：

$$
\boxed{\;\pi(o)=\sum_{i=1}^{E} w_i(o)\,\pi_i(o),\qquad
w(o)=\operatorname{softmax}\big(W_g\,o+b_g\big)\;}
$$

其中 $\pi_i$ 是第 $i$ 个专家（线性或非线性）。

**四种策略类**：

| 类 | 专家 | 门控 |
|---|---|---|
| `LinearPolicy` | — | — （单个线性层） |
| `NonlinearPolicy` | — | — （MLP） |
| `MixtureOfLinearExpertsPolicy` | 线性 | softmax 线性 |
| `MixtureOfNonlinearExpertsPolicy` | MLP | softmax MLP |

**为什么线性专家就够？** 因为 MPC 本身输出的是**时变线性反馈** $u=K(t)x+b(t)$。
在给定模态内，最优控制律近似线性；非线性主要来自模态切换。
混合专家正好把"模态选择"（门控，非线性）与"局部控制律"（专家，线性）分离。

### 14.2.3 门控损失

**文件**：`loss/cross_entropy.py:42`

如果只用哈密顿量损失，门控网络可能学到退化的权重分配
（例如所有输入都用一个专家）。
所以额外用**交叉熵**监督门控：

$$
CE(p_{target}, p_{pred}) = -\sum_{i=1}^{P} p_{target,i}\ln\big(p_{pred,i}+\varepsilon\big)
$$

$p_{target}$ 是 MPC 当前模态的 **one-hot 编码**
（`helper.get_one_hot(data[i].mode, EXPERT_NUM, EXPERT_FOR_MODE)`），
$\varepsilon$ 稳定对数（默认 1e-8）。

**总损失**（`mpcnet.py:291`）：

$$
\mathcal{L}=\mathcal{L}_{experts}+\lambda\,\mathcal{L}_{gating}
$$

```python
empirical_experts_loss = self.experts_loss(x, x, input, u, p, p, dHdxx, dHdux, dHduu, dHdx, dHdu, H)
empirical_gating_loss  = self.gating_loss(x, x, u, u, weights, p, dHdxx, ..., H)
empirical_loss = empirical_experts_loss + self.config.LAMBDA * empirical_gating_loss
```

**`CHEATING` 标志**（`config.CHEATING`）控制是否使用门控损失。
名字里的 "cheating" 指"利用了 MPC 的模态标签"——
这在部署时不可用（策略必须自己判断模态），
但训练时用它加速收敛是合理的。

`EXPERT_FOR_MODE` 是"模态 → 专家索引"的映射，
允许多个模态共享一个专家（例如 trot 的两个相位可以共享，因为对称）。

---

## 14.3 ⭐ 数据生成：DAgger 风格的行为控制器

### 14.3.1 分布漂移问题

纯模仿学习的经典问题：
1. 用专家轨迹训练策略
2. 策略部署后偏离专家分布
3. 遇到训练时没见过的状态 → 输出错误 → 进一步偏离 → 崩溃

**DAgger**（Ross, Gordon, Bagnell 2011）的解法：
**用当前策略生成状态分布，但用专家标注动作**。

### 14.3.2 MPC-Net 的实现

**文件**：`control/MpcnetBehavioralController.cpp:40`

```cpp
vector_t MpcnetBehavioralController::computeInput(scalar_t t, const vector_t& x) {
  if (numerics::almost_eq(alpha_, 0.0)) {
    return learnedControllerPtr_->computeInput(t, x);          // 纯策略
  } else if (numerics::almost_eq(alpha_, 1.0)) {
    return optimalControllerPtr_->computeInput(t, x);          // 纯 MPC
  } else {
    return alpha_ * optimalControllerPtr_->computeInput(t, x)
         + (1 - alpha_) * learnedControllerPtr_->computeInput(t, x);   // 混合
  }
}
```

$$
u = \alpha\,u^{MPC} + (1-\alpha)\,u^{NN}
$$

**$\alpha$ 的调度**（`mpcnet.py:194`）：

```python
alpha = 1.0 - 1.0 * iteration / self.config.LEARNING_ITERATIONS
```

即 $\alpha$ 从 1 **线性退火**到 0：

| 训练阶段 | $\alpha$ | 行为 |
|---|---|---|
| 开始 | 1.0 | 完全由 MPC 驱动 → 数据在专家分布上 |
| 中期 | 0.5 | 混合 → 逐渐探索策略会去的状态 |
| 结束 | 0.0 | 完全由策略驱动 → 数据在策略分布上 |

**这就是 DAgger 的核心**：训练分布逐渐从专家分布过渡到策略分布，
使策略在自己会遇到的状态上也被正确监督。

**注意标注始终来自 MPC**：无论 $\alpha$ 是多少，
`MpcnetDataGeneration` 记录的哈密顿量近似都来自 MPC 求解结果。

### 14.3.3 数据流

```
Python: mpcnet.train()
   │ start_data_generation(policy, alpha)
   │   └─ 把 policy 导出为 ONNX，传给 C++
   ▼
C++: MpcnetRolloutManager (多线程)
   ├─ MpcnetDataGeneration (N 个线程)
   │    每个线程:
   │      1. 随机初始状态 + 随机目标
   │      2. 用 MpcnetBehavioralController(α) 做 rollout
   │      3. 沿轨迹每隔一定间隔:
   │           - 调 MPC 求解
   │           - 记录 (t, x, u, mode, observation, actionTransformation, hamiltonian)
   │      4. 返回 MpcnetData 列表
   │
   └─ MpcnetPolicyEvaluation (M 个线程)
        用 α=0（纯策略）做 rollout，记录 MpcnetMetrics:
          - survivalTime      策略能撑多久不失败
          - incurredHamiltonian  累积的哈密顿量
   ▼
Python: memory.push(...) → 回放缓冲
   │
   ▼
Python: memory.sample(BATCH_SIZE) → 训练一步
```

### 14.3.4 `MpcnetDefinitionBase`：机器人特定的接口

```cpp
class MpcnetDefinitionBase {
  virtual vector_t getObservation(scalar_t t, const vector_t& x) = 0;
  virtual std::pair<matrix_t, vector_t> getActionTransformation(scalar_t t, const vector_t& x) = 0;
  virtual bool isValid(scalar_t t, const vector_t& x) = 0;
};
```

**`getObservation`**：状态 → 网络输入。
通常做归一化、去掉绝对位置（保留相对量）、加入相位信息等。
**这是特征工程的关键**——好的观测能大幅降低网络负担。

**`getActionTransformation`**：网络输出 $a$ → 真实输入：

$$
u = M(t,x)\,a + v(t,x)
$$

用途举例（四足）：
- $v$ = 重力补偿输入（网络只需学"相对于悬停的修正量"）
- $M$ = 缩放矩阵（把网络的 $[-1,1]$ 输出映射到力的物理范围）

**这个变换让网络的输出分布接近标准正态**，训练稳定得多。

**`isValid`**：判断状态是否"还活着"（如四足没摔倒）。
用于 `survivalTime` 指标与提前终止 rollout。

---

## 14.4 部署：ONNX 推理

**文件**：`control/MpcnetOnnxController.h` / `.cpp`

训练完的策略导出为 ONNX：

```python
torch.onnx.export(model=self.policy, args=self.dummy_observation, f=save_path + ".onnx")
```

C++ 侧用 **ONNX Runtime** 加载并推理：

```cpp
class MpcnetOnnxController : public MpcnetControllerBase {
  void loadPolicyModel(const std::string& policyFilePath) override;
  vector_t computeInput(scalar_t t, const vector_t& x) override;
};
```

`computeInput` 的流程：
1. `definitionPtr_->getObservation(t, x)` → 观测向量
2. ONNX Runtime 前向传播 → 动作 $a$
3. `definitionPtr_->getActionTransformation(t, x)` → $(M, v)$
4. 返回 $u = Ma+v$

**为什么用 ONNX 而不是 LibTorch？**
- ONNX Runtime 更轻量（无需完整 PyTorch 运行时）
- 推理延迟更低、更可预测（无 Python GIL、无动态图开销）
- 可部署到更多平台（含嵌入式）

`MpcnetDummyLoopRos` 提供与 `MRT_ROS_Dummy_Loop` 类似的执行循环，
但用 `MpcnetOnnxController` 代替 MPC 策略。

---

## 14.5 训练主循环

**文件**：`python/ocs2_mpcnet_core/mpcnet.py:175`

```python
def train(self):
    # 保存初始策略
    torch.onnx.export(model=self.policy, args=self.dummy_observation, f=save_path + ".onnx")

    # 等第一批数据
    self.start_data_generation(self.policy)
    self.start_policy_evaluation(self.policy)
    while not self.interface.isDataGenerationDone(): time.sleep(1.0)

    for iteration in range(self.config.LEARNING_ITERATIONS):
        alpha = 1.0 - 1.0 * iteration / self.config.LEARNING_ITERATIONS

        # ① 收数据（异步：C++ 在后台生成，Python 不阻塞）
        if self.interface.isDataGenerationDone():
            data = self.interface.getGeneratedData()
            for d in data:
                self.memory.push(d.t, d.x, d.u,
                                 helper.get_one_hot(d.mode, EXPERT_NUM, EXPERT_FOR_MODE),
                                 d.observation, d.actionTransformation, d.hamiltonian)
            self.start_data_generation(self.policy, alpha)     # 立刻启动下一批

        # ② 收评估指标
        if self.interface.isPolicyEvaluationDone():
            metrics = self.interface.getComputedMetrics()
            survival_time = np.mean([m.survivalTime for m in metrics])
            incurred_hamiltonian = np.mean([m.incurredHamiltonian for m in metrics])
            self.start_policy_evaluation(self.policy)

        # ③ 定期保存中间策略
        if iteration % int(0.1 * LEARNING_ITERATIONS) == 0 and iteration > 0:
            torch.onnx.export(...)

        # ④ 采样 + 优化
        (t, x, u, p, observation, M, v, dHdxx, dHdux, dHduu, dHdx, dHdu, H) = \
            self.memory.sample(self.config.BATCH_SIZE)

        def closure():
            self.optimizer.zero_grad()
            action = self.policy(observation)[0]
            input = helper.bmv(M, action) + v                     # 动作变换
            loss = self.experts_loss(x, x, input, u, p, p, dHdxx, dHdux, dHduu, dHdx, dHdu, H)
            loss.backward()
            if GRADIENT_CLIPPING:
                torch.nn.utils.clip_grad_norm_(self.policy.parameters(), GRADIENT_CLIPPING_VALUE)
            return loss

        self.optimizer.step(closure)
```

**三个工程要点**：

1. **数据生成与训练完全异步**：`isDataGenerationDone()` 是非阻塞查询。
   GPU 训练与 CPU 数据生成并行，硬件利用率高。
2. **`optimizer.step(closure)`**：用闭包形式，兼容 L-BFGS 等需要多次求值的优化器。
3. **梯度裁剪**：哈密顿量损失的 Hessian 可能病态（$\mathcal{H}_{uu}$ 量级差异大），
   裁剪防止梯度爆炸。

### 回放缓冲

**文件**：`memory/circular.py`

环形缓冲，固定容量。`push` 满了就覆盖最老的。
`sample(batch_size)` 均匀随机采样。

**为什么用环形而非无限增长？**
因为随着 $\alpha$ 退火，早期数据（专家分布）与后期数据（策略分布）差别很大，
保留太多早期数据会拖慢适应。

---

## 14.6 评估指标

**文件**：`rollout/MpcnetMetrics.h`

```cpp
struct MpcnetMetrics {
  scalar_t survivalTime;          // 策略维持有效状态的时长
  scalar_t incurredHamiltonian;   // 沿轨迹累积的哈密顿量
};
```

- **`survivalTime`**：最直观的指标。四足摔倒、Ballbot 倒下即终止。
- **`incurredHamiltonian`**：控制质量的度量。
  越小说明策略越接近最优。

两者一起看：`survivalTime` 高但 `incurredHamiltonian` 也高
= 策略"能活但很笨"（可能只学会了保守的静止）。

`MpcnetPolicyEvaluation` 用 $\alpha=0$（纯策略）跑 rollout 采集这两个指标，
每次训练迭代都记录到 TensorBoard。

---

## 14.7 两个具体实例

### `ocs2_ballbot_mpcnet`

- 观测：10 维状态直接用（已经归一化得不错）
- 动作变换：恒等
- 策略：`LinearPolicy` 或 `NonlinearPolicy`（无模态，不需要混合专家）
- 训练时间：分钟级

**最好的入门实例**——问题小，训练快，可以完整跑通看效果。

### `ocs2_legged_robot_mpcnet`

- 观测：去掉绝对 xy 位置与 yaw（只保留相对量，实现平移/旋转不变性），
  加入接触相位
- 动作变换：$v$ = 重力补偿输入，$M$ = 对角缩放
- 策略：`MixtureOfNonlinearExpertsPolicy`，专家数 = 独立模态数
- 训练时间：小时级（需要大量 MPC 求解）

**参考**：Reske, Carius, Ma, Farshidian, Hutter (2021),
*Imitation Learning from MPC for Quadrupedal Multi-Gait Control*, ICRA。

---

## 14.8 MPC-Net 的定位与局限

### 优势

| 维度 | MPC | MPC-Net |
|---|---|---|
| 推理延迟 | 10~100 ms | **10~100 μs** |
| 算力需求 | 多核 CPU | 单核 / MCU |
| 内存 | 几十 MB | 几百 KB |
| 训练成本 | 无 | 数小时 GPU + CPU |

### 局限

1. **不保证约束满足**：网络输出可能违反摩擦锥、力矩限。
   部署时通常需要在网络后面加一层安全滤波（如 CBF-QP）。
2. **泛化范围受训练分布限制**：遇到训练时没见过的地形/扰动会失效。
3. **难以改变目标**：MPC 可以随时改变代价函数（换目标位姿），
   网络策略需要把目标编码进观测，且只在训练过的目标范围内有效。
4. **黑盒**：难以调试与验证。

### 与 RL 的关系

| | RL（PPO/SAC） | MPC-Net |
|---|---|---|
| 监督信号 | 标量奖励 | **哈密顿量的二次近似** |
| 样本效率 | 低（$10^7\sim10^9$ 步） | **高**（$10^5\sim10^6$ 步） |
| 需要专家 | 否 | **是**（MPC） |
| 探索 | 需要（熵正则、噪声） | 由 DAgger 的 $\alpha$ 退火提供 |
| 约束处理 | 奖励塑形（难） | 由 MPC 保证（训练数据里） |

**MPC-Net 本质是"用最优控制的一阶/二阶信息代替 RL 的零阶反馈"**，
所以样本效率高出几个数量级——这与"策略梯度只用标量奖励"
vs "监督学习用完整梯度"的差别是同一回事。

---

## 14.9 参考文献

- **Carius, Farshidian, Hutter (2020)**, *MPC-Net: A First Principles Guided Policy Search*,
  IEEE RA-L 5(2):2897–2904 —— **MPC-Net 原始论文**
- **Reske, Carius, Ma, Farshidian, Hutter (2021)**,
  *Imitation Learning from MPC for Quadrupedal Multi-Gait Control*, ICRA
  —— **混合专家 + 四足**
- **Farshidian, Hoeller, Hutter (2020)**, *Deep Value Model Predictive Control*,
  CoRL —— 用网络学值函数（与 MPC-Net 互补的思路）
- **Ross, Gordon, Bagnell (2011)**, *A Reduction of Imitation Learning and Structured
  Prediction to No-Regret Online Learning*, AISTATS —— **DAgger**
- **Jacobs & Jordan et al. (1991)**, *Adaptive Mixtures of Local Experts*, Neural Computation
  —— 混合专家的原始工作
- **Levine & Koltun (2013)**, *Guided Policy Search*, ICML —— 用轨迹优化引导策略学习

**下一章**：数值与工程工具 → [15 工具](15-numerics-tools.md)
