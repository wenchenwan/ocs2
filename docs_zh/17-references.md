# 17 · 参考文献与延伸阅读

本章汇总全文档引用的所有文献，并给出延伸阅读建议。
`ocs2_doc/docs/refs.bib` 里的官方引用**全部包含在内**（标注 †）。

---

## 17.1 OCS2 团队的论文（按主题）

### 核心算法

| 文献 | 内容 | 对应章节 |
|---|---|---|
| † **Farshidian, Kamgarpour, Pardo, Buchli (2017)**, *Sequential Linear Quadratic Optimal Control for Nonlinear Switched Systems*, IFAC 50(1):1463–1469 | **SLQ 原始论文**；切换系统的连续时间 DDP、事件横截条件、切换时刻灵敏度 | [04](04-ocs2-ddp.md), [10](10-integration-rollout.md), [16](16-legacy-modules.md) |
| † **Farshidian, Neunert, Winkler, Rey, Buchli (2017)**, *An Efficient Optimal Planning and Control Framework for Quadrupedal Locomotion*, ICRA, pp. 93–100 | SLQ 在四足上的首个完整应用 | [13](13-robot-examples.md) |
| † **Farshidian, Jelavic, Satapathy, Giftthaler, Buchli (2017)**, *Real-Time Motion Planning of Legged Robots: A Model Predictive Control Approach*, Humanoids, pp. 577–584 | MPC 架构、实时性 | [11](11-mpc-ros.md) |
| † **Sleiman, Farshidian, Hutter (2021)**, *Constraint Handling in Continuous-Time DDP-Based Model Predictive Control*, ICRA | ⭐ **PHR 与光滑 PHR 罚函数**、零空间投影 + 增广拉格朗日 | [04](04-ocs2-ddp.md) §4.4, [08](08-constraints-penalties.md) |
| † **Grandia, Farshidian, Ranftl, Hutter (2019)**, *Feedback MPC for Torque-Controlled Legged Robots*, IROS, pp. 4730–4737 | 高频反馈跟踪、松弛障碍罚函数 | [04](04-ocs2-ddp.md) §4.12, [08](08-constraints-penalties.md), [11](11-mpc-ros.md) |
| † **Grandia, Farshidian, Dosovitskiy, Ranftl, Hutter (2019)**, *Frequency-Aware Model Predictive Control*, IEEE RA-L 4(2):1517–1524 | ⭐ **Loopshaping 的原始论文** | [09](09-loopshaping.md) |

### 机器人应用

| 文献 | 内容 | 章节 |
|---|---|---|
| † **Giftthaler, Farshidian, Sandy, Stadelmann, Buchli (2017)**, *Efficient Kinematic Planning for Mobile Manipulators with Non-Holonomic Constraints Using Optimal Control*, ICRA, pp. 3411–3417 | 移动机械臂、非完整约束 | [13](13-robot-examples.md) §13.5 |
| † **Minniti, Farshidian, Grandia, Hutter (2019)**, *Whole-Body MPC for a Dynamically Stable Mobile Manipulator*, IEEE RA-L 4(4):3687–3694 | **Ballbot** | [13](13-robot-examples.md) §13.4 |
| † **Sleiman, Farshidian, Minniti, Hutter (2021)**, *A Unified MPC Framework for Whole-Body Dynamic Locomotion and Manipulation*, IEEE RA-L 6(3):4688–4695 | 移动操作的统一框架 | [12](12-robot-models.md), [13](13-robot-examples.md) |
| † **Minniti, Grandia, Fäh, Farshidian, Hutter (2021)**, *Model Predictive Robot-Environment Interaction Control for Mobile Manipulation Tasks*, ICRA | 交互力控制 | — |
| † **Gaertner, Bjelonic, Farshidian, Hutter (2021)**, *Collision-Free MPC for Legged Robots in Static and Dynamic Scenes*, ICRA | 球近似 + 距离场避障 | [12](12-robot-models.md) §12.5–12.6, [13](13-robot-examples.md) §13.7 |
| † **Grandia, Taylor, Ames, Hutter (2020)**, *Multi-Layered Safety for Legged Robots via Control Barrier Functions and Model Predictive Control*, arXiv:2011.00032 | CBF + MPC 的分层安全 | — |
| † **Gawel, Blum, Pankert, Krämer, Bartolomei, Ercan, Farshidian, Chli, Gramazio, Siegwart et al. (2019)**, *A Fully-Integrated Sensing and Control System for High-Accuracy Mobile Robotic Building Construction*, IROS, pp. 2300–2307 | 建筑机器人的系统集成 | — |
| † **Mittal, Hoeller, Farshidian, Hutter, Garg (2021)**, *Articulated Object Interaction in Unknown Scenes with Whole-Body Mobile Manipulation*, arXiv:2103.10534 | 铰接物体交互 | — |

### 学习与 MPC 的结合

| 文献 | 内容 | 章节 |
|---|---|---|
| † **Carius, Farshidian, Hutter (2020)**, *MPC-Net: A First Principles Guided Policy Search*, IEEE RA-L 5(2):2897–2904 | ⭐ **MPC-Net 原始论文**，哈密顿量损失 | [14](14-mpcnet.md) |
| † **Reske, Carius, Ma, Farshidian, Hutter (2021)**, *Imitation Learning from MPC for Quadrupedal Multi-Gait Control*, ICRA | 混合专家 + 四足多步态 | [14](14-mpcnet.md) |
| † **Farshidian, Hoeller, Hutter (2020)**, *Deep Value Model Predictive Control*, CoRL, PMLR 100:990–1004 | 学习值函数作为终端代价 | [03](03-ocs2-oc.md) §3.7 |

---

## 17.2 优化算法

### DDP 系列

- **Jacobson, D.H. & Mayne, D.Q. (1970)**, *Differential Dynamic Programming*,
  American Elsevier —— **DDP 原典**
- **Todorov, E. & Li, W. (2005)**, *A Generalized Iterative LQG Method for
  Locally-Optimal Feedback Control of Constrained Nonlinear Stochastic Systems*, ACC
  —— **iLQG/iLQR**，高斯-牛顿近似的来源
- **Tassa, Y., Erez, T., Todorov, E. (2012)**, *Synthesis and Stabilization of Complex
  Behaviors through Online Trajectory Optimization*, IROS
  —— DDP 的正则化、线搜索、改进的 Hessian 处理
- **Tassa, Y., Mansard, N., Todorov, E. (2014)**, *Control-Limited Differential
  Dynamic Programming*, ICRA —— 带箱约束的 DDP
- **Mastalli, C. et al. (2020)**, *Crocoddyl: An Efficient and Versatile Framework for
  Multi-Contact Optimal Control*, ICRA —— 同类库，可对照

### SQP / 内点法 / 多重打靶

- **Bock, H.G. & Plitt, K.J. (1984)**, *A Multiple Shooting Algorithm for Direct
  Solution of Optimal Control Problems*, IFAC World Congress
  —— **多重打靶原典**
- **Diehl, M., Bock, H.G., Schlöder, J.P. (2005)**, *A Real-Time Iteration Scheme for
  Nonlinear Optimization in Optimal Feedback Control*, SIAM J. Control Optim. 43(5)
  —— **RTI（实时迭代）**
- **Wächter, A. & Biegler, L.T. (2006)**, *On the Implementation of an Interior-Point
  Filter Line-Search Algorithm for Large-Scale Nonlinear Programming*,
  Mathematical Programming 106(1):25–57
  —— ⭐ **IPOPT**：滤子线搜索与障碍参数调度，`ocs2_ipm` 与 `FilterLinesearch` 的直接蓝本
  （代码 `SqpSolver.cpp:476` 给出了这个链接）
- **IPOPT 在线文档**：https://coin-or.github.io/Ipopt/OPTIONS.html
  —— `IpmInitialization.h` 的注释直接引用
- **Frison, G. & Diehl, M. (2020)**, *HPIPM: a high-performance quadratic programming
  framework for model predictive control*, IFAC World Congress
- **Frison, G. et al. (2018)**, *BLASFEO: Basic Linear Algebra Subroutines For
  Embedded Optimization*, ACM Trans. Math. Software 44(4)
- **Verschueren, R. et al. (2022)**, *acados — a modular open-source framework for
  fast embedded optimal control*, Math. Prog. Comp. —— 同类框架
- **Katayama, S. & Ohtsuka, T. (2022)**, *Efficient solution method based on inverse
  dynamics for interior point methods of nonlinear MPC*, ICRA
  —— `ocs2_ipm` 的贡献者 Sotaro Katayama 的相关工作（见 `ocs2_ipm/package.xml`）

### 一阶方法

- **Yu, Y., Elango, P., Açıkmeşe, B. (2020)**, *Proportional-Integral Projected
  Gradient Method for Model Predictive Control*, arXiv:2009.06980
  —— ⭐ **PIPG 原始论文**（`PipgSolver.h:50`、`PipgBounds.h:40` 引用）
- **Yu, Y., Elango, P., Açıkmeşe, B., Topcu, U. (2022)**, *Extrapolated
  Proportional-Integral Projected Gradient Method for Conic Optimization*, IEEE L-CSS
- **Ruiz, D. (2001)**, *A Scaling Algorithm to Equilibrate Both Rows and Columns Norms
  in Matrices*, Technical Report RAL-TR-2001-034, Rutherford Appleton Laboratory
  —— ⭐ **Ruiz 均衡**（`Ruzi.h` 注释引用）
- **Stellato, B. et al. (2020)**, *OSQP: An Operator Splitting Solver for Quadratic
  Programs*, Math. Prog. Comp. 12(4) —— 同为一阶 + Ruiz 预条件
- **Chambolle, A. & Pock, T. (2011)**, *A First-Order Primal-Dual Algorithm for Convex
  Problems with Applications to Imaging*, J. Math. Imaging Vis. 40(1) —— PDHG

### 罚函数与增广拉格朗日

- **Hestenes, M.R. (1969)**, *Multiplier and Gradient Methods*, JOTA 4(5):303–320
- **Powell, M.J.D. (1969)**, *A Method for Nonlinear Constraints in Minimization
  Problems*, in *Optimization* (Fletcher, ed.)
- **Rockafellar, R.T. (1974)**, *Augmented Lagrange Multiplier Functions and Duality
  in Nonconvex Programming*, SIAM J. Control 12(2)
  —— ⭐ **PHR 罚函数**（Powell–Hestenes–Rockafellar）
- **Bertsekas, D.P. (1982)**, *Constrained Optimization and Lagrange Multiplier
  Methods*, Academic Press —— 增广拉格朗日的标准教材
- **Conn, A.R., Gould, N.I.M., Toint, P.L. (1991)**, *A Globally Convergent Augmented
  Lagrangian Algorithm for Optimization with General Constraints and Simple Bounds*,
  SIAM J. Numer. Anal. 28(2) —— **罚系数/容差更新策略**（[04](04-ocs2-ddp.md) §4.7.2）
- **Conn, A.R., Gould, N.I.M., Toint, P.L. (1992)**, *LANCELOT: A Fortran Package for
  Large-Scale Nonlinear Optimization*, Springer
- **Feller, C. & Ebenbauer, C. (2017)**, *Relaxed Logarithmic Barrier Function Based
  Model Predictive Control of Linear Systems*, IEEE TAC 62(3)
  —— ⭐ **松弛障碍函数**（[08](08-constraints-penalties.md) §8.2.2）
- **Howell, T., Jackson, B., Manchester, Z. (2019)**, *ALTRO: A Fast Solver for
  Constrained Trajectory Optimization*, IROS —— AL + iLQR，可对照

### 数值优化通用

- **Nocedal, J. & Wright, S.J. (2006)**, *Numerical Optimization*, 2nd ed., Springer
  —— **最重要的通用参考**
  - §3.1 Armijo 条件与回溯线搜索 → [04](04-ocs2-ddp.md) §4.7.3
  - §4.1 信赖域与 $\rho$ 比率 → [04](04-ocs2-ddp.md) §4.7.4
  - §17.2 精确罚函数 → [04](04-ocs2-ddp.md) §4.7.1
  - Ch. 17 罚函数与增广拉格朗日 → [08](08-constraints-penalties.md)
  - Ch. 18 SQP → [05](05-ocs2-sqp.md)
  - Ch. 19 非线性内点法 → [06](06-ocs2-ipm.md)
- **Boyd, S. & Vandenberghe, L. (2004)**, *Convex Optimization*, Cambridge Univ. Press
  —— Ch. 11 障碍法与中心路径
- **Gill, P.E., Murray, W., Wright, M.H. (1981)**, *Practical Optimization*,
  Academic Press —— **修正 Cholesky**（[04](04-ocs2-ddp.md) §4.8）
- **Moré, J.J. (1978)**, *The Levenberg-Marquardt Algorithm: Implementation and Theory*,
  in *Numerical Analysis*, Springer —— LM 的 $\rho$ 阈值 0.25/0.75
- **Frank, M. & Wolfe, P. (1956)**, *An Algorithm for Quadratic Programming*,
  Naval Research Logistics Quarterly 3(1-2):95–110
- **Jaggi, M. (2013)**, *Revisiting Frank-Wolfe: Projection-Free Sparse Convex
  Optimization*, ICML —— `FrankWolfeDescentDirection.h` 的 `\cite jaggi13`

### 最优控制理论

- **Bellman, R. (1957)**, *Dynamic Programming*, Princeton Univ. Press
- **Pontryagin, L.S. et al. (1962)**, *The Mathematical Theory of Optimal Processes*
- **Bryson, A.E. & Ho, Y.-C. (1975)**, *Applied Optimal Control*, Taylor & Francis
  —— 经典教材，Riccati 与 TPBVP 的详细推导
- **Bertsekas, D.P. (2017)**, *Dynamic Programming and Optimal Control*, 4th ed.
- **Rawlings, J.B., Mayne, D.Q., Diehl, M.M. (2017)**, *Model Predictive Control:
  Theory, Computation, and Design*, 2nd ed., Nob Hill —— **MPC 的标准教材**
- **Betts, J.T. (2010)**, *Practical Methods for Optimal Control and Estimation Using
  Nonlinear Programming*, 2nd ed., SIAM

### 切换与混合系统

- **Xu, X. & Antsaklis, P.J. (2004)**, *Optimal Control of Switched Systems Based on
  Parameterization of the Switching Instants*, IEEE TAC 49(1)
  —— **时间尺度变换**（[16](16-legacy-modules.md) §16.1.3）
- **Egerstedt, M., Wardi, Y., Axelsson, H. (2006)**, *Transition-Time Optimization for
  Switched-Mode Dynamical Systems*, IEEE TAC 51(1)
- **Johnson, E.R. & Murphey, T.D. (2011)**, *Second-Order Switching Time Optimization
  for Nonlinear Time-Varying Dynamic Systems*, IEEE TAC 56(8)
- **Van der Schaft, A. & Schumacher, H. (2000)**, *An Introduction to Hybrid Dynamical
  Systems*, Springer —— 混合系统、Zeno 现象

---

## 17.3 机器人动力学与建模

- **Featherstone, R. (2008)**, *Rigid Body Dynamics Algorithms*, Springer
  —— **刚体动力学的标准教材**（ABA、RNEA、CRBA）
- **Orin, D.E. & Goswami, A. (2008)**, *Centroidal Momentum Matrix of a Humanoid Robot:
  Structure and Properties*, IROS —— ⭐ **CMM 的定义**
- **Orin, D.E., Goswami, A., Lee, S.-H. (2013)**, *Centroidal Dynamics of a Humanoid
  Robot*, Autonomous Robots 35(2-3):161–176 —— **质心动力学的系统性阐述**
- **Wensing, P.M. & Orin, D.E. (2016)**, *Improved Computation of the Humanoid
  Centroidal Dynamics and Application for Whole-Body Control*,
  Int. J. Humanoid Robotics 13(1) —— $A_b$ 的结构与高效求逆
- **Dai, H., Valenzuela, A., Tedrake, R. (2014)**, *Whole-body Motion Planning with
  Centroidal Dynamics and Full Kinematics*, Humanoids
- **Carpentier, J., Saurel, G., Buondonno, G., Mirabel, J., Lamiraux, F., Stasse, O.,
  Mansard, N. (2019)**, *The Pinocchio C++ Library*, SII
- **Frigerio, M., Buchli, J., Caldwell, D.G., Semini, C. (2016)**, *RobCoGen: a code
  generator for efficient kinematics and dynamics of articulated robots, based on
  Domain Specific Languages*, J. Software Engineering for Robotics 7(1)
  —— **Ballbot 用的代码生成器**
- **Pan, J., Chitta, S., Manocha, D. (2012)**, *FCL: A General Purpose Library for
  Collision and Proximity Queries*, ICRA
- **Pardo, D., Möller, L., Neunert, M., Winkler, A.W., Buchli, J. (2016)**,
  *Evaluating Direct Transcription and Nonlinear Optimization Methods for Robot
  Motion Planning*, IEEE RA-L 1(2)
- **Winkler, A.W., Bellicoso, C.D., Hutter, M., Buchli, J. (2018)**, *Gait and
  Trajectory Optimization for Legged Systems through Phase-based End-Effector
  Parameterization*, IEEE RA-L 3(3) —— **TOWR**，相位参数化的步态优化
- **Hwangbo, J., Lee, J., Hutter, M. (2018)**, *Per-Contact Iteration Method for
  Solving Contact Dynamics*, IEEE RA-L 3(2) —— **RaiSim** 的接触求解
- **Felzenszwalb, P. & Huttenlocher, D. (2012)**, *Distance Transforms of Sampled
  Functions*, Theory of Computing 8(1) —— EDT 的 $O(n)$ 算法

---

## 17.4 数值方法与软件

- **Griewank, A. & Walther, A. (2008)**, *Evaluating Derivatives: Principles and
  Techniques of Algorithmic Differentiation*, 2nd ed., SIAM
  —— **自动微分的标准教材**
- **Bell, B.M.**, *CppAD: A Package for Differentiation of C++ Algorithms*,
  https://coin-or.github.io/CppAD/
- **Leal, J.R.**, *CppADCodeGen*, https://github.com/joaoleal/CppADCodeGen
- **Guennebaud, G., Jacob, B. et al. (2010)**, *Eigen v3*, https://eigen.tuxfamily.org
- **Ahnert, K. & Mulansky, M. (2011)**, *Odeint – Solving Ordinary Differential
  Equations in C++*, AIP Conf. Proc. 1389
- **Hairer, E., Nørsett, S.P., Wanner, G. (1993)**, *Solving Ordinary Differential
  Equations I: Nonstiff Problems*, 2nd ed., Springer —— Dormand-Prince、自适应步长
- **Higham, N.J. (2002)**, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM
- **Golub, G.H. & Van Loan, C.F. (2013)**, *Matrix Computations*, 4th ed.
  —— QR、Cholesky、特征值算法

### 求根算法（`RootFinder`）

以下四篇在 `ocs2_oc/include/ocs2_oc/rollout/RootFinder.h:51-63` 的注释中被直接引用：

- **Anderson, N. & Björck, Å. (1973)**, *A New High Order Method of Regula Falsi Type
  for Computing a Root of an Equation*, BIT 13(3):253–264
- **Dowell, M. & Jarratt, P. (1971)**, *A Modified Regula Falsi Method for Computing
  the Root of an Equation*, BIT 11(2):168–174 —— **Illinois 方法**
- **Dowell, M. & Jarratt, P. (1972)**, *The "Pegasus" Method for Computing the Root
  of an Equation*, BIT 12(4):503–508
- **Ford, J.A. (1995)**, *Improved Algorithms of Illinois-Type for the Numerical
  Solution of Nonlinear Equations*, Tech. Report, University of Essex

---

## 17.5 控制理论与频域方法

- **Skogestad, S. & Postlethwaite, I. (2005)**, *Multivariable Feedback Control:
  Analysis and Design*, 2nd ed., Wiley —— **loopshaping 的经典教材**
- **Zhou, K., Doyle, J.C., Glover, K. (1996)**, *Robust and Optimal Control*,
  Prentice Hall —— 加权函数设计与鲁棒性
- **Åström, K.J. & Murray, R.M. (2008)**, *Feedback Systems: An Introduction for
  Scientists and Engineers*, Princeton Univ. Press

---

## 17.6 学习方法

- **Ross, S., Gordon, G., Bagnell, D. (2011)**, *A Reduction of Imitation Learning and
  Structured Prediction to No-Regret Online Learning*, AISTATS —— ⭐ **DAgger**
- **Levine, S. & Koltun, V. (2013)**, *Guided Policy Search*, ICML
- **Jacobs, R.A., Jordan, M.I., Nowlan, S.J., Hinton, G.E. (1991)**, *Adaptive Mixtures
  of Local Experts*, Neural Computation 3(1) —— **混合专家**
- **Shazeer, N. et al. (2017)**, *Outrageously Large Neural Networks: The
  Sparsely-Gated Mixture-of-Experts Layer*, ICLR —— 现代 MoE

---

## 17.7 相关开源项目

| 项目 | 链接 | 与 OCS2 的关系 |
|---|---|---|
| **Control Toolbox (CT)** | github.com/ethz-adrl/control-toolbox | 同为 ETH，SLQ 同源 |
| **Crocoddyl** | github.com/loco-3d/crocoddyl | DDP 系，与 Pinocchio 深度耦合 |
| **acados** | github.com/acados/acados | 多重打靶 + HPIPM，思路与 `ocs2_sqp` 一致 |
| **ALTRO / TrajectoryOptimization.jl** | github.com/RoboticExplorationLab | AL + iLQR |
| **TOWR** | github.com/ethz-adrl/towr | 相位参数化的腿足轨迹优化 |
| **HPIPM** | github.com/giaf/hpipm | OCS2 的 QP 求解器 |
| **BLASFEO** | github.com/giaf/blasfeo | HPIPM 的线性代数后端 |
| **Pinocchio** | github.com/stack-of-tasks/pinocchio | 刚体动力学 |
| **HPP-FCL** | github.com/humanoid-path-planner/hpp-fcl | 碰撞检测 |
| **CppAD / CppADCodeGen** | github.com/coin-or/CppAD | 自动微分 |
| **RaiSim** | raisim.com | 物理引擎（需授权） |
| **OSQP** | osqp.org | 一阶 QP，可与 PIPG 对照 |
| **Ipopt** | github.com/coin-or/Ipopt | `ocs2_ipm` 的设计蓝本 |

---

## 17.8 学习路线建议

### 如果你是最优控制新手

1. **Bryson & Ho (1975)** Ch. 1–6：变分法、Riccati、TPBVP
2. **Nocedal & Wright (2006)** Ch. 2–4：无约束优化、线搜索、信赖域
3. 本文档 [04 章](04-ocs2-ddp.md)：把上面两者串起来
4. **Rawlings, Mayne, Diehl (2017)** Ch. 2：MPC 的稳定性

### 如果你要做腿足机器人

1. **Featherstone (2008)** Ch. 1–5：刚体动力学基础
2. **Orin, Goswami, Lee (2013)**：质心动力学
3. 本文档 [12](12-robot-models.md) + [13](13-robot-examples.md)
4. **Farshidian et al. (2017)** ICRA：OCS2 的四足论文
5. **Grandia et al. (2019)** IROS：反馈 MPC

### 如果你要改进求解器

1. **Nocedal & Wright (2006)** Ch. 17–19
2. **Wächter & Biegler (2006)**：滤子线搜索的完整理论
3. 本文档 [05](05-ocs2-sqp.md) + [06](06-ocs2-ipm.md) + [07](07-ocs2-slp.md)
4. **Frison & Diehl (2020)**：HPIPM 的内部结构
5. 读 acados 与 Crocoddyl 的源码做对照

### 如果你要做学习与控制的结合

1. **Carius, Farshidian, Hutter (2020)**：MPC-Net
2. **Ross et al. (2011)**：DAgger
3. 本文档 [14 章](14-mpcnet.md)
4. **Levine & Koltun (2013)**：Guided Policy Search

---

## 17.9 引用 OCS2

仓库里有两份内容相同的 BibTeX 文件（共 17 条，全部是 OCS2 团队的论文）：

- `ocs2_doc/docs/refs.bib`（Sphinx 文档引用）
- `ocs2_doc/tools/sphinx/_static/cite.bib`（供用户下载引用）

**它们不含"引用本工具箱"的通用条目**，
所以引用 OCS2 时的惯例是**引用你实际用到的算法所对应的论文**：

| 你用了 | 引用 |
|---|---|
| SLQ / iLQR | `farshidian2017ocs2` + `farshidian2017slq` |
| 约束处理（AL / 投影） | `sleiman2021constraint` |
| 反馈 MPC | `grandia2019feedback` |
| Loopshaping | `grandia2019frequency` |
| MPC-Net | `carius2020mpcnet`（+ `reske2021imitation` 若用混合专家） |
| 四足运动 | `farshidian2017slq` + `farshidian2017slqmpc` |
| 移动操作 | `sleiman2021locopulation` |
| Ballbot | `minniti2019ballbot` |
| 避障 | `gaertner2021collision` |

外加工具箱本身的链接：https://github.com/leggedrobotics/ocs2

§17.1 的表格给出了这 17 条的完整信息与对应章节。

**下一章**：代码索引 → [18 代码索引](18-code-index.md)
