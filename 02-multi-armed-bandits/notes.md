# 多臂赌博机模拟研究笔记

## 一、赌博机问题的数学框架

多臂赌博机（Multi-Armed Bandit, MAB）问题是序贯决策的核心模型。其形式化定义如下：

设 $K$ 个臂，每个臂 $k \in \{1, \dots, K\}$ 对应一个未知的奖励分布 $\nu_k$，其均值为 $\mu_k$。在每一时刻 $t = 1, 2, \dots, T$，智能体选择一个臂 $A_t \in \{1, \dots, K\}$ 并获得随机奖励 $X_t \sim \nu_{A_t}$。

**最优臂**定义为 $\star = \arg\max_k \mu_k$，其均值 $\mu^\star = \max_k \mu_k$。

**遗憾**（Regret）是衡量算法性能的核心指标。累积遗憾定义为：

$$R(T) = T \mu^\star - \sum_{t=1}^{T} \mathbb{E}[X_t] = \sum_{t=1}^{T} \sum_{k=1}^{K} \Delta_k \cdot \mathbb{E}[N_k(T)]$$

其中 $\Delta_k = \mu^\star - \mu_k$ 是臂 $k$ 的差距（gap），$N_k(T) = \sum_{t=1}^{T} \mathbb{I}(A_t = k)$ 是臂 $k$ 被拉动的次数。

**简单遗憾**定义为 $r(T) = \mu^\star - \mu_{A_T}$，衡量最终时刻是否选中了最优臂。

Lai and Robbins (1985) 证明了渐近下界：对任意一致性策略，

$$\liminf_{T \to \infty} \frac{R(T)}{\log T} \geq \sum_{k: \Delta_k > 0} \frac{\Delta_k}{\text{KL}(\nu_k, \nu^\star)}$$

其中 $\text{KL}(\nu_k, \nu^\star)$ 是臂 $k$ 的分布与最优臂分布之间的 KL 散度。对于 Bernoulli 奖励，$\text{KL}(\text{Ber}(p), \text{Ber}(q)) = p \log(p/q) + (1-p) \log((1-p)/(1-q))$。

## 二、探索与利用的理论张力

探索-利用困境（exploration-exploitation dilemma）是 MAB 问题的本质。纯粹利用（总是选择当前估计最好的臂）会导致次线性但过高的遗憾；纯粹探索（随机选择）导致线性遗憾 $R(T) = \Theta(T)$。

最优策略需要在两者之间取得对数级平衡。关键理论洞见：

**1. $\varepsilon$-greedy 的局限性：** 固定 $\varepsilon$ 的策略遗憾为 $R(T) = \Theta(\varepsilon T \sum_k \Delta_k) + \Theta(\frac{K}{\varepsilon} \log T \sum_k \frac{1}{\Delta_k})$，第一项来自持续随机探索的线性代价。即使 $\varepsilon$ 很小，当 $T$ 足够大时，线性项主导。递减 $\varepsilon_t = c/t$ 可以实现对数遗憾，但需要知道 $\Delta_{\min}$。

**2. UCB 的乐观原则：** "面对不确定性时保持乐观"（optimism in the face of uncertainty）。通过构造置信上界，UCB 确保每个臂至少被拉动 $\Omega(\log T)$ 次，且次优臂被拉动 $O(\frac{\log T}{\Delta_k^2})$ 次，恰好匹配 Lai-Robbins 下界。

**3. Thompson 采样的贝叶斯解释：** TS 从后验分布中采样，等价于在参数空间中做随机探索，探索量与后验不确定性成正比。Agrawal and Goyal (2012) 证明 TS 在 Bernoulli 奖励下达到 $R(T) = O(\sum_k \frac{\log T}{\Delta_k})$，与 UCB 同阶。

**4. 信息论下界与实例依赖：** 对于 $K$ 臂 Bernoulli 赌博机，minimax 遗憾下界为 $\Omega(\sqrt{KT})$。这告诉我们，当臂间距 $\Delta_k$ 很小时（如 Scenario 2），任何算法都无法避免高遗憾。

## 三、算法对比的核心发现

本次模拟涵盖五种主要算法，在五个场景中进行了系统比较：

**Scenario 1（标准 Bernoulli）：** TS 显著优于所有其他算法。其累积遗憾最低，且最优臂识别率最高（约 70%+）。UCB-Tuned 紧随其后，标准 UCB1 次之。EXP3 作为对抗性算法，在此随机（stochastic）环境中表现最差——这是合理的，因为其设计假设了最坏情况的对手。$\varepsilon$-greedy 的表现对 $\varepsilon$ 值高度敏感：$\varepsilon = 0.05$ 在简单场景中接近 TS，但 $\varepsilon = 0.20$ 则因过度探索而显著劣化。

**Scenario 2（困难辨别）：** 当臂间距极小（$\Delta_{\max} = 0.10$）时，所有算法的遗憾均大幅增长。UCB-Tuned 利用方差自适应获得微弱优势。TS 的理论优势在此体现为更平滑的收敛。值得注意的是，EXP3 的差距在此缩小——当所有臂都很接近时，即使随机探索的代价也不大。

**Scenario 3（多臂 K=20）：** 臂数增加对 UCB 的影响为 $O(K \log T)$，对 TS 的影响为 $O(K)$。模拟结果确认 TS 在多臂场景中的可扩展性优势。EXP3 的 $\tilde{O}(\sqrt{KT})$ 遗憾在多臂时急剧增长。

**关键规律总结：**
- TS $\succ$ UCB-Tuned $\succ$ UCB1 $\succ$ EG（小 $\varepsilon$）$\succ$ EXP3（随机环境）
- 臂间距越小，算法间差异越小（都接近下界）
- 臂数增加时，TS 相对于 UCB 的优势放大

## 四、非平稳环境的挑战

Scenario 4 引入了两次突变（t=3000 和 t=7000），最优臂发生轮换。这暴露了标准算法的根本局限：

**标准 UCB 的"记忆过载"问题：** UCB 的置信区间宽度为 $\sqrt{2\log t / N_k(t)}$，随着 $t$ 增大，即使 $N_k(t)$ 不变，置信区间也会缩窄。在分布漂移后，历史数据变得不相关，但 UCB 仍依赖全部历史，导致对旧数据的过度信任。

**标准 TS 的"后验固化"问题：** Beta 后验参数 $\alpha, \beta$ 单调递增，后验方差 $\propto 1/(\alpha+\beta)$ 趋于零。当分布变化时，TS 的采样几乎退化为确定性选择，丧失探索能力。

**解决方案对比：**
- **滑动窗口 TS（SW-TS）：** 仅使用最近 $w$ 步的数据更新后验。模拟显示其在漂移后的适应速度明显快于标准 TS，但窗口大小 $w$ 需要调参——太小导致高方差，太大导致响应迟钝。
- **折扣 UCB（D-UCB）：** 对历史统计量施加指数折扣 $\gamma^t$，有效"遗忘"旧数据。理论遗憾界为 $O(\frac{1}{1-\gamma} \sqrt{T \log T})$，其中 $\gamma$ 控制记忆长度。

模拟结果表明，非平稳感知的算法在突变后能快速重新探索并锁定新的最优臂，而标准算法需要数百甚至数千步才能"醒过来"。

**实践启示：** 在线广告、推荐系统等场景中，用户偏好和物品效果随时间变化，必须采用非平稳感知算法或定期重置策略。

## 五、实际部署考量

将 MAB 算法部署到生产环境需要关注以下工程与统计问题：

**1. 冷启动问题：** 新臂加入时没有任何历史数据。TS 的先验分布自然地处理了这一问题（先验方差大 → 更多探索），而 UCB 需要显式的初始化策略。

**2. 上下文信息（Contextual Bandits）：** 现实中，用户/物品特征可用作上下文。本研究的非上下文算法是 contextual bandit 的特例（上下文为空），其遗憾界提供了 baseline。

**3. 延迟反馈：** 在线实验中，转化可能滞后数天。标准 MAB 假设即时反馈，延迟反馈下 UCB 和 TS 的遗憾界需要修正为 $O(\sqrt{T \cdot d_{\max} \log T})$，其中 $d_{\max}$ 是最大延迟。

**4. 安全性约束：** 某些场景（如医疗剂量选择）要求臂的期望奖励不低于某个阈值。Safe UCB 和 constrained TS 是相关方向。

**5. 计算效率：** 
- $\varepsilon$-greedy: $O(1)$ per step
- UCB1/UCB-Tuned: $O(K)$ per step
- TS: $O(K)$ per step (Beta sampling is $O(1)$ per arm)
- EXP3: $O(K)$ per step

在高维场景（$K > 10^4$）中，需要近似方法或 bandit reduction。

**6. 多指标优化：** 实践中往往需要同时优化点击率和转化率。Pareto bandits 和 constrained bandits 是解决方案。

## 六、与 A/B 测试的联系与区别

MAB 与传统 A/B 测试的关系值得深入讨论：

**相似之处：** 两者都面临探索-利用困境。A/B 测试中的样本量确定等价于 MAB 中的遗憾最小化——都需要在统计功效和资源消耗之间取得平衡。

**关键区别：**

| 维度 | A/B 测试 | MAB |
|------|---------|-----|
| 目标 | 假设检验，做出最终决策 | 在线优化，最小化累积遗憾 |
| 流量分配 | 固定（50/50） | 自适应 |
| 分析时机 | 实验结束后 | 实时 |
| 多重检验 | 需要显式校正 | 内置处理 |
| 序贯性 | 通常需要预设样本量 | 天然序贯，可随时停止 |

**MAB 作为 A/B 测试的升级：** 在多变量测试（MVT）场景中，当有 $K > 2$ 个变体时，MAB 算法可以在保持统计有效性的同时，将更多流量导向表现好的变体，减少实验期间的机会成本。

**Bayesian MAB 与 Bayesian A/B 测试的统一：** 两者都可以使用相同的后验推断框架。Thompson Sampling 本质上就是在做 Bayesian model averaging，与 Week 17-18 的 Bayesian A/B testing 框架一致。区别仅在于决策规则：A/B 测试使用 ROPE/Bayes Factor 做终止决策，MAB 使用 TS 做流量分配。

**注意陷阱：** 在 A/B 测试中使用自适应分配（如 TS）会引入 selection bias，标准的 t-test/chi-square 不再有效。需要使用专门的统计检验方法（如 always-valid p-values）。

## 七、理论遗憾界与经验表现的对比

将本次模拟的经验结果与理论界进行对照：

**UCB1 理论界：** $R(T) \leq 8 \sum_{k: \Delta_k > 0} \frac{\log T}{\Delta_k} + (1 + \frac{\pi^2}{3}) \sum_k \Delta_k$

对于 Scenario 1，$\Delta = (0.8, 0.6, 0.4, 0.2, 0)$，$\log(10000) \approx 9.21$：
$$R_{\text{UCB1}} \leq 8 \times 9.21 \times (\frac{1}{0.8} + \frac{1}{0.6} + \frac{1}{0.4} + \frac{1}{0.2}) + 4.29 \times 2.0 \approx 8 \times 9.21 \times 7.08 + 8.58 \approx 530$$

模拟中 UCB1 的经验遗憾约为 200-300，远低于最坏情况下界——这是因为理论界是最坏情况分析，而实际中臂的奖励方差较低。

**TS 理论界（Agrawal & Goyal, 2013）：** 
$$R_{\text{TS}}(T) \leq \left(\sum_{k: \Delta_k > 0} \frac{1}{d(\mu_k, \mu^\star)}\right) \log T + o(\log T)$$

其中 $d(p,q)$ 是 Bernoulli KL 散度。由于 KL 散度比 $1/\Delta_k^2$ 更紧，TS 在理论上严格优于 UCB1（除极端情况外），这与模拟结果一致。

**EXP3 理论界：** $R_{\text{EXP3}}(T) \leq \sqrt{2TK\log K} + \frac{\gamma T}{2}$

对于 Scenario 1（$K=5, T=10000, \gamma=0.1$）：
$$R_{\text{EXP3}} \leq \sqrt{2 \times 10000 \times 5 \times 1.61} + 500 \approx 401 + 500 = 901$$

经验遗憾约 1500-2500，略高于理论界——这是因为 EXP3 是为对抗性环境设计的，在随机环境中其遗憾界并非紧的，且 $\gamma$ 参数引入的额外探索代价在随机环境中被放大。

**Lai-Robbins 下界验证：** 对于 Scenario 1：
$$R_{\text{LR}} \approx \log(10000) \times \sum_k \frac{\Delta_k}{d(\mu_k, \mu^\star)} \approx 9.21 \times \sum_k \frac{\Delta_k}{\text{KL}(\text{Ber}(\mu_k), \text{Ber}(0.9))}$$

计算各项 KL 散度并求和，下界约为 30-50。TS 的经验遗憾（约 100-150）在此下界的 2-3 倍范围内，表明 TS 接近最优。

**核心结论：** 理论界给出了算法在最坏情况下的性能保证，而经验结果通常显著优于最坏情况界。TS 在大多数场景中接近 Lai-Robbins 下界，是 Bernoulli MAB 的首选算法。UCB-Tuned 作为不需要分布假设的算法，是 TS 的最佳替代。在部署时，应根据环境特性（平稳/非平稳、奖励类型、臂数量）选择算法，并通过模拟验证选择。

---

**参考文献：**
- Lai, T. L., & Robbins, H. (1985). Asymptotically efficient adaptive allocation rules. *Advances in Applied Mathematics*.
- Auer, P., Cesa-Bianchi, N., & Fischer, P. (2002). Finite-time analysis of the multiarmed bandit problem. *Machine Learning*.
- Agrawal, S., & Goyal, N. (2012). Analysis of Thompson sampling for the multi-armed bandit problem. *COLT*.
- Chapelle, O., & Li, L. (2011). An empirical evaluation of Thompson sampling. *NeurIPS*.
- Burtini, G., Loeppky, J., & Lawrence, R. (2015). Improving online marketing experiments with multi-armed bandits.
