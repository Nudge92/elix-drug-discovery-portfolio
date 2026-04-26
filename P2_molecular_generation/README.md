# P2. MCTS-Based Molecular Generation with ChemTSv2

**Tool:** [ChemTSv2](https://github.com/molecule-generator-collection/ChemTSv2) (Elix)  
**Goal:** Design EGFR-selective kinase inhibitors through systematic ablation of reward functions and MCTS policies  
**Key Finding:** Step(1.0→0.1) adaptive scheduling achieves selectivity +0.64 — over 10× the original paper's default

---

## Motivation

ChemTSv2 uses Monte Carlo Tree Search (MCTS) to generate molecules by sequentially selecting SMILES tokens. Two fundamental axes control generation quality: (1) the **reward function** that evaluates generated molecules, and (2) the **MCTS policy** that balances exploration vs. exploitation during search.

The original paper uses fixed defaults (geometric mean reward, UCB1 with c=1.0). This project systematically decomposes both axes to find the optimal combination for a challenging multi-objective task: generating molecules that are selective for EGFR over 8 off-target kinases while maintaining drug-likeness.

---

## Experimental Design

```
ChemTSv2 Ablation Study (8 experiment groups, 100+ runs)
│
├── Reward Function Axis
│   ├── 1.1  Aggregation method (geometric / arithmetic / Tchebycheff)
│   ├── 1.2  Gaussian σ sweep (EGFR / off-target strictness)
│   ├── 1.3  EGFR weight sweep (w = 1–20)
│   ├── Cross experiments (c=0.1 × σ/weight combinations)
│   └── Long-run validation (gen = 30,000)
│
└── MCTS Policy Axis
    ├── 2.1  c_val sweep (0.01 – 2.0)
    ├── 2.2  PUCT (RNN prior integration) ← custom implementation
    └── 2.3  Adaptive c_val (linear / step decay) ← custom implementation
```

**Common conditions:** 18 objectives (9 kinase activities + 4 ADMET + 5 drug-likeness), gen=5,000, 3 seeds (0, 42, 123), ZINC 250K pre-trained RNN.

---

## Results: Top 5 Configurations

| Rank | Configuration | Selectivity | EGFR pIC₅₀ | QED | Unique% |
|---|---|---|---|---|---|
| 1 | **Step(1.0→0.1) adaptive** | **+0.64** | **6.56** | **0.778** | 93.7% |
| 2 | Linear(1.0→0.01) adaptive | +0.56 | 6.51 | 0.749 | 97.4% |
| 3 | UCB1 c=0.1 (fixed) | +0.41 | 6.35 | 0.759 | 91.8% |
| 4 | PUCT cp=0.5 | +0.24 | 6.29 | 0.747 | 95.9% |
| 5 | UCB1 c=1.0 (paper default) | +0.06 | 6.04 | 0.775 | 99.8% |

![Overall Selectivity Ranking](results/figures/overview_selectivity_ranking.png)

---

## Reward Function Axis

### 1.1 Aggregation Method

| Method | Selectivity | EGFR | QED | Valid% | Character |
|---|---|---|---|---|---|
| **Geometric** (baseline) | **+0.06** | 6.04 | 0.775 | 58.9% | Non-compensatory: all properties must be satisfied |
| Arithmetic | −0.05 | 6.04 | 0.772 | 100% | Compensatory: high EGFR can mask poor selectivity |
| Tchebycheff | −0.30 | 5.37 | 0.743 | 50.9% | Ultra-conservative: reward signal collapses |

![Aggregation Comparison](results/figures/reward_1_1_aggregation.png)

**Pharmaceutical Interpretation:**  
Only the geometric mean achieves positive selectivity. Its non-compensatory nature — if *any* property scores zero, the entire reward is zero — forces the optimizer to suppress off-target activity at the formula level. Tchebycheff produces reward values so low (0.023) that MCTS cannot distinguish promising from unpromising branches.

### 1.2 Gaussian σ Sweep

| Condition | EGFR σ | Off-target σ | Selectivity | Reward |
|---|---|---|---|---|
| EGFR strict | 1 | 2 | +0.09 | 0.169 |
| **Baseline** | **2** | **2** | **+0.06** | **0.428** |
| EGFR lenient | 4 | 2 | −0.21 | 0.576 |
| Off-target strict | 2 | 1 | −0.13 | 0.114 |

**Pharmaceutical Interpretation:**  
Tightening off-target σ *decreases* selectivity — counterintuitively. The explanation: with σ=1, most molecules receive near-zero off-target scores, eliminating the gradient signal that MCTS needs to navigate toward selective chemical space. This is the **reward signal collapse** phenomenon.

### 1.3 EGFR Weight Sweep

| Weight | Selectivity | EGFR | QED |
|---|---|---|---|
| w=4 | −0.03 | 5.86 | 0.781 |
| w=8 (baseline) | +0.06 | 6.04 | 0.775 |
| w=16 | +0.06 | 6.20 | 0.759 |

![σ and Weight Ablation](results/figures/reward_1_2_1_3_sigma_weight.png)

**Key Insight:** When combined with exploitation (c=0.1), *lowering* the weight to w=4 actually improves selectivity (+0.96 vs +0.84 for w=16). Heavy exploitation already amplifies the EGFR signal; additional weight creates imbalance.

### Cross Experiments & Long-Run Validation

| Condition | Selectivity | EGFR | QED | Generations |
|---|---|---|---|---|
| c=0.1 baseline | +0.41 | 6.35 | 0.759 | 5K |
| c=0.1 + w=4 | **+0.96** | 6.66 | 0.750 | 5K |
| c=0.1 + w=16 | +0.84 | 6.83 | 0.658 | 5K |
| c=0.1 (long-run) | +1.02 | 7.02 | 0.653 | 30K |
| c=0.1 + EGFR σ=1 (long-run) | **+1.08** | 7.48 | 0.666 | 30K |

![Cross and Long-Run Results](results/figures/cross_longrun.png)

Exploitation effects accumulate over search length. At 30K generations, selectivity exceeds +1.0.

### Reward Axis Lessons

1. **Aggregation method determines optimization direction** — geometric mean is essential for selectivity
2. **σ controls reward signal strength** — too tight causes signal collapse, too loose removes selection pressure
3. **Weight controls inter-property trade-off** — under exploitation, lower weight is paradoxically better

---

## MCTS Policy Axis

### 2.1 c_val Sweep

| c_val | Selectivity | EGFR | QED | Unique% |
|---|---|---|---|---|
| 0.01 | +0.17 | 6.25 | 0.673 | 76.9% |
| **0.1** | **+0.41** | **6.35** | 0.759 | 91.8% |
| 0.5 | +0.13 | 6.08 | 0.775 | 99.6% |
| 1.0 (paper default) | +0.06 | 6.04 | 0.775 | 99.8% |
| 2.0 | +0.03 | 6.03 | 0.764 | 99.9% |

![c_val Sweep](results/figures/mcts_2_1_cval_sweep.png)

**c=0.1 achieves 7× the selectivity of the paper default (c=1.0).** Moderate exploitation causes MCTS to repeatedly explore EGFR-selective chemical subspaces, while excessive exploration (c≥1.0) spreads search too thin.

### 2.2 PUCT — RNN Prior Integration (Custom Implementation)

**Implementation:** New file `policy/puct.py` + patches to `mcts.py` and `utils.py`

```
PUCT = Q(s,a)/N(s,a) + c_puct × P(s,a) × √N(s) / (1 + N(s,a))
```

where P(s,a) is the RNN's predicted probability for token *a* at state *s*.

| Policy | Selectivity | EGFR | QED | Unique% |
|---|---|---|---|---|
| UCB1 c=1.0 | +0.06 | 6.04 | 0.775 | 99.8% |
| **PUCT cp=0.5** | **+0.24** | **6.29** | 0.747 | 95.9% |
| PUCT cp=1.0 | +0.22 | 6.23 | 0.765 | 99.7% |
| PUCT cp=2.0 | +0.02 | 6.09 | 0.761 | 100% |

![PUCT Results](results/figures/mcts_2_2_puct.png)

**Pharmaceutical Interpretation:**  
PUCT (cp=0.5) achieves 4× the selectivity of UCB1 at the same exploration level. The RNN prior steers search toward chemically plausible token sequences, avoiding wasted rollouts on syntactically invalid or pharmacologically irrelevant molecules.

### 2.3 Adaptive c_val — Dynamic Scheduling (Custom Implementation)

**Implementation:** New file `policy/adaptive_ucb1.py` + patch to `mcts.py`

| Policy | Selectivity | EGFR | QED | Unique% |
|---|---|---|---|---|
| UCB1 c=1.0 (fixed) | +0.06 | 6.04 | 0.775 | 99.8% |
| UCB1 c=0.1 (fixed) | +0.41 | 6.35 | 0.759 | 91.8% |
| Linear 1.0→0.1 | +0.23 | 6.20 | 0.781 | 99.0% |
| Linear 1.0→0.01 | +0.56 | 6.51 | 0.749 | 97.4% |
| Linear 0.5→0.1 | +0.14 | 6.14 | 0.797 | 98.6% |
| **Step 1.0→0.1** | **+0.64** | **6.56** | **0.778** | 93.7% |
| Step 1.0→0.01 | +0.79 | 6.84 | 0.676 | 85.7% |

![Adaptive Scheduling Comparison](results/figures/mcts_2_3_adaptive.png)

**Step(1.0→0.1) is the overall best configuration:**
- Selectivity +0.64 (56% higher than fixed c=0.1)
- QED 0.778 (higher than fixed c=0.1)
- Unique% 93.7% (sufficient diversity)

The strategy: **explore broadly for the first 50% of generations** (c=1.0, mapping chemical space), then **exploit aggressively for the remaining 50%** (c=0.1, refining within EGFR-selective regions). This two-phase approach outperforms any fixed policy.

---

## Key Discovery: Reward Signal Collapse

![Reward Signal Collapse](results/figures/reward_signal_collapse.png)

Observed independently in three experiments (Tchebycheff aggregation, off-target σ=1, EGFR σ=1):

> **When reward signal is too weak or uniform, MCTS loses its ability to discriminate between promising and unpromising branches.**

MCTS selects branches based on Q(s,a)/N(s,a) differences in the UCB formula. When most molecules receive similarly low rewards, these differences vanish, and tree search degrades to random sampling. This is a fundamental constraint of MCTS-based molecular generation — distinct from mode collapse in VAEs or training instability in GANs.

**Connection to P3:** This same failure mode reappears when prediction models are inaccurate (P4: kMoL with R²=0.087), confirming that reward signal quality is the critical bottleneck across the entire pipeline.

---

## Custom Code

All custom implementations are included in this repository:

| File | Description |
|---|---|
| [`policy/puct.py`](policy/puct.py) | PUCT policy with RNN prior integration |
| [`policy/adaptive_ucb1.py`](policy/adaptive_ucb1.py) | Adaptive UCB1 with linear and step decay |
| [`reward/dscore_reward_ablation.py`](reward/dscore_reward_ablation.py) | Aggregation method branching |
| [`patches/`](patches/) | Diffs for mcts.py and utils.py modifications |

---

## Reproduction

```bash
conda activate elix
cd /path/to/ChemTSv2

# Baseline (geometric mean, c=1.0)
chemtsv2 -c config/setting_dscore.yaml

# Step adaptive (best configuration)
# Requires adaptive_ucb1.py + mcts.py patch
chemtsv2 -c config/setting_dscore_step.yaml
```

See `configs/` for all experiment configurations.

---

## Limitations

- EGFR activity predictions use LightGBM surrogate models, not experimental assays
- 3 random seeds provide limited variance estimation
- ZINC 250K pre-trained RNN may bias toward certain chemical scaffolds
- Selectivity metric is computed from predicted pIC₅₀ values, not measured Ki ratios
