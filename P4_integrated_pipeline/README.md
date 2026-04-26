# P4. Integrated Pipeline — kMoL → DyRAMO

**Tools:** [kMoL](https://github.com/elix-tech/kmol) + [ChemTSv2](https://github.com/molecule-generator-collection/ChemTSv2) + [DyRAMO](https://github.com/molecule-generator-collection/DyRAMO) (Elix)  
**Goal:** Replace DyRAMO's default LightGBM with kMoL GNN models and measure how prediction accuracy propagates through the generation pipeline  
**Key Finding:** Prediction model accuracy (R²) directly determines reward signal quality, which directly determines generated molecule quality — confirming P2's reward signal collapse at the pipeline level

---

## Motivation

P1 trained kMoL GNN models for ADMET prediction. P2 discovered that weak reward signals cause MCTS to degrade to random search. P3 showed that DyRAMO with LightGBM + AD produces realistic candidates (DSS = 0.369).

**The question:** What happens when we swap LightGBM (R² ≈ 0.7+) for kMoL GNN models (R² = 0.08–0.29)?

This experiment connects all three Elix tools into a single pipeline and tests the weakest link hypothesis: **the entire pipeline is only as good as its least accurate component.**

---

## Pipeline Architecture

```
kMoL (P1)                    ChemTSv2 (P2)              DyRAMO (P3)
┌─────────────┐    reward    ┌──────────────┐   BO      ┌───────────┐
│ GNN Model   │───signal───► │ MCTS Search  │◄─────────►│ Bayesian  │
│ (LEConv/    │              │ (molecule    │  optimize  │ Optimizer │
│  GIN/GAT)   │              │  generation) │  AD levels │           │
└─────────────┘              └──────────────┘            └───────────┘
     R²=0.08~0.29                                        
     vs LightGBM R²≈0.7+                                 
```

---

## Results

### kMoL GNN vs LightGBM

| Model | R² (est.) | Best DSS | Top EGFR pIC₅₀ | Top Reward |
|---|---|---|---|---|
| **LightGBM (P3-A)** | **~0.7+** | **0.369** | **7.97** | **0.957** |
| kMoL LEConv | 0.287 | 0.016 | 3.31 | 0.026 |
| kMoL GIN | 0.227 | 0.013 | 3.16 | 0.026 |
| kMoL GAT | 0.087 | 0.001 | 1.37 | 0.021 |

![kMoL vs LightGBM](results/figures/p4_kmol_vs_lgb.png)

The gap is dramatic: LightGBM achieves DSS 0.369 while the best kMoL model (LEConv) achieves only 0.016 — a **23× difference**. Top EGFR pIC₅₀ drops from 7.97 to 3.31, meaning generated molecules are predicted to be ~4,500× less potent.

### R² → DSS Correlation

| R² | DSS | Ratio to LightGBM |
|---|---|---|
| ~0.7+ | 0.369 | 1.00× |
| 0.287 | 0.016 | 0.04× |
| 0.227 | 0.013 | 0.04× |
| 0.087 | 0.001 | 0.003× |

Even within the kMoL models, the ranking is monotonic: LEConv (R²=0.287) > GIN (R²=0.227) > GAT (R²=0.087), confirming that prediction accuracy is the bottleneck.

---

## Signal Collapse Reconfirmed

This is the same **reward signal collapse** discovered in P2, now at the pipeline level:

- **P2 (reward axis):** Tchebycheff aggregation produces reward ≈ 0.023 → MCTS can't distinguish promising branches
- **P4 (model axis):** kMoL GAT produces reward ≈ 0.021 → same collapse, different cause

The mechanism is identical: when Q(s,a)/N(s,a) values in the UCB formula are uniformly low, the exploration term dominates and MCTS degenerates into random sampling.

---

## Reward Design Spectrum — Complete Picture

![Grand Reward Spectrum](results/figures/grand_reward_spectrum.png)

P4 fills in the left side of the spectrum, confirming that signal collapse can originate from either the reward function design (P2) or the prediction model quality (P4):

```
Signal Collapse ◄────── Reward Design Spectrum ──────► Reward Hacking

P2 Tchebycheff    P4 kMoL GAT    P3-A LGBm+AD     P3-B2 Partial    P3-B AD OFF
(reward=0.023)    (R²=0.087)     ★ OPTIMAL         AD OFF           (DSS=0.907)
```

---

## Why kMoL Underperforms LightGBM

This is **not** a GNN-vs-gradient-boosting comparison. The performance gap comes from training data:

| Factor | kMoL | LightGBM |
|---|---|---|
| Training data | Elix ADMET datasets (hundreds–low thousands) | Task-specific kinase datasets (tens of thousands) |
| Features | Learned GNN embeddings | Morgan fingerprints (2048-bit) |
| Task match | General ADMET → applied to EGFR | Trained directly on EGFR-related data |

With equivalent training data, GNN models can match or exceed fingerprint-based models (as shown in literature). The lesson here is about **data appropriateness**, not model architecture.

---

## Reproduction

```bash
# Step 1: Train kMoL model (P1)
conda activate kmol
cd /path/to/kmol
kmol train data/configs/model/adme/solubility.json

# Step 2: Run DyRAMO with kMoL predictor (P4)
conda activate elix
cd /path/to/DyRAMO
python run.py -c config/setting_dyramo_kmol_le.yaml
```

See `configs/` for LEConv, GIN, and GAT configurations.

---

## Limitations

- kMoL models trained on small ADMET datasets — comparison with LightGBM is confounded by data scale
- Only EGFR prediction was swapped; metabolic stability and permeability still used LightGBM
- No hyperparameter tuning of kMoL models for this specific task
- A fairer comparison would require retraining LightGBM on the same small datasets
