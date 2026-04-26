# Elix Drug Discovery Portfolio

**Systematic exploration of Elix's open-source drug discovery tools — with pharmaceutical interpretation at every step.**

> M.S. in Pharmacy · Applying for Elix AI Engineer  
> Target: EGFR-selective kinase inhibitor design pipeline

---

## The Story: Reward Design Spectrum

Across four projects and 100+ experiments, one unifying insight emerged:

```
Signal Collapse ◄────── Reward Design Spectrum ──────► Reward Hacking

P2 Tchebycheff    P4 kMoL GAT    P3-A LGBm+AD     P3-B2 Partial    P3-B AD OFF
(reward=0.023)    (R²=0.087)     ★ OPTIMAL         AD OFF           (DSS=0.907)
```

| Failure Mode | What Happens | Where Observed |
|---|---|---|
| **Signal Collapse** | Reward too weak or uniform → MCTS loses direction → random search | P2: Tchebycheff aggregation, off-target σ=1; P4: low-accuracy kMoL models |
| **Reward Hacking** | Reward dishonestly high → optimizer exploits extrapolation region → unrealistic molecules | P3: AD disabled, BO exploits out-of-domain predictions |
| **Optimal** | Reward is strong enough, balanced, and honest | P3-A: LightGBM + AD constraint; P2: geometric mean + Step(1.0→0.1) |

**The pharmacist's take:** A drug discovery pipeline is only as good as its reward signal. Too weak, and the generative model wanders. Too generous, and it cheats. The sweet spot requires domain-aware calibration — exactly where pharmaceutical expertise meets ML engineering.

---

## Projects

| # | Project | Tool | Key Result |
|---|---|---|---|
| **P1** | [Property Prediction](P1_property_prediction/) | kMoL | Optimal GNN layer differs by ADMET task: GIN for classification, LEConv for regression, GAT for protein binding |
| **P2** | [Molecular Generation](P2_molecular_generation/) | ChemTSv2 | Step(1.0→0.1) adaptive scheduling achieves selectivity +0.64 — 10× over the original paper's default |
| **P3** | [Multi-objective Optimization](P3_multiobjective/) | DyRAMO | AD ablation reveals reward hacking: DSS inflates from 0.369 to 0.907 when AD is removed |
| **P4** | [Integrated Pipeline](P4_integrated_pipeline/) | kMoL → DyRAMO | Prediction model accuracy directly determines generation quality — signal collapse reappears |
| **P5** | [Molecular Docking](P5_docking/) | AutoDock Vina | Coming Soon |

---

## What Makes This Different

**I'm not just an ML engineer who tunes configs.** I'm a pharmacist who understands *why* GIN's sum aggregation captures solubility-determining functional groups, *why* geometric mean enforces off-target suppression at the formula level, and *why* reward signal collapse is dangerous for real drug discovery pipelines.

Every experiment includes pharmaceutical interpretation that a general ML engineer cannot provide.

### Code Contributions

Beyond configuration, I implemented core algorithmic extensions for ChemTSv2:

| File | Lines | Description |
|---|---|---|
| `policy/puct.py` | 15 | PUCT policy integrating RNN prior into MCTS search |
| `policy/adaptive_ucb1.py` | 20 | Adaptive UCB1 with linear/step decay scheduling |
| `reward/dscore_reward_ablation.py` | 4 | Aggregation method branching (geometric/arithmetic/Tchebycheff) |
| `chemtsv2/mcts.py` patch | 10 | Prior probability storage + progress tracking |
| `chemtsv2/utils.py` patch | 5 | RNN prior probability return |

~50 lines total — but each requires understanding MCTS internals to implement correctly.

---

## Tools Used

| Tool | Developer | Purpose | Repository |
|---|---|---|---|
| [kMoL](https://github.com/elix-tech/kmol) | Elix | GNN-based molecular property prediction | ADMET modeling |
| [ChemTSv2](https://github.com/molecule-generator-collection/ChemTSv2) | Elix | MCTS-based molecular generation | EGFR-selective inhibitor design |
| [DyRAMO](https://github.com/molecule-generator-collection/DyRAMO) | Elix | Dynamic reliability-aware multi-objective optimization | Automated reward calibration |

---

## Repository Structure

```
elix-drug-discovery-portfolio/
├── README.md                          ← You are here
├── P1_property_prediction/            ← kMoL: GNN layer × ADMET task ablation
├── P2_molecular_generation/           ← ChemTSv2: Reward + MCTS policy ablation
├── P3_multiobjective/                 ← DyRAMO: AD ablation (reward hacking)
├── P4_integrated_pipeline/            ← kMoL → DyRAMO: prediction accuracy → generation quality
├── P5_docking/                        ← Coming Soon
└── .gitignore
```

---

## Quick Navigation

- **Best single result:** [P2 — Step adaptive scheduling](P2_molecular_generation/#23-adaptive-c_val) (Selectivity +0.64)
- **Most impactful finding:** [P2 — Reward Signal Collapse](P2_molecular_generation/#key-discovery-reward-signal-collapse)
- **Cross-project insight:** [P3 — Reward Design Spectrum](P3_multiobjective/#reward-design-spectrum)
- **Code I wrote:** [P2 — Custom Policies](P2_molecular_generation/policy/)
