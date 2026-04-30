# Elix Drug Discovery Portfolio

**Systematic exploration of Elix's open-source drug discovery tools — with pharmaceutical interpretation at every step.**

> M.S. in Pharmacy · Applying for Elix AI Engineer  
> Target: EGFR-selective kinase inhibitor design pipeline

---

## The Story: From ML Prediction to Physics-Based Validation

Across five projects and 130+ experiments, two unifying insights emerged:

### 1. Reward Design Spectrum

```
Signal Collapse ◄────── Reward Design Spectrum ──────► Reward Hacking

P2 Tchebycheff    P4 kMoL GAT    P3-A LGBm+AD     P3-B2 Partial    P3-B AD OFF
(reward=0.023)    (R²=0.087)     ★ OPTIMAL         AD OFF           (DSS=0.907)
```

| Failure Mode | What Happens | Where Observed |
|---|---|---|
| **Signal Collapse** | Reward too weak or uniform → MCTS loses direction → random search | P2: Tchebycheff aggregation; P4: low-accuracy kMoL models |
| **Reward Hacking** | Reward dishonestly high → optimizer exploits extrapolation → unrealistic molecules | P3: AD disabled, BO exploits out-of-domain predictions |
| **Optimal** | Reward is strong enough, balanced, and honest | P3-A: LightGBM + AD constraint |

### 2. ML vs Physics-Based Validation Gap

```
ML Prediction (2D patterns)  ──── r = 0.486 ────  GNINA (CNN docking)
                              ──── r = −0.083 ───  Vina (force-field docking)
```

ML-predicted top candidates and physics-based docking top candidates are **different molecules**. Consensus ranking across both tools recovers candidates that neither would find alone — and structural diversity sampling outperforms ML-top selection for docking validation.

---

## Projects

| # | Project | Tool | Experiments | Key Result |
|---|---|---|---|---|
| **P1** | [Property Prediction](P1_property_prediction/) | kMoL | 20+ configs | GNN layer choice depends on ADMET task; FP+XGBoost beats GNN on hepatic clearance |
| **P2** | [Molecular Generation](P2_molecular_generation/) | ChemTSv2 | 33 configs | Step(1.0→0.1) adaptive scheduling: selectivity +0.64, 10× over default |
| **P3+P4** | [Multi-objective Optimization](P3_multiobjective/) | DyRAMO + kMoL | 27 configs | AD threshold cliff at 0.7; Perm AD most critical; kMoL signal collapse reconfirmed |
| **P5** | [Docking Benchmark](P5_docking/) | GNINA + Vina | 4,350 docks | GNINA-Vina consensus identifies mol_0147 as top candidate; diversity sampling finds 6/10 consensus top |

---

## Pipeline Architecture

```
P1. kMoL              P2. ChemTSv2           P3. DyRAMO             P5. Docking
┌────────────┐        ┌──────────────┐        ┌──────────────┐       ┌──────────────┐
│ GNN ADMET  │─reward─│ MCTS molecule│◄──BO──►│ Multi-obj    │─top──►│ GNINA (CNN)  │
│ prediction │ signal │ generation   │        │ optimization │ 200   │ Vina (physics)│
│ (5 layers) │        │ (33 configs) │        │ (AD ablation)│       │ (7 PDB × 2)  │
└────────────┘        └──────────────┘        └──────────────┘       └──────────────┘
                                                                      ↓
                                                                   Consensus
                                                                   ranking
```

---

## What Makes This Different

**I'm not just an ML engineer who tunes configs.** I'm a pharmacist who understands *why* GIN's sum aggregation captures solubility-determining functional groups, *why* Osimertinib's docking score doesn't reflect its true T790M selectivity (covalent binding is not modeled), and *why* permeability models are most vulnerable to extrapolation.

Every experiment includes pharmaceutical interpretation that a general ML engineer cannot provide.

### Code Contributions

| File | Lines | Description |
|---|---|---|
| `policy/puct.py` | 15 | PUCT policy integrating RNN prior into MCTS search |
| `policy/adaptive_ucb1.py` | 20 | Adaptive UCB1 with linear/step decay scheduling |
| `reward/dscore_reward_ablation.py` | 4 | Aggregation method branching (geometric/arithmetic/Tchebycheff) |
| `chemtsv2/mcts.py` patch | 10 | Prior probability storage + progress tracking |
| `chemtsv2/utils.py` patch | 5 | RNN prior probability return |
| `dock_candidates.py` | 350 | 7-PDB GNINA docking pipeline with 2-stage refinement |
| `vina_validation.py` | 400 | Vina cross-validation with GNINA-Vina consensus analysis |

---

## Key Numbers

| Metric | Value |
|---|---|
| Total experiments | 130+ configurations |
| Total docking runs | ~4,350 (GNINA GPU + Vina CPU) |
| PDB structures used | 7 (WT × 4, T790M × 2, allosteric × 1) |
| Reference drugs benchmarked | 7 (1G–4G, multi-seed) |
| Generated molecules docked | 200 (stratified: top50 + random100 + diverse50) |
| Best consensus candidate | mol_0147 (GNINA #2, Vina #1) |
| GNINA-Vina correlation | r = −0.46 (moderate agreement) |

---

## Tools Used

| Tool | Developer | Purpose |
|---|---|---|
| [kMoL](https://github.com/elix-tech/kmol) | Elix | GNN-based molecular property prediction |
| [ChemTSv2](https://github.com/molecule-generator-collection/ChemTSv2) | Elix | MCTS-based molecular generation |
| [DyRAMO](https://github.com/molecule-generator-collection/DyRAMO) | Elix | Dynamic reliability-aware multi-objective optimization |
| [GNINA](https://github.com/gnina/gnina) | Koes Lab | CNN-based molecular docking |
| [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) | Scripps | Physics-based molecular docking |

---

## Repository Structure

```
elix-drug-discovery-portfolio/
├── README.md                          ← You are here
├── P1_property_prediction/            ← kMoL: GNN layer × ADMET task ablation
├── P2_molecular_generation/           ← ChemTSv2: Reward + MCTS policy ablation (33 configs)
├── P3_multiobjective/                 ← DyRAMO + kMoL integration: AD ablation + pipeline test (27 configs)
├── P4_integrated_pipeline/            ← (See P3 Section 4)
└── P5_docking/                        ← GNINA + Vina: 7 PDB × 200 molecules × 2 tools
```

---

## Quick Navigation

- **Best single result:** [P2 — Step adaptive scheduling](P2_molecular_generation/#23-adaptive-c_val) (Selectivity +0.64)
- **Most impactful finding:** [P2 — Reward Signal Collapse](P2_molecular_generation/#key-discovery-reward-signal-collapse)
- **Cross-project insight:** [P3 — Reward Design Spectrum](P3_multiobjective/#summary-four-layers-of-understanding)
- **Physics validation:** [P5 — GNINA vs Vina Consensus](P5_docking/#5-consensus-ranking--top-10)
- **Code I wrote:** [P2 — Custom Policies](P2_molecular_generation/policy/) · [P5 — Docking Scripts](P5_docking/scripts/)
