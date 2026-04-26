# P5. Molecular Docking with AutoDock Vina

**Status:** Coming Soon

---

## Planned Work

Validate top candidates from P3 multi-objective optimization using physics-based docking simulations against the EGFR binding pocket (PDB: 1M17).

### Goals

- Dock top-5 molecules from P3-A (AD ON) against EGFR crystal structure
- Compare predicted pIC₅₀ (LightGBM surrogate) vs docking scores (physics-based)
- Analyze binding poses to verify EGFR selectivity at the structural level
- Provide an independent validation layer for the ML-driven pipeline

### Expected Deliverables

- Docking scores and binding pose visualizations
- Correlation analysis: surrogate model predictions vs docking energy
- Pharmaceutical interpretation of binding interactions
