"""
P3. DyRAMO AD Ablation + P4. kMoL Integration — Analysis & Visualization (v2)
All 5 P3 experiments (A, B, B2, C, D) fully compared
"""

import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams.update({
    'figure.facecolor': '#FAFAFA',
    'axes.facecolor': '#FAFAFA',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'axes.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.3,
})

outdir = '/home/claude/dyramo_figures_v2'
os.makedirs(outdir, exist_ok=True)

# ── Common Data ──
EXP_LABELS = ['A: AD ON\n(Baseline)', 'B: AD OFF\n(All removed)', 'B2: EGFR AD\nOFF (Partial)',
              'C: Fixed\n(No BO)', 'D: EGFR Only\n(Single obj)']
EXP_COLORS = ['#059669', '#DC2626', '#D97706', '#F59E0B', '#7C3AED']
EXP_MARKERS = ['o', 's', 'D', '^', 'v']


# ════════════════════════════════════════════════════════════════
# 1. P3 Overview — DSS Score + Reliability Levels (All 5)
# ════════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('P3: AD Ablation — All 5 Experiments Overview',
             fontsize=16, fontweight='bold', y=1.02)

# DSS comparison
dss_scores = [0.369, 0.907, 0.493, 0.919, 0.458]
bars = ax1.barh(range(len(EXP_LABELS)), dss_scores, color=EXP_COLORS, height=0.6, edgecolor='white')
for bar, val in zip(bars, dss_scores):
    ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
             f'{val:.3f}', va='center', fontsize=11, fontweight='bold')
ax1.set_yticks(range(len(EXP_LABELS)))
ax1.set_yticklabels(EXP_LABELS, fontsize=10)
ax1.set_xlabel('Best DSS Score')
ax1.set_title('DSS Score Comparison')
ax1.invert_yaxis()
ax1.axvline(x=0.5, color='gray', linestyle=':', alpha=0.5)
ax1.text(0.52, 4.5, 'Realistic ← | → Inflated', fontsize=8, color='gray')

# Reliability levels grouped bar
props = ['EGFR', 'Stab', 'Perm']
a_lev = [0.62, 0.46, 0.34]
b_lev = [0.90, 0.84, 0.89]
b2_lev = [0.90, 0.46, 0.45]
c_lev = [0.50, 0.50, 0.50]
d_lev = [0.71, 0.43, 0.43]

x = np.arange(len(props))
w = 0.15
ax2.bar(x - 2*w, a_lev, w, label='A: AD ON', color=EXP_COLORS[0])
ax2.bar(x - w, b_lev, w, label='B: AD OFF', color=EXP_COLORS[1])
ax2.bar(x, b2_lev, w, label='B2: Partial', color=EXP_COLORS[2])
ax2.bar(x + w, c_lev, w, label='C: Fixed', color=EXP_COLORS[3])
ax2.bar(x + 2*w, d_lev, w, label='D: EGFR only', color=EXP_COLORS[4])

ax2.set_xticks(x)
ax2.set_xticklabels(props, fontsize=12)
ax2.set_ylabel('Reliability Level')
ax2.set_title('BO-Selected Reliability Levels')
ax2.set_ylim(0, 1.0)
ax2.legend(fontsize=8, loc='upper right', ncol=2)

plt.tight_layout()
plt.savefig(f'{outdir}/p3_overview_all5.png')
plt.close()
print("✅ P3 overview (all 5) saved")


# ════════════════════════════════════════════════════════════════
# 2. P3 — Top Molecules Comparison (All 5 experiments)
# ════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('P3: Top Molecules — All 5 Experiments Compared',
             fontsize=16, fontweight='bold', y=0.98)

# ── Data: top 5 molecules per experiment ──
# EGFR pIC50
egfr_top = {
    'A': [7.97, 7.94, 7.94, 7.92, 7.87],
    'B': [7.48, 7.39, 7.51, 7.41, 7.35],
    'B2': [6.57, 6.64, 6.76, 6.20, 6.19],
    'C': [7.40, 7.34, 7.35, 7.40, 7.40],
    'D': [7.98, 7.94, 7.92, 7.88, 7.86],
}

# EGFR Tanimoto similarity
sim_top = {
    'A': [0.661, 1.0, 1.0, 0.722, 0.829],
    'B': [0.217, 0.212, 0.509, 0.228, 0.246],
    'B2': [0.395, 0.316, 0.364, 0.468, 0.347],
    'C': [0.633, 0.633, 0.590, 0.755, 0.755],
    'D': [0.911, 1.0, 1.0, 0.872, 0.712],
}

# Reward
reward_top = {
    'A': [0.957, 0.954, 0.954, 0.950, 0.878],
    'B': [0.885, 0.876, 0.872, 0.868, 0.864],
    'B2': [0.754, 0.739, 0.723, 0.720, 0.718],
    'C': [0.899, 0.892, 0.892, 0.876, 0.876],
    'D': [0.878, 0.869, 0.863, 0.852, 0.849],
}

# Metabolic stability
stab_top = {
    'A': [76.1, 73.8, 73.8, 75.4, 43.5],
    'B': [87.4, 76.2, 64.7, 73.0, 78.8],
    'B2': [70.3, 69.6, 58.7, 73.2, 79.1],
    'C': [80.9, 76.3, 73.6, 67.3, 67.3],
    'D': [29.2, 73.8, 46.3, 70.5, 43.5],
}

x_pos = range(5)
exp_keys = ['A', 'B', 'B2', 'C', 'D']
plot_colors = [EXP_COLORS[0], EXP_COLORS[1], EXP_COLORS[2], EXP_COLORS[3], EXP_COLORS[4]]
plot_markers = ['o', 's', 'D', '^', 'v']
plot_labels = ['A: AD ON', 'B: AD OFF', 'B2: EGFR AD OFF', 'C: Fixed (No BO)', 'D: EGFR Only']

# Plot 1: EGFR pIC50
ax = axes[0, 0]
for key, color, marker, label in zip(exp_keys, plot_colors, plot_markers, plot_labels):
    ax.plot(x_pos, egfr_top[key], f'{marker}-', color=color, label=label, linewidth=2, markersize=8)
ax.set_title('EGFR pIC50 (Top 5 Molecules)')
ax.set_ylabel('pIC50')
ax.set_xticks(x_pos)
ax.set_xticklabels(['1st', '2nd', '3rd', '4th', '5th'])
ax.legend(fontsize=9)

# Plot 2: Tanimoto similarity
ax = axes[0, 1]
for key, color, marker, label in zip(exp_keys, plot_colors, plot_markers, plot_labels):
    ax.plot(x_pos, sim_top[key], f'{marker}-', color=color, label=label, linewidth=2, markersize=8)
ax.set_title('EGFR Tanimoto Similarity (vs Training Data)')
ax.set_ylabel('Similarity')
ax.set_xticks(x_pos)
ax.set_xticklabels(['1st', '2nd', '3rd', '4th', '5th'])
ax.axhline(y=0.3, color='gray', linestyle=':', alpha=0.5)
ax.text(4.3, 0.32, 'AD\nthreshold', fontsize=8, color='gray')
ax.legend(fontsize=9)

# Plot 3: Reward
ax = axes[1, 0]
for key, color, marker, label in zip(exp_keys, plot_colors, plot_markers, plot_labels):
    ax.plot(x_pos, reward_top[key], f'{marker}-', color=color, label=label, linewidth=2, markersize=8)
ax.set_title('Reward Score (Top 5)')
ax.set_ylabel('Reward')
ax.set_xticks(x_pos)
ax.set_xticklabels(['1st', '2nd', '3rd', '4th', '5th'])
ax.legend(fontsize=9)

# Plot 4: Metabolic Stability
ax = axes[1, 1]
for key, color, marker, label in zip(exp_keys, plot_colors, plot_markers, plot_labels):
    ax.plot(x_pos, stab_top[key], f'{marker}-', color=color, label=label, linewidth=2, markersize=8)
ax.set_title('Metabolic Stability (Top 5)')
ax.set_ylabel('% Remaining')
ax.set_xticks(x_pos)
ax.set_xticklabels(['1st', '2nd', '3rd', '4th', '5th'])
ax.axhline(y=50, color='gray', linestyle=':', alpha=0.5)
ax.text(4.3, 52, 'Min\ntarget', fontsize=8, color='gray')
ax.legend(fontsize=9)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(f'{outdir}/p3_top_molecules_all5.png')
plt.close()
print("✅ P3 top molecules (all 5) saved")


# ════════════════════════════════════════════════════════════════
# 3. P3 — AD Insight Scatter (All 5 experiments)
# ════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(14, 8))

# Scatter: EGFR pIC50 vs Tanimoto similarity for all experiments
scatter_data = {
    'A: AD ON': {
        'egfr': [7.97, 7.94, 7.94, 7.92, 7.87, 7.88, 7.86, 7.85],
        'sim': [0.661, 1.0, 1.0, 0.722, 0.829, 0.712, 0.736, 0.911],
        'color': '#059669', 'marker': 'o',
    },
    'B: AD OFF': {
        'egfr': [7.48, 7.39, 7.51, 7.41, 7.35, 7.17, 7.41, 7.39],
        'sim': [0.217, 0.212, 0.509, 0.228, 0.245, 0.212, 0.228, 0.212],
        'color': '#DC2626', 'marker': 's',
    },
    'B2: EGFR AD OFF': {
        'egfr': [6.57, 6.64, 6.76, 6.20, 6.19, 6.30, 6.30, 6.30],
        'sim': [0.395, 0.316, 0.364, 0.468, 0.347, 0.413, 0.413, 0.413],
        'color': '#D97706', 'marker': 'D',
    },
    'D: EGFR Only': {
        'egfr': [7.98, 7.94, 7.92, 7.88, 7.86, 7.87, 7.86, 7.85],
        'sim': [0.911, 1.0, 1.0, 0.872, 0.712, 0.736, 0.829, 0.559],
        'color': '#7C3AED', 'marker': 'v',
    },
    'C: Fixed (No BO)': {
        'egfr': [7.40, 7.34, 7.35, 7.40, 7.15, 7.20, 7.12, 7.33],
        'sim': [0.633, 0.633, 0.590, 0.755, 0.691, 0.707, 0.623, 0.688],
        'color': '#F59E0B', 'marker': '^',
    },
}

for label, data in scatter_data.items():
    ax.scatter(data['sim'], data['egfr'], s=120, c=data['color'], marker=data['marker'],
               label=label, edgecolors='#333', linewidth=1, zorder=5, alpha=0.85)

# Regions
ax.axvspan(0, 0.3, alpha=0.06, color='red')
ax.axvspan(0.6, 1.05, alpha=0.06, color='green')

ax.text(0.15, 6.0, 'Out-of-domain\n→ Unreliable\npredictions',
        fontsize=9, ha='center', color='#DC2626', style='italic',
        bbox=dict(boxstyle='round', facecolor='#FEE2E2', alpha=0.7))
ax.text(0.82, 8.05, 'In-domain\n→ Reliable\npredictions',
        fontsize=9, ha='center', color='#059669', style='italic',
        bbox=dict(boxstyle='round', facecolor='#D1FAE5', alpha=0.7))

# Key observations as annotations
ax.annotate('B: High EGFR but\nlow similarity → Hacking',
            xy=(0.22, 7.5), xytext=(0.35, 7.8),
            arrowprops=dict(arrowstyle='->', color='#DC2626', lw=1.5),
            fontsize=9, color='#DC2626', fontweight='bold')

ax.annotate('A & D: High EGFR AND\nhigh similarity → Trustworthy',
            xy=(0.9, 7.95), xytext=(0.55, 8.15),
            arrowprops=dict(arrowstyle='->', color='#059669', lw=1.5),
            fontsize=9, color='#059669', fontweight='bold')

ax.annotate('B2: Mid range — EGFR\ninflated but partially\nconstrained',
            xy=(0.39, 6.6), xytext=(0.55, 6.2),
            arrowprops=dict(arrowstyle='->', color='#D97706', lw=1.5),
            fontsize=9, color='#D97706', fontweight='bold')

ax.annotate('C: No BO but AD intact\n→ moderate EGFR,\nhigh similarity',
            xy=(0.65, 7.35), xytext=(0.38, 7.55),
            arrowprops=dict(arrowstyle='->', color='#F59E0B', lw=1.5),
            fontsize=9, color='#B45309', fontweight='bold')

ax.set_xlabel('Tanimoto Similarity (vs Training Data)', fontsize=12)
ax.set_ylabel('Predicted EGFR pIC50', fontsize=12)
ax.set_title('P3: EGFR pIC50 vs Training Data Similarity — All Experiments\n'
             'AD prevents out-of-domain reward hacking',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=10, loc='lower left')
ax.set_xlim(0.1, 1.05)
ax.set_ylim(5.8, 8.3)

plt.tight_layout()
plt.savefig(f'{outdir}/p3_ad_insight_all5.png')
plt.close()
print("✅ P3 AD insight (all 5) saved")


# ════════════════════════════════════════════════════════════════
# 4. P3 — Experiment D Deep Dive (Single vs Multi-objective)
# ════════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('P3: Single-Objective (D) vs Multi-Objective (A) — The ADMET Trade-off',
             fontsize=15, fontweight='bold', y=1.02)

# Top 5 comparison: A vs D
metrics = ['EGFR\npIC50', 'Metabolic\nStability', 'Permeability', 'EGFR\nSimilarity', 'Reward']

a_vals = [7.94, 74.8, 1.05, 0.88, 0.953]  # averaged top 5
d_vals = [7.93, 44.6, 1.08, 0.92, 0.864]  # averaged top 5

x = np.arange(len(metrics))
w = 0.3
bars_a = ax1.bar(x - w/2, a_vals, w, label='A: Multi-obj + AD', color='#059669')
bars_d = ax1.bar(x + w/2, d_vals, w, label='D: EGFR only', color='#7C3AED')

ax1.set_xticks(x)
ax1.set_xticklabels(metrics, fontsize=9)
ax1.set_ylabel('Score')
ax1.set_title('Average Top 5 Molecules')
ax1.legend(fontsize=10)

# Highlight the key difference
ax1.annotate('D loses\nStab significantly!',
             xy=(1.15, 44.6), xytext=(2, 30),
             arrowprops=dict(arrowstyle='->', color='#DC2626', lw=2),
             fontsize=10, color='#DC2626', fontweight='bold')

# Stab distribution comparison
a_stab_all = [76.1, 73.8, 73.8, 75.4, 70.5, 43.5, 47.2, 38.7]
d_stab_all = [29.2, 73.8, 46.3, 70.5, 43.5, 38.7, 47.2, 29.2]

ax2.hist(a_stab_all, bins=8, alpha=0.7, color='#059669', label='A: Multi-obj', edgecolor='white')
ax2.hist(d_stab_all, bins=8, alpha=0.7, color='#7C3AED', label='D: EGFR only', edgecolor='white')
ax2.axvline(x=50, color='#DC2626', linestyle='--', linewidth=2, label='Min target (50%)')
ax2.set_xlabel('Metabolic Stability (%)')
ax2.set_ylabel('Count')
ax2.set_title('Metabolic Stability Distribution\n(Top Molecules)')
ax2.legend(fontsize=9)

plt.tight_layout()
plt.savefig(f'{outdir}/p3_single_vs_multi.png')
plt.close()
print("✅ P3 single vs multi saved")


# ════════════════════════════════════════════════════════════════
# 5. P4 — kMoL GNN vs LightGBM (unchanged from v1)
# ════════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('P4: kMoL GNN Integration vs LightGBM Baseline',
             fontsize=16, fontweight='bold', y=1.02)

models = ['LightGBM\n(P3-A)', 'kMoL LEConv\n(R²=0.287)', 'kMoL GIN\n(R²=0.227)', 'kMoL GAT\n(R²=0.087)']
dss = [0.369, 0.016, 0.013, 0.001]
model_colors = ['#059669', '#2563EB', '#7C3AED', '#D97706']

bars = ax1.bar(models, dss, color=model_colors, width=0.6, edgecolor='white')
for bar, v in zip(bars, dss):
    label = f'{v:.3f}' if v > 0.05 else f'{v:.4f}'
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
             label, ha='center', fontsize=11, fontweight='bold')
ax1.set_ylabel('Best DSS Score')
ax1.set_title('DSS Score by Prediction Model')

# Top molecule EGFR + Reward
ax2_models = ['LightGBM', 'LEConv', 'GIN', 'GAT']
top_egfr = [7.97, 3.31, 3.16, 1.37]
model_r2 = [None, 0.287, 0.227, 0.087]

bars_e = ax2.bar(range(4), top_egfr, color=model_colors, width=0.6, edgecolor='white')
for i, (bar, v) in enumerate(zip(bars_e, top_egfr)):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
             f'{v:.2f}', ha='center', fontsize=11, fontweight='bold')
    if model_r2[i]:
        ax2.text(bar.get_x() + bar.get_width()/2, 0.3,
                 f'R²={model_r2[i]}', ha='center', fontsize=9, color='white', fontweight='bold')

ax2.set_xticks(range(4))
ax2.set_xticklabels(ax2_models)
ax2.set_ylabel('Best EGFR pIC50')
ax2.set_title('Top Molecule EGFR Quality')

ax2.annotate('Low R² → weak reward\n→ signal collapse',
             xy=(3, 1.37), xytext=(2, 4),
             arrowprops=dict(arrowstyle='->', color='#DC2626', lw=1.5),
             fontsize=10, color='#DC2626', fontweight='bold')

plt.tight_layout()
plt.savefig(f'{outdir}/p4_kmol_vs_lgb.png')
plt.close()
print("✅ P4 kMoL vs LightGBM saved")


# ════════════════════════════════════════════════════════════════
# 6. Grand Story — Reward Design Spectrum (P2→P3→P4)
# ════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(16, 7))

spectrum_labels = [
    'P2: Tchebycheff\n(Signal Collapse)',
    'P4: kMoL GAT\n(Weak Model)',
    'P4: kMoL LEConv\n(Weak Model)',
    'P3-A: LightGBM\n+ AD ON\n★ OPTIMAL',
    'P3-D: EGFR Only\n(Imbalanced)',
    'P3-B2: Partial\nAD OFF',
    'P3-B: All AD OFF\n(Full Hacking)',
]
quality = [0.05, 0.12, 0.18, 0.90, 0.60, 0.45, 0.10]
colors_spectrum = ['#DC2626', '#F59E0B', '#FBBF24', '#059669', '#7C3AED', '#D97706', '#DC2626']

bars = ax.bar(range(len(spectrum_labels)), quality, color=colors_spectrum, width=0.7, edgecolor='white')
bars[3].set_edgecolor('#111')
bars[3].set_linewidth(3)

for bar, v, label in zip(bars, quality, spectrum_labels):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f'{v:.2f}', ha='center', fontsize=10, fontweight='bold')

ax.set_xticks(range(len(spectrum_labels)))
ax.set_xticklabels(spectrum_labels, fontsize=9, ha='center')
ax.set_ylabel('Effective Optimization Quality', fontsize=12)
ax.set_title('The Reward Design Spectrum: Signal Collapse ↔ Optimal ↔ Reward Hacking\n'
             'P2 + P3 + P4 unified: Reward must be strong enough, balanced, AND honest',
             fontsize=14, fontweight='bold')

# Region labels
ax.annotate('', xy=(-0.3, -0.06), xytext=(2.3, -0.06),
            arrowprops=dict(arrowstyle='<->', color='#DC2626', lw=2),
            annotation_clip=False)
ax.text(1, -0.11, 'Reward Too Weak', fontsize=10, color='#DC2626',
        ha='center', fontweight='bold', transform=ax.get_xaxis_transform())

ax.annotate('', xy=(4.7, -0.06), xytext=(6.3, -0.06),
            arrowprops=dict(arrowstyle='<->', color='#DC2626', lw=2),
            annotation_clip=False)
ax.text(5.5, -0.11, 'Reward Hacking', fontsize=10, color='#DC2626',
        ha='center', fontweight='bold', transform=ax.get_xaxis_transform())

ax.text(3, 0.95, '★ Sweet Spot', fontsize=13, color='#059669',
        ha='center', fontweight='bold')

ax.set_ylim(0, 1.05)

plt.tight_layout()
plt.savefig(f'{outdir}/grand_reward_spectrum.png')
plt.close()
print("✅ Grand reward spectrum saved")


print(f"\n📁 All P3/P4 v2 figures saved to {outdir}/")
