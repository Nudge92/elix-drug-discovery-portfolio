"""
P2. ChemTSv2 MCTS Molecular Generation — Ablation Study Visualization
Reward axis (aggregation, σ, weight, cross/long-run) + MCTS policy axis (c_val, PUCT, Adaptive)
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

# ── Style ──
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

outdir = '/home/claude/chemts_figures'
os.makedirs(outdir, exist_ok=True)


# ════════════════════════════════════════════════════════════════
# 1. Grand Overview — All Conditions Selectivity Bar Chart
# ════════════════════════════════════════════════════════════════

conditions = [
    'Step(1→0.1)\nAdaptive',
    'Linear(1→0.01)\nAdaptive',
    'UCB1\nc=0.1',
    'PUCT\ncp=0.5',
    'Linear(1→0.1)\nAdaptive',
    'Linear(0.5→0.1)\nAdaptive',
    'c=0.01\n(fixed)',
    'c=0.5\n(fixed)',
    'UCB1 c=1.0\n(Original)',
    'c=2.0\n(fixed)',
    'PUCT cp=2.0',
]
selectivity = [0.64, 0.56, 0.41, 0.24, 0.23, 0.14, 0.17, 0.13, 0.06, 0.03, 0.02]

# Color by category
cat_colors = {
    'adaptive': '#DC2626',
    'fixed_good': '#2563EB',
    'puct': '#7C3AED',
    'fixed_baseline': '#9CA3AF',
}
colors = [
    cat_colors['adaptive'],    # Step
    cat_colors['adaptive'],    # Linear 1→0.01
    cat_colors['fixed_good'],  # c=0.1
    cat_colors['puct'],        # PUCT 0.5
    cat_colors['adaptive'],    # Linear 1→0.1
    cat_colors['adaptive'],    # Linear 0.5→0.1
    cat_colors['fixed_good'],  # c=0.01
    cat_colors['fixed_baseline'],  # c=0.5
    cat_colors['fixed_baseline'],  # c=1.0 original
    cat_colors['fixed_baseline'],  # c=2.0
    cat_colors['puct'],        # PUCT 2.0
]

fig, ax = plt.subplots(figsize=(16, 7))
bars = ax.barh(range(len(conditions)), selectivity, color=colors, edgecolor='white', height=0.7)

for i, (bar, val) in enumerate(zip(bars, selectivity)):
    ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
            f'+{val:.2f}', va='center', fontsize=10,
            fontweight='bold' if i == 0 else 'normal')

ax.set_yticks(range(len(conditions)))
ax.set_yticklabels(conditions, fontsize=9)
ax.set_xlabel('Selectivity (EGFR pIC50 − mean off-target pIC50)')
ax.set_title('All MCTS Policy Conditions — Selectivity Ranking\n(Higher = Better EGFR Selectivity)')
ax.invert_yaxis()
ax.axvline(x=0.06, color='gray', linestyle=':', alpha=0.7, linewidth=1.5)
ax.text(0.07, len(conditions)-0.5, '← Original\n    (c=1.0)', fontsize=8, color='gray')

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=cat_colors['adaptive'], label='Adaptive (ours)'),
    Patch(facecolor=cat_colors['puct'], label='PUCT (ours)'),
    Patch(facecolor=cat_colors['fixed_good'], label='Fixed c (tuned)'),
    Patch(facecolor=cat_colors['fixed_baseline'], label='Fixed c (baseline)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

plt.tight_layout()
plt.savefig(f'{outdir}/overview_selectivity_ranking.png')
plt.close()
print("✅ Overview saved")


# ════════════════════════════════════════════════════════════════
# 2. Reward Axis — Aggregation Method Comparison
# ════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Reward Axis 1.1: Aggregation Method Comparison',
             fontsize=16, fontweight='bold', y=1.02)

agg_methods = ['Geometric\n(Baseline)', 'Arithmetic', 'Tchebycheff']
agg_colors = ['#059669', '#2563EB', '#DC2626']

# Selectivity
vals = [0.06, -0.05, -0.30]
bars = axes[0].bar(agg_methods, vals, color=agg_colors, width=0.5, edgecolor='white')
axes[0].axhline(y=0, color='black', linewidth=0.8)
axes[0].set_title('Selectivity')
axes[0].set_ylabel('EGFR − Off-target')
for bar, v in zip(bars, vals):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01 * (1 if v >= 0 else -3),
                f'{v:+.2f}', ha='center', fontsize=11, fontweight='bold')

# EGFR pIC50
vals = [6.04, 6.04, 5.37]
axes[1].bar(agg_methods, vals, color=agg_colors, width=0.5, edgecolor='white')
axes[1].set_title('EGFR pIC50')
axes[1].set_ylim(5.0, 6.5)
for i, v in enumerate(vals):
    axes[1].text(i, v + 0.03, f'{v:.2f}', ha='center', fontsize=11)

# QED
vals = [0.775, 0.772, 0.743]
axes[2].bar(agg_methods, vals, color=agg_colors, width=0.5, edgecolor='white')
axes[2].set_title('QED')
axes[2].set_ylim(0.70, 0.80)
for i, v in enumerate(vals):
    axes[2].text(i, v + 0.002, f'{v:.3f}', ha='center', fontsize=11)

plt.tight_layout()
plt.savefig(f'{outdir}/reward_1_1_aggregation.png')
plt.close()
print("✅ Reward 1.1 saved")


# ════════════════════════════════════════════════════════════════
# 3. Reward Axis — σ and Weight Ablation
# ════════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle('Reward Axis: σ Sweep (1.2) & Weight Sweep (1.3)',
             fontsize=16, fontweight='bold', y=1.02)

# σ sweep
sigma_labels = ['EGFR strict\n(σ=1,2)', 'Baseline\n(σ=2,2)', 'EGFR lenient\n(σ=4,2)', 'OT strict\n(σ=2,1)']
sigma_sel = [0.09, 0.06, -0.21, -0.13]
sigma_reward = [0.169, 0.428, 0.576, 0.114]

x = np.arange(len(sigma_labels))
width = 0.35

bars1 = ax1.bar(x - width/2, sigma_sel, width, label='Selectivity', color='#2563EB')
ax1_twin = ax1.twinx()
bars2 = ax1_twin.bar(x + width/2, sigma_reward, width, label='Reward Signal', color='#D97706', alpha=0.7)

ax1.axhline(y=0, color='black', linewidth=0.8)
ax1.set_xticks(x)
ax1.set_xticklabels(sigma_labels, fontsize=9)
ax1.set_ylabel('Selectivity', color='#2563EB')
ax1_twin.set_ylabel('Reward Signal', color='#D97706')
ax1.set_title('σ Sweep — Selectivity vs Reward Signal')

for bar, v in zip(bars1, sigma_sel):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
             f'{v:+.2f}', ha='center', fontsize=9, fontweight='bold')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1_twin.get_legend_handles_labels()
ax1.legend(lines1 + [bars2], labels1 + ['Reward Signal'], fontsize=9, loc='upper left')

# Weight sweep
weight_labels = ['w=4', 'w=8\n(Baseline)', 'w=16']
weight_sel = [-0.03, 0.06, 0.06]
weight_qed = [0.781, 0.775, 0.759]

x2 = np.arange(len(weight_labels))
bars_s = ax2.bar(x2 - width/2, weight_sel, width, label='Selectivity', color='#2563EB')
ax2_twin = ax2.twinx()
bars_q = ax2_twin.bar(x2 + width/2, weight_qed, width, label='QED', color='#059669', alpha=0.7)

ax2.axhline(y=0, color='black', linewidth=0.8)
ax2.set_xticks(x2)
ax2.set_xticklabels(weight_labels, fontsize=10)
ax2.set_ylabel('Selectivity', color='#2563EB')
ax2_twin.set_ylabel('QED', color='#059669')
ax2_twin.set_ylim(0.74, 0.79)
ax2.set_title('EGFR Weight Sweep — Selectivity vs QED')

for bar, v in zip(bars_s, weight_sel):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
             f'{v:+.2f}', ha='center', fontsize=10, fontweight='bold')

lines1, labels1 = ax2.get_legend_handles_labels()
ax2.legend(lines1 + [bars_q], labels1 + ['QED'], fontsize=9, loc='upper left')

plt.tight_layout()
plt.savefig(f'{outdir}/reward_1_2_1_3_sigma_weight.png')
plt.close()
print("✅ Reward 1.2 & 1.3 saved")


# ════════════════════════════════════════════════════════════════
# 4. MCTS Policy — c_val Sweep (Selectivity vs Diversity)
# ════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(10, 6))

c_vals = [0.01, 0.1, 0.5, 1.0, 2.0]
sel = [0.17, 0.41, 0.13, 0.06, 0.03]
unique = [76.9, 91.8, 99.6, 99.8, 99.9]

color_sel = '#2563EB'
color_uni = '#DC2626'

ax.plot(c_vals, sel, 'o-', color=color_sel, linewidth=2.5, markersize=9, label='Selectivity', zorder=5)
for x, y in zip(c_vals, sel):
    ax.annotate(f'+{y:.2f}', (x, y), textcoords='offset points',
                xytext=(0, 12), ha='center', fontsize=9, color=color_sel, fontweight='bold')

ax_twin = ax.twinx()
ax_twin.plot(c_vals, unique, 's--', color=color_uni, linewidth=2, markersize=8, label='Unique%', alpha=0.8)
for x, y in zip(c_vals, unique):
    ax_twin.annotate(f'{y:.1f}%', (x, y), textcoords='offset points',
                     xytext=(0, -15), ha='center', fontsize=8, color=color_uni)

ax.set_xscale('log')
ax.set_xlabel('c_val (log scale)')
ax.set_ylabel('Selectivity', color=color_sel, fontsize=12)
ax_twin.set_ylabel('Unique %', color=color_uni, fontsize=12)
ax.set_title('MCTS Policy 2.1: c_val Sweep — Exploitation vs Exploration Trade-off',
             fontsize=14, fontweight='bold')
ax.set_xticks(c_vals)
ax.set_xticklabels([str(c) for c in c_vals])

ax.axvspan(0.08, 0.12, alpha=0.1, color='green')
ax.text(0.1, 0.35, '← Sweet spot', fontsize=9, color='green', ha='center', style='italic')

lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax_twin.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='center right', fontsize=10)

plt.tight_layout()
plt.savefig(f'{outdir}/mcts_2_1_cval_sweep.png')
plt.close()
print("✅ MCTS 2.1 saved")


# ════════════════════════════════════════════════════════════════
# 5. MCTS Policy — PUCT vs UCB1
# ════════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle('MCTS Policy 2.2: PUCT (RNN Prior) vs UCB1',
             fontsize=16, fontweight='bold', y=1.02)

policies = ['UCB1\nc=1.0', 'PUCT\ncp=0.5', 'PUCT\ncp=1.0', 'PUCT\ncp=2.0']
sel_puct = [0.06, 0.24, 0.22, 0.02]
egfr_puct = [6.04, 6.29, 6.23, 6.09]
policy_colors = ['#9CA3AF', '#7C3AED', '#A78BFA', '#C4B5FD']

bars = ax1.bar(policies, sel_puct, color=policy_colors, width=0.5, edgecolor='white')
ax1.axhline(y=0.06, color='gray', linestyle=':', alpha=0.7)
ax1.set_title('Selectivity')
ax1.set_ylabel('Selectivity')
for bar, v in zip(bars, sel_puct):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
             f'+{v:.2f}', ha='center', fontsize=11, fontweight='bold')

bars2 = ax2.bar(policies, egfr_puct, color=policy_colors, width=0.5, edgecolor='white')
ax2.set_title('EGFR pIC50')
ax2.set_ylim(5.8, 6.5)
for bar, v in zip(bars2, egfr_puct):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
             f'{v:.2f}', ha='center', fontsize=11)

# Annotation
ax1.annotate('4× improvement\nover baseline',
             xy=(1, 0.24), xytext=(2.5, 0.30),
             arrowprops=dict(arrowstyle='->', color='#7C3AED', lw=1.5),
             fontsize=10, color='#7C3AED', fontweight='bold')

plt.tight_layout()
plt.savefig(f'{outdir}/mcts_2_2_puct.png')
plt.close()
print("✅ MCTS 2.2 saved")


# ════════════════════════════════════════════════════════════════
# 6. MCTS Policy — Adaptive c_val Scheduling
# ════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(14, 7))

adaptive_labels = [
    'Step(1→0.1)',
    'Step(1→0.01)',
    'Linear(1→0.01)',
    'UCB1 c=0.1\n(fixed ref)',
    'Linear(1→0.1)',
    'Linear(0.5→0.1)',
    'UCB1 c=1.0\n(original)',
]
adaptive_sel = [0.64, 0.79, 0.56, 0.41, 0.23, 0.14, 0.06]
adaptive_qed = [0.778, 0.676, 0.749, 0.759, 0.781, 0.797, 0.775]
adaptive_uni = [93.7, 85.7, 97.4, 91.8, 99.0, 98.6, 99.8]

# Bubble chart: x=selectivity, y=QED, size=unique%
sizes = [(u/100)*300 for u in adaptive_uni]
scatter_colors = ['#DC2626', '#FF6B6B', '#F59E0B', '#9CA3AF', '#60A5FA', '#93C5FD', '#D1D5DB']

scatter = ax.scatter(adaptive_sel, adaptive_qed, s=sizes, c=scatter_colors,
                     alpha=0.8, edgecolors='#333', linewidth=1.5, zorder=5)

for i, label in enumerate(adaptive_labels):
    offset_x = 0.02 if adaptive_sel[i] < 0.5 else -0.02
    ha = 'left' if adaptive_sel[i] < 0.5 else 'right'
    ax.annotate(label, (adaptive_sel[i], adaptive_qed[i]),
                textcoords='offset points',
                xytext=(15 if ha == 'left' else -15, 8),
                fontsize=9, ha=ha,
                fontweight='bold' if i == 0 else 'normal')

# Highlight best (Step 1→0.1)
ax.annotate('★ BEST: Selectivity +0.64, QED 0.778',
            xy=(0.64, 0.778), xytext=(0.50, 0.80),
            arrowprops=dict(arrowstyle='->', color='#DC2626', lw=2),
            fontsize=11, color='#DC2626', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FEE2E2', edgecolor='#DC2626'))

ax.set_xlabel('Selectivity (higher = better EGFR selectivity)', fontsize=12)
ax.set_ylabel('QED (higher = more drug-like)', fontsize=12)
ax.set_title('MCTS Policy 2.3: Adaptive c_val — Selectivity vs Drug-likeness\n(Bubble size = Unique %)',
             fontsize=14, fontweight='bold')

# Pareto frontier hint
ax.axhline(y=0.75, color='gray', linestyle=':', alpha=0.3)
ax.axvline(x=0.40, color='gray', linestyle=':', alpha=0.3)
ax.text(0.60, 0.66, 'High Selectivity\nLow QED', fontsize=8, color='gray', ha='center', style='italic')
ax.text(0.10, 0.80, 'Low Selectivity\nHigh QED', fontsize=8, color='gray', ha='center', style='italic')

plt.tight_layout()
plt.savefig(f'{outdir}/mcts_2_3_adaptive.png')
plt.close()
print("✅ MCTS 2.3 saved")


# ════════════════════════════════════════════════════════════════
# 7. Cross & Long-run Results
# ════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(12, 6))

cross_labels = ['c=0.1\nbaseline\n(5K)', 'c=0.1\n+w=4\n(5K)', 'c=0.1\n+w=16\n(5K)',
                'c=0.1\nlong\n(30K)', 'c=0.1+σ=1\nlong\n(30K)']
cross_sel = [0.41, 0.96, 0.84, 1.02, 1.08]
cross_egfr = [6.35, 6.66, 6.83, 7.02, 7.48]
cross_qed = [0.759, 0.750, 0.658, 0.653, 0.666]

x = np.arange(len(cross_labels))
width = 0.25

bars1 = ax.bar(x - width, cross_sel, width, label='Selectivity', color='#2563EB')
bars2 = ax.bar(x, [e - 5 for e in cross_egfr], width, label='EGFR pIC50 (−5)', color='#059669')
bars3 = ax.bar(x + width, cross_qed, width, label='QED', color='#D97706')

ax.set_xticks(x)
ax.set_xticklabels(cross_labels, fontsize=9)
ax.set_ylabel('Score')
ax.set_title('Cross & Long-run Experiments — Exploitation Accumulation Effect',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=10)

for bars, vals in [(bars1, cross_sel), (bars3, cross_qed)]:
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{v:.2f}', ha='center', fontsize=8, fontweight='bold')
for bar, v in zip(bars2, cross_egfr):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f'{v:.2f}', ha='center', fontsize=8)

# Arrow showing accumulation
ax.annotate('Exploitation effect\naccumulates over time',
            xy=(3, 1.02), xytext=(1.5, 1.15),
            arrowprops=dict(arrowstyle='->', color='#2563EB', lw=1.5),
            fontsize=10, color='#2563EB', fontweight='bold')

plt.tight_layout()
plt.savefig(f'{outdir}/cross_longrun.png')
plt.close()
print("✅ Cross & Long-run saved")


# ════════════════════════════════════════════════════════════════
# 8. Reward Signal Collapse Illustration
# ════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('"Reward Signal Collapse" — When MCTS Loses Direction',
             fontsize=16, fontweight='bold', y=1.02)

# Simulated distributions
np.random.seed(42)

# Healthy signal
healthy = np.random.beta(2, 5, 500) * 0.8
axes[0].hist(healthy, bins=30, color='#059669', alpha=0.8, edgecolor='white')
axes[0].set_title('Healthy Signal\n(Geometric, σ=2)')
axes[0].set_xlabel('Reward')
axes[0].set_ylabel('Count')
axes[0].axvline(x=np.mean(healthy), color='#DC2626', linestyle='--', linewidth=2, label=f'Mean: {np.mean(healthy):.3f}')
axes[0].legend(fontsize=9)

# Collapsed — too narrow
collapsed = np.random.beta(50, 50, 500) * 0.05
axes[1].hist(collapsed, bins=30, color='#DC2626', alpha=0.8, edgecolor='white')
axes[1].set_title('Collapsed Signal\n(Tchebycheff / OT σ=1)')
axes[1].set_xlabel('Reward')
axes[1].axvline(x=np.mean(collapsed), color='#111', linestyle='--', linewidth=2, label=f'Mean: {np.mean(collapsed):.3f}')
axes[1].legend(fontsize=9)
axes[1].annotate('All molecules get\nsimilar low reward\n→ Random search!',
                 xy=(0.025, 60), fontsize=9, color='#DC2626',
                 fontweight='bold', ha='center',
                 bbox=dict(boxstyle='round', facecolor='#FEE2E2'))

# Too lenient
lenient = np.random.beta(5, 2, 500) * 0.8
axes[2].hist(lenient, bins=30, color='#D97706', alpha=0.8, edgecolor='white')
axes[2].set_title('Too Lenient\n(EGFR σ=4)')
axes[2].set_xlabel('Reward')
axes[2].axvline(x=np.mean(lenient), color='#111', linestyle='--', linewidth=2, label=f'Mean: {np.mean(lenient):.3f}')
axes[2].legend(fontsize=9)
axes[2].annotate('High reward but\nno selection pressure\n→ Low selectivity',
                 xy=(0.5, 60), fontsize=9, color='#D97706',
                 fontweight='bold', ha='center',
                 bbox=dict(boxstyle='round', facecolor='#FEF3C7'))

plt.tight_layout()
plt.savefig(f'{outdir}/reward_signal_collapse.png')
plt.close()
print("✅ Reward Signal Collapse saved")


print(f"\n📁 All P2 figures saved to {outdir}/")
