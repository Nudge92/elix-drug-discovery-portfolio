"""
P1. kMoL ADMET GNN Ablation Study — Visualization
Experiments 1~4: GNN layer comparison, depth/width, PLB encoding, ADME cross-comparison
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

COLORS = {
    'GIN': '#2563EB',
    'GCN': '#7C3AED',
    'GENConv': '#059669',
    'LEConv': '#DC2626',
    'GAT': '#D97706',
}

outdir = '/home/claude/kmol_figures'
os.makedirs(outdir, exist_ok=True)


# ════════════════════════════════════════════════════════════════
# Experiment 1: Solubility — 5 GNN Layers Bar Chart
# ════════════════════════════════════════════════════════════════

layers = ['GIN', 'GENConv', 'GAT', 'GCN', 'LEConv']
metrics = {
    'ROC-AUC':   [0.867, 0.866, 0.861, 0.851, 0.833],
    'PR-AUC':    [0.886, 0.898, 0.878, 0.854, 0.871],
    'Accuracy':  [0.767, 0.728, 0.777, 0.738, 0.728],
    'Precision': [0.720, 0.700, 0.728, 0.694, 0.758],
    'Recall':    [0.983, 0.933, 0.983, 0.983, 0.783],
    "Cohen's κ": [0.483, 0.402, 0.507, 0.413, 0.437],
}

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
fig.suptitle('Experiment 1: Solubility — GNN Layer Comparison (5 Architectures)',
             fontsize=16, fontweight='bold', y=0.98)

for ax, (metric_name, values) in zip(axes.flat, metrics.items()):
    colors = [COLORS[l] for l in layers]
    bars = ax.bar(layers, values, color=colors, width=0.6, edgecolor='white', linewidth=0.5)

    best_idx = np.argmax(values)
    bars[best_idx].set_edgecolor('#111')
    bars[best_idx].set_linewidth(2)

    for i, (bar, v) in enumerate(zip(bars, values)):
        weight = 'bold' if i == best_idx else 'normal'
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f'{v:.3f}', ha='center', va='bottom', fontsize=9, fontweight=weight)

    ax.set_title(metric_name)
    ax.set_ylim(min(values) - 0.08, max(values) + 0.05)
    ax.tick_params(axis='x', rotation=15)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(f'{outdir}/exp1_solubility_layers.png')
plt.close()
print("✅ Exp1 saved")


# ════════════════════════════════════════════════════════════════
# Experiment 2: GIN Depth & Width Line Charts
# ════════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle('Experiment 2: GIN — Depth & Width Ablation',
             fontsize=16, fontweight='bold', y=1.0)

# Depth
depth_x = [2, 4, 6, 8]
depth_metrics = {
    'ROC-AUC':   [0.867, 0.847, 0.810, 0.838],
    'PR-AUC':    [0.886, 0.882, 0.838, 0.875],
    'Accuracy':  [0.767, 0.786, 0.728, 0.777],
    "Cohen's κ": [0.483, 0.533, 0.422, 0.533],
}
depth_colors = ['#2563EB', '#7C3AED', '#059669', '#D97706']
markers = ['o', 's', '^', 'D']

for (name, vals), c, m in zip(depth_metrics.items(), depth_colors, markers):
    ax1.plot(depth_x, vals, marker=m, color=c, label=name, linewidth=2, markersize=7)

ax1.set_xlabel('layers_count')
ax1.set_ylabel('Score')
ax1.set_title('Depth (hidden=128 fixed)')
ax1.set_xticks(depth_x)
ax1.legend(fontsize=9, loc='lower left')
ax1.axvspan(1.5, 2.5, alpha=0.08, color='green', label='Best')

# Width
width_x = [64, 128, 256]
width_metrics = {
    'ROC-AUC':   [0.857, 0.867, 0.842],
    'PR-AUC':    [0.880, 0.886, 0.865],
    'Accuracy':  [0.767, 0.767, 0.718],
    "Cohen's κ": [0.491, 0.483, 0.445],
}

for (name, vals), c, m in zip(width_metrics.items(), width_colors if 'width_colors' in dir() else depth_colors, markers):
    ax2.plot(width_x, vals, marker=m, color=c, label=name, linewidth=2, markersize=7)

ax2.set_xlabel('hidden_features')
ax2.set_ylabel('Score')
ax2.set_title('Width (layers=2 fixed)')
ax2.set_xticks(width_x)
ax2.legend(fontsize=9, loc='lower left')
ax2.axvspan(100, 156, alpha=0.08, color='green')

plt.tight_layout()
plt.savefig(f'{outdir}/exp2_gin_depth_width.png')
plt.close()
print("✅ Exp2 saved")


# ════════════════════════════════════════════════════════════════
# Experiment 3: PLB Affinity — 2×2 Encoding Comparison
# ════════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle('Experiment 3: PLB Affinity — Ligand × Protein Encoding (ROC-AUC)',
             fontsize=16, fontweight='bold', y=1.0)

# 2×2 heatmap
encoding_data = np.array([
    [0.749, 0.721],  # FP row
    [0.698, 0.705],  # Graph row
])

im = ax1.imshow(encoding_data, cmap='YlOrRd', aspect='auto', vmin=0.68, vmax=0.76)
ax1.set_xticks([0, 1])
ax1.set_xticklabels(['BoW', 'CNN'], fontsize=12)
ax1.set_yticks([0, 1])
ax1.set_yticklabels(['Fingerprint', 'Graph'], fontsize=12)
ax1.set_xlabel('Protein Encoding', fontsize=12)
ax1.set_ylabel('Ligand Encoding', fontsize=12)
ax1.set_title('ROC-AUC Heatmap')

for i in range(2):
    for j in range(2):
        val = encoding_data[i, j]
        best = (i == 0 and j == 0)
        ax1.text(j, i, f'{val:.3f}', ha='center', va='center',
                fontsize=16, fontweight='bold' if best else 'normal',
                color='white' if val > 0.73 else 'black')

fig.colorbar(im, ax=ax1, shrink=0.8)

# Bar chart comparison (all metrics)
combos = ['FP+BoW', 'FP+CNN', 'Graph+CNN', 'Graph+BoW']
roc_auc = [0.749, 0.721, 0.705, 0.698]
pr_auc =  [0.807, 0.798, 0.781, 0.778]
accuracy = [0.763, 0.750, 0.731, 0.735]

x = np.arange(len(combos))
w = 0.25

bars1 = ax2.bar(x - w, roc_auc, w, label='ROC-AUC', color='#2563EB')
bars2 = ax2.bar(x, pr_auc, w, label='PR-AUC', color='#7C3AED')
bars3 = ax2.bar(x + w, accuracy, w, label='Accuracy', color='#059669')

ax2.set_xticks(x)
ax2.set_xticklabels(combos, fontsize=10)
ax2.set_ylabel('Score')
ax2.set_title('Multi-Metric Comparison')
ax2.legend(fontsize=9)
ax2.set_ylim(0.65, 0.85)

for bars in [bars1, bars2, bars3]:
    for bar in bars:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.003,
                f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=7.5)

plt.tight_layout()
plt.savefig(f'{outdir}/exp3_plb_encoding.png')
plt.close()
print("✅ Exp3 saved")


# ════════════════════════════════════════════════════════════════
# Experiment 4: ADME 4 Tasks × 5 Layers — Heatmap
# ════════════════════════════════════════════════════════════════

tasks = ['Solubility\n(ROC-AUC)', 'Clint\n(R²)', 'Papp Caco-2\n(R²)', 'Fup Human\n(R²)']
layer_names = ['GIN', 'GCN', 'GENConv', 'LEConv', 'GAT']

heatmap_data = np.array([
    [0.867, 0.485, 0.477, 0.626],  # GIN
    [0.851, 0.494, 0.437, 0.644],  # GCN
    [0.866, 0.504, 0.470, 0.616],  # GENConv
    [0.833, 0.519, 0.493, 0.638],  # LEConv
    [0.861, 0.492, 0.467, 0.649],  # GAT
])

fig, ax = plt.subplots(figsize=(10, 6))
fig.suptitle('Experiment 4: ADME 4 Tasks × 5 GNN Layers',
             fontsize=16, fontweight='bold', y=0.98)

im = ax.imshow(heatmap_data, cmap='RdYlGn', aspect='auto', vmin=0.42, vmax=0.87)

ax.set_xticks(range(len(tasks)))
ax.set_xticklabels(tasks, fontsize=11)
ax.set_yticks(range(len(layer_names)))
ax.set_yticklabels(layer_names, fontsize=12)

# Annotate + highlight best per column
for j in range(len(tasks)):
    col = heatmap_data[:, j]
    best_i = np.argmax(col)
    for i in range(len(layer_names)):
        val = heatmap_data[i, j]
        is_best = (i == best_i)
        ax.text(j, i, f'{val:.3f}', ha='center', va='center',
                fontsize=12, fontweight='bold' if is_best else 'normal',
                color='white' if val < 0.50 else 'black')
        if is_best:
            ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1,
                         fill=False, edgecolor='#111', linewidth=2.5))

fig.colorbar(im, ax=ax, shrink=0.8, label='Score')

# Add subtitle
ax.set_xlabel('\n★ Best layer varies by task: Classification → GIN, Regression → LEConv, Protein binding → GAT',
              fontsize=10, style='italic')

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(f'{outdir}/exp4_adme_heatmap.png')
plt.close()
print("✅ Exp4 saved")

print(f"\n📁 All P1 figures saved to {outdir}/")
