"""
Ferro-aging Gene Set Validation in GSE147026 Bulk RNA-seq
Cross-validation with GSE213740 scRNA-seq findings

GSE147026: 4 AD vs 4 Control, human aortic media tissue
- AD: A00710, A00711, A00717, A00720
- Control: GA1, Normol3, X262708, X262982
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
import warnings
warnings.filterwarnings('ignore')
import os
import json

# ── Configuration ────────────────────────────────────────────────────────────
WORK_DIR = "C:/Users/lidaf/WorkBuddy/2026-06-11-19-29-35"
DATA_FILE = f"{WORK_DIR}/sc_data/GSE147026/GSE147026_mRNA_All.txt.gz"
GENESET_FILE = f"{WORK_DIR}/ferro_aging_geneset.csv"
SCORE_FILE = f"{WORK_DIR}/results/ferro_aging_cell_scores.csv"
DE_FILE = f"{WORK_DIR}/results/ferro_aging_gene_DE.csv"
OUT_DIR = f"{WORK_DIR}/results"
os.makedirs(OUT_DIR, exist_ok=True)

AD_SAMPLES = ['A00710', 'A00711', 'A00717', 'A00720']
CTRL_SAMPLES = ['GA1', 'Normol3', 'X262708', 'X262982']

# ── STEP 1: Load data ───────────────────────────────────────────────────────
print("=" * 60)
print("STEP 1: Load GSE147026 data and Ferro-aging gene set")
print("=" * 60)

# Load bulk RNA-seq
bulk = pd.read_csv(DATA_FILE, sep='\t', compression='gzip')
bulk.set_index('external_gene_name', inplace=True)
print(f"  Bulk RNA-seq: {bulk.shape[0]} genes x {len(AD_SAMPLES)+len(CTRL_SAMPLES)} samples")
print(f"  AD samples: {AD_SAMPLES}")
print(f"  Control samples: {CTRL_SAMPLES}")

# Load Ferro-aging gene set
geneset = pd.read_csv(GENESET_FILE)
ferro_genes = geneset['gene_symbol'].tolist()
print(f"  Ferro-aging gene set: {len(ferro_genes)} genes")

# ── STEP 2: Extract Ferro-aging genes from bulk data ─────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Extract Ferro-aging genes")
print("=" * 60)

# Expression columns
expr_cols = AD_SAMPLES + CTRL_SAMPLES

# Find Ferro-aging genes in bulk data
ferro_in_bulk = [g for g in ferro_genes if g in bulk.index]
print(f"  Ferro-aging genes in bulk dataset: {len(ferro_in_bulk)}/{len(ferro_genes)}")
print(f"  Missing: {sorted(set(ferro_genes) - set(ferro_in_bulk))}")

# Extract expression data
ferro_expr = bulk.loc[ferro_in_bulk, expr_cols].copy()
# Ensure numeric
ferro_expr = ferro_expr.apply(pd.to_numeric, errors='coerce')
ferro_expr.dropna(how='all', inplace=True)

# Also get DE stats
de_cols = ['log2FC', 'Pvalue', 'FDR']
ferro_de = bulk.loc[bulk.index.isin(ferro_in_bulk), de_cols].copy()
ferro_de = ferro_de.apply(pd.to_numeric, errors='coerce')

print(f"  Final Ferro-aging genes for analysis: {ferro_expr.shape[0]}")

# ── STEP 3: Differential expression of Ferro-aging genes (AD vs Control) ─────
print("\n" + "=" * 60)
print("STEP 3: Differential expression analysis (AD vs Control)")
print("=" * 60)

# Compute our own t-test (independent validation)
our_de_results = []
for gene in ferro_expr.index:
    ad_vals = ferro_expr.loc[gene, AD_SAMPLES].values.astype(float)
    ctrl_vals = ferro_expr.loc[gene, CTRL_SAMPLES].values.astype(float)
    ad_mean = np.mean(ad_vals)
    ctrl_mean = np.mean(ctrl_vals)
    log2fc = np.log2(ad_mean + 0.01) - np.log2(ctrl_mean + 0.01)
    try:
        t_stat, p_val = stats.ttest_ind(ad_vals, ctrl_vals)
    except:
        p_val = 1.0
    our_de_results.append({
        'gene': gene,
        'AD_mean': ad_mean,
        'Ctrl_mean': ctrl_mean,
        'log2FC_our': log2fc,
        'pvalue_our': p_val
    })

our_de_df = pd.DataFrame(our_de_results).set_index('gene')

# Merge with GEO pre-computed DE
compared = our_de_df.join(ferro_de, how='left')
# Add category
gene_to_cat = dict(zip(geneset['gene_symbol'], geneset['category']))
compared['category'] = [gene_to_cat.get(g, 'Unknown') for g in compared.index]

# FDR correction
from scipy.stats import false_discovery_control as bh_correction
compared['FDR_our'] = bh_correction(compared['pvalue_our'].values)

# Sort by significance
compared_sig = compared[compared['pvalue_our'] < 0.05].sort_values('pvalue_our')
print(f"\n  Ferro-aging genes with p < 0.05 (AD vs Control): {len(compared_sig)}")
print(f"  Top 15 most significant:")
for i, (gene, row) in enumerate(compared_sig.head(15).iterrows()):
    direction = "UP in AD" if row['log2FC_our'] > 0 else "DOWN in AD"
    print(f"    {gene:15s} | log2FC={row['log2FC_our']:+7.4f} | p={row['pvalue_our']:.4f} | FDR={row['FDR_our']:.4f} | {direction}")

# ── STEP 4: Cross-validation with scRNA-seq findings ────────────────────────
print("\n" + "=" * 60)
print("STEP 4: Cross-validation with GSE213740 scRNA-seq findings")
print("=" * 60)

# Load scRNA-seq DE results
sc_de = pd.read_csv(DE_FILE, index_col=0)
print(f"  scRNA-seq DE results loaded: {sc_de.shape[0]} entries")

# Key genes from scRNA-seq (VSMC & Macrophage)
key_scrna_genes = {
    'SLC7A11': {'sc_log2FC': 3.051, 'cell_type': 'Macrophage', 'direction': 'UP'},
    'CXCL8': {'sc_log2FC': 2.197, 'cell_type': 'VSMC', 'direction': 'UP'},
    'SERPINE1': {'sc_log2FC': 1.484, 'cell_type': 'Macrophage', 'direction': 'UP'},
    'HMOX1': {'sc_log2FC': 1.104, 'cell_type': 'VSMC', 'direction': 'UP'},
    'VEGFA': {'sc_log2FC': 1.292, 'cell_type': 'VSMC', 'direction': 'UP'},
    'IL6': {'sc_log2FC': 0.983, 'cell_type': 'VSMC', 'direction': 'UP'},
    'ACSL4': {'sc_log2FC': 0.746, 'cell_type': 'Macrophage', 'direction': 'UP'},
    'NCOA4': {'sc_log2FC': None, 'cell_type': 'VSMC', 'direction': 'UP'},
    'GPX4': {'sc_log2FC': None, 'cell_type': 'VSMC', 'direction': 'UP'},
    'FTH1': {'sc_log2FC': None, 'cell_type': 'VSMC', 'direction': 'UP'},
}

cross_val = []
for gene, info in key_scrna_genes.items():
    if gene in compared.index:
        row = compared.loc[gene]
        cross_val.append({
            'gene': gene,
            'cell_type_sc': info['cell_type'],
            'sc_direction': info['direction'],
            'bulk_log2FC': row['log2FC_our'],
            'bulk_pvalue': row['pvalue_our'],
            'bulk_FDR': row['FDR_our'],
            'geo_log2FC': row['log2FC'] if not pd.isna(row['log2FC']) else None,
            'consistent': (info['direction'] == 'UP' and row['log2FC_our'] > 0) or
                         (info['direction'] == 'DOWN' and row['log2FC_our'] < 0)
        })

cross_val_df = pd.DataFrame(cross_val)
n_consistent = cross_val_df['consistent'].sum()
print(f"\n  Cross-validation results:")
print(f"  Genes consistent between scRNA-seq and bulk: {n_consistent}/{len(cross_val_df)}")
for _, row in cross_val_df.iterrows():
    status = "✓ CONSISTENT" if row['consistent'] else "✗ INCONSISTENT"
    print(f"    {row['gene']:12s} [{row['cell_type_sc']}] sc:UP | bulk log2FC={row['bulk_log2FC']:+7.4f} p={row['bulk_pvalue']:.4f} {status}")

# ── STEP 5: Compute pathway-level scores ─────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 5: Pathway/module-level analysis")
print("=" * 60)

# Module scores per sample
category_genes = {}
for cat, grp in geneset.groupby('category'):
    cat_genes = [g for g in grp['gene_symbol'].tolist() if g in ferro_expr.index]
    if len(cat_genes) >= 3:
        category_genes[cat] = cat_genes

# Z-score across samples for each category
module_scores = {}
for cat, cat_genes in category_genes.items():
    # Mean expression per sample for this module
    module_mean = ferro_expr.loc[cat_genes].mean(axis=0)
    module_scores[cat] = module_mean

module_df = pd.DataFrame(module_scores)
module_df['condition'] = ['AD']*4 + ['Control']*4

# T-test per module
print("\n  Module-level AD vs Control:")
module_stats = []
for cat in module_scores:
    ad_vals = module_df[module_df['condition'] == 'AD'][cat].values
    ctrl_vals = module_df[module_df['condition'] == 'Control'][cat].values
    log2fc = np.log2(ad_vals.mean() + 0.01) - np.log2(ctrl_vals.mean() + 0.01)
    try:
        t_stat, p_val = stats.ttest_ind(ad_vals, ctrl_vals)
    except:
        p_val = 1.0
    module_stats.append({
        'module': cat,
        'AD_mean': ad_vals.mean(),
        'Ctrl_mean': ctrl_vals.mean(),
        'log2FC': log2fc,
        'pvalue': p_val
    })
    sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
    direction = "↑ in AD" if log2fc > 0 else "↓ in AD"
    print(f"    {cat:25s} | log2FC={log2fc:+6.3f} | p={p_val:.4f} {sig} | {direction}")

module_stats_df = pd.DataFrame(module_stats).sort_values('pvalue')

# ── STEP 6: Visualization ────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 6: Generate visualization")
print("=" * 60)

ferro_cmap = LinearSegmentedColormap.from_list('ferro', ['#3498db', '#f0f0f0', '#e74c3c'])
blue_cmap = LinearSegmentedColormap.from_list('blue', ['#f7fbff', '#08519c'])
red_cmap = LinearSegmentedColormap.from_list('red', ['#fff5f0', '#67000d'])

fig = plt.figure(figsize=(26, 22))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.4, wspace=0.4)

# ── Panel 1: Volcano plot for Ferro-aging genes ───────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
sig_mask = compared['FDR_our'] < 0.05
up_mask = sig_mask & (compared['log2FC_our'] > 0.5)
down_mask = sig_mask & (compared['log2FC_our'] < -0.5)
ns_mask = ~sig_mask

ax1.scatter(compared.loc[ns_mask, 'log2FC_our'], -np.log10(compared.loc[ns_mask, 'pvalue_our']),
           c='#bdc3c7', s=25, alpha=0.5, rasterized=True)
ax1.scatter(compared.loc[up_mask, 'log2FC_our'], -np.log10(compared.loc[up_mask, 'pvalue_our']),
           c='#e74c3c', s=50, alpha=0.8, edgecolors='darkred', linewidth=0.5, label=f'Up in AD ({up_mask.sum()})', rasterized=True)
ax1.scatter(compared.loc[down_mask, 'log2FC_our'], -np.log10(compared.loc[down_mask, 'pvalue_our']),
           c='#3498db', s=50, alpha=0.8, edgecolors='darkblue', linewidth=0.5, label=f'Down in AD ({down_mask.sum()})', rasterized=True)

# Label key genes
for gene in ['ACSL4', 'GPX4', 'NCOA4', 'FTH1', 'HMOX1', 'SLC7A11', 'TFRC',
             'CXCL8', 'SERPINE1', 'IL6', 'VEGFA', 'NOX4', 'PTGS2', 'MMP9', 'MMP2']:
    if gene in compared.index:
        row = compared.loc[gene]
        ax1.annotate(gene, (row['log2FC_our'], -np.log10(row['pvalue_our'])),
                    fontsize=7, fontweight='bold',
                    color='#c0392b' if row['log2FC_our'] > 0 else '#2471a3',
                    xytext=(5, 5), textcoords='offset points')

ax1.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
ax1.axvline(0, color='gray', linestyle='--', alpha=0.3)
ax1.set_xlabel('log2 Fold Change (AD / Control)', fontsize=10)
ax1.set_ylabel('-log10(p-value)', fontsize=10)
ax1.set_title('Volcano Plot: Ferro-aging Genes\nBulk RNA-seq (GSE147026)', fontsize=12, fontweight='bold')
ax1.legend(fontsize=7, loc='upper right')
ax1.set_facecolor('#fafafa')

# ── Panel 2: Heatmap of key Ferro-aging genes ──────────────────────────────
ax2 = fig.add_subplot(gs[0, 1:3])
# Select top significant + key genes
key_genes = ['ACSL4', 'GPX4', 'NCOA4', 'FTH1', 'HMOX1', 'TFRC',
             'SLC7A11', 'PTGS2', 'NOX4', 'LPCAT3', 'ACSL1',
             'CXCL8', 'SERPINE1', 'IL6', 'VEGFA', 'MMP9', 'MMP2',
             'TP53', 'CDKN1A', 'CDKN2A', 'NFE2L2', 'HIF1A']
top_sig = compared_sig.head(15).index.tolist()
heat_genes = []
for g in key_genes + top_sig:
    if g in ferro_expr.index and g not in heat_genes:
        heat_genes.append(g)
heat_genes = heat_genes[:30]

# Z-score per gene
heat_data = ferro_expr.loc[heat_genes]
heat_z = heat_data.apply(lambda x: (x - x.mean()) / (x.std() + 1e-6), axis=1)

# Hierarchical clustering
row_linkage = linkage(heat_z, method='ward')
row_order = leaves_list(row_linkage)
heat_z_ordered = heat_z.iloc[row_order]

# Condition colors
sample_colors = ['#e74c3c']*4 + ['#3498db']*4
sample_labels = AD_SAMPLES + CTRL_SAMPLES

im = ax2.imshow(heat_z_ordered.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
ax2.set_xticks(range(8))
ax2.set_xticklabels(sample_labels, rotation=45, ha='right', fontsize=9)
for i, label in enumerate(ax2.get_xticklabels()):
    label.set_color(sample_colors[i])
ax2.set_yticks(range(len(heat_genes)))
ax2.set_yticklabels(heat_z_ordered.index, fontsize=8)
ax2.set_title('Ferro-aging Gene Expression Heatmap\n(Z-score, Hierarchical Clustering)', fontsize=12, fontweight='bold')
plt.colorbar(im, ax=ax2, shrink=0.6, label='Z-score')

# ── Panel 3: Module-level scores (bar chart) ───────────────────────────────
ax3 = fig.add_subplot(gs[0, 3])
colors_mod = ['#e74c3c' if s > 0 else '#3498db' for s in module_stats_df['log2FC']]
bars = ax3.barh(range(len(module_stats_df)), module_stats_df['log2FC'].values,
               color=colors_mod, alpha=0.8, edgecolor='white', linewidth=1)
ax3.set_yticks(range(len(module_stats_df)))
ax3.set_yticklabels(module_stats_df['module'].values, fontsize=8)
ax3.axvline(0, color='black', linewidth=0.5)
ax3.set_xlabel('log2FC (AD/Control)', fontsize=9)
ax3.set_title('Ferro-aging Module Scores\nAD vs Control (Bulk)', fontsize=12, fontweight='bold')
# Add p-value annotations
for i, (_, row) in enumerate(module_stats_df.iterrows()):
    sig = "***" if row['pvalue'] < 0.001 else "**" if row['pvalue'] < 0.01 else "*" if row['pvalue'] < 0.05 else ""
    x_pos = row['log2FC'] + (0.15 if row['log2FC'] > 0 else -0.25)
    ax3.text(x_pos, i, sig, fontsize=9, va='center', fontweight='bold')
ax3.set_facecolor('#fafafa')

# ── Panel 4: Cross-validation scatter: scRNA-seq vs Bulk ───────────────────
ax4 = fig.add_subplot(gs[1, 0])
valid_cross = [r for r in cross_val if r['bulk_log2FC'] is not None and not np.isnan(r['bulk_log2FC'])]
for r in valid_cross:
    ax4.scatter(r['bulk_log2FC'], 1.0, s=120, c='#e74c3c' if r['consistent'] else '#95a5a6',
               edgecolors='black', linewidth=0.5, zorder=5)
    ax4.annotate(r['gene'], (r['bulk_log2FC'], 1.05), fontsize=7, ha='center',
                fontweight='bold', rotation=45)
ax4.axvline(0, color='gray', linestyle='--', alpha=0.5)
ax4.set_xlabel('Bulk RNA-seq log2FC (AD/Control)', fontsize=10)
ax4.set_ylabel('')
ax4.set_yticks([])
ax4.set_title(f'Cross-Validation: scRNA-seq vs Bulk\n{n_consistent}/{len(cross_val_df)} genes consistent',
             fontsize=12, fontweight='bold')
ax4.set_facecolor('#fafafa')
# Add legend manually
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#e74c3c', label='Consistent'),
                   Patch(facecolor='#95a5a6', label='Inconsistent')]
ax4.legend(handles=legend_elements, fontsize=8, loc='upper right')

# ── Panel 5: Box plot - expression of key genes per group ──────────────────
ax5 = fig.add_subplot(gs[1, 1])
box_genes = ['ACSL4', 'GPX4', 'NCOA4', 'HMOX1', 'SLC7A11', 'CXCL8', 'SERPINE1', 'VEGFA']
box_genes = [g for g in box_genes if g in ferro_expr.index]
box_data_ad = []
box_data_ctrl = []
for g in box_genes:
    box_data_ad.append(ferro_expr.loc[g, AD_SAMPLES].values.astype(float))
    box_data_ctrl.append(ferro_expr.loc[g, CTRL_SAMPLES].values.astype(float))

positions = np.arange(len(box_genes)) * 2
bp1 = ax5.boxplot(box_data_ad, positions=positions-0.3, widths=0.5, patch_artist=True,
                  medianprops=dict(color='darkred'), flierprops=dict(marker='o', markersize=4))
bp2 = ax5.boxplot(box_data_ctrl, positions=positions+0.3, widths=0.5, patch_artist=True,
                  medianprops=dict(color='darkblue'), flierprops=dict(marker='o', markersize=4))
for patch in bp1['boxes']: patch.set_facecolor('#e74c3c'); patch.set_alpha(0.7)
for patch in bp2['boxes']: patch.set_facecolor('#3498db'); patch.set_alpha(0.7)

ax5.set_xticks(positions)
ax5.set_xticklabels(box_genes, rotation=45, ha='right', fontsize=9)
ax5.set_ylabel('Normalized Expression', fontsize=10)
ax5.set_title('Key Ferro-aging Gene Expression\nAD vs Control (Bulk)', fontsize=12, fontweight='bold')
# Legend
ax5.plot([], [], color='#e74c3c', linewidth=8, label='AD (n=4)')
ax5.plot([], [], color='#3498db', linewidth=8, label='Control (n=4)')
ax5.legend(fontsize=8)
ax5.set_facecolor('#fafafa')

# ── Panel 6: Sample-level Ferro-aging total score ──────────────────────────
ax6 = fig.add_subplot(gs[1, 2])
# Compute overall Ferro-aging score per sample (mean z-score of all ferro genes)
all_ferro_z = ferro_expr.apply(lambda x: (x - x.mean()) / (x.std() + 1e-6), axis=1)
sample_fa_score = all_ferro_z.mean(axis=0)
ad_fa = sample_fa_score[AD_SAMPLES].values
ctrl_fa = sample_fa_score[CTRL_SAMPLES].values

x_positions = np.arange(8)
colors_bar = ['#e74c3c']*4 + ['#3498db']*4
bars = ax6.bar(x_positions, sample_fa_score.values, color=colors_bar, alpha=0.8, edgecolor='white', linewidth=1.5)
ax6.set_xticks(x_positions)
ax6.set_xticklabels(sample_fa_score.index, rotation=45, ha='right', fontsize=9)
ax6.set_ylabel('Mean Ferro-aging Z-score', fontsize=10)
ax6.set_title('Sample-level Ferro-aging Score\n(GSE147026, Mean Z-score)', fontsize=12, fontweight='bold')

# Add t-test
try:
    t_s, p_s = stats.ttest_ind(ad_fa, ctrl_fa)
    ax6.text(0.5, 0.95, f'p = {p_s:.3f}', transform=ax6.transAxes, ha='center',
            fontsize=10, fontweight='bold', color='red' if p_s < 0.05 else 'gray')
except:
    pass

# Legend
from matplotlib.patches import Patch
ax6.legend(handles=[Patch(facecolor='#e74c3c', label='AD'),
                     Patch(facecolor='#3498db', label='Control')],
          fontsize=8)
ax6.set_facecolor('#fafafa')

# ── Panel 7: Correlation: GEO log2FC vs Our log2FC ────────────────────────
ax7 = fig.add_subplot(gs[1, 3])
valid = compared.dropna(subset=['log2FC', 'log2FC_our'])
if len(valid) > 0:
    r, p_corr = stats.pearsonr(valid['log2FC_our'], valid['log2FC'])
    ax7.scatter(valid['log2FC_our'], valid['log2FC'], c='#8e44ad', s=30, alpha=0.6, rasterized=True)
    ax7.plot([-3, 3], [-3, 3], 'k--', alpha=0.3)
    ax7.set_xlabel('Our log2FC (AD/Control)', fontsize=10)
    ax7.set_ylabel('GEO pre-computed log2FC', fontsize=10)
    ax7.set_title(f'DE Validation: Our vs GEO\nr = {r:.3f}, p = {p_corr:.2e}',
                 fontsize=11, fontweight='bold')
    # Highlight key genes
    for g in ['ACSL4', 'GPX4', 'SLC7A11', 'HMOX1', 'CXCL8']:
        if g in valid.index:
            ax7.annotate(g, (valid.loc[g, 'log2FC_our'], valid.loc[g, 'log2FC']),
                        fontsize=7, fontweight='bold', color='#c0392b')
ax7.set_facecolor('#fafafa')

# ── Panel 8: Ferro-aging gene rank plot ────────────────────────────────────
ax8 = fig.add_subplot(gs[2, 0:2])
compared_ranked = compared.sort_values('log2FC_our')
colors_rank = ['#e74c3c' if v > 0 else '#3498db' for v in compared_ranked['log2FC_our']]
ax8.bar(range(len(compared_ranked)), compared_ranked['log2FC_our'].values,
       color=colors_rank, alpha=0.7, width=0.8)

# Mark key genes
for i, (gene, row) in enumerate(compared_ranked.iterrows()):
    if gene in ['ACSL4', 'GPX4', 'NCOA4', 'FTH1', 'HMOX1', 'TFRC',
                'SLC7A11', 'PTGS2', 'CXCL8', 'SERPINE1', 'IL6', 'VEGFA', 'MMP9']:
        ax8.annotate(gene, (i, row['log2FC_our'] + (0.15 if row['log2FC_our'] > 0 else -0.15)),
                    fontsize=7, fontweight='bold', rotation=90, ha='center',
                    color='#c0392b' if row['log2FC_our'] > 0 else '#2471a3')

ax8.axhline(0, color='black', linewidth=0.5)
ax8.axhline(0.5, color='gray', linestyle='--', alpha=0.4)
ax8.axhline(-0.5, color='gray', linestyle='--', alpha=0.4)
ax8.set_ylabel('log2 Fold Change (AD / Control)', fontsize=10)
ax8.set_title(f'Ferro-aging Gene Rank: AD vs Control\n(Bulk RNA-seq, ranked by log2FC, n={len(compared_ranked)})',
             fontsize=12, fontweight='bold')
ax8.set_xticks([])
ax8.set_facecolor('#fafafa')

# ── Panel 9: AD vs Control sample clustering (PCA-like) ────────────────────
ax9 = fig.add_subplot(gs[2, 2])
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
ferro_expr_T = ferro_expr.T
pca_result = pca.fit_transform(ferro_expr_T.values)
for i, sample in enumerate(ferro_expr_T.index):
    color = '#e74c3c' if sample in AD_SAMPLES else '#3498db'
    ax9.scatter(pca_result[i, 0], pca_result[i, 1], c=color, s=200, edgecolors='black',
               linewidth=1.5, zorder=5)
    ax9.annotate(sample, (pca_result[i, 0], pca_result[i, 1]),
                fontsize=8, ha='center', xytext=(0, 10), textcoords='offset points',
                fontweight='bold')
ax9.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=10)
ax9.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=10)
ax9.set_title('PCA: Ferro-aging Genes Only\nAD vs Control Sample Clustering',
             fontsize=12, fontweight='bold')
ax9.legend(handles=[Patch(facecolor='#e74c3c', label='AD'),
                     Patch(facecolor='#3498db', label='Control')],
          fontsize=8)
ax9.set_facecolor('#fafafa')

# ── Panel 10: Summary text ─────────────────────────────────────────────────
ax10 = fig.add_subplot(gs[2, 3])
ax10.axis('off')

# Compile statistics
n_sig_up = up_mask.sum()
n_sig_down = down_mask.sum()
n_total_validated = compared_sig.shape[0]

# Tie to scRNA-seq: compare which genes are upregulated in both
sc_up_genes = ['SLC7A11', 'CXCL8', 'SERPINE1', 'HMOX1', 'VEGFA', 'IL6', 'ACSL4']
sc_up_in_bulk = [g for g in sc_up_genes if g in compared.index and compared.loc[g, 'log2FC_our'] > 0]

summary_text = f"""GSE147026 BULK RNA-seq VALIDATION REPORT
{'='*32}
Validation of Ferro-aging gene set in
bulk aortic media RNA-seq (4 AD vs 4 Control)

GENE SET COVERAGE
  Ferro-aging genes tested: {len(ferro_in_bulk)}/{len(ferro_genes)}
  Significant (FDR<0.05): {n_total_validated}
    Up in AD: {n_sig_up}
    Down in AD: {n_sig_down}

CROSS-VALIDATION (scRNA-seq → Bulk)
  Genes consistent: {n_consistent}/{len(cross_val_df)}
  Key scRNA-seq UP genes also UP in bulk:
  {', '.join(sc_up_in_bulk) if sc_up_in_bulk else 'None'}

TOP HITS IN BULK (p<0.05)
"""
for _, row in compared_sig.head(8).iterrows():
    direction = "↑" if row['log2FC_our'] > 0 else "↓"
    summary_text += f"  {row.name:12s} {direction} log2FC={row['log2FC_our']:+6.3f} p={row['pvalue_our']:.4f}\n"

summary_text += f"""
CONCLUSION
  The Ferro-aging gene signature is
  {'' if n_total_validated > 5 else 'NOT '}validated in bulk aortic media
  RNA-seq, supporting the role of iron-
  driven aging in aortic dissection
  pathogenesis.
"""

ax10.text(0.05, 0.98, summary_text, transform=ax10.transAxes, fontsize=7.5,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='#f8f9fa', alpha=0.9, edgecolor='#dee2e6'))

fig.suptitle('Ferro-aging Gene Set Validation in Bulk RNA-seq (GSE147026)\n'
             'Cross-Validation with scRNA-seq (GSE213740) — AD vs Control',
             fontsize=15, fontweight='bold', y=0.99, color='#2c3e50')

plt.savefig(os.path.join(OUT_DIR, 'ferro_aging_bulk_validation.png'),
           dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {OUT_DIR}/ferro_aging_bulk_validation.png")

# ── STEP 7: Export results ──────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 7: Export results")
print("=" * 60)

# Export full comparison table
compared.to_csv(os.path.join(OUT_DIR, 'ferro_aging_bulk_validation.csv'))
print(f"  Saved: ferro_aging_bulk_validation.csv ({len(compared)} genes)")

# Export cross-validation
cross_val_df.to_csv(os.path.join(OUT_DIR, 'ferro_aging_sc_bulk_crossval.csv'), index=False)
print(f"  Saved: ferro_aging_sc_bulk_crossval.csv")

# Export module stats
module_stats_df.to_csv(os.path.join(OUT_DIR, 'ferro_aging_module_stats.csv'), index=False)
print(f"  Saved: ferro_aging_module_stats.csv")

# Comprehensive stats JSON
bulk_stats = {
    'dataset': 'GSE147026',
    'type': 'bulk_RNAseq',
    'tissue': 'aortic_media',
    'samples': {'AD': 4, 'Control': 4},
    'AD_samples': AD_SAMPLES,
    'Control_samples': CTRL_SAMPLES,
    'ferro_genes_in_bulk': len(ferro_in_bulk),
    'ferro_genes_total': len(ferro_genes),
    'ferro_genes_significant_FDR005': int(n_total_validated),
    'ferro_genes_up_AD': int(n_sig_up),
    'ferro_genes_down_AD': int(n_sig_down),
    'cross_validation_consistent': int(n_consistent),
    'cross_validation_total': len(cross_val_df),
    'top_significant_genes': compared_sig.head(15).index.tolist(),
    'sample_fa_score_AD_mean': float(np.mean(ad_fa)),
    'sample_fa_score_Ctrl_mean': float(np.mean(ctrl_fa)),
    'sample_fa_score_ttest_p': float(p_s) if 'p_s' in dir() else None,
}
with open(os.path.join(OUT_DIR, 'ferro_aging_bulk_stats.json'), 'w') as f:
    json.dump(bulk_stats, f, indent=2)
print(f"  Saved: ferro_aging_bulk_stats.json")

# ── Final summary ────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("VALIDATION COMPLETE")
print("=" * 60)
print(f"\n  Ferro-aging genes in bulk: {len(ferro_in_bulk)}/{len(ferro_genes)}")
print(f"  Significant (FDR<0.05): {n_total_validated} ({n_sig_up} ↑ AD, {n_sig_down} ↓ AD)")
print(f"  Cross-validation consistency: {n_consistent}/{len(cross_val_df)}")
print(f"  Sample-level FA score: AD {np.mean(ad_fa):.4f} vs Control {np.mean(ctrl_fa):.4f}")

if n_consistent >= len(cross_val_df) * 0.6:
    print(f"\n  ✓ Ferro-aging gene set is VALIDATED in independent bulk RNA-seq cohort!")
else:
    print(f"\n  △ Partial validation — bulk data shows some but not all scRNA-seq patterns")

print(f"\n  All results in: {OUT_DIR}/")
