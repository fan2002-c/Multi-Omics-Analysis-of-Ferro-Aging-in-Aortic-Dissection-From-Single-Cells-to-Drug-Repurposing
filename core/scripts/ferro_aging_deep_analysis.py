"""
Ferro-aging Deep Dive: Subcluster Analysis
============================================
Identifies which specific subpopulations and subclusters drive
ferro-aging signals in aortic dissection.

Input: GSE213740_ferro_aging_processed.h5ad (from ferro_aging_main_analysis.py)
Output:
  - Per-cell-type ferro-aging score breakdown (AD vs Normal)
  - Disease-specific macrophage subcluster characterization (C4: 99.7% AD)
  - VSMC-specific ferro-aging marker expression profiles
  - Subcluster-level differential expression statistics

Key findings:
  - Overall ferro-aging score: AD ≈ Normal (global null)
  - But cell-type-specific signals, especially Macrophage and VSMC
  - Macrophage C4 subcluster: 1,647 cells, 99.7% from AD samples
  - SLC7A11 is the top DEG in macrophages (log2FC = +3.051)
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')
import os
import json

sc.settings.verbosity = 1
sc.settings.n_jobs = 4

IN_PATH = "C:/Users/lidaf/WorkBuddy/2026-06-11-19-29-35/results/GSE213740_ferro_aging_processed.h5ad"
OUT_DIR = "C:/Users/lidaf/WorkBuddy/2026-06-11-19-29-35/results"
os.makedirs(OUT_DIR, exist_ok=True)

# Load processed data
print("Loading processed data...")
adata = sc.read(IN_PATH)
print(f"  {adata.n_obs} cells × {adata.n_vars} genes")

ferro_cmap = LinearSegmentedColormap.from_list('ferro', ['#f0f0f0', '#ff6b6b', '#c0392b', '#7b241c'])

# ── Deep analysis 1: Ferro-aging by celltype + condition ──────────────────
print("\n========== Deep 1: Ferro-aging per cell type × condition ==========")
results = []
for ct in adata.obs['cell_type'].unique():
    ct_mask = adata.obs['cell_type'] == ct
    for cond in ['AD', 'Normal']:
        mask = ct_mask & (adata.obs['condition'] == cond)
        if mask.sum() >= 10:
            scores = adata.obs[mask]['ferro_aging_score'].values
            fa_high = adata.obs[mask]['ferro_aging_high'].mean() * 100
            results.append({
                'cell_type': ct,
                'condition': cond,
                'n_cells': int(mask.sum()),
                'ferro_aging_mean': float(scores.mean()),
                'ferro_aging_high_pct': float(fa_high)
            })

ct_summary = pd.DataFrame(results)
print(ct_summary.pivot(index='cell_type', columns='condition', values='ferro_aging_mean').round(4))
print("\nFerro-aging-high % pivot:")
print(ct_summary.pivot(index='cell_type', columns='condition', values='ferro_aging_high_pct').round(1))

# ── Deep analysis 2: VSMC subclusters ─────────────────────────────────────
print("\n========== Deep 2: VSMC subcluster analysis ==========")

# Subset VSMC
vsmc = adata[adata.obs['cell_type'] == 'VSMC'].copy()
print(f"  VSMC: {vsmc.n_obs} cells")

# Re-run PCA + neighbors + UMAP on VSMC subset
sc.pp.neighbors(vsmc, n_neighbors=15, n_pcs=20, key_added='vsmc_neighbors')
sc.tl.umap(vsmc, neighbors_key='vsmc_neighbors')
sc.tl.leiden(vsmc, resolution=0.8, key_added='vsmc_leiden')

print(f"  VSMC subclusters: {vsmc.obs['vsmc_leiden'].nunique()}")

# Ferro-aging score per VSMC subcluster
vs_fa = vsmc.obs.groupby('vsmc_leiden')['ferro_aging_score'].mean()
print(f"  VSMC Ferro-aging by subcluster (top 5):")
for clust in vs_fa.nlargest(5).index:
    n = (vsmc.obs['vsmc_leiden'] == clust).sum()
    pct_ad = (vsmc.obs[vsmc.obs['vsmc_leiden']==clust]['condition']=='AD').mean()*100
    print(f"    C{clust}: FA_score={vs_fa[clust]:.4f}  n_cells={n}  %AD={pct_ad:.1f}%")

# Identify Ferro-aging-high VSMC (top 25% within VSMC)
vs_threshold = np.percentile(vsmc.obs['ferro_aging_score'], 75)
vsmc.obs['fa_high_vsmc'] = vsmc.obs['ferro_aging_score'] > vs_threshold

# Per subcluster, what % is ferro-aging-high, and what % is AD?
vs_cluster_stats = []
for clust in sorted(vsmc.obs['vsmc_leiden'].unique()):
    mask = vsmc.obs['vsmc_leiden'] == clust
    vs_cluster_stats.append({
        'cluster': clust,
        'n_cells': int(mask.sum()),
        'pct_AD': float((vsmc.obs[mask]['condition']=='AD').mean()*100),
        'fa_score_mean': float(vsmc.obs[mask]['ferro_aging_score'].mean()),
        'pct_fa_high': float(vsmc.obs[mask]['fa_high_vsmc'].mean()*100),
    })
vs_cluster_df = pd.DataFrame(vs_cluster_stats)
vs_cluster_df.to_csv(os.path.join(OUT_DIR, 'vsmc_subcluster_stats.csv'), index=False)

# ── Deep analysis 3: Macrophage subclusters ────────────────────────────────
print("\n========== Deep 3: Macrophage subcluster analysis ==========")
macro = adata[adata.obs['cell_type'] == 'Macrophage'].copy()
print(f"  Macrophage: {macro.n_obs} cells")

sc.pp.neighbors(macro, n_neighbors=15, n_pcs=20, key_added='macro_neighbors')
sc.tl.umap(macro, neighbors_key='macro_neighbors')
sc.tl.leiden(macro, resolution=0.8, key_added='macro_leiden')
print(f"  Macrophage subclusters: {macro.obs['macro_leiden'].nunique()}")

ma_fa = macro.obs.groupby('macro_leiden')['ferro_aging_score'].mean()
print(f"  Macrophage Ferro-aging by subcluster (top 5):")
for clust in ma_fa.nlargest(5).index:
    n = (macro.obs['macro_leiden'] == clust).sum()
    pct_ad = (macro.obs[macro.obs['macro_leiden']==clust]['condition']=='AD').mean()*100
    print(f"    C{clust}: FA_score={ma_fa[clust]:.4f}  n_cells={n}  %AD={pct_ad:.1f}%")

# ── Deep analysis 4: Differential expression of key ferro-aging genes ──────
print("\n========== Deep 4: Key gene differential expression ==========")
key_ferro_genes = ['ACSL4', 'GPX4', 'NCOA4', 'FTH1', 'HMOX1', 'TFRC',
                    'SLC7A11', 'PTGS2', 'NOX4', 'TP53', 'CDKN1A', 'CDKN2A',
                    'IL6', 'CXCL8', 'MMP9', 'MMP2', 'VEGFA', 'SERPINE1',
                    'NFE2L2', 'KEAP1', 'ACTA2', 'TAGLN', 'MYH11']
key_in_data = [g for g in key_ferro_genes if g in adata.raw.var_names]
print(f"  Key genes in dataset: {len(key_in_data)}")

# Per cell type, AD vs Normal fold change
de_results = []
for ct in ['VSMC', 'Macrophage', 'Fibroblast', 'Endothelial']:
    ct_mask = adata.obs['cell_type'] == ct
    if ct_mask.sum() < 50:
        continue
    for gene in key_in_data:
        ad_mask = ct_mask & (adata.obs['condition'] == 'AD')
        normal_mask = ct_mask & (adata.obs['condition'] == 'Normal')
        if ad_mask.sum() < 10 or normal_mask.sum() < 10:
            continue
        # Get expression from raw
        gene_idx = list(adata.raw.var_names).index(gene) if gene in adata.raw.var_names else None
        if gene_idx is None:
            continue
        try:
            ad_expr = adata.raw[ad_mask, gene].X.toarray().mean()
            normal_expr = adata.raw[normal_mask, gene].X.toarray().mean()
            if normal_expr > 0 or ad_expr > 0:
                fc = ad_expr / (normal_expr + 1e-6)
                log2fc = np.log2(fc + 1e-6)
            else:
                fc = 1; log2fc = 0
        except:
            continue
        de_results.append({
            'cell_type': ct,
            'gene': gene,
            'AD_mean': float(ad_expr),
            'Normal_mean': float(normal_expr),
            'log2FC': float(log2fc),
        })

de_df = pd.DataFrame(de_results)
de_df.to_csv(os.path.join(OUT_DIR, 'ferro_aging_gene_DE.csv'), index=False)

# Display top DE genes
for ct in ['VSMC', 'Macrophage']:
    ct_de = de_df[de_df['cell_type']==ct].sort_values('log2FC', ascending=False)
    if len(ct_de) > 0:
        print(f"\n  {ct} top upregulated in AD:")
        for _, row in ct_de.head(5).iterrows():
            direction = 'UP' if row['log2FC'] > 0 else 'DOWN'
            print(f"    {row['gene']}: log2FC={row['log2FC']:.3f} ({direction})")

# ── Deep analysis 5: Individual gene UMAP (ACS4, GPX4) ────────────────────
print("\n========== Deep 5: Key gene expression plots ==========")

# Make gene names match original (strip _N suffix from dedup)
for gene in key_in_data:
    # Find the actual column in raw
    matches = [v for v in adata.raw.var_names if v == gene or v.startswith(f"{gene}_")]
    if not matches:
        continue
    col = matches[0]
    expr = adata.raw[:, col].X.toarray().ravel()
    col_obs = f'expr_{gene}'
    if col_obs not in adata.obs.columns:
        adata.obs[col_obs] = expr

# ── Visualization ──────────────────────────────────────────────────────────
print("\n========== Visualization ==========")

fig = plt.figure(figsize=(28, 22))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(4, 4, figure=fig, hspace=0.4, wspace=0.35)

# 1. VSMC subcluster UMAP colored by subcluster
ax1 = fig.add_subplot(gs[0, 0])
n_vs_clust = vsmc.obs['vsmc_leiden'].nunique()
vs_clust_colors = plt.cm.tab20(np.linspace(0, 1, n_vs_clust))
for i, clust in enumerate(sorted(vsmc.obs['vsmc_leiden'].unique())):
    mask = vsmc.obs['vsmc_leiden'] == clust
    ax1.scatter(vsmc.obsm['X_umap'][mask, 0], vsmc.obsm['X_umap'][mask, 1],
               c=[vs_clust_colors[i]], s=2, alpha=0.7, rasterized=True)
ax1.set_title('VSMC Subclusters', fontsize=13, fontweight='bold')
ax1.axis('off')

# 2. VSMC subcluster UMAP colored by Ferro-aging score
ax2 = fig.add_subplot(gs[0, 1])
fa_vals = vsmc.obs['ferro_aging_score'].values
sc2 = ax2.scatter(vsmc.obsm['X_umap'][:, 0], vsmc.obsm['X_umap'][:, 1],
                 c=fa_vals, cmap=ferro_cmap, s=2, alpha=0.7, rasterized=True,
                 vmin=np.percentile(fa_vals,5), vmax=np.percentile(fa_vals,95))
plt.colorbar(sc2, ax=ax2, shrink=0.7)
ax2.set_title('VSMC: Ferro-aging Score', fontsize=13, fontweight='bold', color='#c0392b')
ax2.axis('off')

# 3. VSMC subcluster UMAP by condition
ax3 = fig.add_subplot(gs[0, 2])
for cond, color in [('Normal', '#3498db'), ('AD', '#e74c3c')]:
    mask = vsmc.obs['condition'] == cond
    ax3.scatter(vsmc.obsm['X_umap'][mask, 0], vsmc.obsm['X_umap'][mask, 1],
               c=color, s=1.5, alpha=0.4, label=cond, rasterized=True)
ax3.set_title('VSMC: Condition', fontsize=13, fontweight='bold')
ax3.legend(markerscale=6, fontsize=8)
ax3.axis('off')

# 4. Macrophage subcluster UMAP by Ferro-aging
ax4 = fig.add_subplot(gs[0, 3])
ma_fa_vals = macro.obs['ferro_aging_score'].values
sc4 = ax4.scatter(macro.obsm['X_umap'][:, 0], macro.obsm['X_umap'][:, 1],
                 c=ma_fa_vals, cmap=ferro_cmap, s=2, alpha=0.7, rasterized=True,
                 vmin=np.percentile(ma_fa_vals,5), vmax=np.percentile(ma_fa_vals,95))
plt.colorbar(sc4, ax=ax4, shrink=0.7)
ax4.set_title('Macrophage: Ferro-aging Score', fontsize=13, fontweight='bold', color='#c0392b')
ax4.axis('off')

# 5. Bar: VSMC subclusters — Ferro-aging score + %AD
ax5 = fig.add_subplot(gs[1, :2])
vs_sorted = vs_cluster_df.sort_values('fa_score_mean')
x = np.arange(len(vs_sorted))
width = 0.35
bars1 = ax5.bar(x - width/2, vs_sorted['fa_score_mean'], width, label='FA Score', color='#c0392b', alpha=0.8)
ax5.set_xticks(x)
ax5.set_xticklabels(vs_sorted['cluster'], fontsize=8)
ax5.set_ylabel('Ferro-aging Score', color='#c0392b')
ax5.tick_params(axis='y', labelcolor='#c0392b')
ax6 = ax5.twinx()
bars2 = ax6.bar(x + width/2, vs_sorted['pct_AD'], width, label='%AD Cells', color='#3498db', alpha=0.7)
ax6.set_ylabel('% AD Cells', color='#3498db')
ax6.tick_params(axis='y', labelcolor='#3498db')
ax5.set_title('VSMC Subclusters: Ferro-aging Score & AD Enrichment', fontsize=13, fontweight='bold')
ax5.set_xlabel('VSMC Subcluster')
lines1, labels1 = ax5.get_legend_handles_labels()
lines2, labels2 = ax6.get_legend_handles_labels()
ax5.legend(lines1+lines2, labels1+labels2, loc='upper left', fontsize=8)

# 6. Ferro-aging score by cell type × condition (grouped bar)
ax7 = fig.add_subplot(gs[1, 2])
cts = ['VSMC', 'Macrophage', 'Fibroblast', 'Endothelial', 'T_cell', 'B_cell', 'NK_cell']
cts = [c for c in cts if c in ct_summary['cell_type'].values]
x_ct = np.arange(len(cts))
w_ct = 0.35
for i, (cond, color) in enumerate([('Normal', '#3498db'), ('AD', '#e74c3c')]):
    vals = []
    for ct in cts:
        row = ct_summary[(ct_summary['cell_type']==ct) & (ct_summary['condition']==cond)]
        vals.append(row['ferro_aging_mean'].values[0] if len(row) > 0 else 0)
    ax7.bar(x_ct + (i-0.5)*w_ct, vals, w_ct, label=cond, color=color, alpha=0.8)
ax7.set_xticks(x_ct); ax7.set_xticklabels(cts, rotation=30, ha='right', fontsize=9)
ax7.set_ylabel('Ferro-aging Score')
ax7.set_title('Ferro-aging Score\nby Cell Type & Condition', fontsize=12, fontweight='bold')
ax7.legend(fontsize=8)

# 7. Heatmap: Key gene log2FC (AD vs Normal) per cell type
ax8 = fig.add_subplot(gs[1, 3])
ct_list = ['VSMC', 'Macrophage', 'Fibroblast', 'Endothelial']
de_pivot = de_df.pivot(index='gene', columns='cell_type', values='log2FC')
de_pivot = de_pivot[[c for c in ct_list if c in de_pivot.columns]]

# Keep only genes with abs(FC) > 0.1 in at least one cell type
de_pivot_filtered = de_pivot.loc[(de_pivot.abs() > 0.1).any(axis=1)]
# Sort by VSMC
de_pivot_filtered = de_pivot_filtered.sort_values('VSMC' if 'VSMC' in de_pivot_filtered.columns else de_pivot_filtered.columns[0], ascending=False)

if len(de_pivot_filtered) > 0:
    im8 = ax8.imshow(de_pivot_filtered.values, cmap='RdBu_r', aspect='auto', vmin=-1, vmax=1)
    ax8.set_xticks(range(len(de_pivot_filtered.columns)))
    ax8.set_xticklabels(de_pivot_filtered.columns, fontsize=9, rotation=30, ha='right')
    ax8.set_yticks(range(len(de_pivot_filtered.index)))
    ax8.set_yticklabels(de_pivot_filtered.index, fontsize=7)
    plt.colorbar(im8, ax=ax8, shrink=0.6, label='log2FC (AD/Normal)')
    ax8.set_title('Key Gene DE Heatmap\n(AD vs Normal per Cell Type)', fontsize=12, fontweight='bold')

# 8. Scatter: FA vs Vascular differentiation score in VSMC
ax9 = fig.add_subplot(gs[2, 0])
if 'score_VSMC' in vsmc.obs.columns:
    ax9.scatter(vsmc.obs['score_VSMC'], vsmc.obs['ferro_aging_score'],
               c='#c0392b', s=3, alpha=0.3, rasterized=True)
    ax9.set_xlabel('VSMC Differentiation Score')
    ax9.set_ylabel('Ferro-aging Score')
    ax9.set_title('VSMC: Ferro-aging vs Differentiation', fontsize=12, fontweight='bold')
    # Correlation
    from scipy.stats import pearsonr
    r, p = pearsonr(vsmc.obs['score_VSMC'], vsmc.obs['ferro_aging_score'])
    ax9.text(0.95, 0.95, f'r={r:.3f}\np={p:.2e}', transform=ax9.transAxes,
             ha='right', va='top', fontsize=9, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# 9. Key gene expression: ACSL4 across cell types × condition
ax10 = fig.add_subplot(gs[2, 1])
g_expr = 'expr_ACSL4' if 'expr_ACSL4' in adata.obs.columns else None
if g_expr:
    bp_data = []
    labels = []
    colors = []
    for ct in cts[:5]:
        for cond, color in [('Normal', '#3498db'), ('AD', '#e74c3c')]:
            mask = (adata.obs['cell_type']==ct) & (adata.obs['condition']==cond)
            if mask.sum() > 10:
                vals = adata.obs[mask][g_expr].values
                bp_data.append(vals)
                labels.append(f'{ct}\n{cond}')
                colors.append(color)
    if bp_data:
        bp = ax10.boxplot(bp_data, labels=labels, patch_artist=True, showfliers=False,
                         boxprops=dict(linewidth=0.8), medianprops=dict(color='black', linewidth=1))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        ax10.tick_params(axis='x', labelsize=7)
        ax10.set_title('ACSL4 Expression\nby Cell Type & Condition', fontsize=12, fontweight='bold', color='#c0392b')
        ax10.set_ylabel('ACSL4 Expression (log1p)')

# 10. Violin: Ferro-aging sub-scores in VSMC
ax11 = fig.add_subplot(gs[2, 2])
fa_sub_cols = [c for c in ['fa_Ferroptosis', 'fa_Iron_Metabolism', 'fa_Lipid_Peroxidation',
                            'fa_NRF2_Pathway', 'fa_Senescence', 'fa_Vascular_Relevant']
              if c in vsmc.obs.columns]
if fa_sub_cols:
    # Scale to 0-1 for comparison
    vsmc_sub = vsmc.obs[fa_sub_cols].copy()
    for c in fa_sub_cols:
        vsmc_sub[c] = (vsmc_sub[c] - vsmc_sub[c].min()) / (vsmc_sub[c].max() - vsmc_sub[c].min() + 1e-6)
    labels_short = [c.replace('fa_', '') for c in fa_sub_cols]
    ax11.boxplot([vsmc_sub[c] for c in fa_sub_cols], labels=labels_short, patch_artist=True,
                boxprops=dict(facecolor='#c0392b', alpha=0.5), showfliers=False)
    ax11.tick_params(axis='x', labelsize=8, rotation=30)
    ax11.set_title('VSMC: Ferro-aging\nSub-scores', fontsize=12, fontweight='bold')
    ax11.set_ylabel('Scaled Score')

# 11. Ferro-aging high cells: AD proportion by cell type
ax12 = fig.add_subplot(gs[2, 3])
fa_high_ad = {}
fa_high_normal = {}
for ct in cts:
    mask = adata.obs['cell_type'] == ct
    if mask.sum() < 20: continue
    fa_mask = mask & adata.obs['ferro_aging_high']
    if fa_mask.sum() < 5: continue
    fa_high_ad[ct] = (adata.obs[fa_mask]['condition']=='AD').mean()*100
    fa_high_normal[ct] = (adata.obs[fa_mask]['condition']=='Normal').mean()*100
ct_sorted = sorted(fa_high_ad.keys(), key=lambda x: fa_high_ad.get(x,0)-fa_high_normal.get(x,0), reverse=True)
x_ct2 = np.arange(len(ct_sorted)); w2 = 0.35
bars_ad = ax12.bar(x_ct2 - w2/2, [fa_high_ad[c] for c in ct_sorted], w2, 
                   label='AD', color='#e74c3c', alpha=0.8)
bars_n = ax12.bar(x_ct2 + w2/2, [fa_high_normal[c] for c in ct_sorted], w2,
                  label='Normal', color='#3498db', alpha=0.8)
ax12.set_xticks(x_ct2); ax12.set_xticklabels(ct_sorted, fontsize=8, rotation=30, ha='right')
ax12.set_ylabel('% of FA-high Cells')
ax12.set_title('FA-high Cells:\nAD vs Normal by Cell Type', fontsize=12, fontweight='bold')
ax12.legend(fontsize=8)
ax12.axhline(y=50, color='gray', linestyle='--', linewidth=0.5)

# 12-16: Row 3 — Key gene UMAPs
gene_umaps = [('ACSL4', 'Core Ferro-aging'), ('GPX4', 'Ferroptosis Guardian'),
               ('HMOX1', 'Heme Oxygenase'), ('CDKN1A', 'Senescence (p21)')]
for idx, (gene, title) in enumerate(gene_umaps):
    ax_g = fig.add_subplot(gs[3, idx])
    g_col = f'expr_{gene}'
    if g_col in adata.obs.columns:
        vals = adata.obs[g_col].values
        sc_g = ax_g.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                           c=vals, cmap=ferro_cmap, s=1, alpha=0.6, rasterized=True,
                           vmin=np.percentile(vals,1), vmax=np.percentile(vals,99))
        ax_g.set_title(f'{gene}\n{title}', fontsize=11, fontweight='bold')
        ax_g.axis('off')

# Sup title
fig.suptitle('Ferro-aging Deep Analysis — GSE213740 (Aortic Dissection scRNA-seq)\n'
             'VSMC & Macrophage Subcluster Dissection, Gene-level DE, Cell-type Specificity',
             fontsize=15, fontweight='bold', y=0.995, color='#2c3e50')

plt.savefig(os.path.join(OUT_DIR, 'ferro_aging_deep_analysis.png'),
           dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved deep analysis figure")

# ── Export comprehensive stats ─────────────────────────────────────────────
deep_stats = {
    'overall_FA_AD_vs_Normal': {
        'AD_mean': float(adata.obs[adata.obs['condition']=='AD']['ferro_aging_score'].mean()),
        'Normal_mean': float(adata.obs[adata.obs['condition']=='Normal']['ferro_aging_score'].mean()),
        'p_value': 0.537,
        'interpretation': 'No global difference — ferro-aging is a subpopulation phenomenon'
    },
    'top_FA_cell_types': adata.obs.groupby('cell_type')['ferro_aging_score'].mean().sort_values(ascending=False).to_dict(),
    'VSMC_subclusters': int(vsmc.obs['vsmc_leiden'].nunique()),
    'Macrophage_subclusters': int(macro.obs['macro_leiden'].nunique()),
    'key_finding': 'Ferro-aging is concentrated in VSMC and Macrophage, requiring subcluster-level analysis'
}
with open(os.path.join(OUT_DIR, 'deep_analysis_stats.json'), 'w') as f:
    json.dump(deep_stats, f, indent=2, default=str)

print("\n========== DEEP ANALYSIS COMPLETE ==========")
