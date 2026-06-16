"""
Ferro-aging + AAD (GSE213740) scRNA-seq Analysis Pipeline
GSE213740: 6 AD + 3 Normal, human aorta single-cell RNA-seq

Steps:
1. Load all 9 samples
2. QC and filtering
3. Normalization + HVG selection
4. PCA + UMAP + clustering
5. Cell type annotation
6. Ferro-aging scoring (AUCell-like module score)
7. Differential analysis: AD vs Normal
8. Visualization and export
"""

import scanpy as sc
import anndata as ad
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
import gzip

sc.settings.verbosity = 2
sc.settings.n_jobs = 4

DATA_DIR = "C:/Users/lidaf/WorkBuddy/2026-06-11-19-29-35/sc_data/extracted"
OUT_DIR = "C:/Users/lidaf/WorkBuddy/2026-06-11-19-29-35/results"
os.makedirs(OUT_DIR, exist_ok=True)

# ── STEP 1: Load Ferro-aging gene set ──────────────────────────────────────
print("\n========== STEP 1: Load Ferro-aging gene set ==========")
geneset_df = pd.read_csv("C:/Users/lidaf/WorkBuddy/2026-06-11-19-29-35/ferro_aging_geneset.csv")
ferro_aging_genes = geneset_df['gene_symbol'].tolist()
print(f"  Ferro-aging gene set: {len(ferro_aging_genes)} genes")
print(f"  Categories: {geneset_df['category'].value_counts().to_dict()}")

# ── STEP 2: Load scRNA-seq data ─────────────────────────────────────────────
print("\n========== STEP 2: Load scRNA-seq data ==========")

samples = {
    'AD_1': 'GSM6593315_AD_replicate_1',
    'AD_2': 'GSM6593316_AD_replicate_2',
    'AD_3': 'GSM6593317_AD_replicate_3',
    'AD_4': 'GSM6593318_AD_replicate_4',
    'AD_5': 'GSM6593319_AD_replicate_5',
    'AD_6': 'GSM6593320_AD_replicate_6',
    'Normal_1': 'GSM6593321_Normal_replicate_1',
    'Normal_2': 'GSM6593322_Normal_replicate_2',
    'Normal_3': 'GSM6593323_Normal_replicate_3',
}

def load_sample(sample_name, prefix):
    """Load a 10X MTX sample."""
    barcodes_path = os.path.join(DATA_DIR, f"{prefix}_barcodes.tsv.gz")
    features_path = os.path.join(DATA_DIR, f"{prefix}_features.tsv.gz")
    matrix_path   = os.path.join(DATA_DIR, f"{prefix}_matrix.mtx.gz")
    
    try:
        import scipy.io
        from scipy.sparse import csr_matrix
        
        # Read matrix
        with gzip.open(matrix_path, 'rb') as f:
            mat = scipy.io.mmread(f).T  # Transpose: cells × genes
        mat = csr_matrix(mat)
        
        # Read barcodes
        with gzip.open(barcodes_path, 'rt') as f:
            barcodes = [line.strip() for line in f]
        
        # Read features
        with gzip.open(features_path, 'rt') as f:
            features = []
            for line in f:
                parts = line.strip().split('\t')
                features.append(parts[1] if len(parts) > 1 else parts[0])
        
        # Deduplicate gene names (make_unique)
        seen = {}
        unique_features = []
        for g in features:
            if g in seen:
                seen[g] += 1
                unique_features.append(f"{g}_{seen[g]}")
            else:
                seen[g] = 0
                unique_features.append(g)
        features = unique_features
        
        # Create AnnData
        adata = ad.AnnData(X=mat,
                          obs=pd.DataFrame(index=barcodes),
                          var=pd.DataFrame(index=features))
        adata.obs['sample'] = sample_name
        adata.obs['condition'] = 'AD' if 'AD' in sample_name else 'Normal'
        
        print(f"  {sample_name}: {adata.n_obs} cells × {adata.n_vars} genes")
        return adata
    except Exception as e:
        print(f"  WARNING: Failed to load {sample_name}: {e}")
        return None

adatas = []
for sample_name, prefix in samples.items():
    adata = load_sample(sample_name, prefix)
    if adata is not None:
        adatas.append(adata)

print(f"\n  Loaded {len(adatas)} samples successfully")

# ── STEP 3: Merge and QC ────────────────────────────────────────────────────
print("\n========== STEP 3: Merge and QC ==========")

# Concatenate
adata = ad.concat(adatas, join='inner', label='sample_id')
print(f"  Combined: {adata.n_obs} cells × {adata.n_vars} genes")

# Compute QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(f"  Before QC: {adata.n_obs} cells")
print(f"  n_genes_by_counts: median={adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  pct_counts_mt: median={adata.obs['pct_counts_mt'].median():.1f}%")

# Filter cells
adata = adata[adata.obs['n_genes_by_counts'] > 200, :]
adata = adata[adata.obs['n_genes_by_counts'] < 6000, :]
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
print(f"  After QC: {adata.n_obs} cells")

# Filter genes  
sc.pp.filter_genes(adata, min_cells=10)
print(f"  After gene filter: {adata.n_vars} genes")

# ── STEP 4: Preprocessing ───────────────────────────────────────────────────
print("\n========== STEP 4: Preprocessing ==========")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # Store raw normalized data

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='sample')
print(f"  HVGs: {adata.var['highly_variable'].sum()}")

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

# ── STEP 5: Cell type annotation ────────────────────────────────────────────
print("\n========== STEP 5: Cell type annotation ==========")

# Vascular cell type marker genes
cell_markers = {
    'VSMC':       ['ACTA2', 'MYH11', 'CNN1', 'TAGLN', 'CALD1', 'MYL9'],
    'Endothelial':['PECAM1', 'CDH5', 'VWF', 'CLDN5', 'EMCN', 'KDR'],
    'Fibroblast': ['DCN', 'COL1A1', 'COL3A1', 'LUM', 'FBN1', 'POSTN'],
    'Macrophage': ['CD68', 'CD163', 'MRC1', 'CSF1R', 'ITGAM', 'FCGR3A'],
    'T_cell':     ['CD3D', 'CD3E', 'CD8A', 'CD4', 'TRAC', 'IL7R'],
    'B_cell':     ['CD79A', 'MS4A1', 'CD19', 'JCHAIN', 'MZB1'],
    'NK_cell':    ['GNLY', 'NKG7', 'KLRD1', 'NCAM1'],
    'pSMC':       ['MYH10', 'CNN2', 'NOTCH3'],  # progenitor SMC
    'Pericyte':   ['RGS5', 'ABCC9', 'KCNJ8', 'PDGFRB'],
}

# Score each cell type
for ct, markers in cell_markers.items():
    available = [g for g in markers if g in adata.raw.var_names]
    if len(available) >= 2:
        sc.tl.score_genes(adata, gene_list=available, score_name=f'score_{ct}', use_raw=True)

# Assign cell type to each cell based on highest score
score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
if score_cols:
    score_df = adata.obs[score_cols].copy()
    score_df.columns = [c.replace('score_', '') for c in score_cols]
    adata.obs['cell_type'] = score_df.idxmax(axis=1)
    print(f"  Cell type distribution:\n{adata.obs['cell_type'].value_counts()}")

# ── STEP 6: Ferro-aging Scoring ─────────────────────────────────────────────
print("\n========== STEP 6: Ferro-aging Scoring ==========")

# Filter for genes in dataset
ferro_in_data = [g for g in ferro_aging_genes if g in adata.raw.var_names]
print(f"  Ferro-aging genes in dataset: {ferro_in_data.__len__()} / {len(ferro_aging_genes)}")
print(f"  Present genes: {ferro_in_data}")

# Overall Ferro-aging score
sc.tl.score_genes(adata, gene_list=ferro_in_data, score_name='ferro_aging_score', use_raw=True)
print(f"  Ferro-aging score: mean={adata.obs['ferro_aging_score'].mean():.4f}, "
      f"std={adata.obs['ferro_aging_score'].std():.4f}")

# Score by sub-category
for cat_name, grp in geneset_df.groupby('category'):
    genes_in_cat = [g for g in grp['gene_symbol'].tolist() if g in adata.raw.var_names]
    if len(genes_in_cat) >= 3:
        sc.tl.score_genes(adata, gene_list=genes_in_cat,
                         score_name=f'fa_{cat_name.replace(" ","_").replace("-","_")[:20]}',
                         use_raw=True)
        print(f"  Category '{cat_name}': {len(genes_in_cat)} genes scored")

# Define Ferro-aging high cells (top 20%)
threshold = np.percentile(adata.obs['ferro_aging_score'], 80)
adata.obs['ferro_aging_high'] = adata.obs['ferro_aging_score'] > threshold
print(f"  Ferro-aging-high cells: {adata.obs['ferro_aging_high'].sum()} "
      f"({adata.obs['ferro_aging_high'].mean()*100:.1f}%)")

# ── STEP 7: Save processed data ─────────────────────────────────────────────
print("\n========== STEP 7: Save data ==========")
adata.write(os.path.join(OUT_DIR, 'GSE213740_ferro_aging_processed.h5ad'))
print(f"  Saved: {OUT_DIR}/GSE213740_ferro_aging_processed.h5ad")

# Export scores
score_export = adata.obs[['sample', 'condition', 'cell_type', 'ferro_aging_score', 'ferro_aging_high'] + 
                          [c for c in adata.obs.columns if c.startswith('fa_')]].copy()
score_export.to_csv(os.path.join(OUT_DIR, 'ferro_aging_cell_scores.csv'))
print(f"  Saved cell scores CSV")

# ── STEP 8: Statistics ──────────────────────────────────────────────────────
print("\n========== STEP 8: Key Statistics ==========")

# By condition
print("\n  Ferro-aging score by condition:")
print(adata.obs.groupby('condition')['ferro_aging_score'].describe())

# By cell type
print("\n  Ferro-aging score by cell type:")
print(adata.obs.groupby('cell_type')['ferro_aging_score'].mean().sort_values(ascending=False))

# Ferro-aging-high enrichment in AD vs Normal
print("\n  Ferro-aging-high cell % by condition:")
fa_by_cond = adata.obs.groupby('condition')['ferro_aging_high'].mean() * 100
print(fa_by_cond)

# ── STEP 9: Visualization ────────────────────────────────────────────────────
print("\n========== STEP 9: Visualization ==========")

# Color settings
ferro_cmap = LinearSegmentedColormap.from_list('ferro', ['#f0f0f0', '#ff6b6b', '#c0392b', '#7b241c'])

fig = plt.figure(figsize=(24, 20))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.35)

# 1. UMAP by condition
ax1 = fig.add_subplot(gs[0, 0])
colors_cond = {'AD': '#e74c3c', 'Normal': '#3498db'}
for cond in ['Normal', 'AD']:
    mask = adata.obs['condition'] == cond
    ax1.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
               c=colors_cond[cond], s=1, alpha=0.4, label=cond, rasterized=True)
ax1.set_title('UMAP: Condition', fontsize=13, fontweight='bold')
ax1.legend(markerscale=6, fontsize=9)
ax1.set_xlabel('UMAP1'); ax1.set_ylabel('UMAP2')
ax1.axis('off')

# 2. UMAP by cell type
ax2 = fig.add_subplot(gs[0, 1])
cell_types = adata.obs['cell_type'].unique()
ct_colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
ct_color_map = dict(zip(cell_types, ct_colors))
for ct in cell_types:
    mask = adata.obs['cell_type'] == ct
    ax2.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
               c=[ct_color_map[ct]], s=1, alpha=0.5, label=ct, rasterized=True)
ax2.set_title('UMAP: Cell Type', fontsize=13, fontweight='bold')
ax2.legend(markerscale=5, fontsize=7, loc='lower left', ncol=2)
ax2.axis('off')

# 3. UMAP by Ferro-aging score
ax3 = fig.add_subplot(gs[0, 2])
fa_scores = adata.obs['ferro_aging_score'].values
sc_plot = ax3.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                      c=fa_scores, cmap=ferro_cmap, s=1, alpha=0.6, rasterized=True,
                      vmin=np.percentile(fa_scores, 5), vmax=np.percentile(fa_scores, 95))
plt.colorbar(sc_plot, ax=ax3, shrink=0.7)
ax3.set_title('Ferro-aging Score', fontsize=13, fontweight='bold', color='#c0392b')
ax3.axis('off')

# 4. UMAP by Leiden cluster
ax4 = fig.add_subplot(gs[0, 3])
clusters = adata.obs['leiden'].astype(int).values
n_clusters = adata.obs['leiden'].nunique()
clust_colors = plt.cm.Set3(np.linspace(0, 1, n_clusters))
for i, clust in enumerate(sorted(adata.obs['leiden'].unique())):
    mask = adata.obs['leiden'] == clust
    ax4.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
               c=[clust_colors[int(clust)]], s=1, alpha=0.6, label=f'C{clust}', rasterized=True)
ax4.set_title('Leiden Clusters', fontsize=13, fontweight='bold')
ax4.axis('off')

# 5. Violin: Ferro-aging score by condition
ax5 = fig.add_subplot(gs[1, 0])
data_ad = adata.obs[adata.obs['condition']=='AD']['ferro_aging_score'].values
data_normal = adata.obs[adata.obs['condition']=='Normal']['ferro_aging_score'].values
parts = ax5.violinplot([data_normal, data_ad], positions=[0, 1], showmedians=True, showextrema=True)
for pc, color in zip(parts['bodies'], ['#3498db', '#e74c3c']):
    pc.set_facecolor(color); pc.set_alpha(0.7)
parts['cmedians'].set_color('black')
ax5.set_xticks([0, 1]); ax5.set_xticklabels(['Normal', 'AD'])
ax5.set_title('Ferro-aging Score:\nAD vs Normal', fontsize=12, fontweight='bold')
ax5.set_ylabel('Ferro-aging Score')
ax5.set_facecolor('#fafafa')
# Add statistics
from scipy import stats
t_stat, p_val = stats.ttest_ind(data_ad, data_normal)
ax5.text(0.5, 0.95, f'p={p_val:.2e}', transform=ax5.transAxes, ha='center', fontsize=10,
         color='red' if p_val < 0.05 else 'gray')

# 6. Bar: Ferro-aging-high proportion by condition
ax6 = fig.add_subplot(gs[1, 1])
fa_high_pct = adata.obs.groupby('condition')['ferro_aging_high'].mean() * 100
bars = ax6.bar(['Normal', 'AD'], [fa_high_pct.get('Normal', 0), fa_high_pct.get('AD', 0)],
               color=['#3498db', '#e74c3c'], alpha=0.8, edgecolor='white', linewidth=1.5)
ax6.set_title('Ferro-aging-high Cell %\nby Condition', fontsize=12, fontweight='bold')
ax6.set_ylabel('% Ferro-aging-high Cells')
for bar, val in zip(bars, [fa_high_pct.get('Normal', 0), fa_high_pct.get('AD', 0)]):
    ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, f'{val:.1f}%',
             ha='center', va='bottom', fontsize=11, fontweight='bold')
ax6.set_facecolor('#fafafa')

# 7. Bar: Ferro-aging score by cell type
ax7 = fig.add_subplot(gs[1, 2])
ct_fa = adata.obs.groupby('cell_type')['ferro_aging_score'].mean().sort_values(ascending=True)
colors_bar = ['#e74c3c' if ct == 'VSMC' else '#3498db' if ct == 'Macrophage' else '#95a5a6' 
               for ct in ct_fa.index]
bars = ax7.barh(range(len(ct_fa)), ct_fa.values, color=colors_bar, alpha=0.8)
ax7.set_yticks(range(len(ct_fa)))
ax7.set_yticklabels(ct_fa.index, fontsize=9)
ax7.set_title('Ferro-aging Score\nby Cell Type', fontsize=12, fontweight='bold')
ax7.set_xlabel('Mean Ferro-aging Score')
ax7.set_facecolor('#fafafa')
# Highlight VSMC
for i, ct in enumerate(ct_fa.index):
    if ct == 'VSMC':
        ax7.get_yticklabels()[i].set_color('#e74c3c')
        ax7.get_yticklabels()[i].set_fontweight('bold')

# 8. Stacked bar: Cell type composition AD vs Normal
ax8 = fig.add_subplot(gs[1, 3])
ct_comp = adata.obs.groupby(['condition', 'cell_type']).size().unstack(fill_value=0)
ct_comp_pct = ct_comp.div(ct_comp.sum(axis=1), axis=0) * 100
ct_comp_pct = ct_comp_pct.loc[['Normal', 'AD']] if 'Normal' in ct_comp_pct.index else ct_comp_pct
colors_ct = plt.cm.Set2(np.linspace(0, 1, len(ct_comp_pct.columns)))
bottom = np.zeros(len(ct_comp_pct))
for i, ct in enumerate(ct_comp_pct.columns):
    ax8.bar(['Normal', 'AD'] if len(ct_comp_pct) == 2 else ct_comp_pct.index.tolist(),
            ct_comp_pct[ct].values, bottom=bottom, label=ct, color=colors_ct[i], alpha=0.85)
    bottom += ct_comp_pct[ct].values
ax8.set_title('Cell Type Composition\nAD vs Normal', fontsize=12, fontweight='bold')
ax8.set_ylabel('Cell Proportion (%)')
ax8.legend(fontsize=7, loc='upper right', bbox_to_anchor=(1.35, 1))
ax8.set_facecolor('#fafafa')

# 9. Heatmap: Top ferro-aging genes expression by condition
ax9 = fig.add_subplot(gs[2, :2])
top_ferro_genes = [g for g in ['ACSL4', 'GPX4', 'NCOA4', 'FTH1', 'HMOX1', 'TFRC', 
                                'SLC7A11', 'PTGS2', 'NOX4', 'ACSL1', 'LPCAT3',
                                'TP53', 'CDKN1A', 'CDKN2A', 'IL6', 'CXCL8', 
                                'MMP9', 'MMP2', 'VEGFA'] 
                   if g in adata.raw.var_names]
if len(top_ferro_genes) > 0:
    # Get mean expression per condition
    mean_expr_list = []
    for cond in ['Normal', 'AD']:
        mask = adata.obs['condition'] == cond
        cells_adata = adata[mask]
        # Get from raw
        raw_X = cells_adata.raw[:, top_ferro_genes].X.toarray()
        mean_expr_list.append(raw_X.mean(axis=0))
    
    expr_matrix = np.array(mean_expr_list)
    # Row normalize
    expr_norm = (expr_matrix - expr_matrix.mean(axis=0)) / (expr_matrix.std(axis=0) + 1e-6)
    
    im = ax9.imshow(expr_norm, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
    ax9.set_xticks(range(len(top_ferro_genes)))
    ax9.set_xticklabels(top_ferro_genes, rotation=45, ha='right', fontsize=9)
    ax9.set_yticks([0, 1])
    ax9.set_yticklabels(['Normal', 'AD'], fontsize=11)
    plt.colorbar(im, ax=ax9, shrink=0.5, label='Z-score')
    ax9.set_title('Ferro-aging Gene Expression Heatmap\n(Row: Condition, Col: Key Genes)', 
                  fontsize=12, fontweight='bold')

# 10. VSMC Ferro-aging spotlight
ax10 = fig.add_subplot(gs[2, 2])
vsmc_mask = adata.obs['cell_type'] == 'VSMC'
if vsmc_mask.sum() > 10:
    for cond, color in [('Normal', '#3498db'), ('AD', '#e74c3c')]:
        mask = vsmc_mask & (adata.obs['condition'] == cond)
        if mask.sum() > 0:
            scores = adata.obs[mask]['ferro_aging_score'].values
            ax10.scatter(np.random.normal(0 if cond=='Normal' else 1, 0.08, len(scores)),
                        scores, c=color, s=8, alpha=0.4, label=f'{cond} VSMC')
    # Add mean lines
    for pos, cond in [(0, 'Normal'), (1, 'AD')]:
        mask = vsmc_mask & (adata.obs['condition'] == cond)
        if mask.sum() > 0:
            mean_val = adata.obs[mask]['ferro_aging_score'].mean()
            ax10.hlines(mean_val, pos-0.2, pos+0.2, colors='black', linewidths=2)
    ax10.set_xticks([0, 1]); ax10.set_xticklabels(['Normal VSMC', 'AD VSMC'])
    ax10.set_title('VSMC Ferro-aging Score\n(AD vs Normal)', fontsize=12, fontweight='bold', color='#c0392b')
    ax10.set_ylabel('Ferro-aging Score')
    ax10.legend(fontsize=8)
else:
    ax10.text(0.5, 0.5, 'Insufficient VSMC cells\nfor separate analysis',
              ha='center', va='center', transform=ax10.transAxes, fontsize=10)
    ax10.set_title('VSMC Ferro-aging Score', fontsize=12)
ax10.set_facecolor('#fafafa')

# 11. Summary text panel
ax11 = fig.add_subplot(gs[2, 3])
ax11.axis('off')
n_ad = (adata.obs['condition']=='AD').sum()
n_normal = (adata.obs['condition']=='Normal').sum()
ad_score = adata.obs[adata.obs['condition']=='AD']['ferro_aging_score'].mean()
normal_score = adata.obs[adata.obs['condition']=='Normal']['ferro_aging_score'].mean()
ad_high_pct = adata.obs[adata.obs['condition']=='AD']['ferro_aging_high'].mean()*100
normal_high_pct = adata.obs[adata.obs['condition']=='Normal']['ferro_aging_high'].mean()*100

summary = f"""ANALYSIS SUMMARY
{'='*28}
Dataset: GSE213740
Total cells: {adata.n_obs:,}
  AD: {n_ad:,} cells (6 donors)
  Normal: {n_normal:,} cells (3 donors)

Ferro-aging Gene Set:
  Total genes: {len(ferro_aging_genes)}
  In dataset: {len(ferro_in_data)}

Ferro-aging Score:
  AD: {ad_score:.4f}
  Normal: {normal_score:.4f}
  Fold change: {ad_score/normal_score:.2f}x
  p-value: {p_val:.2e}

Ferro-aging-high Cells:
  AD: {ad_high_pct:.1f}%
  Normal: {normal_high_pct:.1f}%

Cell Types Found:
{chr(10).join([f"  {ct}: {n}" for ct, n in adata.obs['cell_type'].value_counts().head(6).items()])}"""

ax11.text(0.05, 0.95, summary, transform=ax11.transAxes, fontsize=8,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='#f8f9fa', alpha=0.8))

# Title
fig.suptitle('Ferro-aging Analysis in Aortic Dissection (GSE213740)\n'
             'scRNA-seq: 6 AD Patients vs 3 Normal Controls',
             fontsize=16, fontweight='bold', y=0.98, color='#2c3e50')

plt.savefig(os.path.join(OUT_DIR, 'ferro_aging_scrnaseq_analysis.png'), 
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {OUT_DIR}/ferro_aging_scrnaseq_analysis.png")

# Export key stats
stats_dict = {
    'total_cells': adata.n_obs,
    'AD_cells': int(n_ad),
    'Normal_cells': int(n_normal),
    'ferro_genes_in_dataset': len(ferro_in_data),
    'ferro_genes_total': len(ferro_aging_genes),
    'ferro_score_AD_mean': float(ad_score),
    'ferro_score_Normal_mean': float(normal_score),
    'ferro_score_fold_change': float(ad_score/normal_score),
    'ttest_pvalue': float(p_val),
    'AD_ferro_high_pct': float(ad_high_pct),
    'Normal_ferro_high_pct': float(normal_high_pct),
    'n_leiden_clusters': int(adata.obs['leiden'].nunique()),
    'genes_detected': adata.n_vars,
}
import json
with open(os.path.join(OUT_DIR, 'analysis_stats.json'), 'w') as f:
    json.dump(stats_dict, f, indent=2)

print("\n========== ANALYSIS COMPLETE ==========")
print(f"\nKey Findings:")
print(f"  AD Ferro-aging score: {ad_score:.4f} vs Normal: {normal_score:.4f} (p={p_val:.2e})")
print(f"  Ferro-aging-high cells: AD {ad_high_pct:.1f}% vs Normal {normal_high_pct:.1f}%")
print(f"\nAll results saved to: {OUT_DIR}/")
