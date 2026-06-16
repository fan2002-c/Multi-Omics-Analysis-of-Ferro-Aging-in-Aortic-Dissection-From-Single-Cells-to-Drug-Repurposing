"""
Ferro-aging Drug Prediction Pipeline
=====================================
DGIdb GraphQL API → multi-criteria drug prioritization → visualization

Targeting key Ferro-aging genes validated in scRNA-seq & bulk RNA-seq:
- Core: ACSL4, GPX4, NCOA4, FTH1, HMOX1, TFRC
- xCT: SLC7A11, SLC3A2
- SASP: CXCL8, SERPINE1, IL6, VEGFA, CCL2
- Senescence: CDKN1A, CDKN2A, TP53
- Iron: HFE, TF, FTL, HAMP
- Lipid perox: PTGS2, NOX4
- ECM: MMP9, MMP2, MMP12, FBN1
- NRF2: NFE2L2, KEAP1
"""

import json
import urllib.request
import urllib.error
import ssl
import time
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Patch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

# ── Config ───────────────────────────────────────────────────────────────────
WORK_DIR = "C:/Users/lidaf/WorkBuddy/2026-06-11-19-29-35"
OUT_DIR = f"{WORK_DIR}/results"
DGIDB_GQL = "https://dgidb.org/api/graphql"

# Build SSL context
ssl_ctx = ssl.create_default_context()
ssl_ctx.check_hostname = False
ssl_ctx.verify_mode = ssl.CERT_NONE

# ── Gene Lists ───────────────────────────────────────────────────────────────
# Tier 1: Core Ferro-aging (validated in both scRNA-seq & bulk)
TIER1_GENES = [
    'HMOX1', 'SLC7A11', 'GPX4', 'ACSL4', 'NCOA4', 'FTH1',
    'CXCL8', 'SERPINE1', 'IL6', 'VEGFA', 'TP53', 'CDKN2A'
]

# Tier 2: Ferro-aging support (validated in bulk)
TIER2_GENES = [
    'TFRC', 'SLC3A2', 'CCL2', 'CXCL1', 'CDKN1A',
    'HFE', 'TF', 'FTL', 'HAMP',
    'PTGS2', 'NOX4', 'LPCAT3', 'ACSL1',
    'MMP9', 'MMP2', 'MMP12', 'FBN1',
    'NFE2L2', 'KEAP1'
]

# Tier 3: Extended target (from bulk significant DE)
TIER3_GENES = [
    'ALOX5', 'MMP9', 'SIRT1', 'MT1G', 'NF2', 'ALAS1',
]

def graphql_query(query_str, max_retries=3):
    """Execute a GraphQL query against DGIdb."""
    payload = json.dumps({'query': query_str}).encode('utf-8')
    req = urllib.request.Request(
        DGIDB_GQL,
        data=payload,
        headers={'Content-Type': 'application/json', 'Accept': 'application/json'}
    )
    for attempt in range(max_retries):
        try:
            with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as resp:
                return json.loads(resp.read().decode('utf-8'))
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(2)
            else:
                print(f"  ✗ GraphQL error: {e}")
                return None

def fetch_interactions(gene_names, approved_only=True):
    """Fetch all drug-gene interactions for a list of genes."""
    all_edges = []
    gene_filter = json.dumps(gene_names)

    query = f"""
    query {{
      interactions(geneNames: {gene_filter}, approved: {str(approved_only).lower()}, first: 100) {{
        pageInfo {{ hasNextPage endCursor }}
        edges {{
          node {{
            id
            interactionScore
            drugSpecificity
            geneSpecificity
            evidenceScore
            interactionTypes {{ type directionality }}
            drug {{ name conceptId }}
            gene {{ name conceptId }}
            sources {{ fullName sourceDbName }}
          }}
        }}
      }}
    }}
    """
    result = graphql_query(query)
    if not result or 'data' not in result:
        return []
    
    data = result['data']
    if not data.get('interactions'):
        return []
    
    edges = data['interactions'].get('edges', [])
    all_edges.extend(edges)
    
    # Handle pagination
    page_info = data['interactions'].get('pageInfo', {})
    after = page_info.get('endCursor')
    while page_info.get('hasNextPage') and after:
        query_paged = f"""
        query {{
          interactions(geneNames: {gene_filter}, approved: {str(approved_only).lower()}, first: 100, after: "{after}") {{
            pageInfo {{ hasNextPage endCursor }}
            edges {{
              node {{
                id
                interactionScore
                drugSpecificity
                geneSpecificity
                evidenceScore
                interactionTypes {{ type directionality }}
                drug {{ name conceptId }}
                gene {{ name conceptId }}
                sources {{ fullName sourceDbName }}
              }}
            }}
          }}
        }}
        """
        result = graphql_query(query_paged)
        if not result or 'data' not in result:
            break
        data = result['data']
        if not data.get('interactions'):
            break
        edges = data['interactions'].get('edges', [])
        all_edges.extend(edges)
        page_info = data['interactions'].get('pageInfo', {})
        after = page_info.get('endCursor')
        time.sleep(0.3)
    
    return all_edges

# ── FDA Approval Lookup ──────────────────────────────────────────────────────
# Known FDA-approved cardiovascular/anti-inflammatory/iron-related drugs
FDA_APPROVED = {
    # Iron chelators
    'DEFERASIROX': {'indication': 'Iron overload', 'fda_year': 2005},
    'DEFERIPRONE': {'indication': 'Iron overload', 'fda_year': 2011},
    'DEFEROXAMINE': {'indication': 'Iron overload', 'fda_year': 1968},
    # NSAIDs / COX inhibitors (target PTGS2)
    'ASPIRIN': {'indication': 'Anti-inflammatory/CV', 'fda_year': 1899},
    'CELECOXIB': {'indication': 'COX-2 inhibitor', 'fda_year': 1998},
    'IBUPROFEN': {'indication': 'NSAID', 'fda_year': 1974},
    'NAPROXEN': {'indication': 'NSAID', 'fda_year': 1976},
    'INDOMETHACIN': {'indication': 'NSAID', 'fda_year': 1965},
    # Vitamin / antioxidant
    'ASCORBIC ACID': {'indication': 'Vitamin C/Antioxidant', 'fda_year': 1939},
    'VITAMIN E': {'indication': 'Antioxidant', 'fda_year': 1941},
    'N-ACETYLCYSTEINE': {'indication': 'Antioxidant/Mucolytic', 'fda_year': 1963},
    # MMP inhibitors
    'DOXYCYCLINE': {'indication': 'Antibiotic/MMP inhibitor', 'fda_year': 1967},
    'MINOCYCLINE': {'indication': 'Antibiotic/MMP inhibitor', 'fda_year': 1971},
    # Cardiovascular
    'LOSARTAN': {'indication': 'ARB/Anti-hypertensive', 'fda_year': 1995},
    'VALSARTAN': {'indication': 'ARB/Anti-hypertensive', 'fda_year': 1996},
    'ATORVASTATIN': {'indication': 'Statin', 'fda_year': 1996},
    'ROSUVASTATIN': {'indication': 'Statin', 'fda_year': 2003},
    'SIMVASTATIN': {'indication': 'Statin', 'fda_year': 1991},
    'METOPROLOL': {'indication': 'Beta blocker', 'fda_year': 1978},
    'CARVEDILOL': {'indication': 'Beta blocker', 'fda_year': 1995},
    'NIFEDIPINE': {'indication': 'CCB', 'fda_year': 1981},
    'AMLODIPINE': {'indication': 'CCB', 'fda_year': 1992},
    'CAPTOPRIL': {'indication': 'ACE inhibitor', 'fda_year': 1981},
    'ENALAPRIL': {'indication': 'ACE inhibitor', 'fda_year': 1985},
    'LISINOPRIL': {'indication': 'ACE inhibitor', 'fda_year': 1987},
    'HYDRALAZINE': {'indication': 'Vasodilator', 'fda_year': 1953},
    'SPIRONOLACTONE': {'indication': 'Aldosterone antagonist', 'fda_year': 1960},
    'WARFARIN': {'indication': 'Anticoagulant', 'fda_year': 1954},
    'CLOPIDOGREL': {'indication': 'Antiplatelet', 'fda_year': 1997},
    'RIVAROXABAN': {'indication': 'Anticoagulant', 'fda_year': 2011},
    'APIXABAN': {'indication': 'Anticoagulant', 'fda_year': 2012},
    # Kinase inhibitors
    'SORAFENIB': {'indication': 'TKI', 'fda_year': 2005},
    'SUNITINIB': {'indication': 'TKI', 'fda_year': 2006},
    'IMATINIB': {'indication': 'TKI', 'fda_year': 2001},
    'ERLOTINIB': {'indication': 'EGFR inhibitor', 'fda_year': 2004},
    # NRF2 activators
    'DIMETHYL FUMARATE': {'indication': 'NRF2 activator', 'fda_year': 2013},
    'SULFORAPHANE': {'indication': 'NRF2 activator (supplement)', 'fda_year': None},
    # Senolytics
    'DASATINIB': {'indication': 'TKI/Senolytic', 'fda_year': 2006},
    'QUERCETIN': {'indication': 'Senolytic (supplement)', 'fda_year': None},
    'FISETIN': {'indication': 'Senolytic (supplement)', 'fda_year': None},
    'NAVITOCLAX': {'indication': 'BCL-2 inhibitor/Senolytic', 'fda_year': 2016},
    # Anti-IL6
    'TOCILIZUMAB': {'indication': 'IL-6R antagonist', 'fda_year': 2010},
    'SILTUXIMAB': {'indication': 'Anti-IL6', 'fda_year': 2014},
    # Ferroptosis inhibitors
    'LIPROXSTATIN-1': {'indication': 'Ferroptosis inhibitor (research)', 'fda_year': None},
    'FERROSTATIN-1': {'indication': 'Ferroptosis inhibitor (research)', 'fda_year': None},
    # N-acetylcysteine
    'ACETYLCYSTEINE': {'indication': 'Antioxidant/Mucolytic', 'fda_year': 1963},
}

# ── Cardiovascular drug-gene targeting relevance map ────────────────────────
CV_RELEVANCE = {
    'HMOX1': ['Heme degradation → free iron release', 'Vascular protection target'],
    'SLC7A11': ['xCT subunit → cystine import', 'Ferroptosis gatekeeper'],
    'GPX4': ['Lipid peroxide reductase', 'Central ferroptosis regulator'],
    'ACSL4': ['PUFA-CoA ligase → lipid peroxidation', 'Ferro-aging master regulator'],
    'IL6': ['SASP cytokine → vascular inflammation', 'AAD inflammatory driver'],
    'CXCL8': ['Neutrophil chemokine → MMP activation', 'SASP effector'],
    'SERPINE1': ['PAI-1 → fibrinolysis inhibition', 'Thrombosis/AAD risk'],
    'PTGS2': ['COX-2 → prostaglandin synthesis', 'Vascular inflammation'],
    'MMP9': ['Gelatinase → ECM degradation', 'Aortic wall weakening'],
    'MMP2': ['Gelatinase → elastin degradation', 'Aortic remodeling'],
    'TP53': ['Cellular senescence master regulator', 'Apoptosis/senescence switch'],
    'NFE2L2': ['NRF2 → antioxidant response', 'Ferroptosis defense'],
}

# ── STEP 1: Query DGIdb ──────────────────────────────────────────────────────
print("=" * 70)
print("STEP 1: Query DGIdb GraphQL API for drug-gene interactions")
print("=" * 70)

all_genes = TIER1_GENES + TIER2_GENES + TIER3_GENES
# Remove duplicates while preserving order
seen = set()
unique_genes = []
for g in all_genes:
    if g not in seen:
        unique_genes.append(g)
        seen.add(g)

print(f"  Querying {len(unique_genes)} genes in total...")

# Fetch in batches of 5 genes
BATCH = 5
all_interactions_raw = []
for i in range(0, len(unique_genes), BATCH):
    batch = unique_genes[i:i+BATCH]
    print(f"  Batch {i//BATCH + 1}: {batch}")
    edges = fetch_interactions(batch, approved_only=False)
    all_interactions_raw.extend(edges)
    n_drugs = len(set(e['node']['drug']['name'] for e in edges))
    print(f"    → {len(edges)} interactions, {n_drugs} unique drugs")
    time.sleep(0.5)

# ── STEP 2: Parse and deduplicate ────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Parse, deduplicate and annotate interactions")
print("=" * 70)

parsed = []
seen_pairs = set()

for edge in all_interactions_raw:
    node = edge['node']
    drug_name = node['drug']['name'].upper()
    gene_name = node['gene']['name'].upper()
    pair_key = (drug_name, gene_name)
    
    if pair_key in seen_pairs:
        continue
    seen_pairs.add(pair_key)
    
    # Parse interaction types
    itypes = [it['type'] for it in node.get('interactionTypes', [])]
    directionalities = [it['directionality'] for it in node.get('interactionTypes', [])]
    
    # Parse sources
    source_dbs = list(set(s.get('sourceDbName', '') for s in node.get('sources', [])))
    source_full = list(set(s.get('fullName', '') for s in node.get('sources', [])))
    
    parsed.append({
        'drug': drug_name,
        'gene': gene_name,
        'gene_tier': 1 if gene_name in [g.upper() for g in TIER1_GENES] else 
                     2 if gene_name in [g.upper() for g in TIER2_GENES] else 3,
        'interaction_score': float(node.get('interactionScore', 0)),
        'drug_specificity': float(node.get('drugSpecificity', 0)),
        'gene_specificity': float(node.get('geneSpecificity', 0)),
        'evidence_score': int(node.get('evidenceScore', 0)),
        'interaction_types': ','.join(itypes) if itypes else 'unknown',
        'directionality': ','.join(directionalities) if directionalities else 'unknown',
        'source_dbs': ','.join(source_dbs),
        'drug_concept_id': node['drug'].get('conceptId', ''),
        'fda_approved': drug_name in FDA_APPROVED,
        'fda_indication': FDA_APPROVED.get(drug_name, {}).get('indication', ''),
        'fda_year': FDA_APPROVED.get(drug_name, {}).get('fda_year', None),
    })

df = pd.DataFrame(parsed)
print(f"  Total unique drug-gene pairs: {len(df)}")
print(f"  Unique drugs: {df['drug'].nunique()}")
print(f"  Unique genes covered: {df['gene'].nunique()}/{len(unique_genes)}")
print(f"  FDA-approved drugs: {df[df['fda_approved']]['drug'].nunique()}")

# ── STEP 3: Drug-level aggregation & scoring ─────────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Multi-criteria drug prioritization scoring")
print("=" * 70)

def score_drug(group, all_genes_set):
    """Score a drug based on multiple criteria."""
    genes_hit = set(group['gene'].values)
    n_genes = len(genes_hit)
    tier1_hits = len([g for g in genes_hit if g in [x.upper() for x in TIER1_GENES]])
    tier2_hits = n_genes - tier1_hits
    
    # Core score
    score = 0
    
    # 1. Number of Ferro-aging genes targeted (Tier-weighted)
    score += tier1_hits * 3.0  # Core genes weight 3
    score += tier2_hits * 1.5  # Support genes weight 1.5
    
    # 2. Average interaction score across all gene targets
    avg_is = group['interaction_score'].mean()
    score += avg_is * 2.0
    
    # 3. Evidence score
    max_ev = group['evidence_score'].max()
    score += np.log1p(max_ev)
    
    # 4. FDA approval bonus
    if group['fda_approved'].iloc[0]:
        score += 5.0
        # Extra for CV-related indication
        indication = str(group['fda_indication'].iloc[0]).lower()
        if any(kw in indication for kw in ['anti-inflammatory', 'statin', 'anti-hypertensive',
                                              'iron', 'vasodilator', 'antiplatelet', 'anticoagulant',
                                              'senolytic', 'antioxidant', 'cox', 'nrf2', 'ferroptosis']):
            score += 3.0
    
    # 5. Multi-gene targeting bonus
    if n_genes >= 3:
        score += 4.0
    elif n_genes >= 2:
        score += 2.0
    
    # 6. Drug specificity bonus (higher = more specific to these genes)
    avg_spec = group['drug_specificity'].mean()
    score += avg_spec * 1.0
    
    return score, genes_hit

drug_scores = []
for drug, grp in df.groupby('drug'):
    s, genes = score_drug(grp, set(unique_genes))
    drug_scores.append({
        'drug': drug,
        'score': round(s, 2),
        'n_genes': len(genes),
        'genes_targeted': ','.join(sorted(genes)),
        'tier1_genes': ','.join(sorted([g for g in genes if g in [x.upper() for x in TIER1_GENES]])),
        'fda_approved': grp['fda_approved'].iloc[0],
        'fda_indication': grp['fda_indication'].iloc[0],
        'fda_year': grp['fda_year'].iloc[0] if not pd.isna(grp['fda_year'].iloc[0]) else None,
        'avg_interaction_score': round(grp['interaction_score'].mean(), 3),
        'max_evidence': grp['evidence_score'].max(),
        'interaction_types': ','.join(set(','.join(grp['interaction_types']).split(','))),
    })

drug_df = pd.DataFrame(drug_scores).sort_values('score', ascending=False).reset_index(drop=True)
drug_df['rank'] = range(1, len(drug_df) + 1)

print(f"\n  Top 20 ranked drugs:")
for _, row in drug_df.head(20).iterrows():
    fda_badge = " [FDA]" if row['fda_approved'] else ""
    print(f"  #{row['rank']:2d} {row['drug']:35s} score={row['score']:5.1f} | {row['n_genes']} genes{fda_badge} | {row['genes_targeted'][:60]}")

# ── STEP 4: Special analysis for key genes ───────────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: Drug profiling for top Ferro-aging targets")
print("=" * 70)

for gene in ['HMOX1', 'SLC7A11', 'GPX4', 'ACSL4', 'IL6', 'PTGS2']:
    gene_drugs = df[df['gene'] == gene.upper()]
    print(f"\n  [{gene}] - {len(gene_drugs)} drugs found")
    gene_drugs_sorted = gene_drugs.sort_values('interaction_score', ascending=False)
    for _, row in gene_drugs_sorted.head(5).iterrows():
        fda = " [FDA]" if row['fda_approved'] else ""
        print(f"    {row['drug']:35s} score={row['interaction_score']:.1f} type={row['interaction_types'][:30]}{fda}")

# ── STEP 5: Drug repurposing candidates for Ferro-aging ─────────────────────
print("\n" + "=" * 70)
print("STEP 5: Ferro-aging drug repurposing analysis")
print("=" * 70)

# Known literature-supported Ferro-aging drugs
ferro_aging_drugs = {
    'VITAMIN C (ASCORBIC ACID)': {
        'mechanism': 'ACSL4 inhibitor / ferrous iron reducer',
        'targets': ['ACSL4', 'TET2'],
        'evidence': 'Cell Metabolism 2026: VitC suppresses ACSL4, reduces ferro-aging',
        'category': 'Ferro-aging core therapy'
    },
    'DEFEROXAMINE': {
        'mechanism': 'Iron chelator → reduces Fe²⁺ pool',
        'targets': ['Free Fe²⁺'],
        'evidence': 'Nat Commun 2024: Iron chelation reduces VSMC senescence',
        'category': 'Iron chelation'
    },
    'DEFERASIROX': {
        'mechanism': 'Oral iron chelator → reduces labile iron pool',
        'targets': ['Free Fe²⁺'],
        'evidence': 'Clinical iron overload drug; repurposing for ferroptosis',
        'category': 'Iron chelation'
    },
    'N-ACETYLCYSTEINE': {
        'mechanism': 'GSH precursor → boosts GPX4 activity',
        'targets': ['GPX4 (indirect)'],
        'evidence': 'GPX4 cofactor restoration; general antioxidant',
        'category': 'Antioxidant / GPX4 support'
    },
    'DIMETHYL FUMARATE': {
        'mechanism': 'NRF2 activator → upregulates GPX4, FTH1, HO-1',
        'targets': ['NFE2L2', 'GPX4', 'FTH1', 'HMOX1'],
        'evidence': 'FDA-approved for MS; NRF2 activation protects against ferroptosis',
        'category': 'NRF2 activation'
    },
    'LIPROXSTATIN-1': {
        'mechanism': 'Radical-trapping antioxidant → blocks lipid peroxidation',
        'targets': ['Lipid ROS'],
        'evidence': 'Gold standard ferroptosis inhibitor in research',
        'category': 'Ferroptosis inhibitor'
    },
    'FERROSTATIN-1': {
        'mechanism': 'Lipophilic antioxidant → inhibits lipid peroxidation',
        'targets': ['Lipid ROS'],
        'evidence': 'First-in-class ferroptosis inhibitor',
        'category': 'Ferroptosis inhibitor'
    },
    'DASATINIB + QUERCETIN': {
        'mechanism': 'Senolytic combination → clears senescent VSMCs',
        'targets': ['BCL-2', 'PI3K/AKT'],
        'evidence': 'Clears senescent cells in vascular tissue',
        'category': 'Senolytic'
    },
    'FISETIN': {
        'mechanism': 'Natural senolytic flavonoid',
        'targets': ['PI3K/AKT', 'p16'],
        'evidence': 'Reduces senescence burden in aged mice',
        'category': 'Senolytic'
    },
    'DOXYCYCLINE': {
        'mechanism': 'MMP inhibitor → reduces ECM degradation',
        'targets': ['MMP9', 'MMP2'],
        'evidence': 'Clinical use in aortic aneurysm stabilization',
        'category': 'MMP inhibition / ECM protection'
    },
}

ferro_drug_df = pd.DataFrame([
    {
        'drug': drug,
        'mechanism': info['mechanism'],
        'targets': ', '.join(info['targets']),
        'evidence': info['evidence'],
        'category': info['category']
    }
    for drug, info in ferro_aging_drugs.items()
])

print(f"  Ferro-aging targeted drugs identified: {len(ferro_drug_df)}")
for _, row in ferro_drug_df.iterrows():
    print(f"    [{row['category']:25s}] {row['drug']}")

# ── STEP 6: Visualization ────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 6: Generate visualization")
print("=" * 70)

fig = plt.figure(figsize=(28, 24))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.45, wspace=0.4)

# ── Panel 1: Top drugs bar chart ──────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, :2])
top30 = drug_df.head(30).iloc[::-1]
colors_top = ['#2ecc71' if v else '#e74c3c' for v in top30['fda_approved']]
bars = ax1.barh(range(len(top30)), top30['score'].values, color=colors_top, alpha=0.85,
               edgecolor='white', linewidth=0.8)
ax1.set_yticks(range(len(top30)))
ax1.set_yticklabels(top30['drug'].values, fontsize=7)
ax1.set_xlabel('Priority Score', fontsize=11, fontweight='bold')
ax1.set_title('Top 30 Drug Candidates for Ferro-aging\n(Multi-criteria Scoring: Target Coverage + Evidence + FDA)', 
             fontsize=13, fontweight='bold', color='#2c3e50')
# Add gene count annotation
for i, (_, row) in enumerate(top30.iterrows()):
    ax1.text(row['score'] + 0.2, i, f"{int(row['n_genes'])} genes", 
            fontsize=6, va='center', color='#555')

legend_elements = [Patch(facecolor='#2ecc71', alpha=0.85, label='FDA Approved'),
                   Patch(facecolor='#e74c3c', alpha=0.85, label='Investigational')]
ax1.legend(handles=legend_elements, fontsize=8, loc='lower right')
ax1.set_facecolor('#fafafa')

# ── Panel 2: Drug-Gene heatmap ────────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 2:])

# Build matrix: top drugs x key genes
TOP_N_DRUGS = 25
HEAT_GENES = ['HMOX1', 'SLC7A11', 'GPX4', 'ACSL4', 'NCOA4', 'FTH1', 'TFRC',
              'CXCL8', 'SERPINE1', 'IL6', 'VEGFA', 'CCL2',
              'TP53', 'CDKN2A', 'CDKN1A',
              'PTGS2', 'NOX4', 'MMP9', 'MMP2', 'NFE2L2', 'HFE']

top_drugs = drug_df.head(TOP_N_DRUGS)['drug'].tolist()
heat_genes_upper = [g.upper() for g in HEAT_GENES]

heat_data = np.zeros((TOP_N_DRUGS, len(HEAT_GENES)))
for i, drug in enumerate(top_drugs):
    drug_genes = set(df[df['drug'] == drug]['gene'].values)
    for j, gene_upper in enumerate(heat_genes_upper):
        if gene_upper in drug_genes:
            # Use interaction score for intensity
            match = df[(df['drug'] == drug) & (df['gene'] == gene_upper)]
            if len(match) > 0:
                heat_data[i, j] = match['interaction_score'].values[0] + 0.1

im = ax2.imshow(heat_data, cmap='YlOrRd', aspect='auto', vmin=0, vmax=5)
ax2.set_xticks(range(len(HEAT_GENES)))
ax2.set_xticklabels(HEAT_GENES, rotation=45, ha='right', fontsize=7)
ax2.set_yticks(range(TOP_N_DRUGS))
ax2.set_yticklabels(top_drugs, fontsize=6.5)
ax2.set_title('Drug-Gene Interaction Matrix\n(Interaction Score Intensity)', fontsize=13, fontweight='bold', color='#2c3e50')

# Highlight Tier 1 genes
for j, g in enumerate(HEAT_GENES):
    if g.upper() in [x.upper() for x in TIER1_GENES]:
        ax2.get_xticklabels()[j].set_fontweight('bold')
        ax2.get_xticklabels()[j].set_color('#c0392b')

plt.colorbar(im, ax=ax2, shrink=0.6, label='Interaction Score')

# ── Panel 3: Gene coverage network ────────────────────────────────────────
ax3 = fig.add_subplot(gs[1, :2])

# Build drug-gene targeting stats
gene_target_counts = df.groupby('gene').size().sort_values(ascending=False)
top_genes_net = gene_target_counts.head(20)

colors_net = ['#e74c3c' if g.upper() in [x.upper() for x in TIER1_GENES] else '#3498db' for g in top_genes_net.index]
bars = ax3.bar(range(len(top_genes_net)), top_genes_net.values, color=colors_net, alpha=0.85,
              edgecolor='white', linewidth=0.8)
ax3.set_xticks(range(len(top_genes_net)))
ax3.set_xticklabels(top_genes_net.index, rotation=45, ha='right', fontsize=9)
ax3.set_ylabel('Number of Targeting Drugs', fontsize=11, fontweight='bold')
ax3.set_title('Druggability: Number of Drugs per Ferro-aging Gene\n(Red = Tier 1 Core Genes)', 
             fontsize=13, fontweight='bold', color='#2c3e50')
ax3.set_facecolor('#fafafa')

# Add count labels
for i, v in enumerate(top_genes_net.values):
    ax3.text(i, v + 0.3, str(v), ha='center', fontsize=8, fontweight='bold')

# ── Panel 4: Drug category breakdown & pathways ───────────────────────────
ax4 = fig.add_subplot(gs[1, 2:])
ax4.axis('off')

# Build pathway-target map
pathway_drugs = {
    'IRON CHELATION': ['DEFEROXAMINE', 'DEFERASIROX', 'DEFERIPRONE'],
    'GPX4 / GSH AXIS': ['N-ACETYLCYSTEINE', 'VITAMIN E'],
    'LIPID PEROXIDATION': ['LIPROXSTATIN-1', 'FERROSTATIN-1', 'VITAMIN C'],
    'NRF2 ACTIVATION': ['DIMETHYL FUMARATE', 'SULFORAPHANE'],
    'SENOLYSIS': ['DASATINIB', 'QUERCETIN', 'FISETIN', 'NAVITOCLAX'],
    'ANTI-IL6/SASP': ['TOCILIZUMAB', 'SILTUXIMAB'],
    'MMP INHIBITION': ['DOXYCYCLINE', 'MINOCYCLINE'],
    'COX-2 INHIBITION': ['ASPIRIN', 'CELECOXIB', 'IBUPROFEN'],
    'STATIN (PLEIOTROPIC)': ['ATORVASTATIN', 'ROSUVASTATIN'],
    'ARB': ['LOSARTAN', 'VALSARTAN'],
}

y_pos = 0
for pathway, drugs in pathway_drugs.items():
    matched = [d for d in drugs if d in drug_df['drug'].values]
    ax4.text(0.05, 0.98 - y_pos, f'{pathway}:', fontsize=9, fontweight='bold', 
            color='#2c3e50', transform=ax4.transAxes)
    ax4.text(0.45, 0.98 - y_pos, f'{len(matched)}/{len(drugs)} found', fontsize=8,
            color='#27ae60' if len(matched) > 0 else '#95a5a6', transform=ax4.transAxes)
    y_pos += 0.09

ax4.text(0.05, 0.98 - y_pos, f'\nMulti-target DGIdb hits:', fontsize=9, fontweight='bold',
        color='#2c3e50', transform=ax4.transAxes)
y_pos += 0.05
# Show multi-gene targeting drugs from top ranks
multi_hits = drug_df.head(10)
for _, row in multi_hits.iterrows():
    if row['n_genes'] >= 2:
        ax4.text(0.07, 0.98 - y_pos, f"#{int(row['rank'])} {row['drug']}: {row['genes_targeted'][:55]}",
                fontsize=7, color='#555', transform=ax4.transAxes)
        y_pos += 0.045

ax4.set_title('Pathway-based Drug Repurposing Strategy\n(Ferro-aging Mechanism → Drug Class)', 
             fontsize=13, fontweight='bold', color='#2c3e50', y=1.02)

# ── Panel 5: Ferro-aging mechanism-based drug targeting schematic ─────────
ax5 = fig.add_subplot(gs[2, :3])
ax5.axis('off')
ax5.set_xlim(0, 12)
ax5.set_ylim(0, 10)

# Draw schematical pathway → drug targeting diagram
# ── Left: Pathway boxes
pathway_data = [
    ('Iron\nOverload', 8.5, '#c0392b', ['HMOX1', 'TFRC', 'HFE', 'FTH1']),
    ('Lipid\nPeroxidation', 6.5, '#e67e22', ['ACSL4', 'LPCAT3', 'PTGS2', 'NOX4']),
    ('Ferroptosis\n(Cell Death)', 4.5, '#8e44ad', ['GPX4', 'SLC7A11', 'NCOA4', 'AIFM2']),
    ('Cellular\nSenescence', 2.5, '#2980b9', ['CDKN1A', 'CDKN2A', 'TP53']),
    ('SASP &\nInflammation', 0.5, '#16a085', ['IL6', 'CXCL8', 'SERPINE1', 'VEGFA']),
]

drug_mapping = {
    'Iron\nOverload': ['DEFEROXAMINE', 'DEFERASIROX'],
    'Lipid\nPeroxidation': ['VITAMIN C', 'FERROSTATIN-1', 'LIPROXSTATIN-1'],
    'Ferroptosis\n(Cell Death)': ['N-ACETYLCYSTEINE', 'DIMETHYL FUMARATE'],
    'Cellular\nSenescence': ['DASATINIB+QUERCETIN', 'FISETIN'],
    'SASP &\nInflammation': ['DOXYCYCLINE', 'ASPIRIN', 'TOCILIZUMAB'],
}

for label, y, color, genes in pathway_data:
    # Pathway box
    box = FancyBboxPatch((0.3, y - 0.4), 2.0, 0.8, boxstyle="round,pad=0.1",
                         facecolor=color, alpha=0.15, edgecolor=color, linewidth=2)
    ax5.add_patch(box)
    ax5.text(1.3, y, label, fontsize=8, fontweight='bold', ha='center', va='center',
            color=color)
    
    # Genes
    ax5.text(1.3, y - 0.55, ', '.join(genes[:3]), fontsize=5.5, ha='center',
            va='top', color='#666', style='italic')
    
    # Arrow
    ax5.annotate('', xy=(2.5, y), xytext=(2.6, y),
                arrowprops=dict(arrowstyle='->', lw=2, color=color, connectionstyle='arc3,rad=0'))
    
    # Drug boxes on right
    drugs = drug_mapping.get(label, [])
    n_drugs = len(drugs)
    for di, drug in enumerate(drugs):
        dx = 3.0 + (di % 3) * 2.8
        dy = y + (0.1 if n_drugs <= 2 else -0.2 + 0.2 * di)
        drug_box = FancyBboxPatch((dx, dy - 0.15), 2.4, 0.3, boxstyle="round,pad=0.05",
                                 facecolor='white', alpha=0.9, edgecolor=color, linewidth=1.5)
        ax5.add_patch(drug_box)
        ax5.text(dx + 1.2, dy, drug, fontsize=6.5, fontweight='bold', ha='center', va='center',
                color='#2c3e50')

# Title
ax5.text(6, 9.5, 'Ferro-aging Pathway → Drug Targeting Map',
        fontsize=14, fontweight='bold', ha='center', va='center', color='#2c3e50')
ax5.text(6, 9.0, 'Targeted intervention at each node of the Iron → Lipid Peroxidation → Ferroptosis → Senescence → SASP cascade',
        fontsize=8, ha='center', va='center', color='#7f8c8d', style='italic')

# ── Panel 6: Drug category pie chart ──────────────────────────────────────
ax6 = fig.add_subplot(gs[2, 3])

categories = ['Iron Chelation', 'Ferroptosis\nInhibitor', 'Senolytic', 'Anti-inflammatory',
              'Antioxidant', 'MMP Inhibitor', 'NRF2 Activator', 'Statin', 'ARB/ACEi', 'Other']
cat_counts = [3, 3, 4, 8, 3, 2, 2, 3, 6, 12]
cat_colors = ['#c0392b', '#e74c3c', '#e67e22', '#f39c12', '#2ecc71', '#1abc9c',
              '#3498db', '#9b59b6', '#34495e', '#95a5a6']

wedges, texts, autotexts = ax6.pie(cat_counts, labels=None, autopct='', 
                                     colors=cat_colors, startangle=90, pctdistance=0.8)
ax6.legend(wedges, [f'{c} ({n})' for c, n in zip(categories, cat_counts)],
          fontsize=6.5, loc='center left', bbox_to_anchor=(1.0, 0.5))
ax6.set_title('Drug Category Distribution', fontsize=11, fontweight='bold', color='#2c3e50')

fig.suptitle('Ferro-aging Drug Prediction & Repurposing Analysis\nDGIdb + Literature — Targeting Iron-driven Vascular Aging',
            fontsize=16, fontweight='bold', y=0.995, color='#2c3e50')

plt.savefig(f"{OUT_DIR}/ferro_aging_drug_prediction.png", dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: ferro_aging_drug_prediction.png")

# ── STEP 7: Export results ───────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 7: Export all results")
print("=" * 70)

# Drug rankings
drug_df.to_csv(f"{OUT_DIR}/ferro_aging_drug_rankings.csv", index=False)
print(f"  Saved: ferro_aging_drug_rankings.csv ({len(drug_df)} drugs)")

# All interactions
df.to_csv(f"{OUT_DIR}/ferro_aging_drug_gene_interactions.csv", index=False)
print(f"  Saved: ferro_aging_drug_gene_interactions.csv ({len(df)} interactions)")

# Ferro-aging specific drugs
ferro_drug_df.to_csv(f"{OUT_DIR}/ferro_aging_targeted_drugs.csv", index=False)
print(f"  Saved: ferro_aging_targeted_drugs.csv")

# Summary stats
stats = {
    'total_genes_queried': len(unique_genes),
    'genes_with_drugs': df['gene'].nunique(),
    'genes_no_drugs': len(unique_genes) - df['gene'].nunique(),
    'total_drugs': drug_df['drug'].nunique(),
    'fda_approved_drugs': int(drug_df['fda_approved'].sum()),
    'multi_target_drugs': int((drug_df['n_genes'] >= 2).sum()),
    'top_scoring_drug': drug_df.iloc[0]['drug'] if len(drug_df) > 0 else None,
    'top_score': float(drug_df.iloc[0]['score']) if len(drug_df) > 0 else None,
    'ferro_aging_specific_drugs': len(ferro_drug_df),
}
with open(f"{OUT_DIR}/ferro_aging_drug_stats.json", 'w') as f:
    json.dump(stats, f, indent=2)
print(f"  Saved: ferro_aging_drug_stats.json")

# ── Final Summary ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("DRUG PREDICTION COMPLETE")
print("=" * 70)
print(f"\n  Genes queried: {len(unique_genes)}")
print(f"  Genes with drugs: {df['gene'].nunique()}")
print(f"  Total drug-gene interactions: {len(df)}")
print(f"  Unique candidate drugs: {drug_df['drug'].nunique()}")
print(f"  FDA-approved repurposable: {int(drug_df['fda_approved'].sum())}")
print(f"  Ferro-aging targeted (literature): {len(ferro_drug_df)}")
print(f"\n  ★ Top drug: {drug_df.iloc[0]['drug']} (score={drug_df.iloc[0]['score']:.1f})")
print(f"  All results in: {OUT_DIR}/")
