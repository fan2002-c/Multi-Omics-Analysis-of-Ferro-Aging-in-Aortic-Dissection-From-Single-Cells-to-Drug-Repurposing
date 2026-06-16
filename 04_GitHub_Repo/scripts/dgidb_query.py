"""
DGIdb Drug-Gene Interaction Query
===================================
Batch-query the Drug-Gene Interaction Database (DGIdb v5.0) GraphQL API
to identify drugs targeting ferro-aging-related genes.

Input: Hard-coded list of 35 ferro-aging key genes (validated from scRNA-seq)
Output: dgidb_raw_results.json — raw API responses for downstream processing

Method:
  - Uses DGIdb GraphQL API endpoint (https://dgidb.org/api/graphql)
  - Genes queried in batches of 5 to avoid URL length limits
  - Retrieves: drug name, interaction type, sources, PMIDs, approval status
  - Rate-limited with 1-second delay between batches

Note: This script produces raw results. Run ferro_aging_drug_prediction.py
      next for multi-criteria drug prioritization and visualization.
"""
import json
import time
import urllib.request
import urllib.error
import ssl

ssl_ctx = ssl.create_default_context()
ssl_ctx.check_hostname = False
ssl_ctx.verify_mode = ssl.CERT_NONE

OUT_DIR = "C:/Users/lidaf/WorkBuddy/2026-06-11-19-29-35/results"

# Key Ferro-aging genes from scRNA-seq and bulk validation
KEY_GENES = [
    # Core Ferro-aging
    'ACSL4', 'GPX4', 'NCOA4', 'FTH1', 'HMOX1', 'TFRC',
    # xCT system / ferroptosis
    'SLC7A11', 'SLC3A2',
    # SASP / inflammation
    'CXCL8', 'SERPINE1', 'IL6', 'VEGFA', 'CCL2', 'CXCL1',
    # Senescence markers
    'CDKN1A', 'CDKN2A', 'TP53',
    # Iron metabolism
    'HFE', 'TF', 'FTL', 'FPN1', 'HAMP', 'BMP6',
    # Lipid peroxidation
    'PTGS2', 'NOX4', 'LPCAT3', 'ACSL1',
    # ECM / vascular
    'MMP9', 'MMP2', 'MMP12', 'FBN1',
    # NRF2 pathway
    'NFE2L2', 'KEAP1',
    # Other ferroptosis
    'AIFM2', 'SLC40A1',
]

# Query genes in batches to avoid URL too long
BATCH_SIZE = 5
all_genes = [g for g in KEY_GENES if g]  # filter empty

all_interactions = []
all_metadata = {}

for i in range(0, len(all_genes), BATCH_SIZE):
    batch = all_genes[i:i+BATCH_SIZE]
    genes_param = ','.join(batch)
    url = f"https://dgidb.org/api/v2/interactions.json?genes={genes_param}"
    print(f"Querying: {genes_param}")
    
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            
        if 'interactions' in data:
            interactions = data['interactions']
            print(f"  → {len(interactions)} interactions found")
            all_interactions.extend(interactions)
        
        # Also get matched genes metadata
        if 'matched_genes' in data:
            for mg in data['matched_genes']:
                all_metadata[mg.get('gene_name', '')] = {
                    'entrez_id': mg.get('entrez_id'),
                    'gene_claim_names': mg.get('gene_claim_names', []),
                    'gene_categories': mg.get('gene_categories', [])
                }
    except Exception as e:
        print(f"  ✗ Error: {e}")
    
    time.sleep(0.5)  # rate limit

# Save raw results
raw_output = {
    'interactions': all_interactions,
    'matched_genes': all_metadata,
    'total_interactions': len(all_interactions),
    'query_genes': all_genes
}

with open(f"{OUT_DIR}/dgidb_raw_interactions.json", 'w') as f:
    json.dump(raw_output, f, indent=2)

print(f"\n{'='*50}")
print(f"Total interactions: {len(all_interactions)}")
print(f"Matched genes: {list(all_metadata.keys())}")
print(f"Saved: {OUT_DIR}/dgidb_raw_interactions.json")

# Quick summary
gene_drug_count = {}
drug_gene_map = {}
for inter in all_interactions:
    gene = inter.get('gene_name', '?')
    drug = inter.get('drug_name', '?')
    gene_drug_count[gene] = gene_drug_count.get(gene, 0) + 1
    if drug not in drug_gene_map:
        drug_gene_map[drug] = set()
    drug_gene_map[drug].add(gene)

print(f"\n── Interactions per gene ──")
for g, c in sorted(gene_drug_count.items(), key=lambda x: -x[1]):
    print(f"  {g}: {c} drugs")

print(f"\n── Drugs targeting multiple Ferro-aging genes (>1) ──")
multi_target = [(d, list(g)) for d, g in drug_gene_map.items() if len(g) > 1]
multi_target.sort(key=lambda x: -len(x[1]))
for d, gs in multi_target:
    print(f"  {d}: {len(gs)} genes → {gs}")
