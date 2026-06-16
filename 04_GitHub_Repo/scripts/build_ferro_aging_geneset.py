#!/usr/bin/env python3
"""
Ferro-aging Gene Set Builder
============================
Build a comprehensive ferro-aging gene set by integrating:
1. Core ferro-aging genes (Cell Metabolism 2026)
2. Iron metabolism genes (KEGG/GO/literature)
3. Ferroptosis regulators (FerrDb V2/V3 literature)
4. Cellular senescence markers (CellAge/literature)
5. Lipid peroxidation genes
6. Nrf2 antioxidant pathway

Output: ferro_aging_geneset.csv with categories, gene symbols, and functional annotations.
"""

import csv
import json
from pathlib import Path

# ============================================================
# 1. CORE FERRO-AGING GENES (Cell Metabolism 2026, Liu et al.)
# ============================================================
CORE_FERRO_AGING = {
    "ACSL4": {
        "full_name": "Acyl-CoA Synthetase Long Chain Family Member 4",
        "role": "Core executor of ferro-aging; catalyzes PUFA esterification; direct target of Vitamin C",
        "category": "Core_FerroAging",
        "subcategory": "Master_Regulator",
        "evidence": "Cell Metabolism 2026 - Identified as master regulator of ferro-aging axis in primates",
        "VC_target": True,
    }
}

# ============================================================
# 2. IRON METABOLISM GENES
# ============================================================
IRON_METABOLISM = {
    # Iron Uptake & Transport
    "TFRC": {
        "full_name": "Transferrin Receptor",
        "role": "Iron uptake; receptor for transferrin-mediated iron import",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Uptake",
    },
    "TFR2": {
        "full_name": "Transferrin Receptor 2",
        "role": "Iron sensing in hepatocytes; hepcidin regulation",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Uptake",
    },
    "SLC11A2": {
        "full_name": "Solute Carrier Family 11 Member 2 (DMT1)",
        "role": "Divalent metal transporter 1; endosomal iron export to cytoplasm",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Transport",
    },
    "SLC40A1": {
        "full_name": "Solute Carrier Family 40 Member 1 (Ferroportin)",
        "role": "Only known iron exporter; exports iron from cells into circulation",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Transport",
    },
    "STEAP1": {
        "full_name": "STEAP Family Member 1",
        "role": "Ferric reductase; reduces Fe3+ to Fe2+ for DMT1 transport",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Reduction",
    },
    "STEAP2": {
        "full_name": "STEAP Family Member 2",
        "role": "Ferric reductase; iron reduction",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Reduction",
    },
    "STEAP3": {
        "full_name": "STEAP Family Member 3",
        "role": "Endosomal ferric reductase; major player in erythroid iron uptake",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Reduction",
    },
    # Iron Storage
    "FTH1": {
        "full_name": "Ferritin Heavy Chain 1",
        "role": "Iron storage; ferroxidase activity; oxidizes Fe2+ to Fe3+ for storage",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Storage",
    },
    "FTL": {
        "full_name": "Ferritin Light Chain",
        "role": "Iron storage; forms ferritin complex with FTH1",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Storage",
    },
    "NCOA4": {
        "full_name": "Nuclear Receptor Coactivator 4",
        "role": "Ferritinophagy receptor; mediates ferritin degradation and iron release",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Storage_Turnover",
    },
    # Iron Regulation
    "HAMP": {
        "full_name": "Hepcidin Antimicrobial Peptide",
        "role": "Master iron regulator; binds ferroportin and triggers its degradation",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Regulation",
    },
    "HFE": {
        "full_name": "Homeostatic Iron Regulator",
        "role": "HFE protein; regulates hepcidin via BMP/SMAD pathway",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Regulation",
    },
    "HJV": {
        "full_name": "Hemojuvelin BMP Co-receptor (HFE2)",
        "role": "BMP co-receptor; activates hepcidin transcription",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Regulation",
    },
    "BMP6": {
        "full_name": "Bone Morphogenetic Protein 6",
        "role": "Key endogenous hepcidin regulator; iron sensing",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Regulation",
    },
    "TMPRSS6": {
        "full_name": "Transmembrane Serine Protease 6 (Matriptase-2)",
        "role": "Cleaves HJV; suppresses hepcidin; promotes iron absorption",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Regulation",
    },
    "B2M": {
        "full_name": "Beta-2-Microglobulin",
        "role": "Component of HFE complex; needed for HFE cell surface expression",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Regulation",
    },
    "CP": {
        "full_name": "Ceruloplasmin",
        "role": "Ferroxidase; oxidizes Fe2+ to Fe3+ for transferrin loading",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Regulation",
    },
    # Iron Response
    "IREB2": {
        "full_name": "Iron Responsive Element Binding Protein 2 (IRP2)",
        "role": "Cytosolic iron sensor; regulates TFRC/FTH1 mRNA via IRE binding",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Sensing",
    },
    "ACO1": {
        "full_name": "Aconitase 1 (IRP1)",
        "role": "Cytosolic iron sensor; dual function aconitase/IRE-binding protein",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Sensing",
    },
    "FBXL5": {
        "full_name": "F-Box And Leucine Rich Repeat Protein 5",
        "role": "Iron-sensing E3 ubiquitin ligase; degrades IRP2 in iron-replete conditions",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Sensing",
    },
    # Iron-Sulfur Cluster
    "ISCU": {
        "full_name": "Iron-Sulfur Cluster Assembly Enzyme",
        "role": "Iron-sulfur cluster scaffold protein",
        "category": "Iron_Metabolism",
        "subcategory": "FeS_Cluster",
    },
    "FXN": {
        "full_name": "Frataxin",
        "role": "Iron chaperone for Fe-S cluster biogenesis; mitochondrial iron homeostasis",
        "category": "Iron_Metabolism",
        "subcategory": "FeS_Cluster",
    },
    "GLRX5": {
        "full_name": "Glutaredoxin 5",
        "role": "Fe-S cluster assembly; mitochondrial iron-sulfur protein maturation",
        "category": "Iron_Metabolism",
        "subcategory": "FeS_Cluster",
    },
    "NFS1": {
        "full_name": "NFS1 Cysteine Desulfurase",
        "role": "Sulfur donor for Fe-S cluster biogenesis",
        "category": "Iron_Metabolism",
        "subcategory": "FeS_Cluster",
    },
    # Heme Metabolism
    "HMOX1": {
        "full_name": "Heme Oxygenase 1",
        "role": "Heme degradation; releases free iron; key antioxidant and ferroptosis player",
        "category": "Iron_Metabolism",
        "subcategory": "Heme_Metabolism",
    },
    "HMOX2": {
        "full_name": "Heme Oxygenase 2",
        "role": "Constitutive heme oxygenase; iron release from heme",
        "category": "Iron_Metabolism",
        "subcategory": "Heme_Metabolism",
    },
    "BLVRA": {
        "full_name": "Biliverdin Reductase A",
        "role": "Converts biliverdin to bilirubin; antioxidant",
        "category": "Iron_Metabolism",
        "subcategory": "Heme_Metabolism",
    },
    "BLVRB": {
        "full_name": "Biliverdin Reductase B",
        "role": "Biliverdin reductase; flavin reductase activity",
        "category": "Iron_Metabolism",
        "subcategory": "Heme_Metabolism",
    },
    "ALAS1": {
        "full_name": "5'-Aminolevulinate Synthase 1",
        "role": "Rate-limiting enzyme in heme biosynthesis (ubiquitous)",
        "category": "Iron_Metabolism",
        "subcategory": "Heme_Biosynthesis",
    },
    "ALAS2": {
        "full_name": "5'-Aminolevulinate Synthase 2",
        "role": "Erythroid-specific heme biosynthesis",
        "category": "Iron_Metabolism",
        "subcategory": "Heme_Biosynthesis",
    },
    "FECH": {
        "full_name": "Ferrochelatase",
        "role": "Final step of heme biosynthesis; inserts Fe2+ into protoporphyrin IX",
        "category": "Iron_Metabolism",
        "subcategory": "Heme_Biosynthesis",
    },
    "SLC48A1": {
        "full_name": "Solute Carrier Family 48 Member 1 (HRG1)",
        "role": "Heme transporter; heme recycling",
        "category": "Iron_Metabolism",
        "subcategory": "Heme_Transport",
    },
    # Lactoferrin / Transferrin
    "TF": {
        "full_name": "Transferrin",
        "role": "Iron carrier protein in blood; delivers iron to cells via TFRC",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Carrier",
    },
    "LTF": {
        "full_name": "Lactotransferrin (Lactoferrin)",
        "role": "Iron-binding glycoprotein; antimicrobial and anti-inflammatory",
        "category": "Iron_Metabolism",
        "subcategory": "Iron_Carrier",
    },
}

# ============================================================
# 3. FERROPTOSIS REGULATORS (from FerrDb literature)
# ============================================================
FERROPTOSIS = {
    # Drivers (promote ferroptosis)
    "LPCAT3": {
        "full_name": "Lysophosphatidylcholine Acyltransferase 3",
        "role": "Ferroptosis driver; incorporates PUFAs into membrane phospholipids",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "ALOX5": {
        "full_name": "Arachidonate 5-Lipoxygenase",
        "role": "Ferroptosis driver; catalyzes lipid peroxidation",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "ALOX12": {
        "full_name": "Arachidonate 12-Lipoxygenase",
        "role": "Ferroptosis driver; lipid peroxidation in p53-mediated ferroptosis",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "ALOX15": {
        "full_name": "Arachidonate 15-Lipoxygenase",
        "role": "Ferroptosis driver; PUFA oxygenation",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "ALOX15B": {
        "full_name": "Arachidonate 15-Lipoxygenase Type B",
        "role": "Ferroptosis driver; lipid peroxidation",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "PTGS2": {
        "full_name": "Prostaglandin-Endoperoxide Synthase 2 (COX-2)",
        "role": "Ferroptosis marker/driver; induced by ferroptotic stimuli",
        "category": "Ferroptosis",
        "subcategory": "Driver_Marker",
    },
    "NOX1": {
        "full_name": "NADPH Oxidase 1",
        "role": "Ferroptosis driver; ROS production",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "NOX4": {
        "full_name": "NADPH Oxidase 4",
        "role": "Ferroptosis driver; constitutive ROS production",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "CYBB": {
        "full_name": "Cytochrome B-245 Beta Chain (NOX2)",
        "role": "Ferroptosis driver; phagocyte ROS production",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "VDAC2": {
        "full_name": "Voltage Dependent Anion Channel 2",
        "role": "Ferroptosis driver; erastin target; mitochondrial outer membrane",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "VDAC3": {
        "full_name": "Voltage Dependent Anion Channel 3",
        "role": "Ferroptosis driver; erastin target",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "SAT1": {
        "full_name": "Spermidine/Spermine N1-Acetyltransferase 1",
        "role": "Ferroptosis driver; polyamine metabolism; p53 target",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "CARS1": {
        "full_name": "Cysteinyl-tRNA Synthetase 1",
        "role": "Ferroptosis driver; transsulfuration pathway",
        "category": "Ferroptosis",
        "subcategory": "Driver",
    },
    "CHAC1": {
        "full_name": "ChaC Glutathione Specific Gamma-Glutamylcyclotransferase 1",
        "role": "Ferroptosis marker; degrades glutathione; ER stress-induced",
        "category": "Ferroptosis",
        "subcategory": "Marker",
    },
    # Suppressors (inhibit ferroptosis)
    "GPX4": {
        "full_name": "Glutathione Peroxidase 4",
        "role": "Master ferroptosis suppressor; reduces lipid hydroperoxides",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "SLC7A11": {
        "full_name": "Solute Carrier Family 7 Member 11 (xCT)",
        "role": "Cystine/glutamate antiporter subunit; critical for GSH synthesis",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "SLC3A2": {
        "full_name": "Solute Carrier Family 3 Member 2 (4F2hc)",
        "role": "Heavy chain of xCT cystine transporter complex",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "AIFM2": {
        "full_name": "Apoptosis Inducing Factor Mitochondria Associated 2 (FSP1)",
        "role": "Ferroptosis suppressor protein 1; CoQ10-mediated lipid radical trapping; parallel to GPX4",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "DHODH": {
        "full_name": "Dihydroorotate Dehydrogenase (Quinone)",
        "role": "Mitochondrial ferroptosis suppressor; CoQ10 regeneration in mitochondria",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "GCH1": {
        "full_name": "GTP Cyclohydrolase 1",
        "role": "Ferroptosis suppressor; BH4 biosynthesis; lipid radical trapping",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "NFE2L2": {
        "full_name": "NFE2 Like BZIP Transcription Factor 2 (Nrf2)",
        "role": "Master antioxidant transcription factor; upregulates ferroptosis suppressors",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "HSPB1": {
        "full_name": "Heat Shock Protein Family B Member 1 (HSP27)",
        "role": "Ferroptosis suppressor; reduces iron uptake via TFRC",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "HSPA5": {
        "full_name": "Heat Shock Protein Family A Member 5 (BiP/GRP78)",
        "role": "ER stress-mediated ferroptosis resistance; binds GPX4",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "PROM2": {
        "full_name": "Prominin 2",
        "role": "Ferroptosis suppressor; promotes ferritin export in exosomes",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "CISD1": {
        "full_name": "CDGSH Iron Sulfur Domain 1 (MitoNEET)",
        "role": "Mitochondrial iron-sulfur protein; inhibits mitochondrial lipid peroxidation",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "CISD2": {
        "full_name": "CDGSH Iron Sulfur Domain 2 (NAF-1)",
        "role": "Iron-sulfur protein; ferroptosis suppressor",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "MT1G": {
        "full_name": "Metallothionein 1G",
        "role": "Metal-binding protein; ferroptosis suppressor",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "OTUB1": {
        "full_name": "OTU Deubiquitinase, Ubiquitin Aldehyde Binding 1",
        "role": "Deubiquitinates SLC7A11; stabilizes xCT; ferroptosis suppressor",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "NF2": {
        "full_name": "NF2, Moesin-Ezrin-Radixin Like (Merlin)",
        "role": "Ferroptosis suppressor via YAP pathway",
        "category": "Ferroptosis",
        "subcategory": "Suppressor",
    },
    "TP53": {
        "full_name": "Tumor Protein P53",
        "role": "Dual role: inhibits SLC7A11 → pro-ferroptosis; but can also suppress ferroptosis",
        "category": "Ferroptosis",
        "subcategory": "Dual_Regulator",
    },
}

# ============================================================
# 4. LIPID PEROXIDATION GENES
# ============================================================
LIPID_PEROXIDATION = {
    "ALOXE3": {
        "full_name": "Arachidonate Lipoxygenase 3",
        "role": "Epidermal-type lipoxygenase; lipid peroxidation",
        "category": "Lipid_Peroxidation",
        "subcategory": "Lipoxygenase",
    },
    "ALOX12B": {
        "full_name": "Arachidonate 12-Lipoxygenase, 12R Type",
        "role": "Lipid peroxidation enzyme",
        "category": "Lipid_Peroxidation",
        "subcategory": "Lipoxygenase",
    },
    "PLA2G6": {
        "full_name": "Phospholipase A2 Group VI",
        "role": "Releases oxidized fatty acids from membranes; phospholipid remodeling",
        "category": "Lipid_Peroxidation",
        "subcategory": "Phospholipase",
    },
    "NOX5": {
        "full_name": "NADPH Oxidase 5",
        "role": "Calcium-dependent ROS production",
        "category": "Lipid_Peroxidation",
        "subcategory": "ROS_Production",
    },
}

# ============================================================
# 5. CELLULAR SENESCENCE MARKERS (CellAge/literature)
# ============================================================
SENESCENCE = {
    # Cell Cycle Arrest
    "CDKN1A": {
        "full_name": "Cyclin Dependent Kinase Inhibitor 1A (p21)",
        "role": "Classical senescence marker; CDK inhibitor; mediates cell cycle arrest",
        "category": "Senescence",
        "subcategory": "Cell_Cycle_Arrest",
    },
    "CDKN2A": {
        "full_name": "Cyclin Dependent Kinase Inhibitor 2A (p16/INK4a)",
        "role": "Classical senescence marker; CDK4/6 inhibitor; aging biomarker",
        "category": "Senescence",
        "subcategory": "Cell_Cycle_Arrest",
    },
    "CDKN2B": {
        "full_name": "Cyclin Dependent Kinase Inhibitor 2B (p15)",
        "role": "CDK inhibitor; TGF-beta induced cell cycle arrest",
        "category": "Senescence",
        "subcategory": "Cell_Cycle_Arrest",
    },
    "RB1": {
        "full_name": "RB Transcriptional Corepressor 1",
        "role": "Retinoblastoma protein; key senescence effector; E2F inhibition",
        "category": "Senescence",
        "subcategory": "Cell_Cycle_Arrest",
    },
    "E2F1": {
        "full_name": "E2F Transcription Factor 1",
        "role": "Cell cycle transcription factor; suppressed by RB in senescence",
        "category": "Senescence",
        "subcategory": "Cell_Cycle",
    },
    # DNA Damage / Epigenetic
    "H2AX": {
        "full_name": "H2A.X Variant Histone (H2AFX)",
        "role": "DNA damage marker (gamma-H2AX); persistent foci in senescence",
        "category": "Senescence",
        "subcategory": "DNA_Damage",
    },
    "HMGA1": {
        "full_name": "High Mobility Group AT-Hook 1",
        "role": "Senescence-associated chromatin remodeling protein",
        "category": "Senescence",
        "subcategory": "Chromatin",
    },
    "HMGA2": {
        "full_name": "High Mobility Group AT-Hook 2",
        "role": "Senescence-associated chromatin remodeling; SASP regulator",
        "category": "Senescence",
        "subcategory": "Chromatin",
    },
    "LMNB1": {
        "full_name": "Lamin B1",
        "role": "Nuclear lamina component; DOWNREGULATED in senescence (negative marker)",
        "category": "Senescence",
        "subcategory": "Nuclear_Lamina",
    },
    # SASP (Senescence-Associated Secretory Phenotype)
    "IL6": {
        "full_name": "Interleukin 6",
        "role": "Classical SASP factor; pro-inflammatory; vascular inflammation",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "CXCL8": {
        "full_name": "C-X-C Motif Chemokine Ligand 8 (IL-8)",
        "role": "SASP chemokine; neutrophil recruitment; angiogenesis",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "IL1A": {
        "full_name": "Interleukin 1 Alpha",
        "role": "SASP factor; upstream regulator of IL6/IL8 in senescence",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "IL1B": {
        "full_name": "Interleukin 1 Beta",
        "role": "SASP factor; inflammasome-activated; pro-inflammatory",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "CCL2": {
        "full_name": "C-C Motif Chemokine Ligand 2 (MCP-1)",
        "role": "SASP chemokine; monocyte recruitment; vascular disease",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "CXCL1": {
        "full_name": "C-X-C Motif Chemokine Ligand 1 (GRO-alpha)",
        "role": "SASP chemokine; neutrophil chemoattractant",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "CXCL2": {
        "full_name": "C-X-C Motif Chemokine Ligand 2 (GRO-beta)",
        "role": "SASP chemokine",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "MMP1": {
        "full_name": "Matrix Metallopeptidase 1 (Collagenase-1)",
        "role": "SASP protease; ECM degradation; tissue remodeling",
        "category": "Senescence",
        "subcategory": "SASP_Protease",
    },
    "MMP3": {
        "full_name": "Matrix Metallopeptidase 3 (Stromelysin-1)",
        "role": "SASP protease; ECM degradation; activates other MMPs",
        "category": "Senescence",
        "subcategory": "SASP_Protease",
    },
    "MMP9": {
        "full_name": "Matrix Metallopeptidase 9 (Gelatinase B)",
        "role": "SASP protease; basement membrane degradation; aneurysm pathogenesis",
        "category": "Senescence",
        "subcategory": "SASP_Protease",
    },
    "MMP12": {
        "full_name": "Matrix Metallopeptidase 12 (Macrophage Elastase)",
        "role": "SASP protease; elastin degradation; aneurysm/emphysema",
        "category": "Senescence",
        "subcategory": "SASP_Protease",
    },
    "SERPINE1": {
        "full_name": "Serpin Family E Member 1 (PAI-1)",
        "role": "SASP factor; key senescence mediator; thrombosis/inflammation link",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "VEGFA": {
        "full_name": "Vascular Endothelial Growth Factor A",
        "role": "SASP angiogenic factor; vascular permeability",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "ICAM1": {
        "full_name": "Intercellular Adhesion Molecule 1",
        "role": "SASP adhesion molecule; leukocyte recruitment to vessels",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    "CSF2": {
        "full_name": "Colony Stimulating Factor 2 (GM-CSF)",
        "role": "SASP factor; myeloid cell differentiation and activation",
        "category": "Senescence",
        "subcategory": "SASP",
    },
    # SA-beta-galactosidase
    "GLB1": {
        "full_name": "Galactosidase Beta 1",
        "role": "SA-beta-galactosidase; classical senescence marker enzyme",
        "category": "Senescence",
        "subcategory": "Marker",
    },
    # Anti-senescence / Longevity
    "SIRT1": {
        "full_name": "Sirtuin 1",
        "role": "Anti-senescence; NAD+-dependent deacetylase; FOXO activation",
        "category": "Senescence",
        "subcategory": "Anti_Senescence",
    },
    "SIRT6": {
        "full_name": "Sirtuin 6",
        "role": "Anti-senescence; chromatin regulation; DNA repair; vascular aging protection",
        "category": "Senescence",
        "subcategory": "Anti_Senescence",
    },
    "TERT": {
        "full_name": "Telomerase Reverse Transcriptase",
        "role": "Telomere maintenance; anti-senescence",
        "category": "Senescence",
        "subcategory": "Anti_Senescence",
    },
    "KL": {
        "full_name": "Klotho",
        "role": "Anti-aging protein; vascular protection; suppresses senescence",
        "category": "Senescence",
        "subcategory": "Anti_Senescence",
    },
    "FOXO3": {
        "full_name": "Forkhead Box O3",
        "role": "Longevity-associated transcription factor; antioxidant defense",
        "category": "Senescence",
        "subcategory": "Anti_Senescence",
    },
}

# ============================================================
# 6. NRF2 ANTIOXIDANT PATHWAY
# ============================================================
NRF2_PATHWAY = {
    "KEAP1": {
        "full_name": "Kelch Like ECH Associated Protein 1",
        "role": "Nrf2 repressor; redox sensor; targets Nrf2 for degradation",
        "category": "NRF2_Pathway",
        "subcategory": "Regulator",
    },
    "NQO1": {
        "full_name": "NAD(P)H Quinone Dehydrogenase 1",
        "role": "Nrf2 target; quinone detoxification; antioxidant",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "GCLC": {
        "full_name": "Glutamate-Cysteine Ligase Catalytic Subunit",
        "role": "Rate-limiting enzyme for GSH synthesis; Nrf2 target",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "GCLM": {
        "full_name": "Glutamate-Cysteine Ligase Modifier Subunit",
        "role": "GSH synthesis; Nrf2 target",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "GSR": {
        "full_name": "Glutathione-Disulfide Reductase",
        "role": "Recycles oxidized glutathione; Nrf2 target",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "TXN": {
        "full_name": "Thioredoxin",
        "role": "Redox protein; Nrf2 target; reduces oxidized proteins",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "TXNRD1": {
        "full_name": "Thioredoxin Reductase 1",
        "role": "Reduces thioredoxin; Nrf2 target; selenium-containing",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "PRDX1": {
        "full_name": "Peroxiredoxin 1",
        "role": "Hydrogen peroxide reduction; Nrf2 target",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "SOD1": {
        "full_name": "Superoxide Dismutase 1",
        "role": "Cytosolic superoxide detoxification; copper/zinc SOD",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "SOD2": {
        "full_name": "Superoxide Dismutase 2 (MnSOD)",
        "role": "Mitochondrial superoxide detoxification",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "CAT": {
        "full_name": "Catalase",
        "role": "H2O2 detoxification; antioxidant enzyme",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "GSTP1": {
        "full_name": "Glutathione S-Transferase Pi 1",
        "role": "Detoxification; Nrf2 target; conjugates GSH to electrophiles",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "GSTM1": {
        "full_name": "Glutathione S-Transferase Mu 1",
        "role": "Detoxification; Nrf2 target",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "AKR1C1": {
        "full_name": "Aldo-Keto Reductase Family 1 Member C1",
        "role": "Nrf2 target; detoxification of lipid peroxidation products",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
    "AKR1B10": {
        "full_name": "Aldo-Keto Reductase Family 1 Member B10",
        "role": "Nrf2 target; detoxification of reactive carbonyls",
        "category": "NRF2_Pathway",
        "subcategory": "Target",
    },
}

# ============================================================
# 7. VASCULAR-SPECIFIC FERRO-AGING RELEVANT GENES
# ============================================================
VASCULAR_RELEVANT = {
    "VCAM1": {
        "full_name": "Vascular Cell Adhesion Molecule 1",
        "role": "Endothelial activation marker; iron-induced expression; vascular inflammation",
        "category": "Vascular_Relevant",
        "subcategory": "Endothelial",
    },
    "SELE": {
        "full_name": "Selectin E (E-Selectin)",
        "role": "Endothelial adhesion molecule; vascular inflammation",
        "category": "Vascular_Relevant",
        "subcategory": "Endothelial",
    },
    "NOS3": {
        "full_name": "Nitric Oxide Synthase 3 (eNOS)",
        "role": "Endothelial function; iron overload impairs eNOS; vascular protection",
        "category": "Vascular_Relevant",
        "subcategory": "Endothelial",
    },
    "EDN1": {
        "full_name": "Endothelin 1",
        "role": "Vasoconstriction; endothelial dysfunction marker; iron-induced",
        "category": "Vascular_Relevant",
        "subcategory": "Endothelial",
    },
    "ACTA2": {
        "full_name": "Actin Alpha 2, Smooth Muscle (alpha-SMA)",
        "role": "VSMC contractile marker; decreases in dedifferentiated/senescent VSMCs",
        "category": "Vascular_Relevant",
        "subcategory": "VSMC",
    },
    "MYH11": {
        "full_name": "Myosin Heavy Chain 11",
        "role": "VSMC-specific contractile marker; aneurysm-related mutations",
        "category": "Vascular_Relevant",
        "subcategory": "VSMC",
    },
    "TAGLN": {
        "full_name": "Transgelin (SM22-alpha)",
        "role": "VSMC marker; decreases in aneurysm",
        "category": "Vascular_Relevant",
        "subcategory": "VSMC",
    },
    "CNN1": {
        "full_name": "Calponin 1",
        "role": "VSMC differentiation marker; downregulated in aneurysm",
        "category": "Vascular_Relevant",
        "subcategory": "VSMC",
    },
    "COL1A1": {
        "full_name": "Collagen Type I Alpha 1 Chain",
        "role": "ECM structural protein; degraded in aneurysm by SASP MMPs",
        "category": "Vascular_Relevant",
        "subcategory": "ECM",
    },
    "COL3A1": {
        "full_name": "Collagen Type III Alpha 1 Chain",
        "role": "Vascular ECM; aneurysm-related mutations",
        "category": "Vascular_Relevant",
        "subcategory": "ECM",
    },
    "ELN": {
        "full_name": "Elastin",
        "role": "Elastic fiber; degraded in aneurysm; target of MMP12 from senescent cells",
        "category": "Vascular_Relevant",
        "subcategory": "ECM",
    },
    "FBN1": {
        "full_name": "Fibrillin 1",
        "role": "Microfibril component; Marfan/aneurysm gene; TGF-beta regulation",
        "category": "Vascular_Relevant",
        "subcategory": "ECM",
    },
    "TGFB1": {
        "full_name": "Transforming Growth Factor Beta 1",
        "role": "Aneurysm pathogenesis; VSMC phenotype; senescence regulator",
        "category": "Vascular_Relevant",
        "subcategory": "Signaling",
    },
    "IL6": {
        "full_name": "Interleukin 6",
        "role": "Vascular inflammation; aneurysm biomarker; SASP factor",
        "category": "Vascular_Relevant",
        "subcategory": "Inflammation",
    },
    "CRP": {
        "full_name": "C-Reactive Protein",
        "role": "Systemic inflammation marker; elevated in aneurysm",
        "category": "Vascular_Relevant",
        "subcategory": "Inflammation",
    },
    "MCP1": {
        "full_name": "Monocyte Chemoattractant Protein 1",
        "role": "Same as CCL2; monocyte recruitment to vessel wall",
        "category": "Vascular_Relevant",
        "subcategory": "Inflammation",
    },
}

# ============================================================
# MERGE all gene sets
# ============================================================
ALL_GENES = {}
for gene_set in [CORE_FERRO_AGING, IRON_METABOLISM, FERROPTOSIS,
                  LIPID_PEROXIDATION, SENESCENCE, NRF2_PATHWAY,
                  VASCULAR_RELEVANT]:
    for gene, info in gene_set.items():
        if gene in ALL_GENES:
            # Merge categories for multi-category genes
            existing = ALL_GENES[gene]
            if "category_list" not in existing:
                existing["category_list"] = [existing["category"]]
            existing["category_list"].append(info["category"])
            existing["category"] = "|".join(existing["category_list"])
            if info.get("subcategory"):
                existing["subcategory"] = f"{existing.get('subcategory','')}|{info['subcategory']}"
        else:
            ALL_GENES[gene] = dict(info)

# Remove category_list helper key before writing
for gene in ALL_GENES:
    ALL_GENES[gene].pop("category_list", None)

# ============================================================
# SAVE to CSV
# ============================================================
output_dir = Path(r"C:\Users\lidaf\WorkBuddy\2026-06-11-19-29-35")
output_file = output_dir / "ferro_aging_geneset.csv"

fieldnames = ["gene_symbol", "full_name", "role", "category", "subcategory", "evidence", "VC_target"]

with open(output_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
    writer.writeheader()
    for gene in sorted(ALL_GENES.keys()):
        row = {"gene_symbol": gene, **ALL_GENES[gene]}
        # Ensure all fields exist
        for field in fieldnames:
            if field not in row:
                row[field] = ""
        writer.writerow(row)

# ============================================================
# STATISTICS
# ============================================================
print("=" * 70)
print("  FERRO-AGING GENE SET - Build Complete")
print("=" * 70)

categories = {}
for gene, info in ALL_GENES.items():
    for cat in info["category"].split("|"):
        cat = cat.strip()
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(gene)

for cat, genes in sorted(categories.items()):
    # deduplicate
    unique_genes = sorted(set(genes))
    print(f"\n  [{cat}]  ({len(unique_genes)} genes)")
    print(f"    {', '.join(unique_genes)}")

total_unique = len(ALL_GENES)
print(f"\n{'=' * 70}")
print(f"  TOTAL: {total_unique} unique genes across {len(categories)} categories")
print(f"  Output: {output_file}")
print(f"{'=' * 70}")

# Save summary JSON
summary = {
    "total_genes": total_unique,
    "categories": {cat: sorted(set(genes)) for cat, genes in categories.items()},
    "gene_count_per_category": {cat: len(set(genes)) for cat, genes in categories.items()},
}

summary_file = output_dir / "ferro_aging_geneset_summary.json"
with open(summary_file, "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print(f"\n  Summary JSON: {summary_file}")
