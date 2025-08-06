import pandas as pd

# Load data
variants = pd.read_csv("results/annotation/thalassemia_gene_variants.csv")
enrichment = pd.read_csv("results/annotation/functional_enrichment_results.csv")

# Step 1: Extract unique gene names
variant_genes = variants['Gene'].dropna().unique()
variant_genes = [g for g in variant_genes if g != "NA"]

# Step 2: Add Gene info to enrichment results (assuming they all came from query_1)
enrichment['Gene'] = ','.join(variant_genes)

# Step 3: Create mapping table
gene_to_pathway = []
for gene in variant_genes:
    for idx, row in enrichment.iterrows():
        gene_to_pathway.append({
            "Gene": gene,
            "Enriched_Term": row['name'],
            "Description": row['description'],
            "p_value": row['p_value']
        })

mapping_df = pd.DataFrame(gene_to_pathway)

# Save result
mapping_df.to_csv("results/annotation/gene_to_pathway_mapping.csv", index=False)
print("✅ Mapped genes to enriched terms ➜ results/annotation/gene_to_pathway_mapping.csv")
