# functional_enrichment.py

import pandas as pd
from gprofiler import GProfiler

# Load your gene annotation CSV
df = pd.read_csv("results/annotation/variant_annotation_summary.csv")

# Filter high-impact genes
high_impact = df[df['Impact'].isin(['HIGH', 'MODERATE'])]

# Drop missing or malformed gene names
genes = high_impact['Gene'].dropna().unique().tolist()

# Initialize gProfiler
gp = GProfiler(return_dataframe=True)

# Run enrichment
results = gp.profile(organism='hsapiens', query=genes)

# Save results
results.to_csv("results/annotation/functional_enrichment_results.csv", index=False)

print("✅ Functional enrichment completed and saved.")
