import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load enrichment results
df = pd.read_csv("results/annotation/functional_enrichment_results.csv")

# Drop missing p-values, if any
df = df[df['p_value'].notna()]

# Compute -log10(p-value) for better visual scaling
df['minus_log10_p'] = -np.log10(df['p_value'])

# Select top 10 terms by p-value
top_terms = df.sort_values('p_value').head(10)

# Plot
plt.figure(figsize=(10, 6))
sns.barplot(data=top_terms, y="name", x="minus_log10_p", color="skyblue")
plt.xlabel("-log10(p-value)")
plt.ylabel("Enriched Term")
plt.title("Top 10 Enriched GO/KEGG Terms")
plt.tight_layout()
plt.savefig("results/annotation/enrichment_plot.png")

print("✅ Enrichment plot saved to results/annotation/enrichment_plot.png")
