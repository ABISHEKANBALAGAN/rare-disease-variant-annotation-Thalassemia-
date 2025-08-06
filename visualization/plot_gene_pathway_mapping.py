import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load mapping data
mapping_df = pd.read_csv("results/annotation/gene_to_pathway_mapping.csv")

# Create graph
G = nx.Graph()

# Add nodes and edges
for _, row in mapping_df.iterrows():
    gene = row["Gene"]
    term = row["Enriched_Term"]
    G.add_node(gene, type="gene")
    G.add_node(term, type="term")
    G.add_edge(gene, term, weight=1)

# Define node colors
node_colors = []
for node in G.nodes(data=True):
    if node[1]["type"] == "gene":
        node_colors.append("skyblue")
    else:
        node_colors.append("lightgreen")

# Draw the graph
plt.figure(figsize=(12, 8))
pos = nx.spring_layout(G, k=0.5, seed=42)  # Layout with fixed randomness

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=1000, alpha=0.9)
nx.draw_networkx_edges(G, pos, edge_color="gray", width=1.5)
nx.draw_networkx_labels(G, pos, font_size=10, font_family="sans-serif")

plt.title("Gene–Pathway Enrichment Network", fontsize=14)
plt.axis("off")
plt.tight_layout()

# Save to file
plt.savefig("results/annotation/gene_pathway_network.png", dpi=300)
plt.show()
