install.packages("ggplot2")
install.packages("dplyr")
install.packages("readr")

# Load required packages
library(ggplot2)
library(dplyr)
library(readr)

#step1: Load variant annotation summary
df <- read_csv("results/annotation/variant_annotation_summary.csv")

# View structure
str(df)
head(df)

#Step2: Plot Variant Impact Distribution

#Count of variants by impact
df %>%
  count(Impact) %>%
  ggplot(aes(x = Impact, y = n, fill = Impact)) +
  geom_bar(stat = "identity") +
  labs(title = "Variant Impact Distribution", y = "Count") +
  theme_minimal()

#save the script
ggsave("results/annotation/variant_impact_distribution.png", width = 6, height = 4)

#step 3: Plot Consequence types (Top15)
df %>%
  count(Consequence, sort = TRUE) %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(Consequence, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 15 Variant Consequences", x = "Consequence", y = "Count") +
  theme_minimal()
#save the script
ggsave("results/annotation/top15_consequences.png", width = 6, height = 4)

#step 4:Gene Summary - Top affected genes
df %>%
  filter(!is.na(Gene) & Gene != "") %>%
  count(Gene, sort = TRUE) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = reorder(Gene, n), y = n)) +
  geom_col(fill = "tomato") +
  coord_flip() +
  labs(title = "Top 10 Genes with Variants", x = "Gene", y = "Variant Count") +
  theme_minimal()
#save the script
ggsave("results/annotation/top10_genes.png", width = 6, height = 4)

#Step 5: Filter variants(Variants in key thalassemia-related genes)
thal_genes <- c("HBB", "HBA1", "HBA2")
thal_variants <- df %>%
  filter(Gene %in% thal_genes)
write_csv(thal_variants, "results/annotation/thalassemia_gene_variants.csv")

#Advance Step: Filter variants by consequence types
# Define important consequences
key_consequences <- c(
  "missense_variant", 
  "stop_gained", 
  "frameshift_variant", 
  "splice_donor_variant", 
  "splice_acceptor_variant", 
  "start_lost"
)

# Filter those variants
filtered_variants <- df %>%
  filter(Consequence %in% key_consequences)

# Save to CSV
write_csv(filtered_variants, "results/annotation/key_functional_variants.csv")

# Print summary
print(paste("✅ Extracted", nrow(filtered_variants), "key functional variants."))
