import subprocess
import os

# Input and output file paths
input_bam = "results/alignment/SRR32809192.sorted.bam"
name_sorted_bam = "results/alignment/SRR32809192.name_sorted.bam"
fixmate_bam = "results/alignment/SRR32809192.fixmate.bam"
coord_sorted_bam = "results/alignment/SRR32809192.coord_sorted.bam"
dedup_bam = "results/alignment/SRR32809192.dedup.bam"

# Step 1: Sort by read name
print("📦 Sorting BAM by read name...")
subprocess.run(f"samtools sort -n -o {name_sorted_bam} {input_bam}", shell=True, check=True)

# Step 2: Fix mate information
print("🔧 Fixing mate information...")
subprocess.run(f"samtools fixmate -m {name_sorted_bam} {fixmate_bam}", shell=True, check=True)

# Step 3: Coordinate sort for markdup
print("📍 Sorting BAM by coordinate...")
subprocess.run(f"samtools sort -o {coord_sorted_bam} {fixmate_bam}", shell=True, check=True)

# Step 4: Mark duplicates
print("🧹 Marking duplicates...")
subprocess.run(f"samtools markdup {coord_sorted_bam} {dedup_bam}", shell=True, check=True)

# Step 5: Index the deduplicated BAM
print("🧭 Indexing deduplicated BAM...")
subprocess.run(f"samtools index {dedup_bam}", shell=True, check=True)

print("✅ Mark duplicates step completed. Output:")
print(f"    ➤ {dedup_bam}")
print(f"    ➤ {dedup_bam}.bai")
