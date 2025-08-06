import os
import subprocess

# Define input/output
vep_vcf = "results/annotation/SRR32809192.vep.vcf"  # adjust path if needed
output_dir = os.path.dirname(vep_vcf)
bgzipped_vcf = vep_vcf + ".gz"
index_file = bgzipped_vcf + ".tbi"

# Compress with bgzip
print(f"🔧 Compressing VCF: {vep_vcf} ➜ {bgzipped_vcf}")
subprocess.run(f"bgzip -c {vep_vcf} > {bgzipped_vcf}", shell=True, check=True)

# Index with tabix
print(f"📌 Indexing VCF for IGV: {bgzipped_vcf}")
subprocess.run(f"tabix -p vcf {bgzipped_vcf}", shell=True, check=True)

# Confirm output
print("\n✅ Output files:")
print(f"- Compressed VCF: {bgzipped_vcf}")
print(f"- Index file    : {index_file}")
