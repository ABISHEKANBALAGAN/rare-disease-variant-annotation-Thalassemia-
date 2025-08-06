import subprocess
import sys

# Install vcfpy if not available
try:
    import vcfpy
except ImportError:
    print("🔧 vcfpy not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "vcfpy"])
    import vcfpy

import pandas as pd

# Input VEP-annotated VCF
reader = vcfpy.Reader.from_path("results/annotation/SRR32809192.vep.vcf")

records = []
for record in reader:
    csq_info = record.INFO.get("CSQ", [])
    for ann in csq_info:
        csq_fields = ann.split('|')
        records.append({
            'CHROM': record.CHROM,
            'POS': record.POS,
            'REF': record.REF,
            'ALT': record.ALT[0].value,
            'Gene': csq_fields[3],
            'Consequence': csq_fields[1],
            'Impact': csq_fields[2],
            'Codons': csq_fields[8],
            'Amino_acids': csq_fields[10],
            'Protein_position': csq_fields[14]
        })

df = pd.DataFrame(records)
df.to_csv('results/annotation/variant_annotation_summary.csv', index=False)
print("✅ Variant annotation summary saved to results/annotation/variant_annotation_summary.csv")
