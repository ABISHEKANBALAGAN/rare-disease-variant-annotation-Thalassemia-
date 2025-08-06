import subprocess
import os

# Paths
bwa_path = "bwa"
reference = "data/reference/chr11.fa"
fastq1 = "data/sub_fastq/SRR32809192_1_sub.fastq"
fastq2 = "data/sub_fastq/SRR32809192_2_sub.fastq"
output_sam = "results/alignment/SRR32809192.sam"

os.makedirs("results/alignment", exist_ok=True)

# Run BWA
cmd = [bwa_path, "mem", reference, fastq1, fastq2]
with open(output_sam, "w") as out:
    result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)

# Check success
if result.returncode == 0:
    print("✅ BWA alignment done.")
else:
    print("❌ Error:")
    print(result.stderr.decode())
