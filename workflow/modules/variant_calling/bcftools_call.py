import subprocess

bam = "results/alignment/SRR32809192.dedup.bam"
ref = "data/reference/chr11.fa"
vcf = "results/variants/SRR32809192.vcf"

# Ensure output folder exists
subprocess.run("mkdir -p results/variants", shell=True)

# Variant calling with bcftools
cmd = f"""
bcftools mpileup -Ou -f {ref} {bam} | \
bcftools call -mv -Ov -o {vcf}
"""

print("🔬 Running bcftools variant calling...")
subprocess.run(cmd, shell=True, check=True)
print("✅ Variant calling done:", vcf)
