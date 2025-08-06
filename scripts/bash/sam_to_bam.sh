# Convert SAM to BAM
samtools view -Sb results/alignment/SRR32809192.sam > results/alignment/SRR32809192.bam

# Sort the BAM file
samtools sort results/alignment/SRR32809192.bam -o results/alignment/SRR32809192.sorted.bam

# Index the sorted BAM
samtools index results/alignment/SRR32809192.sorted.bam
