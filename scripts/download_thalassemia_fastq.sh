#!/usr/bin/env bash

OUT_DIR="data/raw_reads"
mkdir -p "$OUT_DIR"

echo "🔽 Downloading FASTQ files for SRR32809192..."

# Paired-end FASTQ URLs (from your TSV file)
URL1="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR328/092/SRR32809192/SRR32809192_1.fastq.gz"
URL2="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR328/092/SRR32809192/SRR32809192_2.fastq.gz"

# Download using curl
curl -fLo "$OUT_DIR/SRR32809192_1.fastq.gz" "$URL1"
curl -fLo "$OUT_DIR/SRR32809192_2.fastq.gz" "$URL2"

echo "✅ Download complete. Files saved to $OUT_DIR"
