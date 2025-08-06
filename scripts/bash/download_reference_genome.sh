#!/bin/bash

# Set download directory
REF_DIR="data/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

# Download reference genome (GRCh38) from Ensembl via HTTPS
echo "Downloading GRCh38 (hg38) reference genome..."
wget -O GRCh38.fa.gz \
  https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the genome
echo "Unzipping..."
gunzip GRCh38.fa.gz

echo "Download complete. File saved at: $REF_DIR/GRCh38.fa"

