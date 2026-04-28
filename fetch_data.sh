#!/bin/bash

# LocusLift Data Fetcher
# This script downloads necessary liftover chains and sample data.

set -e

echo "📦 Setting up LocusLift data..."

# 1. Create directories
mkdir -p chains
mkdir -p sample-1kg

# 2. Download Liftover Chains (UCSC)
echo "🔗 Downloading UCSC chain files..."
curl -L -o chains/hg19ToHg38.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
curl -L -o chains/hg38ToHg19.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

# 3. Download Sample Data
# NOTE: Replace these URLs with your actual hosted sample data links
echo "🔗 Downloading sample data..."
# curl -L -o sample-1kg/1kg-vcf.vcf.gz "YOUR_URL_HERE"

echo "✅ Done! You can now run LocusLift."
