<p align="center">
  <img width="160" alt="Image" src="https://github.com/user-attachments/assets/c04dc300-5124-47ca-a463-5e2d8ee65b40" />
</p>

<h1 align="center">LocusLift</h1>

<p align="center">
  <b>Genome LiftOver Pipeline (PLINK + VCF)</b><br>
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#usage-reference">Usage Reference</a> •
  <a href="#examples">Examples</a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-1.0-blue" alt="Version" />
  <img src="https://img.shields.io/badge/language-bash-green" alt="Bash" />

A production-ready shell workflow for genome-build coordinate liftover (GRCh37/hg19 ↔ GRCh38/hg38), with native support for PLINK and VCF datasets via UCSC `liftOver` chain mappings.
</p>

---
 Manual genomic liftover is a grind: generating intermediate UCSC `.bed` files, resolving duplicate variant IDs, remapping coordinates, and reconstructing PLINK or VCF outputs. It's tedious, error-prone, and pulls you away from actual research.

**LocusLift** automates the entire pipeline end-to-end. Point it at your dataset and it handles everything, from intermediate file formatting to liftOver execution, sorting, deduplication and cleanup. It is validated against real-world cohort data, including multi-ancestry datasets with duplicate variant IDs, missing positions, and mixed chromosome notations, so edge cases don't become your problem.

---

## Repository Structure

```text
LocusLift/
├── LocusLift.sh
├── README.md
├── chains/
│   ├── hg19ToHg38.over.chain.gz
│   └── hg38ToHg19.over.chain.gz
└── sample-1kg/
    ├── 1kg-plink1.bed
    ├── 1kg-plink1.bim
    ├── 1kg-plink1.fam
    ├── 1kg-plink2.pgen
    ├── 1kg-plink2.psam
    ├── 1kg-plink2.pvar
    ├── 1kg-vcf.vcf.gz
    ├── 1kg-vcf.vcf.gz.csi
    ├── chr16.bed
    ├── chr16.bim
    └── chr16.fam
```

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [Overview](#overview)
- [Features](#features)
- [Supported Formats](#supported-formats)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
- [Usage Reference](#usage-reference)
- [Examples](#examples)
- [Output Directory Structure](#output-directory-structure)
- [Log Files](#log-files)
- [Edge Cases & Caveats](#edge-cases--caveats)
- [Assembly Build Detection Deep Dive](#assembly-build-detection-deep-dive)
- [Force Modes](#force-modes)
- [Troubleshooting](#troubleshooting)
- [Reproducibility Recommendations](#reproducibility-recommendations)
- [Screenshots](#screenshots)
- [References](#references)

---

## Overview

`LocusLift.sh` automates the complex process of converting genomic coordinates between different genome assemblies (like GRCh37 to GRCh38). It manages format conversions, handles duplicate IDs, verifies build compatibility before execution, and securely maps coordinates, producing a clean, ready-to-use dataset with comprehensive logging.

---

## Features

- **Dual Mode Execution**: First-class support for both PLINK and VCF datasets.
- **Smart Build Detection**: Automatically infers the input genome build to prevent accidental liftover mappings.
- **Automated Deduplication**: Removes duplicate IDs in PLINK mode to prevent mapping conflicts.
- **Robust Error Handling**: Traps errors immediately, retains failing logs, and cleans up intermediate files on success.
- **Visual Feedback**: Provides a clean, progress-tracked terminal interface with colored glyphs.
- **Flexible Overrides**: Includes `--force` and `--force-all` flags for expert users bypassing automated safety checks.

---

## Supported Formats

| Mode | Input Formats | Output Formats | Tool Used |
| :--- | :--- | :--- | :--- |
| **PLINK** | `plink1` (`.bed`/`.bim`/`.fam`)<br>`plink2` (`.pgen`/`.pvar`/`.psam`) | `plink1`<br>`plink2` | `plink2` |
| **VCF** | Uncompressed (`.vcf`)<br>Compressed (`.vcf.gz`, `.bcf`) | Uncompressed (`.vcf`)<br>Compressed (`.vcf.gz`) | `bcftools` |

---

## Requirements

The pipeline requires a few standard bioinformatics tools depending on the mode you intend to use.

**Core Dependencies:**
- `bash` (Unix shell)
- `awk`
- UCSC `liftOver` executable
- UCSC Chain file (e.g., `hg38ToHg19.over.chain.gz`)

**Mode-Specific Dependencies:**
- **PLINK Mode**: `plink2`
- **VCF Mode**: `bcftools`

---

## Installation

### 1 — Get the Repository

Choose the method that fits your needs:

**Option A: Lightweight Install (Recommended)**

Downloads only the core script, documentation, and `chains/` directory for the fastest way to get started if you don't need sample datasets.

```bash
mkdir LocusLift && cd LocusLift
git init
git remote add origin https://github.com/bioinumer/LocusLift.git
git sparse-checkout init --cone
git sparse-checkout set LocusLift.sh README.md chains/
git pull origin main
```

**Option B: Full Install**

Clones the entire repository, including `sample-1kg/`. Note that the genomic data in that directory makes this significantly larger and slower to download.

```bash
git clone https://github.com/bioinumer/LocusLift.git
cd LocusLift
```

---

### 2 — Download the UCSC `liftOver` Binary

Select the binary for your OS and architecture:

| OS / Architecture | Download |
| :--- | :--- |
| Linux x86\_64 | [Download](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver) |
| Linux ARM64 | [Download](https://hgdownload.soe.ucsc.edu/admin/exe/linux.aarch64.v492/liftOver) |
| macOS Apple Silicon | [Download](https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/liftOver) |
| macOS Intel | [Download](https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/liftOver) |

Not sure which architecture you have? Run:

```bash
uname -s && uname -m
```

---

### 3 — Install Dependencies

| Tool | Install |
| :--- | :--- |
| **plink2** | Download from [cog-genomics.org](https://www.cog-genomics.org/plink/2.0/) |
| **bcftools** | `brew install bcftools` or build from [htslib.org](https://www.htslib.org/download/) |

---

### 4 — Add Binaries to PATH

Move tools to a directory on your system `PATH` so the script can find them automatically:

```bash
chmod +x liftOver
sudo mv liftOver /usr/local/bin/
```

> **Tip:** If you'd rather not move the binary, you can keep `liftOver` in the same directory as `LocusLift.sh` and the script will find it automatically.

---

### 5 — Verify the Installation

Run the script with the `-h` flag to confirm everything is working:

```bash
./LocusLift.sh -h
```

<img width="803" height="738" alt="LocusLift help menu" src="https://github.com/user-attachments/assets/575a8f7e-a284-4613-8c86-98e91aa2a7fc" />
---

## Quick Start

**PLINK Mode** (GRCh38 → GRCh37):
```bash
./LocusLift.sh \
  -i my_cohort_hg38 \
  -c hg38ToHg19.over.chain.gz \
  -o my_cohort_hg19
```

**VCF Mode** (GRCh38 → hg19, Compressed):
```bash
./LocusLift.sh \
  -I my_cohort_hg38.vcf.gz \
  -c hg38ToHg19.over.chain.gz \
  -O my_cohort_hg19.vcf.gz
```

---

## How It Works

The script follows a rigorous multi-step process to ensure data integrity.

### Step 1: Assembly Build Detection
The script parses the input file (VCF headers or variant coordinates) to infer if the data is GRCh37 or GRCh38. It compares this against the provided chain file to prevent executing a liftover with mismatched builds.
> **Note**: For extra credibility on the build check, you can compare the variants against a reference FASTA file in downstream QC (e.g., using GATK).

### Step 2: Coordinate Extraction
Variants are extracted and converted into the UCSC BED format (0-based start, 1-based end). In PLINK mode, it automatically drops variants mapped to non-primary contigs (unless `--keep-all-chr` is used) and removes duplicate IDs.

### Step 3: UCSC liftOver
The core UCSC `liftOver` binary is executed, mapping coordinates from the source to the target genome build. Unmapped variants are saved separately.

### Step 4: Coordinate Mapping Rewrite
A mapping file is generated linking original variants to their new coordinates.
- **PLINK**: Generates `--update-map` and `--update-chr` tables.
- **VCF**: Prepares an ID rewrite table.

### Step 5: Application & Sorting
- **PLINK**: `plink2` applies the updates, renames variant IDs, and sorts the data into a temporary `.pgen` before final export.
- **VCF**: `awk` rewrites the VCF body in-place. `bcftools` then sorts and indexes the output.

---

## Usage Reference

```text
Usage
  PLINK  ./LocusLift.sh -i <prefix> -c <chain> -o <out_prefix> [opts]
  VCF    ./LocusLift.sh -I <in.vcf[.gz]> -c <chain> -O <out.vcf[.gz]> [opts]

Required
  -c       UCSC chain file  (e.g. hg38ToHg19.over.chain.gz)
  -i / -o  PLINK input / output prefix
  -I / -O  VCF  input / output file

Options
  -f FORMAT        Input  format: plink1|plink2  [auto]
  -t FORMAT        Output format: plink1|plink2  [same as -f]
  -l PATH          liftOver binary               [auto-search]
  -v TMPL          Variant-ID template           [default: chr:pos:ref:alt]
  -k               Keep intermediate files
  --keep-all-chr   Retain non-primary contigs    [default: 1-22/X/Y]
  --force          Skip build-compatibility check for inconclusive inferences
  --force-all      Bypass all safety checks and force liftover despite explicit mismatches
  -h               Show help
```

---

## Examples

<details>
<summary><b>1. PLINK Liftover (hg19 → hg38)</b></summary>

Lifts over a standard PLINK dataset from GRCh37 to GRCh38.
```bash
./LocusLift.sh \
  -i data_hg19 \
  -c hg19ToHg38.over.chain.gz \
  -o data_hg38
```
</details>

<details>
<summary><b>2. Format Conversion During Liftover</b></summary>

Takes `plink1` input and outputs `plink2` format, keeping all alternative chromosomes.
```bash
./LocusLift.sh \
  -i data_hg19 -c hg19ToHg38.over.chain.gz -o lifted_data \
  -f plink1 -t plink2 --keep-all-chr
```
</details>

<details>
<summary><b>3. VCF Compressed Input/Output</b></summary>

Lifts a compressed VCF and automatically generates a `.csi` index.
```bash
./LocusLift.sh \
  -I cohort.hg38.vcf.gz \
  -c hg38ToHg19.over.chain.gz \
  -O cohort.hg19.vcf.gz
```
</details>

<details>
<summary><b>4. Custom liftOver Binary Path</b></summary>

Specify exactly where the binary is located if it's not in your PATH.
```bash
./LocusLift.sh \
  -i cohort -c chain.gz -o lifted \
  -l /opt/software/ucsc/liftOver
```
</details>

<details>
<summary><b>5. Forcing Execution (--force-all)</b></summary>

Bypass a known mismatch (e.g., lifting GRCh38 data when the script thinks it's GRCh37). *Use with extreme caution.*
```bash
./LocusLift.sh \
  -i cohort -c hg38ToHg19.over.chain.gz -o lifted \
  --force-all
```
</details>

---

## Output Directory Structure

The script neatly organizes outputs into directories:

```text
output_prefix/
├── logs/
│   ├── LocusLift.log               # Summary of the successful run
│   └── LocusLift.rmdup.log         # (PLINK only) Duplicate removal log
├── lifted-data/
│   ├── my_cohort.bed/.pgen/.vcf    # Final output data
│   ├── my_cohort.bim/.pvar
│   ├── my_cohort.fam/.psam
│   └── my_cohort.unmapped.bed      # Variants that failed to lift
```

---

## Log Files

- **`LocusLift.log`**: Generated on success. Contains timestamps, variant counts, drop rates, and CLI parameters.
- **`LocusLift.error.log`**: Generated on failure. Contains detailed tracebacks, error messages, and detected builds.
- **`LocusLift.failed-step.log`**: Retained on failure. Contains the raw `stderr`/`stdout` of the specific tool (plink/bcftools) that crashed.
- **`LocusLift.rmdup.log`**: Generated if duplicate IDs were found and removed in PLINK mode.

---

## Edge Cases & Caveats

1. **VCF Caveat (Coordinate Only)**: VCF mode performs coordinate-only liftover. It does *not* reverse-complement REF/ALT alleles if the strand swaps, nor does it check against a target reference FASTA. For rigorous allele-aware VCF liftover, use Picard `LiftoverVcf`.
2. **chr16-Only Datasets**: The size difference between GRCh37 and GRCh38 for chr16 is too small to confidently infer the build automatically. The script will warn you and prompt for continuation.
3. **Non-Primary Contigs**: By default, variants on `GL*`, `KI*`, or other decoy contigs are dropped. Use `--keep-all-chr` to retain them.
4. **Duplicate Variant IDs**: In PLINK mode, `plink2 --rm-dup exclude-mismatch` is automatically invoked if duplicate IDs are detected, as UCSC `liftOver` requires unique IDs to map back correctly.
5. **Variant ID Update Behaviour**:
   - **VCF mode** checks whether each existing ID matches the old positional format (e.g., `chr:pos:ref:alt` or `chr_pos_ref_alt`) before overwriting it. IDs that match are updated to the new lifted coordinates; rsIDs (e.g., `rs123456`) and missing IDs (`.`) are preserved.
   - **PLINK mode** uses `--set-all-var-ids`, which unconditionally overwrites every variant ID with the new coordinate format (default: `chr:pos:ref:alt`). ⚠️ This permanently erases any existing rsIDs or custom annotations.

---

## Assembly Build Detection Deep Dive

To prevent disastrous wrong-direction liftovers, the script implements a two-tier inference engine:
1. **Header Analysis (VCF only)**: Scans for explicit build tags (`GRCh37`, `hg38`) and specific sequence markers (like `NC_000001.10` vs `NC_000001.11`).
2. **Coordinate Analysis**: Scans the terminal coordinates (maximum position) for each primary chromosome and compares them against the known exact lengths of GRCh37 and GRCh38. For example, if a variant on chr1 is at position `249,200,000`, it *must* be GRCh37, because GRCh38 chr1 is only `248,956,422` bases long.

---

## Force Modes

Sometimes the automatic detection fails or is overly cautious.

- **`--force`**: Use this when the script returns `inconclusive` or `unknown` for the build detection. It bypasses the "unknown" block but will still fail if an explicit *mismatch* is detected.
- **`--force-all`**: The nuclear option. Bypasses all safety checks. It will execute the liftover even if the script is 100% certain that the input is GRCh38 and the chain file expects GRCh37. Use entirely at your own risk.

---

## Troubleshooting

- **`liftOver binary appears incompatible`**: You downloaded the wrong binary for your OS. For Apple Silicon Macs, ensure you download `macOSX.arm64`.
- **`liftOver mapped 0 variants`**: You likely specified the wrong chain file direction (e.g., trying to map hg19 data using an hg38ToHg19 chain).
- **`awk: syntax error`**: Ensure you are using a standard modern `awk` (GNU awk or macOS default).

---

## Reproducibility Recommendations

For academic and clinical settings:
1. Record exact versions of `liftOver`, `plink2`, and `bcftools`.
2. Archive the `.run.log` alongside your final datasets.
3. For VCFs, perform a downstream sanity check by calculating allele frequencies and validating REF/ALT consistency against the target reference FASTA.

---

## Screenshots

<!-- screenshot: Build Detection Failure -->
### Build Detection Failure
<img width="827" height="406" alt="Image" src="https://github.com/user-attachments/assets/b082cb25-1310-4f5f-8187-1c55420a950c" />

<!-- screenshot: Successful Vcf Run -->
### Successful Vcf Run
<img width="826" height="468" alt="Image" src="https://github.com/user-attachments/assets/abfc9de9-42f6-485f-9c5f-25880b4dae56" />

<!-- screenshot: LocusLift log file -->
### LocusLift Log File
<img width="882" height="636" alt="Image" src="https://github.com/user-attachments/assets/2c47696c-ebdf-4c5f-bc8a-447d722315ea" />

---

## References

[1] Monti, R. (2023). HapMap3-1KG: 1000 Genomes processed genotypes in PLINK (.bed/bim/fam) format. figshare. https://doi.org/10.6084/m9.figshare.20802700.v1

[2] UCSC BED format FAQ (0-based `chromStart`, half-open interval semantics): [UCSC FAQ](https://genome.ucsc.edu/FAQ/FAQformat.html)

[3] VCF specification (POS is 1-based, fixed fields): [HTS Specs](https://samtools.github.io/hts-specs/VCFv4.3.pdf)

[4] UCSC command-line liftOver usage and caveats: [UCSC FAQ Downloads](https://genome.ucsc.edu/FAQ/FAQdownloads)

[5] bcftools manual (VCF/BCF support, query/view/sort/index): [bcftools](https://samtools.github.io/bcftools/bcftools.html)

[6] Picard/GATK LiftoverVcf documentation (allele-aware VCF liftover behavior): [Broad Institute](https://gatk.broadinstitute.org/hc/en-us/articles/360036363632-LiftoverVcf-Picard)

---
*Last updated: 2026-04-26*
