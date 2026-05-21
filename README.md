<p align="center">
  <img width="160" alt="LocusLift logo" src="https://github.com/user-attachments/assets/c04dc300-5124-47ca-a463-5e2d8ee65b40" />
</p>

<h1 align="center">LocusLift</h1>

<p align="center">
  <b>End-to-end genome coordinate liftover for PLINK and VCF datasets</b>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-1.0-blue" alt="Version 1.0" />
  <img src="https://img.shields.io/badge/language-bash-green" alt="Bash" />
</p>

<p align="center">
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#usage-reference">Usage Reference</a> •
  <a href="#examples">Examples</a> •
  <a href="#troubleshooting">Troubleshooting</a>
</p>

---

Manual genomic liftover is tedious: generate UCSC `.bed` intermediates, resolve duplicate variant IDs, remap coordinates, rebuild PLINK or VCF outputs, all before your actual analysis begins. **LocusLift** automates the entire pipeline. Point it at your dataset and it handles everything, including build detection, deduplication, liftover, sorting, and cleanup.

Validated against real-world cohort data with duplicate IDs, missing positions, and mixed chromosome notations.

---

## Table of Contents

- [Features](#features)
- [Supported Formats](#supported-formats)
- [Requirements](#requirements)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
- [Usage Reference](#usage-reference)
- [Examples](#examples)
- [Output Structure](#output-structure)
- [Log Files](#log-files)
- [Edge Cases & Caveats](#edge-cases--caveats)
- [Assembly Build Detection](#assembly-build-detection)
- [Force Modes](#force-modes)
- [Troubleshooting](#troubleshooting)
- [Reproducibility Recommendations](#reproducibility-recommendations)
- [Screenshots](#screenshots)
- [References](#references)

---

## Features

| | |
|---|---|
| **Dual-mode execution** | First-class support for both PLINK and VCF datasets |
| **Smart build detection** | Infers input assembly automatically to prevent wrong-direction liftovers |
| **Automatic deduplication** | Removes duplicate IDs in PLINK mode before mapping |
| **Robust error handling** | Traps errors immediately, retains diagnostic logs, cleans up on success |
| **Visual progress feedback** | Color-coded terminal output with step-by-step status |
| **Flexible overrides** | `--force` and `--force-all` flags for expert users |

---

## Supported Formats

| Mode | Input | Output | Tool |
|:---|:---|:---|:---|
| **PLINK** | `.bed/.bim/.fam` (plink1), `.pgen/.pvar/.psam` (plink2) | plink1 or plink2 | `plink2` |
| **VCF** | `.vcf`, `.vcf.gz`, `.bcf` | `.vcf`, `.vcf.gz` | `bcftools` |

---

## Requirements

**Core (all modes):**
- `bash`
- `awk`
- UCSC `liftOver` binary
- UCSC chain file (e.g., `hg38ToHg19.over.chain.gz`)

**PLINK mode:** `plink2`

**VCF mode:** `bcftools`

---

## Repository Structure

```
LocusLift/
├── LocusLift.sh                    # Main pipeline script
├── README.md
└── chains/                         # Pre-downloaded UCSC chain files
    ├── hg19ToHg38.over.chain.gz
    └── hg38ToHg19.over.chain.gz
```

> A subset of the 1000 Genomes dataset  is available in `sample-1kg.zip`, included in [GitHub Releases](https://github.com/bioinumer/LocusLift/releases) for testing purposes.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/bioinumer/LocusLift.git
cd LocusLift
```

### 2. Download the UCSC `liftOver` binary

Choose the binary for your OS and architecture:

| Platform | Download |
|:---|:---|
| Linux x86_64 | [Download](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver) |
| Linux ARM64 | [Download](https://hgdownload.soe.ucsc.edu/admin/exe/linux.aarch64.v492/liftOver) |
| macOS Apple Silicon | [Download](https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/liftOver) |
| macOS Intel | [Download](https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/liftOver) |

Not sure which architecture you have? Run `uname -s && uname -m`.

Once downloaded, make it executable and move it to your `PATH`:

```bash
chmod +x liftOver
sudo mv liftOver /usr/local/bin/
```

Alternatively, keep it in the LocusLift directory — the script will find it automatically.

### 3. Install mode-specific dependencies

- **plink2:** [cog-genomics.org/plink/2.0](https://www.cog-genomics.org/plink/2.0/)
- **bcftools:** `brew install bcftools` or from [htslib.org](https://www.htslib.org/download/)

### 4. Verify installation

```bash
./LocusLift.sh -h
```

<img width="803" height="738" alt="LocusLift help menu" src="https://github.com/user-attachments/assets/575a8f7e-a284-4613-8c86-98e91aa2a7fc" />

### 5. Download the sample dataset (optional)

A subset of 1000 Genomes Project data [[1](#references)] is available in PLINK 1.9, PLINK 2.0, and VCF formats.

```bash
curl -L -O https://github.com/bioinumer/LocusLift/releases/download/v1.0/sample-1kg.zip
unzip sample-1kg.zip
```

---

## Quick Start

**PLINK mode** (GRCh38 → GRCh37):
```bash
./LocusLift.sh \
  -i my_cohort_hg38 \
  -c hg38ToHg19.over.chain.gz \
  -o my_cohort_hg19
```

**VCF mode** (GRCh38 → GRCh37, compressed):
```bash
./LocusLift.sh \
  -I my_cohort_hg38.vcf.gz \
  -c hg38ToHg19.over.chain.gz \
  -O my_cohort_hg19.vcf.gz
```

---

## How It Works

### Step 1: Assembly build detection

The script parses the input file (VCF headers or variant coordinates) to infer whether the data is GRCh37 or GRCh38, then checks this against the provided chain file. If a mismatch is detected, execution halts before any data is modified.

### Step 2: Coordinate extraction

Variants are converted to UCSC BED format (0-based start, 1-based end). Variants on non-primary contigs (i.e., NOT chr 1-22,X,Y) are dropped by default (override with `--keep-all-chr`) and duplicate IDs are removed automatically.

### Step 3: UCSC liftOver

The `liftOver` binary maps coordinates from the source to the target assembly. Unmapped variants are saved to a separate `.unmapped.bed` file.

### Step 4: Coordinate rewrite

A mapping file links original variants to their new coordinates:
- **PLINK:** generates `--update-map` and `--update-chr` tables
- **VCF:** generates an ID rewrite table

### Step 5: Application and sorting

- **PLINK:** `plink2` applies coordinate updates, renames variant IDs, and exports the final dataset
- **VCF:** `awk` rewrites coordinates in-place; `bcftools` sorts and indexes the output

---

## Usage Reference

```
Usage
  PLINK  ./LocusLift.sh -i <prefix> -c <chain> -o <out_prefix> [opts]
  VCF    ./LocusLift.sh -I <in.vcf[.gz]> -c <chain> -O <out.vcf[.gz]> [opts]

Required
  -c         UCSC chain file (e.g., hg38ToHg19.over.chain.gz)
  -i / -o    PLINK input / output prefix
  -I / -O    VCF input / output file

Options
  -f FORMAT        Input format: plink1|plink2  [auto-detected]
  -t FORMAT        Output format: plink1|plink2  [same as -f]
  -l PATH          liftOver binary path  [auto-searched]
  -v TMPL          Variant ID template  [default: chr:pos:ref:alt]
  -k               Keep intermediate files
  --keep-all-chr   Retain non-primary contigs  [default: chr1–22, X, Y]
  --force          Skip build check when inference is inconclusive
  --force-all      Bypass all safety checks, including explicit mismatches
  -h               Show this help
```

---

## Examples

<details>
<summary><b>1. PLINK liftover (hg19 → hg38)</b></summary>

```bash
./LocusLift.sh \
  -i data_hg19 \
  -c hg19ToHg38.over.chain.gz \
  -o data_hg38
```
</details>

<details>
<summary><b>2. Format conversion during liftover</b></summary>

Convert from plink1 input to plink2 output while retaining all alternative chromosomes:

```bash
./LocusLift.sh \
  -i data_hg19 -c hg19ToHg38.over.chain.gz -o lifted_data \
  -f plink1 -t plink2 --keep-all-chr
```
</details>

<details>
<summary><b>3. Compressed VCF with automatic indexing</b></summary>

Lifts a compressed VCF and generates a `.csi` index automatically:

```bash
./LocusLift.sh \
  -I cohort.hg38.vcf.gz \
  -c hg38ToHg19.over.chain.gz \
  -O cohort.hg19.vcf.gz
```
</details>

<details>
<summary><b>4. Custom liftOver binary path</b></summary>

```bash
./LocusLift.sh \
  -i cohort -c chain.gz -o lifted \
  -l /opt/software/ucsc/liftOver
```
</details>

<details>
<summary><b>5. Custom PLINK variant ID template</b></summary>

Changes the PLINK output variant IDs with `-v`. Quote the template so the shell does not expand `$r` and `$a`.

```bash
./LocusLift.sh \
  -i data_hg19 -c hg19ToHg38.over.chain.gz -o data_hg38 \
  -v 'chr@:#:$r:$a'
```
</details>

<details>
<summary><b>6. Bypassing safety checks (use with caution)</b></summary>

Forces execution despite a detected build mismatch. Only use when you are certain the detection is wrong.

```bash
./LocusLift.sh \
  -i cohort -c hg38ToHg19.over.chain.gz -o lifted \
  --force-all
```
</details>

---

## Output Structure

```
output_prefix/
├── logs/
│   ├── LocusLift.log           # Run summary (timestamps, variant counts, drop rates)
│   └── LocusLift.rmdup.log     # (PLINK only) Duplicate removal details
└── lifted-data/
    ├── my_cohort.bed           # Final output data (format depends on -t)
    ├── my_cohort.bim
    ├── my_cohort.fam
    └── my_cohort.unmapped.bed  # Variants that failed to lift
```

---

## Log Files

| File | Generated when | Contents |
|:---|:---|:---|
| `LocusLift.log` | Run succeeds | Timestamps, variant counts, drop rates, CLI parameters |
| `LocusLift.error.log` | Run fails | Tracebacks, error messages, detected builds |
| `LocusLift.failed-step.log` | Run fails | Raw stderr/stdout of the failing tool |
| `LocusLift.rmdup.log` | Duplicates found (PLINK) | List of removed duplicate IDs |

---

## Edge Cases & Caveats

**VCF liftover is coordinate-only.** REF/ALT alleles are not reverse-complemented if the strand flips, and the output is not validated against a target reference FASTA. For allele-aware VCF liftover, use Picard `LiftoverVcf` instead [[6](#references)].

**Chromosome 16 datasets.** The size difference between GRCh37 and GRCh38 for chr16 is too small to resolve build inference confidently. The script will warn you and prompt before continuing.

**Non-primary contigs.** Variants on `GL*`, `KI*`, or other decoy contigs are dropped by default. Use `--keep-all-chr` to retain them.

**Duplicate variant IDs.** In PLINK mode, `plink2 --rm-dup exclude-mismatch` is invoked automatically if duplicates are detected. UCSC `liftOver` requires unique IDs to map variants back correctly.

**Variant ID update behavior.** VCF mode only overwrites IDs that match a recognized positional format (e.g., `chr:pos:ref:alt` or `chr_pos_ref_alt`). rsIDs and `.` entries are preserved. PLINK mode uses `--set-all-var-ids`, which **unconditionally overwrites every variant ID**—including rsIDs and custom annotations—with the new coordinate format.

---

## Assembly Build Detection

To prevent wrong-direction liftovers, the script uses a two-tier detection engine:

**Tier 1 — Header analysis (VCF only).** Scans for explicit build tags (`GRCh37`, `hg38`) and sequence length markers (e.g., `NC_000001.10` vs. `NC_000001.11`).

**Tier 2 — Coordinate analysis.** Scans the maximum observed position per chromosome and compares it against the known lengths of GRCh37 and GRCh38. For example, a variant at position 249,200,000 on chr1 cannot be GRCh38 (which has a chr1 length of 248,956,422 bp), so the data is inferred as GRCh37.

> For additional confidence downstream, compare lifted variants against a reference FASTA (e.g., with GATK).

---

## Force Modes

| Flag | When to use |
|:---|:---|
| `--force` | Build detection returned `inconclusive` or `unknown`. Bypasses the block but still fails on an explicit mismatch. |
| `--force-all` | Overrides all safety checks, including confirmed mismatches. Use only when you are certain the automated detection is wrong. |

---

## Troubleshooting

**`liftOver binary appears incompatible`**
You downloaded the wrong binary for your platform. For Apple Silicon Macs, use the `macOSX.arm64` binary, not `macOSX.x86_64`.

**`liftOver mapped 0 variants`**
The chain file direction likely does not match your input build. Verify you are using `hg19ToHg38` for GRCh37 input, or `hg38ToHg19` for GRCh38 input.

**`awk: syntax error`**
Ensure you are using a standard, modern `awk` (GNU awk or macOS default). Some minimal containers ship without it.

---

## Reproducibility Recommendations

1. Record exact versions of `liftOver`, `plink2`, and `bcftools` alongside your results.
2. Archive the `LocusLift.log` file with your final datasets.
3. For VCFs, validate REF/ALT consistency against the target reference FASTA by computing allele frequencies and checking for strand flips.

---

## Screenshots

### Build detection failure

<img width="827" height="406" alt="Build detection failure" src="https://github.com/user-attachments/assets/b082cb25-1310-4f5f-8187-1c55420a950c" />

### Successful VCF run

<img width="826" height="468" alt="Successful VCF run" src="https://github.com/user-attachments/assets/abfc9de9-42f6-485f-9c5f-25880b4dae56" />

### Log file output

<img width="882" height="636" alt="LocusLift log file" src="https://github.com/user-attachments/assets/2c47696c-ebdf-4c5f-bc8a-447d722315ea" />

---

## References

[1] Monti, R. (2023). HapMap3-1KG: 1000 Genomes processed genotypes in PLINK format. figshare. https://doi.org/10.6084/m9.figshare.20802700.v1

[2] UCSC BED format (0-based `chromStart`, half-open intervals): https://genome.ucsc.edu/FAQ/FAQformat.html

[3] VCF specification v4.3 (1-based POS): https://samtools.github.io/hts-specs/VCFv4.3.pdf

[4] UCSC liftOver usage and caveats: https://genome.ucsc.edu/FAQ/FAQdownloads

[5] bcftools manual: https://samtools.github.io/bcftools/bcftools.html

[6] Picard LiftoverVcf (allele-aware VCF liftover): https://gatk.broadinstitute.org/hc/en-us/articles/360036363632-LiftoverVcf-Picard
