# HCMV Amplicon QIIME2 Nextflow Workflow

## Description

This Nextflow DSL2 workflow processes paired-end amplicon sequencing data using QIIME2. It imports demultiplexed reads from a QIIME2 manifest, merges paired reads, performs quality filtering and dereplication, filters low-frequency features, and clusters representative sequences against a reference sequence database using open-reference clustering.

The workflow is intended for targeted amplicon sequencing data, such as HCMV multi-locus amplicon data.

## Requirements

- Nextflow
- Docker or Singularity
- QIIME2 Amplicon container: `quay.io/qiime2/amplicon:2024.5`

Available profiles:

- `docker`
- `singularity`

## Input files

| Input | Description | Default |
|---|---|---|
| `--manifest` | QIIME2 paired-end FASTQ manifest file using `PairedEndFastqManifestPhred33V2` format | `manifest.txt` |
| `--ref_seqs` | QIIME2 reference sequence artifact (`.qza`) for open-reference clustering | `closed_picking_ref_07-30-2024.qza` |

The FASTQ paths listed in the manifest must be accessible from the execution environment and visible inside the selected container.

## Basic usage

```bash
nextflow run main_open.nf \
  -profile docker \
  --manifest manifest.txt \
  --ref_seqs hcmv_reference.qza \
  --outdir results
```

For Singularity, use:

```bash
nextflow run main_open.nf \
  -profile singularity \
  --manifest manifest.txt \
  --ref_seqs hcmv_reference.qza \
  --outdir results
```

If the FASTQ files are outside the project directory, provide a bind path:

```bash
nextflow run main_open.nf \
  -profile singularity \
  --manifest manifest.txt \
  --ref_seqs hcmv_reference.qza \
  --bind_path /path/to/fastq_parent_directory \
  --outdir results
```

## Common parameters

| Parameter | Description | Default |
|---|---|---|
| `--outdir` | Output directory | `results` |
| `--min_frequency` | Minimum total feature frequency retained after dereplication | `10` |
| `--min_mergelen` | Minimum merged-read length; also used as minimum dereplicated sequence length | `350` |
| `--cluster_perc` | Percent identity for open-reference clustering | `0.99` |
| `--threads` | CPU threads for read merging and clustering | `8` |
| `--bind_path` | Optional Singularity bind path for external input data | `null` |

## Example command

```bash
nextflow run main_open.nf \
  -profile docker \
  --manifest data/manifest.txt \
  --ref_seqs refs/hcmv_amplicon_reference.qza \
  --min_frequency 10 \
  --min_mergelen 350 \
  --cluster_perc 0.99 \
  --threads 8 \
  --outdir results_hcmv
```

## Expected outputs

The workflow writes QIIME2 artifacts and visualizations into the output directory:

| Directory | Contents |
|---|---|
| `01_demux/` | Imported demultiplexed paired-end reads and demux summary |
| `02_merged/` | Merged and unmerged paired reads |
| `03_quality/` | Quality-filtered reads and filtering statistics |
| `04_derep/` | Dereplicated feature table and representative sequences |
| `05_filtered_features/` | Feature table and representative sequences after low-frequency filtering |
| `06_clustering/` | Open-reference clustered feature table, clustered sequences, updated reference sequences, and QIIME2 visualizations |

Key output file types:

- `.qza`: QIIME2 artifacts for downstream analysis
- `.qzv`: QIIME2 visualizations viewable with QIIME2 View

## Resuming a run

To resume from cached completed steps after interruption or parameter-compatible changes:

```bash
nextflow run main_open.nf \
  -profile docker \
  --manifest manifest.txt \
  --ref_seqs hcmv_reference.qza \
  --outdir results \
  -resume
```
