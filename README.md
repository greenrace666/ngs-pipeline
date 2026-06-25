# NGS Pipeline

NGS Pipeline is a small Windows-friendly sequencing analysis app with both a CLI and a Streamlit UI. It performs basic read QC, lightweight alignment, variant calling, optional annotation, and report generation.

## Features

- FASTQ quality control with read and base summaries
- FASTA reference loading and in-memory k-mer indexing
- Seeded local alignment on both forward and reverse-complement strands
- Variant calling from aligned reads
- Optional real variant annotation via MyVariant.info
- Markdown report generation
- Streamlit web interface for submitting and reviewing jobs
- CLI workflow for batch runs
- Designed to run on Windows with standard Python tooling

## Install

This project is intended to be installed with `uv`.

```powershell
uv sync
```

If you only want the runtime dependencies without creating a full project environment:

```powershell
uv pip install --system biopython click levenshtein myvariant pyyaml streamlit vcfpy
```

## Quickstart

Fastest way to install and launch the latest Windows build from GitHub Releases:

```powershell
irm https://github.com/greenrace666/ngs-pipeline/releases/latest/download/installer.ps1 | iex
```

Safer two-step option:

```powershell
irm https://github.com/greenrace666/ngs-pipeline/releases/latest/download/installer.ps1 -OutFile installer.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File .\installer.ps1
```

## Input

The pipeline expects:

- `FASTQ` or `FQ` reads
- `FASTA` or `FA` reference genome
- YAML config for CLI execution

### YAML Config Example

```yaml
fastq_input: C:\data\reads.fastq
reference_genome: C:\data\reference.fasta
threads: 4
enable_annotation: true
annotation_assembly: hg38
```

Required keys:

- `fastq_input`
- `reference_genome`

Optional keys:

- `threads` default: `4`
- `enable_annotation` default: `false`
- `annotation_assembly` default: `hg38`; supported values depend on MyVariant.info, commonly `hg38` and `hg19`

## Output

Each run writes a results directory containing:

- `alignments.tsv` - alignment summary table
- `variants.vcf.gz` - called variants
- `variants_annotated.vcf.gz` - annotated VCF when annotation is enabled
- `report.md` - Markdown analysis report

The Streamlit job view also stores:

- `job_info.json` - job status and progress
- `*_results.zip` - downloadable archive of completed outputs

## CLI

Show help:

```powershell
python main.py --help
python main.py run --help
```

Run the pipeline:

```powershell
python main.py run --config config.yml --output results
```

## Streamlit UI

Launch the web app:

```powershell
python main.py
```

or, if you prefer the standard Streamlit entrypoint:

```powershell
streamlit run main.py
```

## Behavior

- Reads are aligned against the reference using a simple seed-and-compare strategy.
- Alignments are stored as TSV rather than BAM, which keeps the implementation compact and Windows-safe.
- Variant calls are written to compressed VCF output.
- Annotation uses the `myvariant` Python client to query MyVariant.info and adds real annotation fields such as `ANN`, `MV_GENE`, `MV_CLNSIG`, `MV_CADD_PHRED`, and `MV_AF` when available. This requires internet access when annotation is enabled.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
