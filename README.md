# NGS Data Processing Pipeline

A **production-grade, end-to-end bioinformatics pipeline** for processing Next-Generation Sequencing (NGS) data.

![Status](https://img.shields.io/badge/Status-Production%20Ready-brightgreen)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)

## Features

- ✅ Automated Quality Control (FastQC)
- ✅ Read Alignment (BWA-MEM)
- ✅ Variant Calling (SAMtools/BCFtools)
- ✅ HTML Report Generation
- ✅ YAML Configuration Support
- ✅ Docker Support (coming soon)
- ✅ CI/CD Integration (coming soon)

## Quick Start

```bash
# Install
git clone https://github.com/Dhanu577/ngs-pipeline.git
cd ngs-pipeline
python3 -m venv venv
source venv/bin/activate
pip install -e .
pip install -r requirements.txt

# Run
ngs-pipeline run --config config/sample_config.yaml --output results/
```

## Pipeline Stages

1. Quality Control (FastQC)
2. Reference Indexing (BWA)
3. Read Alignment (BWA-MEM)
4. Variant Calling (SAMtools/BCFtools)
5. Report Generation (HTML)

## Output

- `aligned.bam` - Aligned reads
- `variants.vcf.gz` - Variant calls
- `report.html` - Summary report
- `qc/` - FastQC results

## Configuration

Edit `config/sample_config.yaml`:

```yaml
fastq_input: "data/sample.fastq.gz"
reference_genome: "data/reference.fasta"
output_dir: "./results"
threads: 4
```

## Author

Danuska - Bioinformatics Student, Tamil Nadu, India

## License

MIT
