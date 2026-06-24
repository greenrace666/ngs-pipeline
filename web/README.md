# NGS Pipeline Web UI

A Flask-based web interface for the NGS pipeline.

## Features

- Submit jobs through web interface
- Upload FASTQ and reference files
- Monitor job progress in real-time
- Download results and reports
- Beautiful responsive UI

## Installation

```bash
pip install -r web/requirements.txt
```

## Running

```bash
python app.py
```

Then open: http://localhost:5000

## Usage

1. Go to "Submit Job"
2. Upload your FASTQ and reference files
3. Configure parameters (threads, annotation)
4. Click "Submit Job"
5. Monitor progress in "My Jobs"
6. Download results when complete
