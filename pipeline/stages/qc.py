import subprocess
import json
from pathlib import Path
from pipeline.utils.logger import setup_logger

logger = setup_logger(__name__)

def run_fastqc(fastq_file, output_dir):
    """Run FastQC on input FASTQ file"""
    logger.info(f"Running FastQC on {fastq_file}...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        'fastqc',
        '-o', str(output_dir),
        '--quiet',
        fastq_file
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"FastQC failed: {result.stderr}")
        raise RuntimeError(f"FastQC error: {result.stderr}")
    
    logger.info(f"? FastQC completed. Results in {output_dir}")
    return output_dir

def parse_fastq_basic(fastq_file):
    """Parse FASTQ file and return basic statistics"""
    read_count = 0
    total_bases = 0
    quality_scores = []
    
    with open(fastq_file, 'r') as f:
        lines = f.readlines()
        
    for i in range(1, len(lines), 4):
        if i + 2 < len(lines):
            seq = lines[i].strip()
            qual = lines[i + 2].strip()
            
            read_count += 1
            total_bases += len(seq)
            quality_scores.extend([ord(q) - 33 for q in qual])
    
    avg_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
    
    return {
        'reads': read_count,
        'total_bases': total_bases,
        'avg_quality': round(avg_quality, 2),
        'reads_per_mb': round(read_count / (total_bases / 1e6), 0) if total_bases > 0 else 0,
    }
