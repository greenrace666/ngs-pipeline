from pathlib import Path
import subprocess
from pipeline.stages.qc import run_fastqc, parse_fastq_basic
from pipeline.stages.align import index_reference, align_reads
from pipeline.stages.variant_call import call_variants
from pipeline.stages.report import generate_html_report
from pipeline.utils.logger import setup_logger

logger = setup_logger(__name__)

def run_full_pipeline(config, output_dir):
    """Execute full NGS pipeline"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    fastq_input = config.get('fastq_input')
    reference = config.get('reference_genome')
    threads = config.get('threads', 4)
    
    logger.info("=" * 60)
    logger.info("NGS PIPELINE STARTED")
    logger.info("=" * 60)
    
    # Stage 1: QC
    logger.info("\n[STAGE 1] Quality Control")
    qc_dir = output_dir / 'qc'
    run_fastqc(fastq_input, qc_dir)
    qc_stats = parse_fastq_basic(fastq_input)
    logger.info(f"QC Stats: {qc_stats}")
    
    # Stage 2: Index Reference
    logger.info("\n[STAGE 2] Reference Indexing")
    index_reference(reference)
    
    # Stage 3: Alignment
    logger.info("\n[STAGE 3] Alignment")
    bam_file = output_dir / 'aligned.bam'
    align_reads(fastq_input, reference, str(bam_file), threads=threads)
    
    # Stage 4: Variant Calling
    logger.info("\n[STAGE 4] Variant Calling")
    vcf_file = output_dir / 'variants.vcf.gz'
    call_variants(str(bam_file), reference, str(vcf_file))
    
    # Count variants
    try:
        result = subprocess.run(['bcftools', 'view', '-H', str(vcf_file)], 
                              capture_output=True, text=True)
        variant_count = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
    except:
        variant_count = 0
    
    logger.info(f"Variants found: {variant_count}")
    
    # Stage 5: Report
    logger.info("\n[STAGE 5] Report Generation")
    generate_html_report(output_dir, qc_stats, variant_count)
    
    logger.info("\n" + "=" * 60)
    logger.info("? PIPELINE COMPLETED SUCCESSFULLY")
    logger.info("=" * 60)
    logger.info(f"Results saved to: {output_dir}")
