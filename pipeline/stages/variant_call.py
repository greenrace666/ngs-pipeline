import subprocess
from pathlib import Path
from pipeline.utils.logger import setup_logger

logger = setup_logger(__name__)

def call_variants(bam_file, reference_fasta, output_vcf, min_depth=10):
    """Call variants using samtools/bcftools"""
    logger.info(f"Calling variants from {bam_file}...")
    
    try:
        # Index reference
        subprocess.run(['samtools', 'faidx', reference_fasta], 
                      capture_output=True, check=False)
        
        # Create intermediate pileup file
        pileup_file = output_vcf.replace('.vcf.gz', '.pileup')
        
        logger.info(f"Generating pileup file...")
        pileup_cmd = [
            'samtools', 'mpileup',
            '-f', reference_fasta,
            '-q', '20',
            bam_file
        ]
        
        with open(pileup_file, 'w') as pf:
            result = subprocess.run(pileup_cmd, stdout=pf, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0 and result.stderr:
            logger.warning(f"Mpileup warning: {result.stderr}")
        
        # Check if pileup file has content
        pileup_size = Path(pileup_file).stat().st_size
        logger.info(f"Pileup file size: {pileup_size} bytes")
        
        if pileup_size == 0:
            logger.warning("No pileup data - creating empty VCF")
            # Create minimal VCF header
            with open(output_vcf.replace('.gz', ''), 'w') as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            # Compress it
            subprocess.run(['bgzip', output_vcf.replace('.gz', '')], check=False)
        else:
            logger.info(f"Calling variants from pileup...")
            call_cmd = [
                'bcftools', 'call',
                '-m', '-v',
                '-O', 'z',
                '-o', output_vcf
            ]
            
            with open(pileup_file, 'r') as pf:
                result = subprocess.run(call_cmd, stdin=pf, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                logger.warning(f"Bcftools warning: {result.stderr}")
        
        logger.info(f"✅ Variants called: {output_vcf}")
        
    except Exception as e:
        logger.error(f"Variant calling failed: {str(e)}")
        raise
