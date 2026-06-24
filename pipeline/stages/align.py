import subprocess
from pathlib import Path
from pipeline.utils.logger import setup_logger

logger = setup_logger(__name__)

def index_reference(reference_fasta):
    """Index reference genome with BWA"""
    logger.info(f"Indexing reference genome: {reference_fasta}")
    
    cmd = ['bwa', 'index', reference_fasta]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"BWA indexing failed: {result.stderr}")
        raise RuntimeError(f"BWA indexing error: {result.stderr}")
    
    logger.info("? Reference indexed")

def align_reads(fastq_file, reference_fasta, output_bam, threads=4):
    """Align reads to reference using BWA"""
    logger.info(f"Aligning {fastq_file} to {reference_fasta}...")
    
    cmd = [
        'bwa', 'mem',
        '-t', str(threads),
        reference_fasta,
        fastq_file
    ]
    
    logger.info(f"Running: {' '.join(cmd)}")
    
    try:
        with open(output_bam, 'w') as out_f:
            bwa_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            sort_cmd = ['samtools', 'sort', '-o', output_bam, '-']
            sort_proc = subprocess.Popen(sort_cmd, stdin=bwa_proc.stdout, 
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            bwa_proc.stdout.close()
            
            stdout, stderr = sort_proc.communicate()
            
            if sort_proc.returncode != 0:
                logger.error(f"Samtools sort failed: {stderr.decode()}")
                raise RuntimeError(f"Sort error: {stderr.decode()}")
        
        subprocess.run(['samtools', 'index', output_bam], check=True)
        logger.info(f"? Alignment complete: {output_bam}")
        
    except Exception as e:
        logger.error(f"Alignment failed: {str(e)}")
        raise
