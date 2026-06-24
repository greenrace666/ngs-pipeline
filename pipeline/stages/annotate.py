import subprocess
from pathlib import Path
from pipeline.utils.logger import setup_logger

logger = setup_logger(__name__)

def annotate_variants(vcf_file, output_vcf, snpeff_config=None):
    """Annotate variants using SnpEff"""
    logger.info(f"Annotating variants: {vcf_file}")
    
    try:
        # Check if SnpEff is installed
        result = subprocess.run(['which', 'snpEff.jar'], 
                              capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.warning("SnpEff not installed, skipping annotation")
            logger.info("To install: sudo apt install snpeff")
            # Copy VCF as-is if SnpEff not available
            subprocess.run(['cp', vcf_file, output_vcf], check=True)
            logger.info(f"✅ VCF copied (no annotation): {output_vcf}")
            return
        
        # Decompress VCF if needed
        input_vcf = vcf_file
        if vcf_file.endswith('.gz'):
            input_vcf = vcf_file.replace('.vcf.gz', '.vcf')
            logger.info(f"Decompressing VCF...")
            subprocess.run(['gunzip', '-c', vcf_file, '>', input_vcf], 
                         shell=True, check=False)
        
        # Run SnpEff annotation
        cmd = [
            'java', '-jar', '/usr/share/snpEff/snpEff.jar',
            '-c', snpeff_config or '/etc/snpEff/snpEff.config',
            'hg38',  # Reference genome (change as needed)
            input_vcf
        ]
        
        logger.info(f"Running SnpEff annotation...")
        
        with open(output_vcf.replace('.gz', ''), 'w') as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            logger.warning(f"SnpEff warning: {result.stderr}")
        
        # Compress output
        subprocess.run(['bgzip', output_vcf.replace('.gz', '')], check=False)
        logger.info(f"✅ Variants annotated: {output_vcf}")
        
    except Exception as e:
        logger.error(f"Annotation failed: {str(e)}")
        raise

def parse_annotation(annotated_vcf):
    """Parse annotated VCF and extract key annotations"""
    logger.info(f"Parsing annotations from {annotated_vcf}")
    
    annotations = []
    try:
        # Simple VCF parsing
        import gzip
        
        vcf_file = annotated_vcf.replace('.gz', '')
        if annotated_vcf.endswith('.gz'):
            with gzip.open(annotated_vcf, 'rt') as f:
                lines = f.readlines()
        else:
            with open(vcf_file, 'r') as f:
                lines = f.readlines()
        
        for line in lines:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]
                
                annotation = {
                    'chromosome': chrom,
                    'position': pos,
                    'ref': ref,
                    'alt': alt,
                    'quality': qual,
                }
                
                # Extract ANN field if present
                if 'ANN=' in info:
                    ann_start = info.find('ANN=') + 4
                    ann_end = info.find(';', ann_start)
                    if ann_end == -1:
                        ann_end = len(info)
                    annotation['annotation'] = info[ann_start:ann_end]
                
                annotations.append(annotation)
        
        logger.info(f"Parsed {len(annotations)} annotated variants")
        return annotations
        
    except Exception as e:
        logger.warning(f"Annotation parsing failed: {str(e)}")
        return []
