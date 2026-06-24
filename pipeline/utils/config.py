import yaml
from pathlib import Path

def load_config(config_path):
    """Load YAML configuration file"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    required_fields = ['fastq_input', 'reference_genome', 'output_dir']
    for field in required_fields:
        if field not in config:
            raise ValueError(f"Missing required config field: {field}")
    
    return config

def get_default_config():
    """Return default configuration"""
    return {
        'fastq_input': '',
        'reference_genome': '',
        'output_dir': './results',
        'threads': 4,
        'align_tool': 'bwa',
        'variant_caller': 'samtools',
        'min_depth': 10,
    }
