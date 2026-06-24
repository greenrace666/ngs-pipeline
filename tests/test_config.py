import pytest
from pipeline.utils.config import load_config, get_default_config

def test_get_default_config():
    """Test default configuration is valid"""
    config = get_default_config()
    
    assert 'fastq_input' in config
    assert 'reference_genome' in config
    assert 'output_dir' in config
    assert config['threads'] == 4

def test_config_has_required_fields():
    """Test that config has all required fields"""
    config = get_default_config()
    required = ['fastq_input', 'reference_genome', 'output_dir']
    
    for field in required:
        assert field in config
