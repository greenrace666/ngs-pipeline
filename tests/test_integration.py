import pytest
import subprocess
from pathlib import Path

def test_pipeline_help():
    """Test that CLI help works"""
    result = subprocess.run(
        ['ngs-pipeline', 'run', '--help'],
        capture_output=True,
        text=True
    )
    
    assert result.returncode == 0
    assert 'Run the full NGS pipeline' in result.stdout

def test_pipeline_requires_config():
    """Test that pipeline requires config file"""
    result = subprocess.run(
        ['ngs-pipeline', 'run'],
        capture_output=True,
        text=True
    )
    
    # Should fail without required arguments
    assert result.returncode != 0
