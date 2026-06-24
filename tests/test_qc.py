import pytest
from pathlib import Path
from pipeline.stages.qc import parse_fastq_basic

def test_parse_fastq():
    """Test FASTQ parsing returns correct stats"""
    test_fastq = "data/test.fastq"
    
    stats = parse_fastq_basic(test_fastq)
    
    # Verify all required keys are present
    assert 'reads' in stats
    assert 'total_bases' in stats
    assert 'avg_quality' in stats
    
    # Verify values are correct for test data (2 reads, 40 bases each)
    assert stats['reads'] == 2
    assert stats['total_bases'] == 80
    assert stats['avg_quality'] == 40.0

def test_parse_fastq_empty_file(tmp_path):
    """Test parsing empty FASTQ file"""
    empty_fastq = tmp_path / "empty.fastq"
    empty_fastq.write_text("")
    
    stats = parse_fastq_basic(str(empty_fastq))
    
    assert stats['reads'] == 0
    assert stats['total_bases'] == 0
