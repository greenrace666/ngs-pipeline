import pytest
from pipeline.utils.logger import setup_logger

def test_logger_creation():
    """Test logger setup"""
    logger = setup_logger("test_logger")
    
    assert logger is not None
    assert logger.name == "test_logger"
    assert len(logger.handlers) > 0
