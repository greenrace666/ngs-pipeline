from setuptools import setup, find_packages

setup(
    name="ngs-pipeline",
    version="0.1.0",
    description="End-to-end NGS data processing pipeline",
    author="Your Name",
    packages=find_packages(),
    install_requires=[
        'click>=8.0.0',
        'pyyaml>=5.4.0',
        'pandas>=1.3.0',
        'jinja2>=3.0.0',
        'plotly>=5.0.0',
        'matplotlib>=3.4.0',
    ],
    entry_points={
        'console_scripts': [
            'ngs-pipeline=pipeline.main:cli',
        ],
    },
    python_requires='>=3.8',
)
