import click
import yaml
from pathlib import Path
from pipeline.workflows.full_pipeline import run_full_pipeline
from pipeline.utils.logger import setup_logger

logger = setup_logger(__name__)

@click.group()
def cli():
    """NGS Data Processing Pipeline"""
    pass

@cli.command()
@click.option('--config', required=True, type=click.Path(exists=True), 
              help='Path to config YAML file')
@click.option('--output', required=True, type=click.Path(), 
              help='Output directory')
def run(config, output):
    """Run the full NGS pipeline"""
    click.echo(f"?? Starting NGS Pipeline")
    click.echo(f"?? Config: {config}")
    click.echo(f"?? Output: {output}")
    
    with open(config, 'r') as f:
        config_data = yaml.safe_load(f)
    
    Path(output).mkdir(parents=True, exist_ok=True)
    
    try:
        run_full_pipeline(config_data, output)
        click.echo("? Pipeline completed successfully!")
    except Exception as e:
        click.echo(f"? Error: {str(e)}", err=True)

if __name__ == '__main__':
    cli()
