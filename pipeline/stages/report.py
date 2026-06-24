import json
from pathlib import Path
from jinja2 import Template
from pipeline.utils.logger import setup_logger

logger = setup_logger(__name__)

HTML_TEMPLATE = """
<html>
<head>
    <title>NGS Pipeline Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; }
        h1 { color: #333; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
        .success { color: green; }
    </style>
</head>
<body>
    <h1>NGS Pipeline Analysis Report</h1>
    
    <div class="section">
        <h2>Summary</h2>
        <table>
            <tr>
                <th>Metric</th>
                <th>Value</th>
            </tr>
            {% for key, value in summary.items() %}
            <tr>
                <td>{{ key }}</td>
                <td>{{ value }}</td>
            </tr>
            {% endfor %}
        </table>
    </div>
    
    <div class="section">
        <h2>Quality Control</h2>
        <pre>{{ qc_data }}</pre>
    </div>
    
    <div class="section">
        <h2 class="success">? Pipeline Complete</h2>
    </div>
</body>
</html>
"""

def generate_html_report(output_dir, qc_stats, variant_count):
    """Generate HTML report"""
    logger.info(f"Generating HTML report...")
    
    template = Template(HTML_TEMPLATE)
    
    summary = {
        'Total Reads': qc_stats.get('reads', 'N/A'),
        'Total Bases': qc_stats.get('total_bases', 'N/A'),
        'Avg Quality': qc_stats.get('avg_quality', 'N/A'),
        'Variants Called': variant_count,
    }
    
    html_content = template.render(
        summary=summary,
        qc_data=json.dumps(qc_stats, indent=2)
    )
    
    report_path = Path(output_dir) / 'report.html'
    with open(report_path, 'w') as f:
        f.write(html_content)
    
    logger.info(f"? Report generated: {report_path}")
    return report_path
