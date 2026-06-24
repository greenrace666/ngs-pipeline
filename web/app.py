import os
import json
import uuid
from datetime import datetime
from pathlib import Path
from flask import Flask, render_template, request, jsonify, send_file
import subprocess
import threading
import yaml

app = Flask(__name__, template_folder='templates')
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB max upload

# Job tracking
JOBS = {}
UPLOAD_FOLDER = Path('./uploads')
RESULTS_FOLDER = Path('./web_results')
UPLOAD_FOLDER.mkdir(exist_ok=True)
RESULTS_FOLDER.mkdir(exist_ok=True)

class PipelineJob:
    def __init__(self, job_id):
        self.job_id = job_id
        self.status = 'pending'  # pending, running, completed, failed
        self.progress = 0
        self.output_dir = RESULTS_FOLDER / job_id
        self.config = {}
        self.created_at = datetime.now()
        self.started_at = None
        self.completed_at = None
        self.error = None
        
    def to_dict(self):
        return {
            'job_id': self.job_id,
            'status': self.status,
            'progress': self.progress,
            'created_at': self.created_at.isoformat(),
            'started_at': self.started_at.isoformat() if self.started_at else None,
            'completed_at': self.completed_at.isoformat() if self.completed_at else None,
            'error': self.error,
        }

def run_pipeline_job(job_id, fastq_file, reference_file, config):
    """Run pipeline in background"""
    job = JOBS[job_id]
    job.status = 'running'
    job.started_at = datetime.now()
    job.progress = 10
    
    try:
        # Prepare output directory
        job.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy files to job directory
        import shutil
        job_fastq = job.output_dir / Path(fastq_file).name
        job_ref = job.output_dir / Path(reference_file).name
        shutil.copy(fastq_file, job_fastq)
        shutil.copy(reference_file, job_ref)
        
        job.progress = 20
        
        # Update config with paths
        config['fastq_input'] = str(job_fastq)
        config['reference_genome'] = str(job_ref)
        config['output_dir'] = str(job.output_dir)
        
        # Write config file
        config_file = job.output_dir / 'pipeline_config.yaml'
        with open(config_file, 'w') as f:
            yaml.dump(config, f)
        
        job.progress = 30
        
        # Run pipeline
        cmd = ['ngs-pipeline', 'run', '--config', str(config_file), '--output', str(job.output_dir)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        
        job.progress = 90
        
        if result.returncode == 0:
            job.status = 'completed'
            job.progress = 100
        else:
            job.status = 'failed'
            job.error = result.stderr[-500:] if result.stderr else 'Unknown error'
        
        job.completed_at = datetime.now()
        
    except subprocess.TimeoutExpired:
        job.status = 'failed'
        job.error = 'Pipeline execution timed out'
        job.completed_at = datetime.now()
    except Exception as e:
        job.status = 'failed'
        job.error = str(e)
        job.completed_at = datetime.now()

@app.route('/')
def index():
    """Home page"""
    return render_template('index.html')

@app.route('/api/job', methods=['POST'])
def create_job():
    """Create a new pipeline job"""
    try:
        job_id = str(uuid.uuid4())[:8]
        
        # Get uploaded files
        if 'fastq' not in request.files or 'reference' not in request.files:
            return jsonify({'error': 'Missing files'}), 400
        
        fastq_file = request.files['fastq']
        reference_file = request.files['reference']
        
        if fastq_file.filename == '' or reference_file.filename == '':
            return jsonify({'error': 'Empty files'}), 400
        
        # Save uploads
        fastq_path = UPLOAD_FOLDER / f"{job_id}_{fastq_file.filename}"
        ref_path = UPLOAD_FOLDER / f"{job_id}_{reference_file.filename}"
        fastq_file.save(fastq_path)
        reference_file.save(ref_path)
        
        # Get config
        config = {
            'threads': int(request.form.get('threads', 4)),
            'enable_annotation': request.form.get('enable_annotation') == 'true',
            'fastqc_enabled': True,
        }
        
        # Create job
        job = PipelineJob(job_id)
        job.config = config
        JOBS[job_id] = job
        
        # Run pipeline in background
        thread = threading.Thread(target=run_pipeline_job
thread.daemon = True
        thread.start()
        
        return jsonify({'job_id': job_id}), 201
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/job/<job_id>', methods=['GET'])
def get_job(job_id):
    """Get job status"""
    if job_id not in JOBS:
        return jsonify({'error': 'Job not found'}), 404
    
    job = JOBS[job_id]
    return jsonify(job.to_dict())

@app.route('/api/jobs', methods=['GET'])
def list_jobs():
    """List all jobs"""
    jobs = [job.to_dict() for job in JOBS.values()]
    return jsonify(jobs)

@app.route('/api/job/<job_id>/report', methods=['GET'])
def get_report(job_id):
    """Download report"""
    if job_id not in JOBS:
        return jsonify({'error': 'Job not found'}), 404
    
    job = JOBS[job_id]
    report_path = job.output_dir / 'report.html'
    
    if not report_path.exists():
        return jsonify({'error': 'Report not found'}), 404
    
    return send_file(report_path, as_attachment=True)

@app.route('/api/job/<job_id>/download', methods=['GET'])
def download_results(job_id):
    """Download all results as ZIP"""
    if job_id not in JOBS:
        return jsonify({'error': 'Job not found'}), 404
    
    job = JOBS[job_id]
    
    if not job.output_dir.exists():
        return jsonify({'error': 'Results not found'}), 404
    
    # Create ZIP file
    import zipfile
    zip_path = RESULTS_FOLDER / f"{job_id}_results.zip"
    
    with zipfile.ZipFile(zip_path, 'w') as zf:
        for file in job.output_dir.rglob('*'):
            if file.is_file():
                zf.write(file, arcname=file.relative_to(job.output_dir))
    
    return send_file(zip_path, as_attachment=True)

@app.route('/jobs')
def jobs_page():
    """Jobs page"""
    return render_template('jobs.html')

@app.route('/submit')
def submit_page():
    """Submit job page"""
    return render_template('submit.html')

@app.errorhandler(404)
def not_found(error):
    return jsonify({'error': 'Not found'}), 404

@app.errorhandler(500)
def internal_error(error):
    return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
