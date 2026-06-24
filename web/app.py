import uuid
import yaml
import zipfile
import subprocess
import threading
from datetime import datetime
from pathlib import Path
from flask import Flask, render_template, request, jsonify, send_file

app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024

JOBS = {}
# Already defined above with APP_DIR
import os
# Get the directory where app.py is located
APP_DIR = Path(__file__).parent.parent  # Go up one level from web/ to root
RESULTS_FOLDER = APP_DIR / 'web_results'
UPLOAD_FOLDER = APP_DIR / 'uploads'
UPLOAD_FOLDER.mkdir(exist_ok=True)
RESULTS_FOLDER.mkdir(exist_ok=True)

class PipelineJob:
    def __init__(self, job_id):
        self.job_id = job_id
        self.status = 'pending'
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
    job = JOBS[job_id]
    job.status = 'running'
    job.started_at = datetime.now()
    job.progress = 10
    try:
        job.output_dir.mkdir(parents=True, exist_ok=True)
        import shutil
        job_fastq = job.output_dir / Path(fastq_file).name
        job_ref = job.output_dir / Path(reference_file).name
        shutil.copy(fastq_file, job_fastq)
        shutil.copy(reference_file, job_ref)
        job.progress = 20
        config['fastq_input'] = str(job_fastq)
        config['reference_genome'] = str(job_ref)
        config['output_dir'] = str(job.output_dir)
        config_file = job.output_dir / 'pipeline_config.yaml'
        with open(config_file, 'w') as f:
            yaml.dump(config, f)
        job.progress = 30
        cmd = ['ngs-pipeline', 'run', '--config', str(config_file), '--output', str(job.output_dir)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        job.progress = 90
        if result.returncode == 0:
            job.status = 'completed'
            job.progress = 100
        else:
            job.status = 'failed'
            job.error = result.stderr[-500:] if result.stderr else 'Unknown error'
    except subprocess.TimeoutExpired:
        job.status = 'failed'
        job.error = 'Pipeline timed out'
    except Exception as e:
        job.status = 'failed'
        job.error = str(e)
    finally:
        job.completed_at = datetime.now()

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/job', methods=['POST'])
def create_job():
    try:
        if 'fastq' not in request.files or 'reference' not in request.files:
            return jsonify({'error': 'Missing files'}), 400
        fastq_file = request.files['fastq']
        reference_file = request.files['reference']
        if not fastq_file.filename or not reference_file.filename:
            return jsonify({'error': 'Empty filename'}), 400
        job_id = str(uuid.uuid4())[:8]
        fastq_path = UPLOAD_FOLDER / f"{job_id}_{fastq_file.filename}"
        ref_path = UPLOAD_FOLDER / f"{job_id}_{reference_file.filename}"
        fastq_file.save(fastq_path)
        reference_file.save(ref_path)
        config = {
            'threads': int(request.form.get('threads', 4)),
            'enable_annotation': request.form.get('enable_annotation') == 'on',
            'fastqc_enabled': True,
        }
        job = PipelineJob(job_id)
        job.config = config
        JOBS[job_id] = job
        t = threading.Thread(target=run_pipeline_job, args=(job_id, str(fastq_path), str(ref_path), config))
        t.daemon = True
        t.start()
        return jsonify({'job_id': job_id}), 201
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/jobs', methods=['GET'])
def list_jobs():
    jobs = sorted(JOBS.values(), key=lambda j: j.created_at, reverse=True)
    return jsonify([j.to_dict() for j in jobs])

@app.route('/api/job/<job_id>', methods=['GET'])
def get_job(job_id):
    if job_id not in JOBS:
        return jsonify({'error': 'Job not found'}), 404
    return jsonify(JOBS[job_id].to_dict())

@app.route('/api/job/<job_id>/error', methods=['GET'])
def get_job_error(job_id):
    if job_id not in JOBS:
        return jsonify({'error': 'Job not found'}), 404
    job = JOBS[job_id]
    return jsonify({
        'job_id': job_id,
        'status': job.status,
        'error': job.error,
    })


@app.route('/api/job/<job_id>/report')
def get_report(job_id):
    if job_id not in JOBS:
        return jsonify({'error': 'Job not found'}), 404
    job = JOBS[job_id]
    if not job.output_dir.exists():
        return jsonify({'error': 'No output folder'}), 404
    report = job.output_dir / 'report.html'
    if not report.exists():
        return jsonify({'error': 'Report not generated'}), 404
    return send_file(report, as_attachment=True)

@app.route('/api/job/<job_id>/download')
def download_results(job_id):
    if job_id not in JOBS:
        return jsonify({'error': 'Job not found'}), 404
    job = JOBS[job_id]
    if not job.output_dir.exists():
        return jsonify({'error': 'No results folder'}), 404
    zip_path = RESULTS_FOLDER / f"{job_id}_results.zip"
    with zipfile.ZipFile(zip_path, 'w') as zf:
        files_added = 0
        for f in job.output_dir.rglob('*'):
            if f.is_file():
                zf.write(f, arcname=f.relative_to(job.output_dir))
                files_added += 1
    if files_added == 0:
        zip_path.unlink(missing_ok=True)
        return jsonify({'error': 'No output files generated'}), 404
    return send_file(zip_path, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
