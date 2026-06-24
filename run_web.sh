#!/bin/bash

# Activate virtual environment
source venv/bin/activate

# Install web dependencies
pip install -r web/requirements.txt

# Run Flask app
python web/app.py
