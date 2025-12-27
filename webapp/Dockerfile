# Dockerfile for Kinase Classification Web App

FROM python:3.12-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    hmmer \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY app.py predictor.py motif_highlighter.py ./

# Copy necessary data files from parent directory
# Note: In production, these should be mounted or downloaded
COPY ../supervised_results/logistic_regression_model.joblib supervised_results/
COPY ../kinases_domains_e0.01_layers_mid kinases_domains_e0.01_layers_mid/
COPY ../kinases_domains_e0.01.csv .
COPY ../PF00069.hmm .

# Expose port
EXPOSE 7860

# Set environment variables
ENV GRADIO_SERVER_NAME="0.0.0.0"
ENV GRADIO_SERVER_PORT=7860

# Run the app
CMD ["python", "app.py"]

