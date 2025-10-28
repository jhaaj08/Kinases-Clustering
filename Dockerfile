# Production Dockerfile for Kinase Classification Web App
# Optimized for Cloud Run / App Runner / any container platform

FROM python:3.12-slim as base

# Install system dependencies
RUN apt-get update && apt-get install -y \
    hmmer \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy requirements first (for layer caching)
COPY webapp/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY webapp/app.py webapp/predictor.py webapp/motif_highlighter.py ./
COPY extract_motif_features.py ./

# Copy model artifacts and data
# In production, consider downloading these from GCS/S3 instead
COPY supervised_results/logistic_regression_model.joblib supervised_results/logistic_regression_model.joblib
COPY kinases_domains_e0.01_layers_mid/ kinases_domains_e0.01_layers_mid/
COPY data/processed/kinases_domains_e0.01.csv data/processed/kinases_domains_e0.01.csv
COPY data/hmm_profiles/PF00069.hmm data/hmm_profiles/PF00069.hmm

# Create non-root user for security
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

# Environment variables
ENV PYTHONUNBUFFERED=1
ENV GRADIO_SERVER_NAME="0.0.0.0"
ENV GRADIO_SERVER_PORT=7860

# Expose port
EXPOSE 7860

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:7860/ || exit 1

# Run the application
CMD ["python", "app.py"]

