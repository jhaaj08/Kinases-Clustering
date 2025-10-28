# Deployment Guide: Kinase Classifier Web App

This guide covers deploying the Kinase Classifier as part of a larger computational biology website with multiple tools.

## Architecture Overview

```
yourdomain.com
├── /                    → Main website (Vercel/Netlify)
├── /tools/kinase        → This app (Cloud Run/App Runner)
├── /tools/other-tool    → Future tools (Cloud Run/App Runner)
└── /docs                → Documentation
```

---

## Quick Start: Local Development

### 1. Using Docker Compose (Recommended)

```bash
# Build and run
docker-compose up --build

# Access at http://localhost:7860
```

### 2. Direct Python

```bash
cd /path/to/Kinases-Clustering
python webapp/app.py
```

---

## Production Deployment Options

### Option A: Google Cloud Run (Recommended)

**Why Cloud Run?**
- Autoscaling to zero (pay per use)
- Built-in HTTPS and custom domains
- Easy subpath routing with Load Balancer
- Minimal ops, great for multiple services

#### Prerequisites
- Google Cloud account
- `gcloud` CLI installed
- Docker installed

#### Step 1: Build and Push Image

```bash
# Set your project
export PROJECT_ID=your-project-id
export REGION=us-central1

# Configure Docker for GCR
gcloud auth configure-docker

# Build and push
docker build -t gcr.io/$PROJECT_ID/kinase-classifier:latest .
docker push gcr.io/$PROJECT_ID/kinase-classifier:latest
```

#### Step 2: Deploy to Cloud Run

```bash
# Deploy service
gcloud run deploy kinase-classifier \
  --image gcr.io/$PROJECT_ID/kinase-classifier:latest \
  --platform managed \
  --region $REGION \
  --allow-unauthenticated \
  --port 7860 \
  --memory 8Gi \
  --cpu 4 \
  --timeout 300 \
  --min-instances 0 \
  --max-instances 10

# Get the service URL
gcloud run services describe kinase-classifier --region $REGION --format 'value(status.url)'
```

#### Step 3: Set Up Custom Domain with Subpath

For `yourdomain.com/tools/kinase`:

1. **Create HTTPS Load Balancer**:
```bash
# Create serverless NEG
gcloud compute network-endpoint-groups create kinase-neg \
  --region=$REGION \
  --network-endpoint-type=serverless \
  --cloud-run-service=kinase-classifier

# Create backend service
gcloud compute backend-services create kinase-backend \
  --global

# Add NEG to backend
gcloud compute backend-services add-backend kinase-backend \
  --global \
  --network-endpoint-group=kinase-neg \
  --network-endpoint-group-region=$REGION

# Create URL map with path matcher
gcloud compute url-maps create compbio-lb \
  --default-service kinase-backend

# Add path rule for /tools/kinase
gcloud compute url-maps add-path-matcher compbio-lb \
  --path-matcher-name=tools \
  --default-service=kinase-backend \
  --path-rules="/tools/kinase/*=kinase-backend"

# Create HTTPS proxy and forwarding rule
gcloud compute ssl-certificates create compbio-cert \
  --domains=yourdomain.com

gcloud compute target-https-proxies create compbio-proxy \
  --url-map=compbio-lb \
  --ssl-certificates=compbio-cert

gcloud compute forwarding-rules create compbio-https-rule \
  --global \
  --target-https-proxy=compbio-proxy \
  --ports=443
```

2. **Update Cloud Run service with ROOT_PATH**:
```bash
gcloud run services update kinase-classifier \
  --region $REGION \
  --set-env-vars ROOT_PATH=/tools/kinase
```

3. **Point your domain** to the load balancer IP

#### Alternative: Subdomain Approach (Simpler)

For `kinase.yourdomain.com`:

```bash
# Map custom domain directly to Cloud Run
gcloud run services add-iam-policy-binding kinase-classifier \
  --region=$REGION \
  --member="allUsers" \
  --role="roles/run.invoker"

# Map domain
gcloud run domain-mappings create \
  --service kinase-classifier \
  --domain kinase.yourdomain.com \
  --region $REGION
```

Then add DNS A/AAAA records pointing to the Cloud Run IP.

---

### Option B: AWS App Runner

#### Prerequisites
- AWS account
- AWS CLI configured
- ECR repository created

#### Deploy

```bash
# Login to ECR
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com

# Build and push
docker build -t YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/kinase-classifier:latest .
docker push YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/kinase-classifier:latest

# Create App Runner service (via console or CLI)
aws apprunner create-service \
  --service-name kinase-classifier \
  --source-configuration '{
    "ImageRepository": {
      "ImageIdentifier": "YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/kinase-classifier:latest",
      "ImageRepositoryType": "ECR",
      "ImageConfiguration": {
        "Port": "7860",
        "RuntimeEnvironmentVariables": {
          "ROOT_PATH": "/tools/kinase"
        }
      }
    },
    "AutoDeploymentsEnabled": true
  }' \
  --instance-configuration '{
    "Cpu": "2 vCPU",
    "Memory": "8 GB"
  }'
```

---

### Option C: Simple Platforms (Render/Railway/Fly.io)

#### Render

1. Connect your GitHub repo
2. Create new "Web Service"
3. Select "Docker" as environment
4. Set port to `7860`
5. Add environment variable: `ROOT_PATH=/tools/kinase` (if using subpath)
6. Deploy

#### Railway

```bash
# Install Railway CLI
npm i -g @railway/cli

# Login and init
railway login
railway init

# Deploy
railway up
```

---

## CI/CD with GitHub Actions

The included `.github/workflows/deploy.yml` automatically:
- Builds Docker image on push to `main`
- Pushes to GitHub Container Registry (ghcr.io)
- Tags with version on git tags

### Enable GitHub Container Registry

1. Go to repo Settings → Secrets and variables → Actions
2. No secrets needed for GHCR (uses `GITHUB_TOKEN`)

### Optional: Auto-deploy to Cloud Run

Uncomment the Cloud Run deployment step in `.github/workflows/deploy.yml` and add these secrets:

```bash
# Create service account
gcloud iam service-accounts create github-actions \
  --display-name="GitHub Actions"

# Grant permissions
gcloud projects add-iam-policy-binding $PROJECT_ID \
  --member="serviceAccount:github-actions@$PROJECT_ID.iam.gserviceaccount.com" \
  --role="roles/run.admin"

# Create key
gcloud iam service-accounts keys create key.json \
  --iam-account=github-actions@$PROJECT_ID.iam.gserviceaccount.com

# Add to GitHub Secrets:
# GCP_PROJECT_ID
# GCP_SA_KEY (contents of key.json)
```

---

## Scaling to Multiple Tools

### 1. Standardize Structure

```
your-compbio-site/
├── tools/
│   ├── kinase-classifier/
│   │   ├── Dockerfile
│   │   ├── app.py
│   │   └── requirements.txt
│   ├── protein-retrieval/
│   │   ├── Dockerfile
│   │   ├── app.py
│   │   └── requirements.txt
│   └── ...
├── frontend/              # Main website
└── infrastructure/        # Terraform/Pulumi configs
```

### 2. Shared Base Image

Create a base image with common dependencies:

```dockerfile
# base.Dockerfile
FROM python:3.12-slim
RUN apt-get update && apt-get install -y hmmer curl && rm -rf /var/lib/apt/lists/*
RUN pip install torch fair-esm scikit-learn gradio
```

Then each tool extends it:
```dockerfile
FROM your-registry/compbio-base:latest
COPY requirements.txt .
RUN pip install -r requirements.txt
# ... rest of Dockerfile
```

### 3. Infrastructure as Code (Optional)

Use Terraform or Pulumi to manage all services:

```hcl
# terraform/main.tf
module "kinase_classifier" {
  source = "./modules/cloud-run-service"
  name = "kinase-classifier"
  image = "gcr.io/PROJECT/kinase-classifier:latest"
  path = "/tools/kinase"
}

module "other_tool" {
  source = "./modules/cloud-run-service"
  name = "other-tool"
  image = "gcr.io/PROJECT/other-tool:latest"
  path = "/tools/other"
}
```

---

## Monitoring & Operations

### Health Checks

Built into Dockerfile:
```dockerfile
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:7860/ || exit 1
```

### Logging

Cloud Run/App Runner automatically collect logs. View with:
```bash
# Cloud Run
gcloud run services logs read kinase-classifier --region $REGION

# AWS
aws logs tail /aws/apprunner/kinase-classifier --follow
```

### Monitoring

- **Cloud Run**: Use Google Cloud Monitoring (automatic)
- **App Runner**: Use CloudWatch (automatic)
- **Custom**: Add Sentry for error tracking

### Alerts

Set up uptime monitoring:
- Google Cloud Monitoring: Uptime checks
- AWS: CloudWatch Alarms
- External: UptimeRobot, Pingdom

---

## Security Best Practices

### 1. Authentication (Optional)

Add basic auth:
```python
# In app.py
app.launch(
    auth=("username", "password"),  # Or list of tuples for multiple users
    # ... other params
)
```

Or use environment variable:
```bash
export GRADIO_AUTH="user1:pass1,user2:pass2"
```

### 2. Rate Limiting

- **Cloud Run**: Use Cloud Armor
- **App Runner**: Use AWS WAF
- **Application-level**: Add middleware

### 3. HTTPS

All platforms provide automatic HTTPS. Enforce it:
```python
# In app.py (if behind proxy)
app.launch(
    ssl_verify=True,
    # ...
)
```

---

## Cost Optimization

### Cloud Run
- **Scale to zero**: Set `--min-instances 0`
- **Right-size**: Start with 2 CPU / 4GB, adjust based on metrics
- **Concurrency**: Set `--concurrency 10` (Gradio handles this well)

### Estimated Costs (Cloud Run)
- Idle: $0/month (scales to zero)
- Light use (100 requests/day): ~$5-10/month
- Moderate use (1000 requests/day): ~$20-40/month

---

## Troubleshooting

### Container won't start
```bash
# Test locally first
docker run -p 7860:7860 kinase-classifier

# Check logs
docker logs <container-id>
```

### Out of memory
- Increase memory limit: `--memory 8Gi` (Cloud Run)
- Check ESM-2 model loading (largest memory consumer)

### Slow cold starts
- Set `--min-instances 1` to keep one instance warm
- Or use Cloud Run "always allocated CPU" tier

### Subpath routing issues
- Ensure `ROOT_PATH` environment variable is set
- Check load balancer path rules
- Test with: `curl https://yourdomain.com/tools/kinase/`

---

## Next Steps

1. **Test locally** with Docker Compose
2. **Deploy to staging** (Cloud Run with test domain)
3. **Set up CI/CD** (GitHub Actions)
4. **Configure custom domain** (subpath or subdomain)
5. **Add monitoring** (uptime checks, error tracking)
6. **Repeat for additional tools**

---

## Support

- **Issues**: GitHub Issues
- **Documentation**: See README.md
- **Cloud Run docs**: https://cloud.google.com/run/docs
- **App Runner docs**: https://docs.aws.amazon.com/apprunner/

---

**Generated**: 2025-01-28  
**Maintainer**: See git commit history

