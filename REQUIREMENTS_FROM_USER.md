# What I Need From You to Deploy

## Immediate Actions (Choose Your Platform)

### Option 1: Google Cloud Run (Recommended - Best for Multiple Tools)

**What you need:**
1. **Google Cloud account** with billing enabled
2. **Install gcloud CLI**: https://cloud.google.com/sdk/docs/install
3. **Your domain name** (if you want custom domain)

**Steps:**
```bash
# 1. Login and set project
gcloud auth login
gcloud config set project YOUR_PROJECT_ID

# 2. Enable required APIs
gcloud services enable run.googleapis.com
gcloud services enable containerregistry.googleapis.com

# 3. Test locally first
cd /Users/ajikumar/codefiles/Kinases-Clustering
docker-compose up --build
# Visit http://localhost:7860 to verify

# 4. Build and deploy
docker build -t gcr.io/YOUR_PROJECT_ID/kinase-classifier:latest .
docker push gcr.io/YOUR_PROJECT_ID/kinase-classifier:latest

gcloud run deploy kinase-classifier \
  --image gcr.io/YOUR_PROJECT_ID/kinase-classifier:latest \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --port 7860 \
  --memory 8Gi \
  --cpu 4

# You'll get a URL like: https://kinase-classifier-xxx-uc.a.run.app
```

**For custom domain** (e.g., `kinase.yourdomain.com`):
```bash
gcloud run domain-mappings create \
  --service kinase-classifier \
  --domain kinase.yourdomain.com \
  --region us-central1
```
Then add the DNS records it tells you about.

---

### Option 2: AWS App Runner (Alternative)

**What you need:**
1. **AWS account**
2. **AWS CLI installed and configured**
3. **Create ECR repository** (via AWS Console)

**Steps:**
```bash
# 1. Login to ECR
aws ecr get-login-password --region us-east-1 | \
  docker login --username AWS --password-stdin \
  YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com

# 2. Build and push
docker build -t YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/kinase-classifier:latest .
docker push YOUR_ACCOUNT.dkr.ecr.us-east-1.amazonaws.com/kinase-classifier:latest

# 3. Create service via AWS Console:
# - Go to App Runner
# - Create service from ECR
# - Select your image
# - Set port to 7860
# - Set CPU: 2 vCPU, Memory: 8 GB
```

---

### Option 3: Render (Easiest - No CLI Required)

**What you need:**
1. **Render account** (free tier available)
2. **GitHub repo** (push this code)

**Steps:**
1. Go to https://render.com
2. Click "New +" → "Web Service"
3. Connect your GitHub repo
4. Settings:
   - Environment: Docker
   - Region: Choose closest
   - Instance Type: Standard (4GB RAM minimum)
   - Port: 7860
5. Click "Create Web Service"

You'll get a URL like: `https://kinase-classifier.onrender.com`

---

## GitHub Actions CI/CD Setup (Optional but Recommended)

**What you need:**
1. Push code to GitHub
2. Enable GitHub Actions (it's automatic)

**What it does:**
- Automatically builds Docker image on every push to `main`
- Pushes to GitHub Container Registry
- Tags versions when you create git tags

**No additional setup needed!** Just push to GitHub.

---

## For Your Future Multi-Tool Website

### Decision: Subpath vs Subdomain

**Subpath** (e.g., `yourdomain.com/tools/kinase`):
- ✅ Cleaner URL structure
- ✅ Single SSL cert
- ❌ Requires load balancer setup (more complex)
- **Best for**: Professional multi-tool site

**Subdomain** (e.g., `kinase.yourdomain.com`):
- ✅ Simpler setup
- ✅ Each tool independent
- ❌ Multiple DNS records
- **Best for**: Getting started quickly

### Recommendation
Start with **subdomains** now, migrate to **subpaths** later when you have 5+ tools.

---

## Cost Estimates

### Google Cloud Run
- **Idle**: $0/month (scales to zero)
- **Light use** (100 predictions/day): ~$5-10/month
- **Moderate use** (1000 predictions/day): ~$20-40/month

### AWS App Runner
- **Always running**: ~$25-50/month minimum
- **Higher traffic**: Scales automatically

### Render
- **Free tier**: Available but limited (512MB RAM - too small for this)
- **Starter**: $7/month (1GB RAM - might work for light use)
- **Standard**: $25/month (4GB RAM - recommended)

---

## What I've Set Up For You

✅ **Production-ready Dockerfile** - Optimized, secure, with health checks  
✅ **Docker Compose** - For local testing  
✅ **GitHub Actions CI/CD** - Automatic builds and deployments  
✅ **Subpath support** - Ready for `yourdomain.com/tools/kinase`  
✅ **Environment configuration** - Easy to customize  
✅ **Complete deployment guide** - Step-by-step for all platforms  

---

## Quick Start (Test Locally Right Now)

```bash
cd /Users/ajikumar/codefiles/Kinases-Clustering

# Option 1: Docker Compose (easiest)
docker-compose up --build

# Option 2: Direct Docker
docker build -t kinase-classifier .
docker run -p 7860:7860 kinase-classifier

# Visit http://localhost:7860
```

---

## Next Steps

1. **Test locally** (5 minutes)
   ```bash
   docker-compose up --build
   ```

2. **Choose platform** (see options above)
   - Easiest: Render (no CLI)
   - Best for future: Google Cloud Run
   - AWS users: App Runner

3. **Deploy** (15-30 minutes)
   - Follow platform-specific steps above

4. **Get your URL** and test

5. **Optional: Set up custom domain**

---

## Questions to Answer

Please let me know:

1. **Which platform do you prefer?**
   - [ ] Google Cloud Run (recommended)
   - [ ] AWS App Runner
   - [ ] Render
   - [ ] Other: ___________

2. **Do you have a domain name?**
   - [ ] Yes: ___________
   - [ ] No (will use platform URL for now)

3. **Subpath or subdomain preference?**
   - [ ] Subpath: `yourdomain.com/tools/kinase`
   - [ ] Subdomain: `kinase.yourdomain.com`
   - [ ] Don't care / decide later

4. **GitHub repo?**
   - [ ] Yes, it's at: ___________
   - [ ] No, need to create one

---

## I Can Help You With

- Setting up the cloud account
- Running the deployment commands
- Configuring custom domains
- Setting up CI/CD
- Troubleshooting any issues

Just let me know which platform you choose and I'll guide you through it! 🚀

