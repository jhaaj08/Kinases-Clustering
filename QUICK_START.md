# Quick Start: Deploy in 15 Minutes

## Step 1: Test Locally (2 minutes)

```bash
cd /Users/ajikumar/codefiles/Kinases-Clustering
docker-compose up --build
```

Open http://localhost:7860 - if it works, you're ready to deploy!

Press `Ctrl+C` to stop.

---

## Step 2: Choose Your Path

### Path A: Google Cloud Run (Recommended)

**Prerequisites:**
```bash
# Install gcloud CLI (if not installed)
# macOS:
brew install google-cloud-sdk

# Login
gcloud auth login
gcloud config set project YOUR_PROJECT_ID
```

**Deploy (5 commands):**
```bash
cd /Users/ajikumar/codefiles/Kinases-Clustering

# Build
docker build -t gcr.io/YOUR_PROJECT_ID/kinase-classifier:latest .

# Push
docker push gcr.io/YOUR_PROJECT_ID/kinase-classifier:latest

# Deploy
gcloud run deploy kinase-classifier \
  --image gcr.io/YOUR_PROJECT_ID/kinase-classifier:latest \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --port 7860 \
  --memory 8Gi \
  --cpu 4 \
  --timeout 300

# Get URL
gcloud run services describe kinase-classifier --region us-central1 --format 'value(status.url)'
```

**Done!** You'll get a URL like: `https://kinase-classifier-xxx-uc.a.run.app`

---

### Path B: Render (No CLI Required)

1. Push code to GitHub
2. Go to https://render.com → New → Web Service
3. Connect GitHub repo
4. Settings:
   - Environment: **Docker**
   - Instance Type: **Standard (4GB)**
   - Port: **7860**
5. Click "Create Web Service"

**Done!** You'll get a URL like: `https://kinase-classifier.onrender.com`

---

## Step 3: Custom Domain (Optional)

### For Cloud Run:
```bash
gcloud run domain-mappings create \
  --service kinase-classifier \
  --domain kinase.yourdomain.com \
  --region us-central1
```

Then add the DNS records it shows you.

### For Render:
1. Go to service settings → Custom Domains
2. Add your domain
3. Update DNS as instructed

---

## Troubleshooting

**"Out of memory"**
→ Increase memory: `--memory 8Gi` (Cloud Run) or upgrade instance (Render)

**"Container failed to start"**
→ Test locally first: `docker-compose up --build`

**"Too slow"**
→ Keep instance warm: `--min-instances 1` (Cloud Run)

---

## What's Next?

- ✅ Your app is live!
- 📊 Monitor at: Cloud Console (GCP) or Render Dashboard
- 🔄 Auto-deploys on git push (if using GitHub Actions)
- 🚀 Add more tools following the same pattern

---

**Need help?** See `DEPLOYMENT.md` for detailed guides or `REQUIREMENTS_FROM_USER.md` for what I need from you.

