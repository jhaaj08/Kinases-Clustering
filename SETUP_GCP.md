# Google Cloud Run Setup - Step by Step

## ✅ What's Already Done

- [x] Google Cloud SDK installed
- [x] Deployment script created (`deploy-to-gcp.sh`)
- [x] Docker configuration ready
- [x] App configured for Cloud Run

## 🎯 What You Need to Do (5 Steps)

### Step 1: Authenticate with Google Cloud (2 minutes)

Open your terminal and run:

```bash
gcloud auth login
```

This will:
- Open your browser
- Ask you to sign in with your Google account
- Grant permissions to gcloud CLI

**Note**: Use the Google account you want to use for billing.

---

### Step 2: Create a Google Cloud Project (2 minutes)

#### Option A: Create New Project

```bash
# Replace YOUR-PROJECT-ID with something unique (e.g., kinase-classifier-2025)
gcloud projects create YOUR-PROJECT-ID --name="Kinase Classifier"

# Set it as active
gcloud config set project YOUR-PROJECT-ID
```

#### Option B: Use Existing Project

```bash
# List your projects
gcloud projects list

# Set the one you want
gcloud config set project YOUR-EXISTING-PROJECT-ID
```

**Important**: You'll need to enable billing for this project:
1. Go to: https://console.cloud.google.com/billing
2. Link a billing account to your project
3. New users get $300 free credits!

---

### Step 3: Deploy! (10 minutes)

```bash
cd /Users/ajikumar/codefiles/Kinases-Clustering
./deploy-to-gcp.sh
```

The script will:
- ✅ Enable required APIs (Cloud Run, Cloud Build)
- ✅ Build your Docker image in the cloud (no local Docker needed!)
- ✅ Deploy to Cloud Run
- ✅ Give you a public URL

**What to expect:**
- First deployment: 5-10 minutes (building Docker image)
- Future deployments: 2-3 minutes (cached layers)

---

### Step 4: Test Your App (1 minute)

After deployment, you'll get a URL like:
```
https://kinase-classifier-xxxxx-uc.a.run.app
```

Visit it in your browser! 🎉

---

### Step 5: (Optional) Custom Domain (5 minutes)

If you have a domain (e.g., `yourdomain.com`):

```bash
# For subdomain: kinase.yourdomain.com
gcloud run domain-mappings create \
  --service kinase-classifier \
  --domain kinase.yourdomain.com \
  --region us-central1

# Follow the instructions to add DNS records
```

For subpath routing (`yourdomain.com/tools/kinase`), see `DEPLOYMENT.md` for Load Balancer setup.

---

## 🔧 Configuration Options

You can customize the deployment by setting environment variables before running the script:

```bash
# Use different region
export REGION=us-east1

# Increase memory/CPU
export MEMORY=16Gi
export CPU=8

# Keep instance warm (faster response, costs more)
export MIN_INSTANCES=1

# Deploy
./deploy-to-gcp.sh
```

---

## 💰 Cost Estimates

With default settings (scales to zero):

| Usage | Monthly Cost |
|-------|-------------|
| Idle (no traffic) | $0 |
| 100 predictions/day | $5-10 |
| 1000 predictions/day | $20-40 |
| 10,000 predictions/day | $80-150 |

**Free tier**: First 2 million requests/month are free!

---

## 📊 Monitoring & Management

### View Logs
```bash
gcloud run services logs read kinase-classifier --region us-central1
```

### View in Console
https://console.cloud.google.com/run

### Update Deployment
Just run the script again:
```bash
./deploy-to-gcp.sh
```

### Delete Service
```bash
gcloud run services delete kinase-classifier --region us-central1
```

---

## 🚨 Troubleshooting

### "APIs not enabled"
The script enables them automatically, but if you see errors:
```bash
gcloud services enable run.googleapis.com
gcloud services enable cloudbuild.googleapis.com
```

### "Billing not enabled"
Go to: https://console.cloud.google.com/billing
Link a billing account to your project.

### "Permission denied"
Make sure you're the owner of the project:
```bash
gcloud projects get-iam-policy YOUR-PROJECT-ID
```

### "Build failed"
Check the build logs:
```bash
gcloud builds list --limit 5
gcloud builds log BUILD_ID
```

Common issues:
- Large files in repo (check `.dockerignore`)
- Missing dependencies (check `Dockerfile`)

---

## 🎓 Learning Resources

- **Cloud Run Docs**: https://cloud.google.com/run/docs
- **Pricing Calculator**: https://cloud.google.com/products/calculator
- **Quotas & Limits**: https://cloud.google.com/run/quotas

---

## 🔄 For Future Tools

To add more tools to your website:

1. Create new directory: `tools/tool-name/`
2. Add Dockerfile and app code
3. Deploy:
```bash
cd tools/tool-name
gcloud run deploy tool-name --source . --region us-central1 ...
```

Each tool is independent, scales separately, and has its own URL!

---

## ✅ Quick Checklist

- [ ] Run `gcloud auth login`
- [ ] Create or select project
- [ ] Enable billing
- [ ] Run `./deploy-to-gcp.sh`
- [ ] Visit the URL and test
- [ ] (Optional) Set up custom domain

---

**Need help?** Just ask! I can guide you through any step.

