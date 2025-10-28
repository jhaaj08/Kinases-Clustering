# 🚀 START HERE: Deploy Your Kinase Classifier

## What's Ready

✅ **Everything is set up!** Your app is ready to deploy to Google Cloud Run.

✅ **Google Cloud SDK** is installed on your Mac

✅ **Deployment script** is ready (`deploy-to-gcp.sh`)

---

## Next: 3 Simple Steps (10 minutes total)

### Step 1: Login to Google Cloud (2 min)

Open Terminal and run:

```bash
gcloud auth login
```

- Browser will open
- Sign in with your Google account
- Click "Allow"

---

### Step 2: Create a Project (2 min)

```bash
# Create a new project (choose a unique ID)
gcloud projects create kinase-classifier-2025 --name="Kinase Classifier"

# Set it as active
gcloud config set project kinase-classifier-2025

# Enable billing (required)
# Go to: https://console.cloud.google.com/billing
# Link your project to a billing account
# (New users get $300 free credits!)
```

---

### Step 3: Deploy! (5-10 min)

```bash
cd /Users/ajikumar/codefiles/Kinases-Clustering
./deploy-to-gcp.sh
```

**That's it!** The script will:
- Enable required APIs
- Build your Docker image in the cloud
- Deploy to Cloud Run
- Give you a public URL

---

## Your URL

After deployment, you'll get something like:

```
https://kinase-classifier-xxxxx-uc.a.run.app
```

**Share this URL** - your app is live! 🎉

---

## Cost

- **Idle**: $0/month (scales to zero)
- **Light use**: $5-10/month
- **Free tier**: First 2M requests/month free

---

## Need Help?

See `SETUP_GCP.md` for detailed instructions and troubleshooting.

---

## Quick Commands

```bash
# View logs
gcloud run services logs read kinase-classifier --region us-central1

# Update app (after code changes)
./deploy-to-gcp.sh

# Delete service
gcloud run services delete kinase-classifier --region us-central1
```

---

**Ready? Run these 3 commands:**

```bash
# 1. Login
gcloud auth login

# 2. Create project
gcloud projects create kinase-classifier-2025 --name="Kinase Classifier"
gcloud config set project kinase-classifier-2025

# 3. Deploy
cd /Users/ajikumar/codefiles/Kinases-Clustering
./deploy-to-gcp.sh
```

🎯 **That's all you need to do!**

