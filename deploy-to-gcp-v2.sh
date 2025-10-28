#!/bin/bash
# Deploy Kinase Classifier to Google Cloud Run (Version 2 - with better error handling)
# Usage: ./deploy-to-gcp-v2.sh

set -e

echo "=========================================="
echo "Kinase Classifier - Cloud Run Deploy v2"
echo "=========================================="
echo ""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check gcloud
if ! command -v gcloud &> /dev/null; then
    echo -e "${RED}Error: gcloud CLI not found${NC}"
    exit 1
fi

# Get project
PROJECT=$(gcloud config get-value project 2>/dev/null)
if [ -z "$PROJECT" ]; then
    echo -e "${RED}No project selected${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Using project: $PROJECT${NC}"

# Variables
SERVICE_NAME="kinase-classifier"
REGION="${REGION:-us-central1}"
IMAGE_NAME="gcr.io/$PROJECT/$SERVICE_NAME"

echo ""
echo "Configuration:"
echo "  Service: $SERVICE_NAME"
echo "  Region: $REGION"
echo "  Image: $IMAGE_NAME"
echo ""

# Enable APIs
echo -e "${YELLOW}Enabling APIs...${NC}"
gcloud services enable run.googleapis.com --quiet
gcloud services enable cloudbuild.googleapis.com --quiet
gcloud services enable containerregistry.googleapis.com --quiet
echo -e "${GREEN}✓ APIs enabled${NC}"
echo ""

# Build using Cloud Build with explicit timeout
echo -e "${YELLOW}Building Docker image (this may take 10-15 minutes)...${NC}"
echo "Building in the cloud with Cloud Build..."
echo ""

gcloud builds submit \
  --tag "$IMAGE_NAME" \
  --timeout=20m \
  --machine-type=e2-highcpu-8 \
  .

if [ $? -ne 0 ]; then
    echo -e "${RED}Build failed!${NC}"
    echo ""
    echo "Debugging steps:"
    echo "1. Check build logs:"
    echo "   gcloud builds list --limit=1"
    echo "   gcloud builds log BUILD_ID"
    echo ""
    echo "2. Check if files are too large:"
    echo "   du -sh ."
    echo ""
    echo "3. Try viewing logs in console:"
    echo "   https://console.cloud.google.com/cloud-build/builds?project=$PROJECT"
    exit 1
fi

echo -e "${GREEN}✓ Build complete${NC}"
echo ""

# Deploy to Cloud Run
echo -e "${YELLOW}Deploying to Cloud Run...${NC}"

gcloud run deploy $SERVICE_NAME \
  --image "$IMAGE_NAME" \
  --platform managed \
  --region $REGION \
  --allow-unauthenticated \
  --port 7860 \
  --memory 8Gi \
  --cpu 4 \
  --timeout 300 \
  --min-instances 0 \
  --max-instances 10 \
  --set-env-vars "GRADIO_SERVER_NAME=0.0.0.0,GRADIO_SERVER_PORT=7860"

echo ""
echo -e "${GREEN}=========================================="
echo "✓ Deployment Complete!"
echo "==========================================${NC}"
echo ""

# Get URL
SERVICE_URL=$(gcloud run services describe $SERVICE_NAME --region $REGION --format 'value(status.url)')

echo -e "${GREEN}Your app is live at:${NC}"
echo -e "${GREEN}$SERVICE_URL${NC}"
echo ""
echo "Test it: curl $SERVICE_URL"
echo ""

