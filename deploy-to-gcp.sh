#!/bin/bash
# Deploy Kinase Classifier to Google Cloud Run
# Usage: ./deploy-to-gcp.sh

set -e

echo "=========================================="
echo "Kinase Classifier - Google Cloud Run Deploy"
echo "=========================================="
echo ""

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Check if gcloud is installed
if ! command -v gcloud &> /dev/null; then
    echo -e "${RED}Error: gcloud CLI not found${NC}"
    echo "Please run: brew install google-cloud-sdk"
    exit 1
fi

# Check if authenticated
if ! gcloud auth list --filter=status:ACTIVE --format="value(account)" &> /dev/null; then
    echo -e "${YELLOW}Not authenticated. Running gcloud auth login...${NC}"
    gcloud auth login
fi

# Get current project
CURRENT_PROJECT=$(gcloud config get-value project 2>/dev/null)

if [ -z "$CURRENT_PROJECT" ]; then
    echo -e "${YELLOW}No project selected.${NC}"
    echo "Please select or create a project:"
    echo ""
    echo "Option 1: List existing projects"
    echo "  gcloud projects list"
    echo ""
    echo "Option 2: Create new project"
    echo "  gcloud projects create YOUR-PROJECT-ID --name='Kinase Classifier'"
    echo ""
    echo "Then set the project:"
    echo "  gcloud config set project YOUR-PROJECT-ID"
    echo ""
    exit 1
fi

echo -e "${GREEN}✓ Using project: $CURRENT_PROJECT${NC}"
echo ""

# Set variables
SERVICE_NAME="kinase-classifier"
REGION="${REGION:-us-central1}"
MEMORY="${MEMORY:-8Gi}"
CPU="${CPU:-4}"
TIMEOUT="${TIMEOUT:-300}"
MAX_INSTANCES="${MAX_INSTANCES:-10}"
MIN_INSTANCES="${MIN_INSTANCES:-0}"

echo "Deployment Configuration:"
echo "  Service Name: $SERVICE_NAME"
echo "  Region: $REGION"
echo "  Memory: $MEMORY"
echo "  CPU: $CPU"
echo "  Timeout: ${TIMEOUT}s"
echo "  Min Instances: $MIN_INSTANCES"
echo "  Max Instances: $MAX_INSTANCES"
echo ""

# Enable required APIs
echo -e "${YELLOW}Enabling required APIs...${NC}"
gcloud services enable run.googleapis.com --quiet
gcloud services enable cloudbuild.googleapis.com --quiet
echo -e "${GREEN}✓ APIs enabled${NC}"
echo ""

# Deploy using Cloud Build (no local Docker needed!)
echo -e "${YELLOW}Deploying to Cloud Run...${NC}"
echo "This will build the Docker image in the cloud and deploy it."
echo "This may take 5-10 minutes..."
echo ""

gcloud run deploy $SERVICE_NAME \
  --source . \
  --platform managed \
  --region $REGION \
  --allow-unauthenticated \
  --port 7860 \
  --memory $MEMORY \
  --cpu $CPU \
  --timeout $TIMEOUT \
  --min-instances $MIN_INSTANCES \
  --max-instances $MAX_INSTANCES \
  --set-env-vars "GRADIO_SERVER_NAME=0.0.0.0,GRADIO_SERVER_PORT=7860"

echo ""
echo -e "${GREEN}=========================================="
echo "✓ Deployment Complete!"
echo "==========================================${NC}"
echo ""

# Get the service URL
SERVICE_URL=$(gcloud run services describe $SERVICE_NAME --region $REGION --format 'value(status.url)')

echo -e "${GREEN}Your app is live at:${NC}"
echo -e "${GREEN}$SERVICE_URL${NC}"
echo ""
echo "Test it:"
echo "  curl $SERVICE_URL"
echo ""
echo "View logs:"
echo "  gcloud run services logs read $SERVICE_NAME --region $REGION"
echo ""
echo "Update deployment:"
echo "  ./deploy-to-gcp.sh"
echo ""
echo -e "${GREEN}Next steps:${NC}"
echo "1. Visit the URL above to test your app"
echo "2. (Optional) Set up custom domain:"
echo "   gcloud run domain-mappings create --service $SERVICE_NAME --domain kinase.yourdomain.com --region $REGION"
echo ""

