# ============================================================================
# Kinase Classification Pipeline - Makefile
# ============================================================================
# 
# This Makefile provides dependency-tracked builds for reproducible experiments.
# Each run creates an isolated directory under runs/<RUN_ID>/.
#
# Usage:
#   make all                    # Fresh run with auto-generated timestamp
#   make all RUN_ID=exp_v2      # Named run
#   make all RUN_ID=exp_v2 FORCE=1  # Force overwrite existing run
#   make verify                 # Verify current run
#   make zip                    # Create Zenodo package
#   make clean                  # Remove current run
#
# ============================================================================

# Configuration
RUN_ID ?= $(shell date +%Y-%m-%d_%H%M%S)
RUN_DIR := runs/$(RUN_ID)
PYTHON := python3
CURRENT_LINK := runs/current

# Reproducibility settings (ensures deterministic results)
export PYTHONHASHSEED := 0
export OMP_NUM_THREADS := 1
export MKL_NUM_THREADS := 1
export OPENBLAS_NUM_THREADS := 1

# Source data (project-level, not per-run)
SRC_DOMAINS_FASTA := data/domains/domains_E001.fasta
SRC_DOMAINS_COORDS := data/domains/domain_coords_E001.tsv
SRC_LABELS := data/processed/labels.csv
SRC_HMM_PROFILES := data/hmm_profiles

# ============================================================================
# PHONY TARGETS
# ============================================================================
.PHONY: all clean verify zip help init check-fresh \
        manifests splits embeddings clustering supervised \
        calibration baselines retrieval tables figures

# Default target
all: $(RUN_DIR)/MANIFEST.txt
	@echo ""
	@echo "============================================================"
	@echo "BUILD COMPLETE: $(RUN_DIR)"
	@echo "============================================================"
	@echo "Run 'make verify' to validate the package"
	@echo "Run 'make zip' to create Zenodo package"

# Help target
help:
	@echo "Kinase Classification Pipeline"
	@echo ""
	@echo "Usage:"
	@echo "  make all                    Fresh run (auto-generated ID)"
	@echo "  make all RUN_ID=name        Named run"
	@echo "  make all FORCE=1            Force overwrite"
	@echo "  make verify                 Verify current run"
	@echo "  make zip                    Create Zenodo ZIP"
	@echo "  make clean                  Remove current run"
	@echo "  make list                   List all runs"
	@echo ""
	@echo "Individual steps:"
	@echo "  make manifests              Step 6: Create dataset manifests"
	@echo "  make embeddings             Step 8: Link embeddings"
	@echo "  make splits                 Step 10: Create homology-aware splits"
	@echo "  make clustering             Step 9: Run k-means clustering"
	@echo "  make supervised             Step 11: Train classifiers"
	@echo "  make calibration            Step 12: Calibration"
	@echo "  make baselines              Step 13: Baseline comparisons"
	@echo "  make retrieval              Step 14: Retrieval experiment"
	@echo "  make tables                 Step 15: Generate manuscript tables"
	@echo "  make figures                Step 16: Generate manuscript figures"
	@echo ""
	@echo "Current RUN_ID: $(RUN_ID)"
	@echo "Current RUN_DIR: $(RUN_DIR)"

# ============================================================================
# INITIALIZATION
# ============================================================================

# Check for fresh run (abort if exists without FORCE)
check-fresh:
ifndef FORCE
	@if [ -d "$(RUN_DIR)" ]; then \
		echo "ERROR: Run directory '$(RUN_DIR)' already exists."; \
		echo "Options:"; \
		echo "  1. Use 'make FORCE=1' to overwrite"; \
		echo "  2. Specify a new RUN_ID"; \
		echo "  3. Delete the existing directory"; \
		exit 1; \
	fi
endif

# Initialize run directory
init: check-fresh
	@mkdir -p $(RUN_DIR)/data/manifests
	@mkdir -p $(RUN_DIR)/data/splits
	@mkdir -p $(RUN_DIR)/data/domains
	@mkdir -p $(RUN_DIR)/embeddings/esm2_t33_650M
	@mkdir -p $(RUN_DIR)/results/clustering
	@mkdir -p $(RUN_DIR)/results/supervised
	@mkdir -p $(RUN_DIR)/results/calibration
	@mkdir -p $(RUN_DIR)/results/baselines
	@mkdir -p $(RUN_DIR)/results/retrieval
	@mkdir -p $(RUN_DIR)/results/layer_comparison
	@mkdir -p $(RUN_DIR)/tables
	@mkdir -p $(RUN_DIR)/figures
	@echo '{"run_id": "$(RUN_ID)", "created_at": "$(shell date -Iseconds)"}' > $(RUN_DIR)/run_config.json
	@rm -f $(CURRENT_LINK)
	@ln -s $(RUN_ID) $(CURRENT_LINK)
	@echo "[init] Created run directory: $(RUN_DIR)"

# ============================================================================
# STEP 6: MANIFESTS (Dataset Membership)
# ============================================================================
$(RUN_DIR)/data/manifests/supervised_eligible.txt: $(SRC_DOMAINS_COORDS) $(SRC_LABELS) | init
	@echo ""
	@echo "[Step 6] Creating manifests..."
	$(PYTHON) pipeline/step_06_manifests.py --run-dir $(RUN_DIR)
	@echo "[Step 6] ✓ Manifests created"

manifests: $(RUN_DIR)/data/manifests/supervised_eligible.txt

# ============================================================================
# STEP 8: EMBEDDINGS
# ============================================================================
$(RUN_DIR)/embeddings/esm2_t33_650M/ids.txt: $(SRC_DOMAINS_FASTA) | init
	@echo ""
	@echo "[Step 8] Copying/linking embeddings..."
	$(PYTHON) pipeline/step_08_embeddings.py --run-dir $(RUN_DIR)
	@echo "[Step 8] ✓ Embeddings ready"

embeddings: $(RUN_DIR)/embeddings/esm2_t33_650M/ids.txt

# ============================================================================
# STEP 10: SPLITS (Homology-aware)
# ============================================================================
$(RUN_DIR)/data/splits/split40_train.txt: $(RUN_DIR)/data/manifests/supervised_eligible.txt
	@echo ""
	@echo "[Step 10] Creating homology-aware splits..."
	$(PYTHON) pipeline/step_10_splits.py --run-dir $(RUN_DIR)
	@echo "[Step 10] ✓ Splits created"

splits: $(RUN_DIR)/data/splits/split40_train.txt

# ============================================================================
# STEP 9: CLUSTERING
# ============================================================================
$(RUN_DIR)/results/clustering/clustering_registry.json: $(RUN_DIR)/embeddings/esm2_t33_650M/ids.txt $(RUN_DIR)/data/manifests/supervised_eligible.txt
	@echo ""
	@echo "[Step 9] Running clustering experiments..."
	$(PYTHON) pipeline/step_09_clustering.py --run-dir $(RUN_DIR)
	@echo "[Step 9] ✓ Clustering complete"

clustering: $(RUN_DIR)/results/clustering/clustering_registry.json

# ============================================================================
# STEP 11: SUPERVISED LEARNING
# ============================================================================
$(RUN_DIR)/results/supervised/supervised_registry.json: $(RUN_DIR)/data/splits/split40_train.txt $(RUN_DIR)/embeddings/esm2_t33_650M/ids.txt
	@echo ""
	@echo "[Step 11] Running supervised learning..."
	$(PYTHON) pipeline/step_11_supervised.py --run-dir $(RUN_DIR)
	@echo "[Step 11] ✓ Supervised learning complete"

supervised: $(RUN_DIR)/results/supervised/supervised_registry.json

# ============================================================================
# STEP 12: CALIBRATION
# ============================================================================
$(RUN_DIR)/results/calibration/split40_calibration.json: $(RUN_DIR)/results/supervised/supervised_registry.json
	@echo ""
	@echo "[Step 12] Running calibration..."
	$(PYTHON) pipeline/step_12_calibration.py --run-dir $(RUN_DIR)
	@echo "[Step 12] ✓ Calibration complete"

calibration: $(RUN_DIR)/results/calibration/split40_calibration.json

# ============================================================================
# STEP 13: BASELINES
# ============================================================================
$(RUN_DIR)/results/baselines/baselines_split40.csv: $(RUN_DIR)/data/splits/split40_train.txt $(RUN_DIR)/embeddings/esm2_t33_650M/ids.txt
	@echo ""
	@echo "[Step 13] Running baselines..."
	$(PYTHON) pipeline/step_13_baselines.py --run-dir $(RUN_DIR)
	@echo "[Step 13] ✓ Baselines complete"

baselines: $(RUN_DIR)/results/baselines/baselines_split40.csv

# ============================================================================
# STEP 14: RETRIEVAL
# ============================================================================
$(RUN_DIR)/results/retrieval/split40_retrieval.json: $(RUN_DIR)/data/splits/split40_train.txt $(RUN_DIR)/embeddings/esm2_t33_650M/ids.txt
	@echo ""
	@echo "[Step 14] Running retrieval experiment..."
	$(PYTHON) pipeline/step_14_retrieval.py --run-dir $(RUN_DIR)
	@echo "[Step 14] ✓ Retrieval complete"

retrieval: $(RUN_DIR)/results/retrieval/split40_retrieval.json

# ============================================================================
# STEP 15: BUILD MANUSCRIPT NUMBERS
# ============================================================================
$(RUN_DIR)/results/manuscript_numbers.json: \
		$(RUN_DIR)/results/clustering/clustering_registry.json \
		$(RUN_DIR)/results/supervised/supervised_registry.json \
		$(RUN_DIR)/results/calibration/split40_calibration.json \
		$(RUN_DIR)/results/baselines/baselines_split40.csv \
		$(RUN_DIR)/results/retrieval/split40_retrieval.json
	@echo ""
	@echo "[Step 15] Building manuscript numbers and tables..."
	$(PYTHON) pipeline/step_15_build_numbers.py --run-dir $(RUN_DIR)
	@echo "[Step 15] ✓ Manuscript numbers built"

tables: $(RUN_DIR)/results/manuscript_numbers.json

# ============================================================================
# STEP 16: FIGURES
# ============================================================================
$(RUN_DIR)/figures/figure_registry.json: $(RUN_DIR)/results/manuscript_numbers.json
	@echo ""
	@echo "[Step 16] Generating manuscript figures..."
	$(PYTHON) pipeline/step_16_figures.py --run-dir $(RUN_DIR)
	@echo "[Step 16] ✓ Figures generated"

figures: $(RUN_DIR)/figures/figure_registry.json

# ============================================================================
# MANIFEST (SHA256 hashes)
# ============================================================================
$(RUN_DIR)/MANIFEST.txt: $(RUN_DIR)/figures/figure_registry.json
	@echo ""
	@echo "[Final] Generating MANIFEST.txt with SHA256 hashes..."
	$(PYTHON) pipeline/generate_manifest.py --run-dir $(RUN_DIR)
	@echo "[Final] ✓ MANIFEST.txt generated"

# ============================================================================
# VERIFICATION (uses runs/current by default)
# ============================================================================
verify:
	@if [ ! -e "runs/current" ] && [ ! -d "$(RUN_DIR)" ]; then \
		echo "Error: No run to verify."; \
		echo "  Either run 'make all' first, or specify: make verify RUN_ID=xxx"; \
		exit 1; \
	fi
	@VERIFY_PATH="runs/current"; \
	if [ -n "$(filter RUN_ID=%,$(MAKEFLAGS))" ] || [ "$(origin RUN_ID)" = "command line" ]; then \
		VERIFY_PATH="$(RUN_DIR)"; \
	fi; \
	echo ""; \
	echo "============================================================"; \
	echo "VERIFYING PACKAGE: $$VERIFY_PATH"; \
	echo "============================================================"; \
	$(PYTHON) scripts/verify_package.py $$VERIFY_PATH

# ============================================================================
# PACKAGE FOR ZENODO (uses runs/current by default)
# ============================================================================
zip:
	@if [ ! -e "runs/current" ]; then \
		echo "Error: No run to package. Run 'make all' first."; \
		exit 1; \
	fi
	@ZIP_DIR="runs/current"; \
	ZIP_NAME="kinase_data_$$(basename $$(readlink runs/current)).zip"; \
	echo ""; \
	echo "Creating Zenodo package..."; \
	cd $$ZIP_DIR && zip -r ../../$$ZIP_NAME . -x "*.DS_Store" -x "__MACOSX/*"; \
	echo "✓ Created: $$ZIP_NAME"

# ============================================================================
# UTILITIES
# ============================================================================

# List all runs
list:
	@echo "Available runs:"
	@ls -lt runs/ 2>/dev/null | grep -v "^total" | grep -v "current" | head -10 || echo "  (none)"
	@echo ""
	@if [ -L "$(CURRENT_LINK)" ]; then \
		echo "Current: $$(readlink $(CURRENT_LINK))"; \
	fi

# Clean current run
clean:
ifdef RUN_ID
	@echo "Removing: $(RUN_DIR)"
	rm -rf $(RUN_DIR)
	@if [ -L "$(CURRENT_LINK)" ] && [ "$$(readlink $(CURRENT_LINK))" = "$(RUN_ID)" ]; then \
		rm -f $(CURRENT_LINK); \
		echo "Removed symlink: $(CURRENT_LINK)"; \
	fi
else
	@echo "Specify RUN_ID to clean, e.g.: make clean RUN_ID=2025-01-01_000000"
endif

# Clean all runs (dangerous!)
clean-all:
	@echo "WARNING: This will delete ALL runs!"
	@read -p "Are you sure? [y/N] " confirm && [ "$$confirm" = "y" ] && rm -rf runs/* || echo "Aborted."

# Create .gitkeep
$(shell mkdir -p runs && touch runs/.gitkeep)

