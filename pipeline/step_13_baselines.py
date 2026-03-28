#!/usr/bin/env python3
"""
Step 13: Baseline Comparisons

Runs baseline methods (k-NN, Random, MLP, Motifs-only LR) for comparison
with the main logistic regression model.

k-NN uses layers20_30_mean embeddings (cosine distance) as per manuscript.
Motifs-only LR uses 30 handcrafted regex features from domain sequences.

Usage:
    python pipeline/step_13_baselines.py --run-dir runs/2025-01-01_000000/

Outputs:
    - results/baselines/baselines_split40.csv
    - results/baselines/knn_split40.json
    - results/baselines/random_split40.json
    - results/baselines/mlp_split40.json
    - results/baselines/motifs_split40.json
"""

import argparse
import json
import re
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.dummy import DummyClassifier
from sklearn.metrics import accuracy_score, f1_score, log_loss

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.membership import load_split

# Project root
PROJECT_ROOT = Path(__file__).parent.parent

# Configuration
RANDOM_STATE = 42

# --- Kinase motif patterns (30 features) ---
# Each entry: (feature_name, regex_pattern)
MOTIF_PATTERNS = [
    # Catalytic loop / DFG motif
    ("DFG_present",       r"DFG"),
    ("DFG_flanking_A",    r".DFG"),          # residue before DFG
    ("DFG_flanking_B",    r"DFG."),          # residue after DFG
    ("DFG_in",            r"DFG[IL]"),        # DFG-in conformation indicator
    ("DFG_out",           r"DFG[^IL]"),       # DFG-out
    # HRD / catalytic aspartate
    ("HRD_present",       r"HRD"),
    ("HRD_K",             r"HRD[LIV]K"),     # HRD...K motif
    # APE motif (activation loop end)
    ("APE_present",       r"APE"),
    ("APE_L",             r"APEL"),
    # P-loop / glycine-rich loop
    ("GxGxxG",            r"G.G..G"),
    ("GXGXFG",            r"G.G.FG"),
    ("VAIK",              r"VA[IL]K"),        # conserved Lys in beta3
    ("VAIK_alt",          r"V[AI][IL]K"),
    # Activation loop
    ("DFGxxx",            r"DFG..."),         # short DFG context
    ("activation_ELK",    r"E[LIV]K"),
    # Gatekeeper residue context (T/M/L common)
    ("gatekeeper_T",      r"[LI]T[IV]"),
    ("gatekeeper_M",      r"[LI]M[IV]"),
    # C-helix glutamate (E..K salt bridge)
    ("ExxxK",             r"E...[LI]K"),
    ("ExK",               r"E.K"),
    # Hydrophobic spine residues context
    ("HxH_spine",         r"[HY][^ACDEFGHIKLMNPQRSTVWY][HY]"),
    # WxGxGxxG (some kinases)
    ("WxGxG",             r"W.G.G"),
    # RD pocket
    ("RD_kinase",         r"R[^P]D[LIVF]"),
    # Autophosphorylation site (TxY or YxxM)
    ("TxY_motif",         r"T.Y"),
    ("YxxM_motif",        r"Y..M"),
    # Src-family SH2 binding
    ("pY_motif",          r"Y[DE][DE]"),
    # AGC-specific (hydrophobic motif)
    ("Fxx[FY]",           r"F..[FY]"),
    # CK1 specific (pS/pT-x-x-S)
    ("CK1_motif",         r"S..S"),
    # CMGC-specific insert
    ("KxK_CMGC",          r"K.K"),
    # TK-specific NPXY
    ("NPXY",              r"NP.Y"),
    # Generic acidic cluster
    ("acidic_cluster",    r"[DE]{3,}"),
]


def extract_motif_features(sequences):
    """
    Extract 30 binary motif features from a list of amino-acid sequences.
    Returns (N, 30) float32 array.
    """
    n = len(sequences)
    X = np.zeros((n, len(MOTIF_PATTERNS)), dtype=np.float32)
    for i, seq in enumerate(sequences):
        seq_upper = seq.upper()
        for j, (_, pattern) in enumerate(MOTIF_PATTERNS):
            X[i, j] = 1.0 if re.search(pattern, seq_upper) else 0.0
    return X


def load_domain_sequences(run_dir):
    """Load domain FASTA sequences, return {uniprot_id: sequence}."""
    fasta_file = PROJECT_ROOT / "data" / "domains" / "domains_E001.fasta"
    if not fasta_file.exists():
        return {}

    id_to_seq = {}
    current_id = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    id_to_seq[current_id] = "".join(current_seq)
                # Header may be ">sp|UniProtID|..." or just ">UniProtID"
                header = line[1:].split()[0]
                parts = header.split("|")
                current_id = parts[1] if len(parts) >= 2 else parts[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        id_to_seq[current_id] = "".join(current_seq)
    return id_to_seq


def main():
    parser = argparse.ArgumentParser(description="Run baseline comparisons")
    parser.add_argument("--run-dir", type=str, required=True)
    args = parser.parse_args()

    run_dir     = Path(args.run_dir)
    results_dir = run_dir / "results" / "baselines"
    results_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Step 13: Baseline Comparisons")
    print("=" * 60)

    # Load embedding IDs
    ids_file = run_dir / "embeddings" / "esm2_t33_650M" / "ids.txt"
    with open(ids_file) as f:
        embedding_ids = [line.strip() for line in f if line.strip()]
    id_to_idx = {uid: i for i, uid in enumerate(embedding_ids)}

    # Load labels
    labels_file = PROJECT_ROOT / "data" / "processed" / "labels.csv"
    labels_df   = pd.read_csv(labels_file)
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))

    # Load embeddings
    emb_layer33 = np.load(
        run_dir / "embeddings" / "esm2_t33_650M" / "domain_E001_layer33_mean.npy"
    )
    emb_mid_file = run_dir / "embeddings" / "esm2_t33_650M" / "domain_E001_layers20_30_mean.npy"
    emb_mid = np.load(emb_mid_file) if (emb_mid_file.exists() or emb_mid_file.is_symlink()) else emb_layer33

    split_name = "split40"
    train_ids  = load_split(split_name, "train", run_dir)
    test_ids   = load_split(split_name, "test",  run_dir)
    print(f"\nUsing {split_name}: {len(train_ids)} train, {len(test_ids)} test")

    train_idx    = [id_to_idx[uid] for uid in train_ids if uid in id_to_idx]
    test_idx     = [id_to_idx[uid] for uid in test_ids  if uid in id_to_idx]
    train_labels = [id_to_label[uid] for uid in train_ids if uid in id_to_label]
    test_labels  = [id_to_label[uid] for uid in test_ids  if uid in id_to_label]

    X_train_33  = emb_layer33[train_idx]
    X_test_33   = emb_layer33[test_idx]
    X_train_mid = emb_mid[train_idx]
    X_test_mid  = emb_mid[test_idx]

    classes = sorted(set(train_labels))

    summary_rows = []

    def run_and_save(key, name, model, X_tr, X_te, top3=False):
        print(f"\n{name}...")
        model.fit(X_tr, train_labels)
        y_pred = model.predict(X_te)
        try:
            y_proba = model.predict_proba(X_te)
            ll = float(log_loss(test_labels, y_proba, labels=model.classes_))
        except AttributeError:
            y_proba = None
            ll = None

        acc = float(accuracy_score(test_labels, y_pred))
        f1  = float(f1_score(test_labels, y_pred, average='macro', zero_division=0))
        print(f"  Accuracy: {acc:.4f}  Macro-F1: {f1:.4f}")

        # Top-3 accuracy
        top3_acc = None
        if top3 and y_proba is not None:
            top3_preds = np.argsort(y_proba, axis=1)[:, -3:]
            cls_list   = list(model.classes_)
            correct    = sum(
                cls_list.index(tl) in list(tp)
                for tl, tp in zip(test_labels, top3_preds)
                if tl in cls_list
            )
            top3_acc = float(correct / len(test_labels))
            print(f"  Top-3 Acc: {top3_acc:.4f}")

        result = {"name": name, "accuracy": acc, "macro_f1": f1, "log_loss": ll}
        if top3_acc is not None:
            result["top3_accuracy"] = top3_acc

        result_file = results_dir / f"{key}_{split_name}.json"
        with open(result_file, 'w') as fh:
            json.dump(result, fh, indent=2)
        print(f"  ✓ Saved: {result_file.name}")

        summary_rows.append({
            "Method":    name,
            "Accuracy":  f"{acc:.4f}",
            "Macro-F1":  f"{f1:.4f}",
            "Log-loss":  f"{ll:.4f}" if ll is not None else "N/A",
            "Top-3 Acc": f"{top3_acc:.4f}" if top3_acc is not None else "—",
        })

    # k-NN (k=5, cosine) — uses layers20_30_mean as per manuscript
    run_and_save(
        "knn", "k-NN (k=5, cosine, layers20-30)",
        KNeighborsClassifier(n_neighbors=5, metric='cosine'),
        X_train_mid, X_test_mid, top3=True
    )

    # Random stratified
    run_and_save(
        "random", "Random (stratified)",
        DummyClassifier(strategy='stratified', random_state=RANDOM_STATE),
        X_train_33, X_test_33, top3=True
    )

    # MLP 256→64 (two hidden layers as per manuscript)
    run_and_save(
        "mlp", "MLP (256→64)",
        MLPClassifier(
            hidden_layer_sizes=(256, 64),
            max_iter=500,
            random_state=RANDOM_STATE
        ),
        X_train_33, X_test_33, top3=True
    )

    # Motifs-only LR
    print("\nMotifs-only LR (30 regex features)...")
    id_to_seq = load_domain_sequences(run_dir)
    if id_to_seq:
        train_seqs = [id_to_seq.get(uid, "") for uid in train_ids if uid in id_to_label]
        test_seqs  = [id_to_seq.get(uid, "") for uid in test_ids  if uid in id_to_label]
        X_train_motif = extract_motif_features(train_seqs)
        X_test_motif  = extract_motif_features(test_seqs)
        run_and_save(
            "motifs", "Motifs-only LR (30 features)",
            LogisticRegression(
                max_iter=1000,
                random_state=RANDOM_STATE,
                class_weight='balanced',
                solver='lbfgs',
                multi_class='multinomial'
            ),
            X_train_motif, X_test_motif
        )
    else:
        print("  ⚠ Domain FASTA not found — skipping motifs baseline")

    # Summary CSV
    summary_df   = pd.DataFrame(summary_rows)
    summary_file = results_dir / f"baselines_{split_name}.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"\n✓ Saved: {summary_file}")

    metadata = {
        "step": 13,
        "name": "Baseline Comparisons",
        "timestamp": datetime.now().isoformat(),
        "split": split_name,
        "n_train": len(train_idx),
        "n_test":  len(test_idx),
        "n_classes": len(classes)
    }
    with open(results_dir / "baselines_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    print("\n" + "=" * 60)
    print("Step 13 COMPLETE")
    print("=" * 60)
    print(f"\n{'Method':<35} {'Accuracy':<12} {'Macro-F1':<12}")
    print("-" * 60)
    for r in summary_rows:
        print(f"{r['Method']:<35} {r['Accuracy']:<12} {r['Macro-F1']:<12}")
    print("-" * 60)


if __name__ == "__main__":
    main()
