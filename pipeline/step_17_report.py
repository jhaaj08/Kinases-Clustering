#!/usr/bin/env python3
"""
Step 17: Generate Reviewer HTML Report

Produces a single self-contained REPORT.html in the run directory that embeds
all figures (base64), renders all tables, and displays all key manuscript
numbers.  No external dependencies — opens offline in any browser.

Usage:
    python pipeline/step_17_report.py --run-dir runs/run_current_data
"""

import argparse
import base64
import csv
import json
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_json(path: Path):
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def load_csv_rows(path: Path):
    """Return (headers, rows) from a CSV file."""
    if not path.exists():
        return [], []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        rows = [r for r in reader if any(cell.strip() for cell in r)]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def encode_image(path: Path) -> str:
    """Return base64-encoded PNG data-URI string, or empty string if missing."""
    if not path.exists():
        return ""
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return f"data:image/png;base64,{data}"


# ---------------------------------------------------------------------------
# HTML building blocks
# ---------------------------------------------------------------------------

CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
    font-size: 14px; line-height: 1.55; color: #222;
    background: #f5f5f5; padding: 24px;
}
.page { max-width: 960px; margin: 0 auto; background: #fff;
        padding: 40px 48px; border-radius: 6px;
        box-shadow: 0 2px 12px rgba(0,0,0,.08); }
h1 { font-size: 22px; font-weight: 700; color: #1a1a2e; margin-bottom: 4px; }
h2 { font-size: 16px; font-weight: 700; color: #1a1a2e;
     margin: 36px 0 12px; border-bottom: 2px solid #e0e0e0; padding-bottom: 6px; }
h3 { font-size: 14px; font-weight: 600; color: #333; margin: 24px 0 8px; }
.meta { font-size: 12px; color: #666; margin-bottom: 28px; }
.kpi-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px;
            margin-bottom: 8px; }
.kpi { background: #f0f4ff; border-left: 4px solid #3b5bdb;
       border-radius: 4px; padding: 12px 16px; }
.kpi .label { font-size: 11px; color: #555; text-transform: uppercase;
              letter-spacing: .5px; }
.kpi .value { font-size: 22px; font-weight: 700; color: #1a1a2e; }
.kpi .sub   { font-size: 11px; color: #777; margin-top: 2px; }
table { border-collapse: collapse; width: 100%; margin-bottom: 8px; font-size: 13px; }
th { background: #f0f4ff; color: #1a1a2e; font-weight: 600;
     text-align: left; padding: 8px 10px;
     border: 1px solid #d0d8f0; }
td { padding: 7px 10px; border: 1px solid #e0e0e0; }
tr:nth-child(even) td { background: #fafbff; }
tr.primary td { background: #e8f0fe; font-weight: 600; }
tr.baseline td { color: #555; }
.figure-block { margin: 28px 0; }
.figure-block img { width: 100%; border: 1px solid #e0e0e0;
                    border-radius: 4px; display: block; }
.figure-block .caption { font-size: 12px; color: #555; margin-top: 6px;
                          font-style: italic; }
.check { color: #2d9d2d; font-weight: 600; }
.warn  { color: #b45309; font-weight: 600; }
.section-note { font-size: 12px; color: #666; margin-top: 4px; }
@media print {
    body { background: white; padding: 0; }
    .page { box-shadow: none; padding: 20px; }
}
"""


def html_table(headers, rows, row_classes=None):
    th = "".join(f"<th>{h}</th>" for h in headers)
    body = ""
    for i, row in enumerate(rows):
        cls = (row_classes[i] if row_classes and i < len(row_classes) else "")
        tds = "".join(f"<td>{cell}</td>" for cell in row)
        body += f'<tr class="{cls}">{tds}</tr>\n'
    return f"<table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table>"


def kpi_card(label, value, sub=""):
    return (f'<div class="kpi">'
            f'<div class="label">{label}</div>'
            f'<div class="value">{value}</div>'
            f'<div class="sub">{sub}</div>'
            f'</div>')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_report(run_dir: Path) -> str:
    numbers = load_json(run_dir / "results" / "manuscript_numbers.json")
    fig_registry = load_json(run_dir / "figures" / "figure_registry.json")
    run_config = load_json(run_dir / "run_config.json")

    run_id = run_config.get("run_id", run_dir.name)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    numbers_ts = numbers.get("timestamp", "")[:19].replace("T", " ")

    # ---- Key metrics ----
    clust = numbers.get("clustering", {})
    cal = numbers.get("calibration", {})
    ret = numbers.get("retrieval", {})
    ds = numbers.get("dataset", {})

    best_ari = clust.get("best_ARI", 0)
    improvement = clust.get("improvement_percent", 0)
    cal_acc = cal.get("layer33_mean", {}).get("calibrated_accuracy", 0)
    cal_ece = cal.get("layer33_mean", {}).get("calibrated_ece", 0)
    p1 = ret.get("P@1", 0)
    mrr = ret.get("MRR", 0)

    kpis = (
        kpi_card("Best Clustering ARI", f"{best_ari:.3f}",
                 f"+{improvement:.1f}% vs. Layer 33 baseline") +
        kpi_card("Calibrated Accuracy", f"{cal_acc*100:.1f}%",
                 "Layer 33 LR, Platt scaling, split40") +
        kpi_card("Retrieval P@1", f"{p1:.3f}",
                 f"MRR = {mrr:.3f}") +
        kpi_card("ECE (calibrated)", f"{cal_ece:.3f}",
                 "Layer 33 after Platt scaling") +
        kpi_card("Supervised N", f"{ds.get('supervised_eligible_n', '?')}",
                 f"{ds.get('supervised_eligible_classes', '?')} classes, 40% split") +
        kpi_card("Clustering N", f"{clust.get('n_sequences', '?')}",
                 f"k={clust.get('k', '?')} families")
    )

    # ---- Dataset table ----
    split_info = numbers.get("splits", {})
    ds_headers = ["Stage", "N sequences", "N classes"]
    ds_rows = [
        ["Domains (E &lt; 0.01, excl. Other)",
         ds.get("domain_E001_n", ""),
         ds.get("domain_E001_classes", "")],
        ["Supervised-eligible",
         ds.get("supervised_eligible_n", ""),
         ds.get("supervised_eligible_classes", "")],
    ]
    for sname in ["split40", "split50", "split70"]:
        si = split_info.get(sname, {})
        ds_rows.append([f"{sname} train", si.get("n_train", ""), ds.get("supervised_eligible_classes", "")])
        ds_rows.append([f"{sname} test",  si.get("n_test", ""),  ds.get("supervised_eligible_classes", "")])

    # ---- Table 1: Supervised Classification ----
    t1s_headers, t1s_rows_raw = load_csv_rows(run_dir / "tables" / "Table1_supervised.csv")
    t1s_classes = []
    for row in t1s_rows_raw:
        type_col = row[4] if len(row) > 4 else ""
        if "calibrated)" in (row[0] if row else "") and "uncalibrated" not in row[0]:
            if "Layer 33" in row[0]:
                t1s_classes.append("primary")
            else:
                t1s_classes.append("")
        elif type_col == "Baseline":
            t1s_classes.append("baseline")
        else:
            t1s_classes.append("")
    # Format accuracy as percentage for display
    t1s_rows_display = []
    for row in t1s_rows_raw:
        display = list(row)
        try:
            display[1] = f"{float(row[1])*100:.1f}%"
        except (ValueError, IndexError):
            pass
        t1s_rows_display.append(display)

    # ---- Table S1: Clustering ----
    ts1_headers, ts1_rows = load_csv_rows(run_dir / "tables" / "TableS1.csv")

    # ---- Table S2: Baselines ----
    ts2_headers, ts2_rows = load_csv_rows(run_dir / "tables" / "TableS2.csv")

    # ---- Figures ----
    figures_html = ""
    fig_order = [
        "Figure1_clustering_ari.png",
        "Figure2_confusion_matrix.png",
        "Figure3_homology_classification.png",
        "Figure4_pooling_comparison.png",
        "Figure5_calibration.png",
        "Figure6_retrieval_pr.png",
    ]
    fig_meta = fig_registry.get("figures", {})
    for i, fname in enumerate(fig_order, 1):
        src = encode_image(run_dir / "figures" / fname)
        meta = fig_meta.get(fname, {})
        title = meta.get("title", fname)
        cited = meta.get("cited_in", "")
        if src:
            figures_html += (
                f'<div class="figure-block">'
                f'<img src="{src}" alt="Figure {i}: {title}">'
                f'<div class="caption"><strong>Figure {i}.</strong> {title}'
                + (f" ({cited})" if cited else "") +
                f'</div></div>\n'
            )
        else:
            figures_html += f'<p class="warn">Figure {i} ({fname}) not found.</p>\n'

    # ---- Verification checklist ----
    checks = []

    # Dataset hierarchy
    dom_n = ds.get("domain_E001_n", 0)
    sup_n = ds.get("supervised_eligible_n", 0)
    ok = dom_n >= sup_n > 0
    checks.append((ok, f"Domain N ({dom_n}) &ge; Supervised N ({sup_n})"))

    # Split integrity (split40)
    s40 = split_info.get("split40", {})
    total = s40.get("n_train", 0) + s40.get("n_test", 0)
    ok = total == sup_n
    checks.append((ok, f"split40 train+test ({total}) = supervised_eligible ({sup_n})"))

    # Calibrated accuracy sanity
    ok = 0.7 < cal_acc < 1.0
    checks.append((ok, f"Calibrated accuracy ({cal_acc*100:.1f}%) in expected range (70–100%)"))

    # ARI improvement
    ok = improvement > 100
    checks.append((ok, f"ARI improvement ({improvement:.1f}%) &gt; 100%"))

    # Figures present
    n_figs = sum(1 for f in fig_order if (run_dir / "figures" / f).exists())
    ok = n_figs == 6
    checks.append((ok, f"All 6 figures present ({n_figs}/6)"))

    check_html = "<ul style='list-style:none;padding:0'>"
    for ok, msg in checks:
        icon = '<span class="check">&#10003;</span>' if ok else '<span class="warn">&#10007;</span>'
        check_html += f"<li style='margin:4px 0'>{icon} {msg}</li>"
    check_html += "</ul>"

    # ---- Assemble HTML ----
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Kinase Classification — Reviewer Report ({run_id})</title>
<style>{CSS}</style>
</head>
<body>
<div class="page">

<h1>Kinase Functional Classification with ESM-2 Layer Selection</h1>
<div class="meta">
  Run ID: <strong>{run_id}</strong> &nbsp;|&nbsp;
  Results timestamp: {numbers_ts} &nbsp;|&nbsp;
  Report generated: {generated_at}
</div>

<h2>Key Results</h2>
<div class="kpi-grid">{kpis}</div>
<p class="section-note">Primary claim: mid-layer averaging (layers 20–30) improves unsupervised clustering ARI
by +{improvement:.1f}% over the final layer, while the final layer excels at supervised classification.</p>

<h2>Dataset Summary</h2>
{html_table(ds_headers, ds_rows)}

<h2>Table 1. Supervised Classification Performance (40% Identity Split)</h2>
<p class="section-note">Highlighted rows = ESM-2 LR calibrated (primary method). Grey rows = baselines.</p>
{html_table(t1s_headers, t1s_rows_display, t1s_classes) if t1s_rows_display else "<p>Table1_supervised.csv not found — re-run step 15.</p>"}

<h2>Table S1. Clustering Performance Across ESM-2 Layer Configurations</h2>
{html_table(ts1_headers, ts1_rows) if ts1_rows else "<p>TableS1.csv not found.</p>"}

<h2>Table S2. Baseline Comparisons (40% Identity Split)</h2>
{html_table(ts2_headers, ts2_rows) if ts2_rows else "<p>TableS2.csv not found.</p>"}

<h2>Figures</h2>
{figures_html}

<h2>Verification Checklist</h2>
{check_html}
<p class="section-note" style="margin-top:8px">
  Full integrity check: <code>make verify RUN_ID={run_id}</code>
</p>

</div>
</body>
</html>"""
    return html


def main():
    parser = argparse.ArgumentParser(description="Generate reviewer HTML report")
    parser.add_argument("--run-dir", required=True, help="Run directory path")
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise SystemExit(f"Run directory not found: {run_dir}")

    print("=" * 60)
    print("Step 17: Generate Reviewer HTML Report")
    print("=" * 60)

    html = build_report(run_dir)

    out = run_dir / "REPORT.html"
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)

    size_kb = out.stat().st_size // 1024
    print(f"\n✓ Saved: {out}  ({size_kb} KB, self-contained)")
    print(f"\nOpen in browser:")
    print(f"  open {out}   # macOS")
    print(f"  xdg-open {out}   # Linux")


if __name__ == "__main__":
    main()
