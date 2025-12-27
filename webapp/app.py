#!/usr/bin/env python3
"""
Kinase Classification Web Application
======================================

Gradio-based web interface for kinase sequence classification.

Features:
- Sequence input (FASTA or raw)
- Automatic domain extraction (HMMER)
- ESM-2 embedding generation (mid-layer averaging)
- Top-3 family predictions with calibrated probabilities
- Motif analysis and visualization
- Nearest training exemplars
- Confidence-based flagging
- JSON export

Usage:
    python webapp/app.py

Then open: http://localhost:7860
"""

import sys
import os
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import gradio as gr
import numpy as np
import pandas as pd
import json
from predictor import KinasePredictor
from motif_highlighter import MotifHighlighter, highlight_sequence_html
import plotly.graph_objects as go


# Initialize predictor (loads models on startup)
print("Loading models...")
predictor = KinasePredictor(
    model_path='supervised_results/logistic_regression_model.joblib',
    embeddings_path='kinases_domains_e0.01_layers_mid/esm2_embeddings.npy',
    index_path='kinases_domains_e0.01_layers_mid/esm2_index.csv',
    labels_path='kinases_domains_e0.01.csv',
    hmm_path='PF00069.hmm'
)
print("✅ Models loaded!")

motif_highlighter = MotifHighlighter()


def parse_fasta(text):
    """Parse FASTA or raw sequence."""
    text = text.strip()
    
    if text.startswith('>'):
        # FASTA format
        lines = text.split('\n')
        header = lines[0][1:].strip()
        sequence = ''.join(lines[1:]).replace(' ', '').upper()
        return header, sequence
    else:
        # Raw sequence
        sequence = text.replace(' ', '').replace('\n', '').upper()
        return "Input_Sequence", sequence


def predict_kinase(sequence_input):
    """Main prediction function."""
    try:
        # Parse input
        header, sequence = parse_fasta(sequence_input)
        
        if len(sequence) < 50:
            return ("❌ Error: Sequence too short (minimum 50 amino acids)", 
                    "", "", "", "", "")
        
        # Validate sequence
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        if not all(aa in valid_aa for aa in sequence):
            return ("❌ Error: Invalid amino acids in sequence", 
                    "", "", "", "", "")
        
        # Run prediction
        result = predictor.predict(sequence, header)
        
        if result['status'] == 'error':
            return (f"❌ Error: {result['message']}", "", "", "", "", "")
        
        # Format results
        summary = format_summary(result)
        predictions = format_predictions(result)
        motifs = format_motifs(result)
        exemplars = format_exemplars(result)
        sequence_viz = highlight_sequence_html(
            result['domain_sequence'], 
            result['motifs']
        )
        json_export = json.dumps(result, indent=2)
        
        return summary, predictions, motifs, exemplars, sequence_viz, json_export
        
    except Exception as e:
        return (f"❌ Error: {str(e)}", "", "", "", "", "")


def format_summary(result):
    """Format prediction summary."""
    domain_info = f"**Domain**: {result['domain_start']}-{result['domain_end']} ({result['domain_length']} aa)"
    
    top_pred = result['top_predictions'][0]
    confidence = "🟢 High" if result['confidence_flag'] == 'high' else "🟡 Needs Review"
    
    summary = f"""
# Prediction Summary

{domain_info}

## Top Prediction
**Family**: {top_pred['family']}  
**Probability**: {top_pred['probability']:.1%}  
**Confidence**: {confidence}

{result['recommendation']}
"""
    return summary


def format_predictions(result):
    """Format top-3 predictions with bar chart."""
    predictions = result['top_predictions']
    
    # Create plotly bar chart
    families = [p['family'] for p in predictions]
    probs = [p['probability'] for p in predictions]
    
    fig = go.Figure(data=[
        go.Bar(
            x=probs,
            y=families,
            orientation='h',
            marker=dict(
                color=['#2ecc71' if i == 0 else '#3498db' for i in range(len(families))]
            ),
            text=[f"{p:.1%}" for p in probs],
            textposition='auto',
        )
    ])
    
    fig.update_layout(
        title="Top-3 Predictions",
        xaxis_title="Probability",
        yaxis_title="",
        height=250,
        margin=dict(l=120, r=20, t=40, b=40),
        xaxis=dict(range=[0, 1])
    )
    
    return fig


def format_motifs(result):
    """Format motif information."""
    motifs = result['motifs']
    
    motif_text = "# Kinase Motifs\n\n"
    
    core_motifs = ['VAIK', 'HRD', 'DFG', 'APE', 'P-loop']
    
    for name in core_motifs:
        if name in motifs and motifs[name]['found']:
            m = motifs[name]
            motif_text += f"✅ **{name}**: {m['sequence']} (position {m['position']})\n"
        else:
            motif_text += f"❌ **{name}**: Not found\n"
    
    motif_text += f"\n**Gatekeeper**: {motifs.get('gatekeeper', {}).get('residue', 'N/A')}\n"
    
    if 'K_E_distance' in motifs:
        motif_text += f"**K-E Salt Bridge Distance**: {motifs['K_E_distance']} residues\n"
    
    if 'motif_integrity_score' in motifs:
        score = motifs['motif_integrity_score']
        motif_text += f"\n**Integrity Score**: {score:.2f}/1.0 "
        if score > 0.7:
            motif_text += "🟢 Good"
        elif score > 0.5:
            motif_text += "🟡 Acceptable"
        else:
            motif_text += "🔴 Low"
    
    return motif_text


def format_exemplars(result):
    """Format nearest training exemplars."""
    if 'nearest_exemplars' not in result or not result['nearest_exemplars']:
        return "No exemplars available"
    
    exemplar_text = "# Nearest Training Sequences\n\n"
    
    for i, ex in enumerate(result['nearest_exemplars'][:5], 1):
        exemplar_text += f"**{i}. {ex['uniprot_id']}** ({ex['family']})\n"
        exemplar_text += f"   - Similarity: {ex['similarity']:.3f}\n"
        exemplar_text += f"   - Protein: {ex['protein_name'][:60]}...\n\n"
    
    return exemplar_text


# Create Gradio interface
with gr.Blocks(title="Kinase Family Classifier", theme=gr.themes.Soft()) as app:
    gr.Markdown("""
    # 🧬 Kinase Family Classifier
    
    Predict kinase family from protein sequence using ESM-2 embeddings with mid-layer averaging.
    
    **Supported families**: AGC, CAMK, CK1, CMGC, STE, TK, TKL, Atypical
    
    **Citation**: Layer Selection in Protein Language Models for Kinase Classification (2025)
    """)
    
    with gr.Row():
        with gr.Column(scale=1):
            sequence_input = gr.Textbox(
                label="Protein Sequence",
                placeholder="Paste FASTA format or raw sequence...\n\nExample:\n>MyKinase\nMENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPS...",
                lines=10
            )
            
            predict_btn = gr.Button("🔬 Predict Family", variant="primary", size="lg")
            
            gr.Markdown("""
            ### Example Sequences
            Click to load:
            """)
            
            with gr.Row():
                example_cmgc = gr.Button("CDK2 (CMGC)", size="sm")
                example_tk = gr.Button("EGFR (TK)", size="sm")
                example_camk = gr.Button("CaMK2 (CAMK)", size="sm")
    
    with gr.Row():
        summary_output = gr.Markdown(label="Summary")
    
    with gr.Tabs():
        with gr.Tab("📊 Predictions"):
            predictions_plot = gr.Plot(label="Top-3 Predictions")
        
        with gr.Tab("🎯 Motifs"):
            with gr.Row():
                with gr.Column():
                    motifs_output = gr.Markdown(label="Motif Analysis")
                with gr.Column():
                    sequence_viz = gr.HTML(label="Sequence Visualization")
        
        with gr.Tab("🔍 Similar Sequences"):
            exemplars_output = gr.Markdown(label="Nearest Training Exemplars")
        
        with gr.Tab("💾 Export"):
            json_output = gr.Code(label="Full Results (JSON)", language="json")
    
    # Wire up predict button
    predict_btn.click(
        fn=predict_kinase,
        inputs=sequence_input,
        outputs=[summary_output, predictions_plot, motifs_output, 
                 exemplars_output, sequence_viz, json_output]
    )
    
    # Example sequences
    def load_example_cdk2():
        return ">CDK2_HUMAN\nMENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTTCDLKICDFGLARIADPEHDHTGFLTEYVATRWYRAPEIMLNSKGYTKSIDIWS"
    
    def load_example_egfr():
        return ">EGFR_HUMAN\nFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWEL"
    
    def load_example_camk2():
        return ">CAMK2_HUMAN\nSTTITSKEKKDGKAVVFKGVNVATANKATGAKVKALIDKFLKQKWDNESKLKSKKAKGAAKGAAILSLLDVKFKKLDRDGNTTQIVKTSKELGGNILSVMEYNPHEGITLYINAKHLETHVLGHQFLTGEDQSMLVEKGVTGSAFGTVYKGQPEGEVVKIPVAIKVLREATTPPKANREILDE"
    
    example_cmgc.click(fn=load_example_cdk2, outputs=sequence_input)
    example_tk.click(fn=load_example_egfr, outputs=sequence_input)
    example_camk.click(fn=load_example_camk2, outputs=sequence_input)
    
    gr.Markdown("""
    ---
    
    **Model**: ESM-2 650M with mid-layer averaging (layers 20-33)  
    **Domain Extraction**: HMMER with Pfam PF00069 (E-value 0.01)  
    **Classifier**: Calibrated Logistic Regression  
    **Training**: Homology-aware splits (40% identity threshold)
    
    **Disclaimer**: Predictions are for research purposes only. Experimental validation recommended.
    
    [GitHub](https://github.com/jhaaj08/Kinases-Clustering) | [Paper](link-to-paper) | [Dataset](link-to-zenodo)
    """)


if __name__ == "__main__":
    print("="*80)
    print("Starting Kinase Classification Web App")
    print("="*80)
    print("\nServer will start at: http://localhost:7860")
    print("Press Ctrl+C to stop\n")
    
    app.launch(
        server_name="0.0.0.0",
        server_port=7860,
        share=False,
        show_error=True
    )

