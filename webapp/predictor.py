"""
Kinase Prediction Module
========================

Handles the complete prediction pipeline:
1. Domain extraction (HMMER)
2. Embedding generation (ESM-2)
3. Classification (trained model)
4. Exemplar retrieval
5. Motif analysis
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import joblib
import torch
import esm
from sklearn.metrics.pairwise import cosine_similarity
from extract_motif_features import extract_all_features


class KinasePredictor:
    def __init__(self, model_path, embeddings_path, index_path, labels_path, hmm_path):
        """Initialize predictor with pre-trained models and data."""
        print("Loading classification model...")
        model_data = joblib.load(model_path)
        self.model = model_data['model']
        self.scaler = model_data['scaler']
        self.label_encoder = model_data['label_encoder']
        self.classes = self.label_encoder.classes_
        
        print("Loading ESM-2 model...")
        self.esm_model, self.alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.esm_model.eval()
        
        # Use CPU for web app (more accessible)
        self.device = 'cpu'
        self.esm_model = self.esm_model.to(self.device)
        
        print("Loading reference embeddings...")
        self.train_embeddings = np.load(embeddings_path)
        self.train_index = pd.read_csv(index_path)
        self.train_labels = pd.read_csv(labels_path)
        
        # Merge for exemplar retrieval
        self.train_data = self.train_index.merge(
            self.train_labels[['uniprot_id', 'kinome_group_major', 'protein_name']],
            on='uniprot_id',
            how='left'
        )
        
        self.hmm_path = hmm_path
        
        # Check HMMER
        try:
            subprocess.run(['hmmsearch', '-h'], capture_output=True, timeout=5)
            self.hmmer_available = True
        except:
            self.hmmer_available = False
            print("⚠️  HMMER not available - domain extraction disabled")
    
    def extract_domain(self, sequence, header="query"):
        """Extract kinase domain using HMMER."""
        if not self.hmmer_available:
            return None, "HMMER not available"
        
        # Create temp FASTA
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            fasta_file = f.name
            f.write(f">{header}\n{sequence}\n")
        
        # Create temp output
        with tempfile.NamedTemporaryFile(mode='w', suffix='.domtblout', delete=False) as f:
            output_file = f.name
        
        # Run HMMER
        try:
            cmd = ['hmmsearch', '--domtblout', output_file, '-E', '0.01', 
                   self.hmm_path, fasta_file]
            subprocess.run(cmd, capture_output=True, timeout=30)
            
            # Parse results
            domains = []
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.split()
                    if len(fields) >= 23:
                        domains.append({
                            'evalue': float(fields[12]),
                            'score': float(fields[13]),
                            'env_from': int(fields[19]),
                            'env_to': int(fields[20]),
                        })
            
            # Clean up
            os.unlink(fasta_file)
            os.unlink(output_file)
            
            if not domains:
                return None, "No kinase domain found (E-value > 0.01)"
            
            # Get best domain
            best = sorted(domains, key=lambda x: (x['evalue'], -x['score']))[0]
            domain_seq = sequence[best['env_from']-1:best['env_to']]
            
            return {
                'sequence': domain_seq,
                'start': best['env_from'],
                'end': best['env_to'],
                'evalue': best['evalue'],
                'score': best['score']
            }, None
            
        except subprocess.TimeoutExpired:
            return None, "HMMER timeout"
        except Exception as e:
            return None, f"HMMER error: {str(e)}"
    
    def generate_embedding(self, sequence):
        """Generate ESM-2 embedding (mid-layer averaging)."""
        layers = list(range(20, 34))  # Layers 20-33
        
        batch_converter = self.alphabet.get_batch_converter()
        data = [("seq", sequence)]
        labels, strs, tokens = batch_converter(data)
        tokens = tokens.to(self.device)
        
        with torch.no_grad():
            results = self.esm_model(tokens, repr_layers=layers, return_contacts=False)
        
        # Average across layers
        layer_reps = [results["representations"][layer][0] for layer in layers]
        reps = torch.stack(layer_reps).mean(dim=0)  # (seq_len, dim)
        
        # Mean pooling (exclude special tokens)
        pad_idx = self.alphabet.padding_idx
        cls_idx = self.alphabet.cls_idx
        eos_idx = self.alphabet.eos_idx
        
        tok = tokens[0]
        mask = (tok != pad_idx) & (tok != cls_idx) & (tok != eos_idx)
        residue_reps = reps[mask].cpu().numpy()
        
        # Sequence-level embedding
        seq_embedding = residue_reps.mean(axis=0)
        
        return seq_embedding
    
    def classify(self, embedding):
        """Classify using trained model."""
        # Standardize
        embedding_scaled = self.scaler.transform([embedding])
        
        # Predict
        proba = self.model.predict_proba(embedding_scaled)[0]
        
        # Top 3
        top_indices = np.argsort(proba)[::-1][:3]
        top_predictions = [
            {
                'family': self.classes[i],
                'probability': float(proba[i])
            }
            for i in top_indices
        ]
        
        # Confidence flag
        max_prob = proba[top_indices[0]]
        confidence = 'high' if max_prob >= 0.7 else 'low'
        
        return top_predictions, confidence
    
    def find_exemplars(self, embedding, n=5):
        """Find nearest training sequences."""
        # Compute similarities
        similarities = cosine_similarity([embedding], self.train_embeddings)[0]
        
        # Top n
        top_indices = np.argsort(similarities)[::-1][:n]
        
        exemplars = []
        for idx in top_indices:
            exemplars.append({
                'uniprot_id': self.train_data.iloc[idx]['uniprot_id'],
                'family': self.train_data.iloc[idx]['kinome_group_major'],
                'protein_name': self.train_data.iloc[idx]['protein_name'],
                'similarity': float(similarities[idx])
            })
        
        return exemplars
    
    def analyze_motifs(self, domain_sequence):
        """Extract and analyze motifs."""
        # Create temporary dataframe for motif extraction
        temp_df = pd.DataFrame({
            'uniprot_id': ['query'],
            'sequence': [domain_sequence]
        })
        
        # Extract features
        try:
            result_df = extract_all_features(temp_df)
            motifs = {}
            
            # Extract key motifs
            for col in result_df.columns:
                if col.startswith('motif_') or col.startswith('has_') or col == 'gatekeeper_residue':
                    motifs[col] = result_df.iloc[0][col]
            
            return motifs
        except:
            return {}
    
    def predict(self, sequence, header="query"):
        """Complete prediction pipeline."""
        result = {
            'status': 'success',
            'header': header,
            'sequence_length': len(sequence)
        }
        
        try:
            # Step 1: Extract domain
            domain_result, error = self.extract_domain(sequence, header)
            
            if error:
                result['status'] = 'error'
                result['message'] = error
                return result
            
            result['domain_sequence'] = domain_result['sequence']
            result['domain_start'] = domain_result['start']
            result['domain_end'] = domain_result['end']
            result['domain_length'] = len(domain_result['sequence'])
            result['domain_evalue'] = domain_result['evalue']
            
            # Step 2: Generate embedding
            embedding = self.generate_embedding(domain_result['sequence'])
            result['embedding_dim'] = len(embedding)
            
            # Step 3: Classify
            top_predictions, confidence = self.classify(embedding)
            result['top_predictions'] = top_predictions
            result['confidence_flag'] = confidence
            
            # Recommendation
            if confidence == 'high':
                result['recommendation'] = "✅ High confidence prediction - can be used directly"
            else:
                result['recommendation'] = "⚠️  Low confidence - recommend manual review or experimental validation"
            
            # Step 4: Find exemplars
            exemplars = self.find_exemplars(embedding, n=5)
            result['nearest_exemplars'] = exemplars
            
            # Step 5: Motif analysis
            motifs = self.analyze_motifs(domain_result['sequence'])
            result['motifs'] = motifs
            
            return result
            
        except Exception as e:
            result['status'] = 'error'
            result['message'] = str(e)
            return result

