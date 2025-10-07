#!/usr/bin/env python3
"""
Provenance tracking utilities for reproducible research.

Captures:
- Data source versions (UniProt, Pfam, InterPro)
- Tool versions (HMMER, CD-HIT, Python packages)
- Processing steps and parameters
- Checksums of key files
"""

import os
import json
import hashlib
import subprocess
import datetime
from pathlib import Path
from typing import Dict, Any, Optional


class ProvenanceTracker:
    """Track data provenance for reproducibility."""
    
    def __init__(self, output_dir="data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.provenance_file = self.output_dir / "provenance.json"
        
        # Load existing provenance or create new
        if self.provenance_file.exists():
            with open(self.provenance_file, 'r') as f:
                self.data = json.load(f)
        else:
            self.data = {
                "created_at": datetime.datetime.now().isoformat(),
                "updated_at": None,
                "uniprot": {},
                "pfam_interpro": {},
                "hmmer": {},
                "cdhit": {},
                "environment": {},
                "processing_steps": [],
                "files": {},
                "splits": {},
                "inclusion_exclusion_rules": {}
            }
    
    def add_uniprot_info(self, query: str, response_headers: Dict = None, 
                        count: int = None):
        """Record UniProt download information."""
        self.data["uniprot"] = {
            "query": query,
            "downloaded_at": datetime.datetime.now().isoformat(),
            "count": count,
            "source": "UniProt SwissProt (reviewed only)",
            "api_endpoint": "https://rest.uniprot.org/uniprotkb/search",
        }
        
        if response_headers:
            # Extract release date from headers if available
            if 'date' in response_headers:
                self.data["uniprot"]["response_date"] = response_headers['date']
            if 'x-uniprot-release' in response_headers:
                self.data["uniprot"]["release"] = response_headers['x-uniprot-release']
        
        self.save()
    
    def add_pfam_info(self, pfam_ids: list, urls: Dict, response_metadata: Dict = None):
        """Record Pfam/InterPro HMM download information."""
        self.data["pfam_interpro"] = {
            "pfam_ids": pfam_ids,
            "pfam_names": {
                "PF00069": "Protein kinase domain (Ser/Thr/Tyr)",
                "PF07714": "Protein tyrosine kinase"
            },
            "endpoints": urls,
            "downloaded_at": datetime.datetime.now().isoformat(),
            "source": "InterPro/Pfam via EBI API",
        }
        
        if response_metadata:
            self.data["pfam_interpro"]["response_metadata"] = response_metadata
        
        self.save()
    
    def add_hmmer_info(self, version: str = None, command: str = None, 
                       params: Dict = None):
        """Record HMMER version and parameters."""
        # Get HMMER version
        if version is None:
            try:
                result = subprocess.run(['hmmsearch', '-h'], 
                                      capture_output=True, text=True, timeout=5)
                # Parse version from output
                for line in result.stdout.split('\n'):
                    if 'HMMER' in line and ('3.3' in line or '3.4' in line):
                        version = line.strip()
                        break
            except:
                version = "Unknown"
        
        self.data["hmmer"] = {
            "version": version,
            "command": command or "hmmsearch --domtblout OUTPUT -E THRESHOLD HMM FASTA",
            "parameters": params or {},
            "coordinate_system": "1-based (HMMER) → 0-based Python slice via start-1:end",
            "boundary_type": "envelope (env_from, env_to)",
        }
        
        self.save()
    
    def add_cdhit_info(self, version: str = None, thresholds: Dict = None):
        """Record CD-HIT version and thresholds."""
        # Get CD-HIT version
        if version is None:
            try:
                result = subprocess.run(['cd-hit', '-h'], 
                                      capture_output=True, text=True, timeout=5)
                # Parse version
                for line in result.stdout.split('\n')[:5]:
                    if 'CD-HIT' in line:
                        version = line.strip()
                        break
            except:
                version = "Unknown"
        
        self.data["cdhit"] = {
            "version": version,
            "thresholds": thresholds or {},
        }
        
        self.save()
    
    def add_environment_info(self):
        """Record Python environment and packages."""
        import sys
        import platform
        
        # Get package versions
        packages = {}
        try:
            import torch
            packages['torch'] = torch.__version__
        except:
            pass
        
        try:
            import esm
            packages['fair-esm'] = esm.__version__ if hasattr(esm, '__version__') else '2.0.0'
        except:
            pass
        
        try:
            import sklearn
            packages['scikit-learn'] = sklearn.__version__
        except:
            pass
        
        try:
            import pandas as pd
            packages['pandas'] = pd.__version__
        except:
            pass
        
        try:
            import numpy as np
            packages['numpy'] = np.__version__
        except:
            pass
        
        # Git commit hash
        git_commit = "Unknown"
        try:
            result = subprocess.run(['git', 'rev-parse', 'HEAD'], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                git_commit = result.stdout.strip()[:8]
        except:
            pass
        
        self.data["environment"] = {
            "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            "platform": platform.platform(),
            "os": platform.system(),
            "packages": packages,
            "git_commit": git_commit,
        }
        
        self.save()
    
    def add_processing_step(self, step_name: str, params: Dict, 
                           input_count: int, output_count: int, 
                           description: str = ""):
        """Record a processing step."""
        step = {
            "name": step_name,
            "timestamp": datetime.datetime.now().isoformat(),
            "description": description,
            "parameters": params,
            "input_count": input_count,
            "output_count": output_count,
            "removed_count": input_count - output_count,
        }
        
        self.data["processing_steps"].append(step)
        self.save()
    
    def add_file_info(self, file_id: str, file_path: str, 
                     description: str = "", compute_hash: bool = True):
        """Record file metadata and checksum."""
        file_info = {
            "path": str(file_path),
            "description": description,
            "timestamp": datetime.datetime.now().isoformat(),
        }
        
        if compute_hash and os.path.exists(file_path):
            file_info["sha256"] = self._compute_sha256(file_path)
            file_info["size_mb"] = os.path.getsize(file_path) / (1024 * 1024)
        
        self.data["files"][file_id] = file_info
        self.save()
    
    def add_inclusion_exclusion_rules(self, rules: Dict):
        """Record data inclusion/exclusion criteria."""
        self.data["inclusion_exclusion_rules"] = {
            "updated_at": datetime.datetime.now().isoformat(),
            **rules
        }
        self.save()
    
    def add_split_info(self, split_info: Dict):
        """Record train/test split information."""
        self.data["splits"] = {
            "updated_at": datetime.datetime.now().isoformat(),
            **split_info
        }
        self.save()
    
    def _compute_sha256(self, file_path: str, chunk_size: int = 8192) -> str:
        """Compute SHA256 hash of file."""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(chunk_size), b""):
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()
    
    def save(self):
        """Save provenance data to JSON."""
        self.data["updated_at"] = datetime.datetime.now().isoformat()
        
        with open(self.provenance_file, 'w') as f:
            json.dump(self.data, f, indent=2)
    
    def get_summary(self) -> str:
        """Get human-readable summary."""
        lines = []
        lines.append("="*80)
        lines.append("DATA PROVENANCE SUMMARY")
        lines.append("="*80)
        lines.append("")
        
        if self.data.get("uniprot"):
            lines.append("UniProt:")
            up = self.data["uniprot"]
            lines.append(f"  Downloaded: {up.get('downloaded_at', 'N/A')}")
            lines.append(f"  Count: {up.get('count', 'N/A')}")
            lines.append(f"  Query: {up.get('query', 'N/A')}")
            lines.append("")
        
        if self.data.get("pfam_interpro"):
            lines.append("Pfam/InterPro:")
            pf = self.data["pfam_interpro"]
            lines.append(f"  Profiles: {', '.join(pf.get('pfam_ids', []))}")
            lines.append(f"  Downloaded: {pf.get('downloaded_at', 'N/A')}")
            lines.append("")
        
        if self.data.get("hmmer"):
            lines.append("HMMER:")
            hm = self.data["hmmer"]
            lines.append(f"  Version: {hm.get('version', 'N/A')}")
            lines.append(f"  Parameters: {hm.get('parameters', {})}")
            lines.append("")
        
        if self.data.get("environment"):
            lines.append("Environment:")
            env = self.data["environment"]
            lines.append(f"  Python: {env.get('python_version', 'N/A')}")
            lines.append(f"  OS: {env.get('os', 'N/A')}")
            lines.append(f"  Git commit: {env.get('git_commit', 'N/A')}")
            lines.append("")
        
        if self.data.get("processing_steps"):
            lines.append(f"Processing Steps: {len(self.data['processing_steps'])}")
            for step in self.data["processing_steps"]:
                lines.append(f"  - {step['name']}: {step['input_count']} → {step['output_count']}")
            lines.append("")
        
        lines.append(f"Full details: {self.provenance_file}")
        lines.append("="*80)
        
        return "\n".join(lines)


def get_tool_version(tool_name: str) -> Optional[str]:
    """Get version of a command-line tool."""
    try:
        # Try common version flags
        for flag in ['-h', '--version', '-v', 'version']:
            try:
                result = subprocess.run([tool_name, flag], 
                                      capture_output=True, text=True, timeout=5)
                output = result.stdout + result.stderr
                
                # Parse version from output
                if tool_name in output.lower():
                    for line in output.split('\n')[:10]:
                        if any(v in line.lower() for v in ['version', tool_name.lower()]):
                            return line.strip()
                
                return output.split('\n')[0].strip()
            except:
                continue
    except:
        pass
    
    return None
