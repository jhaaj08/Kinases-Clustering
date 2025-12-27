"""
Motif Highlighter Module
=========================

Highlights kinase motifs in sequences with color-coded HTML.
"""

import re


class MotifHighlighter:
    """Highlight kinase motifs in sequence."""
    
    MOTIF_PATTERNS = {
        'VAIK': (r'[VIL][AG][IV]K', '#e74c3c', 'β3-Lys (ATP binding)'),
        'HRD': (r'HRD', '#3498db', 'Catalytic loop'),
        'DFG': (r'DFG', '#2ecc71', 'Activation loop'),
        'APE': (r'APE', '#f39c12', 'Activation loop'),
        'P-loop': (r'G.G..G', '#9b59b6', 'ATP phosphate coordination'),
    }
    
    def find_motifs(self, sequence):
        """Find all motif positions in sequence."""
        motifs = {}
        
        for name, (pattern, color, description) in self.MOTIF_PATTERNS.items():
            matches = list(re.finditer(pattern, sequence))
            if matches:
                motifs[name] = {
                    'found': True,
                    'sequence': matches[0].group(),
                    'position': matches[0].start() + 1,  # 1-based
                    'color': color,
                    'description': description
                }
            else:
                motifs[name] = {
                    'found': False,
                    'color': color,
                    'description': description
                }
        
        return motifs
    
    def get_gatekeeper(self, sequence, motifs):
        """Identify gatekeeper residue (DFG-15 position)."""
        if 'DFG' in motifs and motifs['DFG']['found']:
            dfg_pos = motifs['DFG']['position'] - 1  # 0-based
            gatekeeper_pos = dfg_pos - 15
            
            if 0 <= gatekeeper_pos < len(sequence):
                return {
                    'residue': sequence[gatekeeper_pos],
                    'position': gatekeeper_pos + 1,  # 1-based
                    'color': '#e67e22'
                }
        
        return None


def highlight_sequence_html(sequence, motifs_dict):
    """
    Generate HTML with highlighted motifs.
    
    Args:
        sequence: Protein sequence
        motifs_dict: Motif information from extract_all_features
    
    Returns:
        HTML string with colored motifs
    """
    highlighter = MotifHighlighter()
    motifs = highlighter.find_motifs(sequence)
    gatekeeper = highlighter.get_gatekeeper(sequence, motifs)
    
    # Create position-to-motif mapping
    positions = {}  # position -> (motif_name, color)
    
    for name, info in motifs.items():
        if info['found']:
            start = info['position'] - 1  # 0-based
            length = len(info['sequence'])
            for i in range(start, start + length):
                positions[i] = (name, info['color'])
    
    # Add gatekeeper
    if gatekeeper:
        positions[gatekeeper['position'] - 1] = ('Gatekeeper', gatekeeper['color'])
    
    # Build HTML
    html = '<div style="font-family: monospace; font-size: 14px; line-height: 2.0; background: #f8f9fa; padding: 20px; border-radius: 8px;">'
    
    # Legend
    html += '<div style="margin-bottom: 15px; padding: 10px; background: white; border-radius: 5px;">'
    html += '<strong>Motif Legend:</strong><br>'
    for name, info in motifs.items():
        if info['found']:
            html += f'<span style="background: {info["color"]}; color: white; padding: 2px 6px; margin: 2px; border-radius: 3px;">{name}</span> {info["description"]}<br>'
    if gatekeeper:
        html += f'<span style="background: {gatekeeper["color"]}; color: white; padding: 2px 6px; margin: 2px; border-radius: 3px;">GK</span> Gatekeeper (DFG-15)<br>'
    html += '</div>'
    
    # Sequence with highlights
    html += '<div style="word-wrap: break-word;">'
    
    for i, aa in enumerate(sequence):
        if i % 60 == 0 and i > 0:
            html += '<br>'
        
        if i in positions:
            motif_name, color = positions[i]
            html += f'<span style="background: {color}; color: white; padding: 1px 2px; border-radius: 2px;" title="{motif_name}">{aa}</span>'
        else:
            html += f'<span style="color: #666;">{aa}</span>'
    
    html += '</div>'
    html += '</div>'
    
    return html

