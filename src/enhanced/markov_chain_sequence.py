"""
Markov Chain Model for Protein Sequence Prediction

This module implements a Markov Chain model for generating protein sequences
with metal-binding motifs and topology constraints based on the context-free grammar.
"""

import numpy as np
import pandas as pd
from collections import defaultdict, Counter
import re
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class MarkovChainSequenceGenerator:
    """Markov Chain model for protein sequence generation with metal-binding motifs."""
    
    def __init__(self, order: int = 2):
        """
        Initialize Markov Chain sequence generator.
        
        Parameters:
        -----------
        order : int
            Order of the Markov chain (number of previous residues to consider)
        """
        self.order = order
        self.transition_matrix = defaultdict(Counter)
        self.metal_binding_motifs = {
            'Zn2+': ['C', 'H', 'E', 'D'],  # Cysteine, Histidine, Glutamate, Aspartate
            'Cu2+': ['C', 'H', 'E', 'D', 'M'],  # + Methionine
            'Fe2+': ['C', 'H', 'E', 'D', 'Y'],  # + Tyrosine
            'Mg2+': ['E', 'D', 'N', 'Q'],  # Glutamate, Aspartate, Asparagine, Glutamine
            'Ca2+': ['E', 'D', 'N', 'Q', 'S', 'T']  # + Serine, Threonine
        }
        
        # Amino acid properties for topology constraints
        self.aa_properties = {
            'helix_preference': {
                'A': 1.45, 'R': 0.79, 'N': 0.73, 'D': 0.98, 'C': 0.77,
                'E': 1.53, 'Q': 1.17, 'G': 0.53, 'H': 1.24, 'I': 1.00,
                'L': 1.34, 'K': 1.07, 'M': 1.20, 'F': 1.12, 'P': 0.59,
                'S': 0.79, 'T': 0.82, 'W': 1.14, 'Y': 0.61, 'V': 1.14
            },
            'strand_preference': {
                'A': 0.97, 'R': 0.90, 'N': 0.65, 'D': 0.80, 'C': 1.30,
                'E': 0.26, 'Q': 1.23, 'G': 0.81, 'H': 0.71, 'I': 1.60,
                'L': 1.22, 'K': 0.74, 'M': 1.67, 'F': 1.28, 'P': 0.62,
                'S': 0.72, 'T': 1.20, 'W': 1.19, 'Y': 1.29, 'V': 1.65
            }
        }
        
    def train_on_metalloprotein_data(self, sequences: List[str], metal_ions: List[str] = None):
        """
        Train the Markov chain on metalloprotein sequences.
        
        Parameters:
        -----------
        sequences : List[str]
            List of protein sequences
        metal_ions : List[str], optional
            List of metal ions associated with each sequence
        """
        print(f"Training Markov chain (order {self.order}) on {len(sequences)} sequences...")
        
        # Build transition matrix
        for i, sequence in enumerate(sequences):
            metal_ion = metal_ions[i] if metal_ions else 'Zn2+'
            self._update_transition_matrix(sequence, metal_ion)
        
        # Normalize transition probabilities
        self._normalize_transition_matrix()
        
        print(f"Training completed. Transition matrix built with {len(self.transition_matrix)} states.")
    
    def _update_transition_matrix(self, sequence: str, metal_ion: str):
        """Update transition matrix with a single sequence."""
        # Get preferred residues for this metal ion
        preferred_residues = self.metal_binding_motifs.get(metal_ion, ['C', 'H', 'E', 'D'])
        
        # Process sequence with sliding window
        for i in range(len(sequence) - self.order):
            # Current state (previous residues)
            current_state = sequence[i:i + self.order]
            next_residue = sequence[i + self.order]
            
            # Boost probability for metal-binding residues
            if next_residue in preferred_residues:
                self.transition_matrix[current_state][next_residue] += 2.0
            else:
                self.transition_matrix[current_state][next_residue] += 1.0
    
    def _normalize_transition_matrix(self):
        """Normalize transition probabilities."""
        for state in self.transition_matrix:
            total = sum(self.transition_matrix[state].values())
            for residue in self.transition_matrix[state]:
                self.transition_matrix[state][residue] /= total
    
    def generate_sequence_with_topology(self, topology_string: str, target_length: int, 
                                      metal_ion: str = 'Zn2+', binding_sites: List[int] = None) -> str:
        """
        Generate a protein sequence that matches the given topology string.
        
        Parameters:
        -----------
        topology_string : str
            Topology string from the grammar (e.g., '-C+0+B+0-B-1+C-1')
        target_length : int
            Target length of the protein sequence
        metal_ion : str
            Target metal ion for binding
        binding_sites : List[int], optional
            Positions where metal-binding motifs should be placed
        
        Returns:
        --------
        str : Generated protein sequence
        """
        print(f"Generating sequence for topology: {topology_string}")
        print(f"Target length: {target_length}, Metal ion: {metal_ion}")
        
        # Parse topology string to get structural elements
        structural_elements = self._parse_topology_string(topology_string)
        
        # Calculate segment lengths
        segment_lengths = self._calculate_segment_lengths(structural_elements, target_length)
        
        # Generate sequence segment by segment
        sequence = ""
        current_pos = 0
        
        for i, (element, length) in enumerate(zip(structural_elements, segment_lengths)):
            print(f"  Generating segment {i+1}: {element} (length {length})")
            
            # Generate sequence for this structural element
            segment = self._generate_structural_segment(
                element, length, metal_ion, 
                binding_sites, current_pos
            )
            
            sequence += segment
            current_pos += length
        
        # Ensure the sequence has the target length
        if len(sequence) < target_length:
            # Add flexible linker
            remaining = target_length - len(sequence)
            linker = self._generate_linker_sequence(remaining, metal_ion)
            sequence += linker
        elif len(sequence) > target_length:
            # Truncate to target length
            sequence = sequence[:target_length]
        
        print(f"Generated sequence length: {len(sequence)}")
        return sequence
    
    def _parse_topology_string(self, topology_string: str) -> List[Dict]:
        """
        Parse topology string into structural elements.
        
        Example: '-C+0+B+0-B-1+C-1' -> [{'layer': 'C', 'orientation': '-', 'position': '0'}, ...]
        """
        elements = []
        
        # Split by layer boundaries (A, B, C, D)
        pattern = r'([+-][ABCD][+-]?\d*)'
        matches = re.findall(pattern, topology_string)
        
        for match in matches:
            if len(match) >= 2:
                orientation = match[0]
                layer = match[1]
                position = match[2:] if len(match) > 2 else '0'
                
                elements.append({
                    'layer': layer,
                    'orientation': orientation,
                    'position': position
                })
        
        return elements
    
    def _calculate_segment_lengths(self, structural_elements: List[Dict], target_length: int) -> List[int]:
        """Calculate appropriate lengths for each structural segment."""
        n_segments = len(structural_elements)
        
        if n_segments == 0:
            return [target_length]
        
        # Base lengths for different structural elements
        base_lengths = {
            'A': 15,  # Alpha helix
            'B': 8,   # Beta strand
            'C': 12,  # Mixed
            'D': 10   # Other
        }
        
        # Calculate total base length
        total_base_length = sum(base_lengths.get(elem['layer'], 10) for elem in structural_elements)
        
        # Scale lengths to match target
        scale_factor = target_length / total_base_length
        
        lengths = []
        for elem in structural_elements:
            base_length = base_lengths.get(elem['layer'], 10)
            scaled_length = max(5, int(base_length * scale_factor))
            lengths.append(scaled_length)
        
        # Adjust to exactly match target length
        total_generated = sum(lengths)
        if total_generated != target_length:
            # Distribute difference proportionally
            diff = target_length - total_generated
            for i in range(len(lengths)):
                if diff == 0:
                    break
                if diff > 0:
                    lengths[i] += 1
                    diff -= 1
                else:
                    if lengths[i] > 5:
                        lengths[i] -= 1
                        diff += 1
        
        return lengths
    
    def _generate_structural_segment(self, element: Dict, length: int, metal_ion: str,
                                   binding_sites: List[int], current_pos: int) -> str:
        """Generate sequence for a specific structural element."""
        layer = element['layer']
        orientation = element['orientation']
        
        # Get preferred residues for this structural element
        if layer == 'A':  # Alpha helix
            preferred_residues = self._get_helix_preferred_residues()
        elif layer == 'B':  # Beta strand
            preferred_residues = self._get_strand_preferred_residues()
        else:  # Mixed or other
            preferred_residues = self._get_mixed_preferred_residues()
        
        # Check if this segment should contain a binding site
        segment_binding_sites = []
        if binding_sites:
            for site in binding_sites:
                if current_pos <= site < current_pos + length:
                    segment_binding_sites.append(site - current_pos)
        
        # Generate sequence
        sequence = ""
        for i in range(length):
            if i in segment_binding_sites:
                # Place metal-binding residue
                residue = self._select_metal_binding_residue(metal_ion, layer)
            else:
                # Use Markov chain to select residue
                if len(sequence) >= self.order:
                    current_state = sequence[-self.order:]
                    residue = self._select_next_residue(current_state, preferred_residues)
                else:
                    # Start with preferred residue for this structure
                    residue = np.random.choice(preferred_residues)
            
            sequence += residue
        
        return sequence
    
    def _get_helix_preferred_residues(self) -> List[str]:
        """Get residues preferred in alpha helices."""
        helix_residues = []
        for aa, score in self.aa_properties['helix_preference'].items():
            if score > 1.0:
                helix_residues.extend([aa] * int(score * 10))
        return helix_residues if helix_residues else ['A', 'L', 'E', 'K']
    
    def _get_strand_preferred_residues(self) -> List[str]:
        """Get residues preferred in beta strands."""
        strand_residues = []
        for aa, score in self.aa_properties['strand_preference'].items():
            if score > 1.0:
                strand_residues.extend([aa] * int(score * 10))
        return strand_residues if strand_residues else ['V', 'I', 'F', 'Y']
    
    def _get_mixed_preferred_residues(self) -> List[str]:
        """Get residues for mixed structural elements."""
        return ['A', 'G', 'S', 'T', 'N', 'Q']
    
    def _select_metal_binding_residue(self, metal_ion: str, layer: str) -> str:
        """Select appropriate metal-binding residue based on metal ion and structure."""
        preferred = self.metal_binding_motifs.get(metal_ion, ['C', 'H', 'E', 'D'])
        
        # Adjust based on structural context
        if layer == 'A':  # Alpha helix
            # Prefer H, E, D in helices
            helix_compatible = [r for r in preferred if r in ['H', 'E', 'D']]
            return np.random.choice(helix_compatible) if helix_compatible else 'H'
        elif layer == 'B':  # Beta strand
            # Prefer C, E, D in strands
            strand_compatible = [r for r in preferred if r in ['C', 'E', 'D']]
            return np.random.choice(strand_compatible) if strand_compatible else 'C'
        else:
            return np.random.choice(preferred)
    
    def _select_next_residue(self, current_state: str, preferred_residues: List[str]) -> str:
        """Select next residue using Markov chain with preference weighting."""
        if current_state in self.transition_matrix:
            # Get transition probabilities
            transitions = self.transition_matrix[current_state]
            
            # Weight by structural preferences
            weighted_transitions = {}
            for residue, prob in transitions.items():
                weight = 1.0
                if residue in preferred_residues:
                    weight = 2.0  # Boost preferred residues
                weighted_transitions[residue] = prob * weight
            
            # Normalize
            total = sum(weighted_transitions.values())
            if total > 0:
                for residue in weighted_transitions:
                    weighted_transitions[residue] /= total
                
                # Sample from weighted distribution
                residues = list(weighted_transitions.keys())
                probs = list(weighted_transitions.values())
                return np.random.choice(residues, p=probs)
        
        # Fallback to preferred residues
        return np.random.choice(preferred_residues)
    
    def _generate_linker_sequence(self, length: int, metal_ion: str) -> str:
        """Generate flexible linker sequence."""
        linker_residues = ['G', 'S', 'T', 'A']
        sequence = ""
        
        for i in range(length):
            if len(sequence) >= self.order:
                current_state = sequence[-self.order:]
                residue = self._select_next_residue(current_state, linker_residues)
            else:
                residue = np.random.choice(linker_residues)
            sequence += residue
        
        return sequence
    
    def validate_sequence(self, sequence: str, topology_string: str) -> Dict:
        """
        Validate generated sequence against topology constraints.
        
        Returns:
        --------
        dict : Validation results
        """
        validation = {
            'length_ok': len(sequence) > 0,
            'metal_binding_motifs': self._count_metal_binding_motifs(sequence),
            'structural_compatibility': self._check_structural_compatibility(sequence, topology_string),
            'overall_score': 0.0
        }
        
        # Calculate overall score
        score = 0.0
        if validation['length_ok']:
            score += 0.3
        if validation['metal_binding_motifs'] > 0:
            score += 0.4
        if validation['structural_compatibility']:
            score += 0.3
        
        validation['overall_score'] = score
        return validation
    
    def _count_metal_binding_motifs(self, sequence: str) -> int:
        """Count metal-binding motifs in the sequence."""
        motifs = 0
        for metal_ion, residues in self.metal_binding_motifs.items():
            for residue in residues:
                motifs += sequence.count(residue)
        return motifs
    
    def _check_structural_compatibility(self, sequence: str, topology_string: str) -> bool:
        """Check if sequence is compatible with topology."""
        # Simple check: ensure sequence has reasonable length
        if len(sequence) < 10:
            return False
        
        # Check for basic structural elements
        has_helix_residues = any(r in sequence for r in ['A', 'L', 'E', 'K'])
        has_strand_residues = any(r in sequence for r in ['V', 'I', 'F', 'Y'])
        
        return has_helix_residues or has_strand_residues

def load_metalloprotein_training_data() -> Tuple[List[str], List[str]]:
    """
    Load training data for metalloprotein sequences.
    This is a simplified version - in practice, you'd load from a database.
    
    Returns:
    --------
    Tuple[List[str], List[str]] : (sequences, metal_ions)
    """
    # Example metalloprotein sequences (simplified)
    sequences = [
        "MHHHHHHSSGGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFPTSRE",
        "MHHHHHHSSGGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFPTSRE",
        "MHHHHHHSSGGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFPTSRE"
    ]
    
    metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
    
    return sequences, metal_ions

if __name__ == "__main__":
    # Example usage
    generator = MarkovChainSequenceGenerator(order=2)
    
    # Load training data
    sequences, metal_ions = load_metalloprotein_training_data()
    
    # Train the model
    generator.train_on_metalloprotein_data(sequences, metal_ions)
    
    # Generate sequence with topology
    topology = "-C+0+B+0-B-1+C-1"
    sequence = generator.generate_sequence_with_topology(
        topology, target_length=150, metal_ion='Zn2+'
    )
    
    print(f"Generated sequence: {sequence}")
    
    # Validate sequence
    validation = generator.validate_sequence(sequence, topology)
    print(f"Validation results: {validation}") 