"""
Topology Grammar Parser for Protein Structure Analysis

This module implements the context-free grammar for parsing protein topology strings
as described in the variational encoder documentation.
"""

import re
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
import warnings
warnings.filterwarnings('ignore')

@dataclass
class TopologyElement:
    """Represents a single element in a topology string."""
    orientation: str  # '+' or '-'
    layer: str       # 'A', 'B', 'C', or 'D'
    position: str    # Position identifier (e.g., '0', '+1', '-2')
    
    def __str__(self):
        return f"{self.orientation}{self.layer}{self.position}"

class TopologyGrammarParser:
    """Parser for protein topology strings using context-free grammar."""
    
    def __init__(self):
        """Initialize the topology grammar parser."""
        
        # Define the grammar components as per the documentation
        self.non_terminals = {
            '<topology>',
            '<element>',
            '<orientation>',
            '<layer>',
            '<position>'
        }
        
        self.terminals = {
            '+', '-',
            'A', 'B', 'C', 'D',
            '+0', '+1', '+2', '+3', '+4', '+5', '+6',
            '-1', '-2', '-3', '-4', '-5', '-6'
        }
        
        # Production rules
        self.production_rules = {
            '<topology>': [
                ['<element>', '<topology>'],
                ['<element>']
            ],
            '<element>': [
                ['<orientation>', '<layer>', '<position>']
            ],
            '<orientation>': [
                ['+'],
                ['-']
            ],
            '<layer>': [
                ['A'], ['B'], ['C'], ['D']
            ],
            '<position>': [
                ['+0'], ['+1'], ['+2'], ['+3'], ['+4'], ['+5'], ['+6'],
                ['-1'], ['-2'], ['-3'], ['-4'], ['-5'], ['-6']
            ]
        }
        
        # Start symbol
        self.start_symbol = '<topology>'
        
        # Structural element definitions
        self.structural_elements = {
            'A': {
                'name': 'Alpha Helix',
                'description': 'Alpha helical structure',
                'typical_length': 15,
                'preferred_residues': ['A', 'L', 'E', 'K', 'R', 'Q']
            },
            'B': {
                'name': 'Beta Strand',
                'description': 'Beta sheet strand',
                'typical_length': 8,
                'preferred_residues': ['V', 'I', 'F', 'Y', 'W', 'T']
            },
            'C': {
                'name': 'Mixed Structure',
                'description': 'Mixed alpha-beta structure',
                'typical_length': 12,
                'preferred_residues': ['A', 'G', 'S', 'T', 'N', 'Q']
            },
            'D': {
                'name': 'Other Structure',
                'description': 'Other structural elements',
                'typical_length': 10,
                'preferred_residues': ['G', 'P', 'S', 'T', 'A']
            }
        }
    
    def parse_topology_string(self, topology_string: str) -> List[TopologyElement]:
        """
        Parse a topology string into structured elements.
        
        Parameters:
        -----------
        topology_string : str
            Topology string to parse (e.g., '-C+0+B+0-B-1+C-1')
        
        Returns:
        --------
        List[TopologyElement] : Parsed topology elements
        
        Raises:
        -------
        ValueError : If topology string is invalid
        """
        if not topology_string:
            raise ValueError("Topology string cannot be empty")
        
        # Validate the topology string format
        if not self._validate_topology_format(topology_string):
            raise ValueError(f"Invalid topology string format: {topology_string}")
        
        # Parse the string into elements
        elements = []
        pattern = r'([+-][ABCD][+-]?\d*)'
        matches = re.findall(pattern, topology_string)
        
        for match in matches:
            if len(match) >= 2:
                orientation = match[0]
                layer = match[1]
                position = match[2:] if len(match) > 2 else '0'
                
                # Validate individual components
                if orientation not in ['+', '-']:
                    raise ValueError(f"Invalid orientation: {orientation}")
                if layer not in ['A', 'B', 'C', 'D']:
                    raise ValueError(f"Invalid layer: {layer}")
                if not self._validate_position(position):
                    raise ValueError(f"Invalid position: {position}")
                
                element = TopologyElement(orientation, layer, position)
                elements.append(element)
        
        return elements
    
    def _validate_topology_format(self, topology_string: str) -> bool:
        """Validate the overall format of a topology string."""
        # Check for valid characters
        valid_chars = set('+-ABCD0123456')
        if not all(c in valid_chars for c in topology_string):
            return False
        
        # Check for proper element structure
        pattern = r'([+-][ABCD][+-]?\d*)'
        matches = re.findall(pattern, topology_string)
        
        # Reconstruct string from matches
        reconstructed = ''.join(matches)
        
        return reconstructed == topology_string
    
    def _validate_position(self, position: str) -> bool:
        """Validate position identifier."""
        valid_positions = {
            '+0', '+1', '+2', '+3', '+4', '+5', '+6',
            '-1', '-2', '-3', '-4', '-5', '-6'
        }
        return position in valid_positions
    
    def validate_sequence_compatibility(self, sequence: str, topology_elements: List[TopologyElement]) -> Dict:
        """
        Validate if a protein sequence is compatible with the topology.
        
        Parameters:
        -----------
        sequence : str
            Protein sequence to validate
        topology_elements : List[TopologyElement]
            Parsed topology elements
        
        Returns:
        --------
        dict : Validation results
        """
        validation = {
            'compatible': True,
            'issues': [],
            'warnings': [],
            'structural_analysis': {}
        }
        
        # Calculate expected sequence length
        expected_length = self._calculate_expected_length(topology_elements)
        
        if len(sequence) < expected_length * 0.5:
            validation['compatible'] = False
            validation['issues'].append(f"Sequence too short: {len(sequence)} < {expected_length * 0.5}")
        
        if len(sequence) > expected_length * 2.0:
            validation['warnings'].append(f"Sequence longer than expected: {len(sequence)} > {expected_length * 2.0}")
        
        # Analyze structural compatibility
        structural_analysis = self._analyze_structural_compatibility(sequence, topology_elements)
        validation['structural_analysis'] = structural_analysis
        
        # Check for structural element preferences
        for element in topology_elements:
            layer_info = self.structural_elements[element.layer]
            preferred_residues = layer_info['preferred_residues']
            
            # Count preferred residues in sequence
            preferred_count = sum(1 for residue in sequence if residue in preferred_residues)
            preferred_ratio = preferred_count / len(sequence)
            
            if preferred_ratio < 0.1:  # Less than 10% preferred residues
                validation['warnings'].append(
                    f"Low proportion of {layer_info['name']} preferred residues: {preferred_ratio:.2f}"
                )
        
        return validation
    
    def _calculate_expected_length(self, topology_elements: List[TopologyElement]) -> int:
        """Calculate expected sequence length based on topology elements."""
        total_length = 0
        
        for element in topology_elements:
            layer_info = self.structural_elements[element.layer]
            total_length += layer_info['typical_length']
        
        return total_length
    
    def _analyze_structural_compatibility(self, sequence: str, topology_elements: List[TopologyElement]) -> Dict:
        """Analyze structural compatibility between sequence and topology."""
        analysis = {
            'helix_content': 0.0,
            'strand_content': 0.0,
            'mixed_content': 0.0,
            'other_content': 0.0,
            'structural_balance': 'good'
        }
        
        # Count residues by structural preference
        helix_residues = ['A', 'L', 'E', 'K', 'R', 'Q']
        strand_residues = ['V', 'I', 'F', 'Y', 'W', 'T']
        mixed_residues = ['A', 'G', 'S', 'T', 'N', 'Q']
        other_residues = ['G', 'P', 'S', 'T', 'A']
        
        helix_count = sum(1 for r in sequence if r in helix_residues)
        strand_count = sum(1 for r in sequence if r in strand_residues)
        mixed_count = sum(1 for r in sequence if r in mixed_residues)
        other_count = sum(1 for r in sequence if r in other_residues)
        
        total = len(sequence)
        analysis['helix_content'] = helix_count / total
        analysis['strand_content'] = strand_count / total
        analysis['mixed_content'] = mixed_count / total
        analysis['other_content'] = other_count / total
        
        # Check structural balance
        max_content = max(analysis['helix_content'], analysis['strand_content'], 
                         analysis['mixed_content'], analysis['other_content'])
        
        if max_content > 0.6:
            analysis['structural_balance'] = 'biased'
        elif max_content < 0.1:
            analysis['structural_balance'] = 'unbalanced'
        
        return analysis
    
    def generate_topology_summary(self, topology_elements: List[TopologyElement]) -> Dict:
        """
        Generate a summary of the topology structure.
        
        Parameters:
        -----------
        topology_elements : List[TopologyElement]
            Parsed topology elements
        
        Returns:
        --------
        dict : Topology summary
        """
        summary = {
            'total_elements': len(topology_elements),
            'element_types': {},
            'orientation_balance': {'positive': 0, 'negative': 0},
            'estimated_length': 0,
            'structural_complexity': 'simple'
        }
        
        # Count element types
        for element in topology_elements:
            layer = element.layer
            if layer not in summary['element_types']:
                summary['element_types'][layer] = 0
            summary['element_types'][layer] += 1
            
            # Count orientations
            if element.orientation == '+':
                summary['orientation_balance']['positive'] += 1
            else:
                summary['orientation_balance']['negative'] += 1
        
        # Calculate estimated length
        summary['estimated_length'] = self._calculate_expected_length(topology_elements)
        
        # Determine structural complexity
        if len(topology_elements) <= 2:
            summary['structural_complexity'] = 'simple'
        elif len(topology_elements) <= 5:
            summary['structural_complexity'] = 'moderate'
        else:
            summary['structural_complexity'] = 'complex'
        
        return summary
    
    def suggest_sequence_modifications(self, sequence: str, topology_elements: List[TopologyElement]) -> List[str]:
        """
        Suggest modifications to improve sequence-topology compatibility.
        
        Parameters:
        -----------
        sequence : str
            Current protein sequence
        topology_elements : List[TopologyElement]
            Target topology elements
        
        Returns:
        --------
        List[str] : List of modification suggestions
        """
        suggestions = []
        
        # Analyze current sequence
        structural_analysis = self._analyze_structural_compatibility(sequence, topology_elements)
        
        # Check for structural imbalances
        if structural_analysis['helix_content'] < 0.1:
            suggestions.append("Consider adding more helix-favoring residues (A, L, E, K, R, Q)")
        
        if structural_analysis['strand_content'] < 0.1:
            suggestions.append("Consider adding more strand-favoring residues (V, I, F, Y, W, T)")
        
        if structural_analysis['structural_balance'] == 'biased':
            suggestions.append("Sequence is biased toward one structural type - consider balancing")
        
        # Check length compatibility
        expected_length = self._calculate_expected_length(topology_elements)
        if len(sequence) < expected_length * 0.7:
            suggestions.append(f"Sequence may be too short for topology (current: {len(sequence)}, expected: ~{expected_length})")
        
        if len(sequence) > expected_length * 1.5:
            suggestions.append(f"Sequence may be too long for topology (current: {len(sequence)}, expected: ~{expected_length})")
        
        return suggestions
    
    def get_structural_guidance(self, topology_elements: List[TopologyElement]) -> Dict:
        """
        Get structural guidance for sequence design based on topology.
        
        Parameters:
        -----------
        topology_elements : List[TopologyElement]
            Parsed topology elements
        
        Returns:
        --------
        dict : Structural guidance information
        """
        guidance = {
            'segment_guidance': [],
            'overall_recommendations': [],
            'metal_binding_considerations': []
        }
        
        for i, element in enumerate(topology_elements):
            layer_info = self.structural_elements[element.layer]
            
            segment_guide = {
                'segment_index': i + 1,
                'layer_type': element.layer,
                'layer_name': layer_info['name'],
                'preferred_residues': layer_info['preferred_residues'],
                'typical_length': layer_info['typical_length'],
                'orientation': element.orientation,
                'position': element.position
            }
            
            guidance['segment_guidance'].append(segment_guide)
        
        # Overall recommendations
        if len(topology_elements) > 3:
            guidance['overall_recommendations'].append("Complex topology - ensure proper transitions between segments")
        
        # Metal binding considerations
        guidance['metal_binding_considerations'].extend([
            "Consider placing metal-binding residues in accessible regions",
            "Avoid placing metal-binding sites in buried hydrophobic regions",
            "Ensure proper coordination geometry for target metal ion"
        ])
        
        return guidance

def example_topology_analysis():
    """Example usage of the topology grammar parser."""
    parser = TopologyGrammarParser()
    
    # Example topology strings from the documentation
    example_topologies = [
        "-C+0+B+0-B-1+C-1",  # Bacterial Protein (PDB ID: 2CU6)
        "+A+0",              # Simple Alpha Helix
        "+B+0-B+1+B+2"       # Beta Sheet
    ]
    
    for topology_string in example_topologies:
        print(f"\nAnalyzing topology: {topology_string}")
        
        try:
            # Parse topology
            elements = parser.parse_topology_string(topology_string)
            print(f"  Parsed {len(elements)} elements:")
            
            for i, element in enumerate(elements):
                layer_info = parser.structural_elements[element.layer]
                print(f"    {i+1}. {element} - {layer_info['name']}")
            
            # Generate summary
            summary = parser.generate_topology_summary(elements)
            print(f"  Summary: {summary['structural_complexity']} structure, "
                  f"estimated length: {summary['estimated_length']} residues")
            
            # Get structural guidance
            guidance = parser.get_structural_guidance(elements)
            print(f"  Guidance: {len(guidance['segment_guidance'])} segments identified")
            
        except ValueError as e:
            print(f"  Error: {e}")

if __name__ == "__main__":
    example_topology_analysis() 