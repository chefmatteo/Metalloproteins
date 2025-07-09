"""
Enhanced Metalloprotein Pipeline with Markov Chain and GROMACS Integration

This module integrates:
1. Markov Chain sequence generation with topology constraints
2. GROMACS structural optimization
3. Topology grammar parsing and validation
4. Enhanced binding efficiency prediction
"""

import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import our enhanced modules
from .markov_chain_sequence import MarkovChainSequenceGenerator
from .gromacs_optimizer import GROMACSOptimizer
from .topology_grammar_parser import TopologyGrammarParser
from .enhanced_binding_kinetics import EnhancedBindingKinetics
from .enhanced_binding_site_identification import EnhancedBindingSiteIdentifier
from ..pdb_processor import PDBProcessor
from ..brownian_simulation import BrownianSimulation

class EnhancedMetalloproteinPipelineWithMarkovGromacs:
    """Enhanced pipeline with Markov Chain and GROMACS integration."""
    
    def __init__(self, config_path="config/enhanced/enhanced_config.yaml"):
        """
        Initialize enhanced pipeline with all components.
        
        Parameters:
        -----------
        config_path : str
            Path to configuration file
        """
        self.config = self._load_config(config_path)
        
        # Initialize all components
        self.markov_generator = MarkovChainSequenceGenerator(order=2)
        self.gromacs_optimizer = GROMACSOptimizer()
        self.topology_parser = TopologyGrammarParser()
        self.enhanced_kinetics = EnhancedBindingKinetics(self.config)
        self.binding_site_identifier = EnhancedBindingSiteIdentifier(self.config)
        self.pdb_processor = PDBProcessor()
        self.brownian_sim = BrownianSimulation()
        
        # Create output directories
        self.output_dir = Path("output/enhanced_with_markov_gromacs")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        (self.output_dir / "sequences").mkdir(exist_ok=True)
        (self.output_dir / "structures").mkdir(exist_ok=True)
        (self.output_dir / "gromacs").mkdir(exist_ok=True)
        (self.output_dir / "analysis").mkdir(exist_ok=True)
        (self.output_dir / "plots").mkdir(exist_ok=True)
        
        # Train Markov model if training data is available
        self._train_markov_model()
    
    def _load_config(self, config_path):
        """Load configuration from YAML file."""
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    
    def _train_markov_model(self):
        """Train the Markov Chain model on metalloprotein data."""
        try:
            # Load training data (in practice, this would come from a database)
            sequences, metal_ions = self._load_training_data()
            if sequences:
                print("Training Markov Chain model on metalloprotein sequences...")
                self.markov_generator.train_on_metalloprotein_data(sequences, metal_ions)
            else:
                print("No training data available - using default Markov model")
        except Exception as e:
            print(f"Warning: Could not train Markov model: {e}")
    
    def _load_training_data(self):
        """Load training data for Markov model."""
        # In practice, this would load from a database or file
        # For now, we'll use some example sequences
        sequences = [
            "MHHHHHHSSGGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFPTSRE",
            "MHHHHHHSSGGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFPTSRE"
        ]
        metal_ions = ['Zn2+', 'Cu2+']
        return sequences, metal_ions
    
    def design_and_optimize_protein(self, topology_string: str, target_metal_ion: str,
                                  target_length: int = 150, binding_sites: List[int] = None,
                                  environmental_conditions: Dict = None) -> Dict:
        """
        Design and optimize a protein using the complete pipeline.
        
        Parameters:
        -----------
        topology_string : str
            Target topology string (e.g., '-C+0+B+0-B-1+C-1')
        target_metal_ion : str
            Target metal ion for binding
        target_length : int
            Target protein length
        binding_sites : List[int], optional
            Specific positions for metal binding sites
        environmental_conditions : Dict, optional
            Environmental conditions for binding analysis
        
        Returns:
        --------
        dict : Complete design and analysis results
        """
        print("="*80)
        print("ENHANCED METALLOPROTEIN DESIGN AND OPTIMIZATION PIPELINE")
        print("="*80)
        
        # Step 1: Parse and validate topology
        print("\n1. Parsing and validating topology...")
        topology_analysis = self._analyze_topology(topology_string)
        
        # Step 2: Generate sequence using Markov Chain
        print("\n2. Generating protein sequence using Markov Chain...")
        sequence_generation = self._generate_sequence_with_topology(
            topology_string, target_length, target_metal_ion, binding_sites
        )
        
        # Step 3: Validate sequence-topology compatibility
        print("\n3. Validating sequence-topology compatibility...")
        sequence_validation = self._validate_sequence_compatibility(
            sequence_generation['sequence'], topology_analysis['elements']
        )
        
        # Step 4: Optimize structure using GROMACS
        print("\n4. Optimizing protein structure using GROMACS...")
        structure_optimization = self._optimize_protein_structure(
            sequence_generation['sequence'], topology_string, [target_metal_ion]
        )
        
        # Step 5: Analyze binding efficiency
        print("\n5. Analyzing binding efficiency...")
        binding_analysis = self._analyze_binding_efficiency(
            structure_optimization['final_structure'], target_metal_ion, environmental_conditions
        )
        
        # Step 6: Generate comprehensive report
        print("\n6. Generating comprehensive report...")
        comprehensive_report = self._generate_comprehensive_report(
            topology_analysis, sequence_generation, sequence_validation,
            structure_optimization, binding_analysis, target_metal_ion
        )
        
        # Save all results
        self._save_complete_results(
            topology_analysis, sequence_generation, sequence_validation,
            structure_optimization, binding_analysis, comprehensive_report
        )
        
        return {
            'topology_analysis': topology_analysis,
            'sequence_generation': sequence_generation,
            'sequence_validation': sequence_validation,
            'structure_optimization': structure_optimization,
            'binding_analysis': binding_analysis,
            'comprehensive_report': comprehensive_report
        }
    
    def _analyze_topology(self, topology_string: str) -> Dict:
        """Analyze topology string using grammar parser."""
        try:
            # Parse topology
            elements = self.topology_parser.parse_topology_string(topology_string)
            
            # Generate summary
            summary = self.topology_parser.generate_topology_summary(elements)
            
            # Get structural guidance
            guidance = self.topology_parser.get_structural_guidance(elements)
            
            return {
                'topology_string': topology_string,
                'elements': elements,
                'summary': summary,
                'guidance': guidance,
                'valid': True
            }
        
        except ValueError as e:
            return {
                'topology_string': topology_string,
                'elements': [],
                'summary': {},
                'guidance': {},
                'valid': False,
                'error': str(e)
            }
    
    def _generate_sequence_with_topology(self, topology_string: str, target_length: int,
                                       target_metal_ion: str, binding_sites: List[int]) -> Dict:
        """Generate protein sequence using Markov Chain with topology constraints."""
        
        # Generate sequence
        sequence = self.markov_generator.generate_sequence_with_topology(
            topology_string, target_length, target_metal_ion, binding_sites
        )
        
        # Validate generated sequence
        validation = self.markov_generator.validate_sequence(sequence, topology_string)
        
        # Save sequence to file
        sequence_file = self.output_dir / "sequences" / "generated_sequence.fasta"
        with open(sequence_file, 'w') as f:
            f.write(f">Generated_Protein topology={topology_string} metal={target_metal_ion}\n")
            f.write(sequence + "\n")
        
        return {
            'sequence': sequence,
            'length': len(sequence),
            'validation': validation,
            'sequence_file': str(sequence_file),
            'topology_string': topology_string,
            'target_metal_ion': target_metal_ion
        }
    
    def _validate_sequence_compatibility(self, sequence: str, topology_elements) -> Dict:
        """Validate sequence compatibility with topology."""
        
        # Use topology parser to validate
        validation = self.topology_parser.validate_sequence_compatibility(sequence, topology_elements)
        
        # Get modification suggestions
        suggestions = self.topology_parser.suggest_sequence_modifications(sequence, topology_elements)
        validation['suggestions'] = suggestions
        
        return validation
    
    def _optimize_protein_structure(self, sequence: str, topology_string: str,
                                  metal_ions: List[str]) -> Dict:
        """Optimize protein structure using GROMACS."""
        
        # Run GROMACS optimization
        optimization_results = self.gromacs_optimizer.optimize_protein_structure(
            sequence=sequence,
            topology_string=topology_string,
            metal_ions=metal_ions,
            output_dir=str(self.output_dir / "gromacs")
        )
        
        # Validate binding site geometry
        if optimization_results['success']:
            geometry_validation = self.gromacs_optimizer.validate_binding_site_geometry(
                optimization_results['final_structure'], metal_ions[0]
            )
            optimization_results['geometry_validation'] = geometry_validation
        
        return optimization_results
    
    def _analyze_binding_efficiency(self, structure_file: str, target_metal_ion: str,
                                  environmental_conditions: Dict = None) -> Dict:
        """Analyze binding efficiency using enhanced pipeline."""
        
        # Load the optimized structure
        self.pdb_processor.load_pdb(structure_file)
        protein_structure = self.pdb_processor.structure
        
        # Identify binding sites
        binding_site_results = self.binding_site_identifier.identify_binding_sites(
            protein_structure, ""
        )
        
        # Set up environmental conditions
        if environmental_conditions is None:
            environmental_conditions = {
                'temperature': 298.15,
                'pH': 7.0,
                'pressure': 1.0,
                'redox_potential': 0.0
            }
        
        # Run binding kinetics analysis
        if binding_site_results['binding_sites']:
            kinetics_solution = self.enhanced_kinetics.solve_enhanced_kinetics(
                binding_site_results['binding_sites'],
                [target_metal_ion],
                {target_metal_ion: 1e-6},
                time_span=(0, 1000)
            )
            
            efficiency_results = self.enhanced_kinetics.calculate_enhanced_binding_efficiency(
                kinetics_solution,
                binding_site_results['binding_sites'],
                [target_metal_ion],
                {target_metal_ion: 1e-6}
            )
        else:
            kinetics_solution = None
            efficiency_results = None
        
        return {
            'binding_sites': binding_site_results,
            'kinetics_solution': kinetics_solution,
            'efficiency_results': efficiency_results,
            'environmental_conditions': environmental_conditions
        }
    
    def _generate_comprehensive_report(self, topology_analysis: Dict, sequence_generation: Dict,
                                     sequence_validation: Dict, structure_optimization: Dict,
                                     binding_analysis: Dict, target_metal_ion: str) -> str:
        """Generate comprehensive analysis report."""
        
        report = f"""
# ENHANCED METALLOPROTEIN DESIGN AND OPTIMIZATION REPORT

## Executive Summary
This report presents the results of a comprehensive protein design and optimization pipeline
that integrates Markov Chain sequence generation, GROMACS structural optimization, and
enhanced binding efficiency prediction.

## Design Parameters
- Target Topology: {topology_analysis['topology_string']}
- Target Metal Ion: {target_metal_ion}
- Generated Sequence Length: {sequence_generation['length']}

## 1. Topology Analysis
"""
        
        if topology_analysis['valid']:
            summary = topology_analysis['summary']
            report += f"""
- Structural Complexity: {summary['structural_complexity']}
- Total Elements: {summary['total_elements']}
- Estimated Length: {summary['estimated_length']} residues
- Orientation Balance: {summary['orientation_balance']['positive']} positive, {summary['orientation_balance']['negative']} negative
"""
        else:
            report += f"- Error: {topology_analysis.get('error', 'Unknown error')}\n"
        
        report += f"""
## 2. Sequence Generation
- Generated Sequence: {sequence_generation['sequence'][:50]}...
- Sequence Length: {sequence_generation['length']}
- Validation Score: {sequence_generation['validation']['overall_score']:.3f}
- Metal Binding Motifs: {sequence_generation['validation']['metal_binding_motifs']}

## 3. Sequence-Topology Compatibility
- Compatible: {sequence_validation['compatible']}
- Issues: {len(sequence_validation['issues'])}
- Warnings: {len(sequence_validation['warnings'])}
"""
        
        if sequence_validation['suggestions']:
            report += "\n### Suggestions for Improvement:\n"
            for suggestion in sequence_validation['suggestions']:
                report += f"- {suggestion}\n"
        
        report += f"""
## 4. Structure Optimization
- Optimization Success: {structure_optimization['success']}
- Final Structure: {structure_optimization['final_structure']}
"""
        
        if 'energy_minimization' in structure_optimization:
            em = structure_optimization['energy_minimization']
            report += f"""
### Energy Minimization Results:
- Initial Energy: {em.get('initial_energy', 'N/A'):.2f} kJ/mol
- Final Energy: {em.get('final_energy', 'N/A'):.2f} kJ/mol
- Energy Reduction: {em.get('energy_reduction', 'N/A'):.1f}%
"""
        
        report += f"""
## 5. Binding Efficiency Analysis
"""
        
        if binding_analysis['efficiency_results']:
            eff = binding_analysis['efficiency_results']
            report += f"""
- Overall Binding Efficiency: {eff['overall_efficiency']:.3f}
- Binding Sites Identified: {binding_analysis['binding_sites']['total_sites']}
- Average Consensus Score: {binding_analysis['binding_sites']['average_consensus_score']:.3f}
"""
        else:
            report += "- No binding sites identified or analysis failed\n"
        
        report += f"""
## 6. Environmental Conditions
- Temperature: {binding_analysis['environmental_conditions']['temperature']} K
- pH: {binding_analysis['environmental_conditions']['pH']}
- Pressure: {binding_analysis['environmental_conditions']['pressure']} atm
- Redox Potential: {binding_analysis['environmental_conditions']['redox_potential']} V

## 7. Recommendations
"""
        
        # Generate recommendations based on results
        recommendations = []
        
        if sequence_validation['overall_score'] < 0.5:
            recommendations.append("Consider regenerating sequence with improved topology compatibility")
        
        if binding_analysis['binding_sites']['total_sites'] == 0:
            recommendations.append("No binding sites identified - consider modifying sequence to include metal-binding motifs")
        
        if structure_optimization.get('energy_reduction', 0) < 50:
            recommendations.append("Structure optimization may need more steps or different parameters")
        
        for rec in recommendations:
            report += f"- {rec}\n"
        
        report += f"""
## Technical Details
- Markov Chain Order: {self.markov_generator.order}
- GROMACS Force Field: {self.gromacs_optimizer.force_field}
- Topology Grammar: Context-free grammar with {len(self.topology_parser.terminals)} terminal symbols
- Pipeline Version: Enhanced with Markov Chain and GROMACS integration

---
Report generated by Enhanced Metalloprotein Pipeline with Markov Chain and GROMACS Integration
"""
        
        return report
    
    def _save_complete_results(self, topology_analysis: Dict, sequence_generation: Dict,
                              sequence_validation: Dict, structure_optimization: Dict,
                              binding_analysis: Dict, comprehensive_report: str):
        """Save all results to files."""
        import json
        
        # Save topology analysis
        with open(self.output_dir / "analysis" / "topology_analysis.json", 'w') as f:
            json.dump(topology_analysis, f, indent=2, default=str)
        
        # Save sequence generation results
        with open(self.output_dir / "analysis" / "sequence_generation.json", 'w') as f:
            json.dump(sequence_generation, f, indent=2, default=str)
        
        # Save sequence validation
        with open(self.output_dir / "analysis" / "sequence_validation.json", 'w') as f:
            json.dump(sequence_validation, f, indent=2, default=str)
        
        # Save structure optimization results
        with open(self.output_dir / "analysis" / "structure_optimization.json", 'w') as f:
            json.dump(structure_optimization, f, indent=2, default=str)
        
        # Save binding analysis
        with open(self.output_dir / "analysis" / "binding_analysis.json", 'w') as f:
            json.dump(binding_analysis, f, indent=2, default=str)
        
        # Save comprehensive report
        with open(self.output_dir / "comprehensive_report.md", 'w') as f:
            f.write(comprehensive_report)
        
        print(f"\nAll results saved to: {self.output_dir}")
    
    def run_design_experiment(self, experiment_config: Dict) -> Dict:
        """
        Run a complete design experiment with multiple topologies and metal ions.
        
        Parameters:
        -----------
        experiment_config : Dict
            Configuration for the experiment
        
        Returns:
        --------
        dict : Experiment results
        """
        print("="*80)
        print("RUNNING DESIGN EXPERIMENT")
        print("="*80)
        
        results = {}
        
        for topology in experiment_config['topologies']:
            for metal_ion in experiment_config['metal_ions']:
                print(f"\nDesigning protein for topology: {topology}, metal: {metal_ion}")
                
                try:
                    result = self.design_and_optimize_protein(
                        topology_string=topology,
                        target_metal_ion=metal_ion,
                        target_length=experiment_config.get('target_length', 150),
                        environmental_conditions=experiment_config.get('environmental_conditions')
                    )
                    
                    results[f"{topology}_{metal_ion}"] = result
                    
                except Exception as e:
                    print(f"Error in design for {topology}_{metal_ion}: {e}")
                    results[f"{topology}_{metal_ion}"] = {'error': str(e)}
        
        # Generate experiment summary
        self._generate_experiment_summary(results, experiment_config)
        
        return results
    
    def _generate_experiment_summary(self, results: Dict, experiment_config: Dict):
        """Generate summary of design experiment."""
        
        summary = f"""
# DESIGN EXPERIMENT SUMMARY

## Experiment Configuration
- Topologies: {experiment_config['topologies']}
- Metal Ions: {experiment_config['metal_ions']}
- Target Length: {experiment_config.get('target_length', 150)}

## Results Summary
"""
        
        successful_designs = 0
        total_designs = len(results)
        
        for design_key, result in results.items():
            if 'error' not in result:
                successful_designs += 1
                if result['binding_analysis']['efficiency_results']:
                    efficiency = result['binding_analysis']['efficiency_results']['overall_efficiency']
                    summary += f"- {design_key}: Efficiency = {efficiency:.3f}\n"
                else:
                    summary += f"- {design_key}: No binding sites identified\n"
            else:
                summary += f"- {design_key}: Failed - {result['error']}\n"
        
        summary += f"""
## Statistics
- Total Designs: {total_designs}
- Successful Designs: {successful_designs}
- Success Rate: {successful_designs/total_designs*100:.1f}%
"""
        
        # Save experiment summary
        with open(self.output_dir / "experiment_summary.md", 'w') as f:
            f.write(summary)
        
        print(f"Experiment summary saved to: {self.output_dir / 'experiment_summary.md'}")

def main():
    """Main function to demonstrate the enhanced pipeline."""
    
    # Initialize enhanced pipeline
    pipeline = EnhancedMetalloproteinPipelineWithMarkovGromacs()
    
    # Example design experiment
    experiment_config = {
        'topologies': [
            "-C+0+B+0-B-1+C-1",  # Bacterial protein topology
            "+A+0",              # Simple alpha helix
            "+B+0-B+1+B+2"       # Beta sheet
        ],
        'metal_ions': ['Zn2+', 'Cu2+'],
        'target_length': 150,
        'environmental_conditions': {
            'temperature': 298.15,
            'pH': 7.0,
            'pressure': 1.0,
            'redox_potential': 0.0
        }
    }
    
    # Run design experiment
    results = pipeline.run_design_experiment(experiment_config)
    
    print("\n" + "="*80)
    print("DESIGN EXPERIMENT COMPLETED!")
    print("="*80)
    print(f"Results saved to: {pipeline.output_dir}")
    print("Check the output directory for detailed analysis and reports.")

if __name__ == "__main__":
    main() 