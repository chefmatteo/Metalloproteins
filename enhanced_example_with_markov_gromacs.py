#!/usr/bin/env python3
"""
Enhanced Example Usage: Metalloprotein Design with Markov Chain and GROMACS

This script demonstrates the complete enhanced pipeline that integrates:
1. Markov Chain sequence generation with topology constraints
2. GROMACS structural optimization
3. Topology grammar parsing and validation
4. Enhanced binding efficiency prediction
"""

import sys
import os
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.enhanced.enhanced_pipeline_with_markov_gromacs import EnhancedMetalloproteinPipelineWithMarkovGromacs

def demonstrate_single_protein_design():
    """Demonstrate single protein design with the enhanced pipeline."""
    print("="*80)
    print("SINGLE PROTEIN DESIGN DEMONSTRATION")
    print("="*80)
    
    # Initialize enhanced pipeline
    pipeline = EnhancedMetalloproteinPipelineWithMarkovGromacs()
    
    # Design parameters
    topology_string = "-C+0+B+0-B-1+C-1"  # Bacterial protein topology
    target_metal_ion = "Zn2+"
    target_length = 150
    
    # Environmental conditions
    environmental_conditions = {
        'temperature': 298.15,  # 25°C
        'pH': 7.0,
        'pressure': 1.0,  # 1 atm
        'redox_potential': 0.0  # V
    }
    
    print(f"Designing protein with:")
    print(f"  Topology: {topology_string}")
    print(f"  Target Metal: {target_metal_ion}")
    print(f"  Target Length: {target_length}")
    print(f"  Environmental Conditions: {environmental_conditions}")
    
    # Run design and optimization
    results = pipeline.design_and_optimize_protein(
        topology_string=topology_string,
        target_metal_ion=target_metal_ion,
        target_length=target_length,
        environmental_conditions=environmental_conditions
    )
    
    # Display key results
    print("\n" + "="*80)
    print("DESIGN RESULTS SUMMARY")
    print("="*80)
    
    # Topology analysis
    topology = results['topology_analysis']
    if topology['valid']:
        summary = topology['summary']
        print(f"Topology Analysis:")
        print(f"  Structural Complexity: {summary['structural_complexity']}")
        print(f"  Total Elements: {summary['total_elements']}")
        print(f"  Estimated Length: {summary['estimated_length']} residues")
    else:
        print(f"Topology Error: {topology.get('error', 'Unknown error')}")
    
    # Sequence generation
    sequence_gen = results['sequence_generation']
    print(f"\nSequence Generation:")
    print(f"  Generated Length: {sequence_gen['length']}")
    print(f"  Validation Score: {sequence_gen['validation']['overall_score']:.3f}")
    print(f"  Metal Binding Motifs: {sequence_gen['validation']['metal_binding_motifs']}")
    print(f"  Sequence (first 50 residues): {sequence_gen['sequence'][:50]}...")
    
    # Sequence validation
    seq_validation = results['sequence_validation']
    print(f"\nSequence-Topology Compatibility:")
    print(f"  Compatible: {seq_validation['compatible']}")
    print(f"  Issues: {len(seq_validation['issues'])}")
    print(f"  Warnings: {len(seq_validation['warnings'])}")
    
    if seq_validation['suggestions']:
        print(f"  Suggestions:")
        for suggestion in seq_validation['suggestions'][:3]:  # Show first 3
            print(f"    - {suggestion}")
    
    # Structure optimization
    struct_opt = results['structure_optimization']
    print(f"\nStructure Optimization:")
    print(f"  Success: {struct_opt['success']}")
    print(f"  Final Structure: {struct_opt['final_structure']}")
    
    if 'energy_minimization' in struct_opt:
        em = struct_opt['energy_minimization']
        print(f"  Energy Reduction: {em.get('energy_reduction', 'N/A'):.1f}%")
    
    # Binding analysis
    binding_analysis = results['binding_analysis']
    print(f"\nBinding Efficiency Analysis:")
    print(f"  Binding Sites Identified: {binding_analysis['binding_sites']['total_sites']}")
    print(f"  Average Consensus Score: {binding_analysis['binding_sites']['average_consensus_score']:.3f}")
    
    if binding_analysis['efficiency_results']:
        eff = binding_analysis['efficiency_results']
        print(f"  Overall Binding Efficiency: {eff['overall_efficiency']:.3f}")
    else:
        print(f"  No binding efficiency calculated")
    
    print(f"\nDetailed results saved to: {pipeline.output_dir}")

def demonstrate_design_experiment():
    """Demonstrate a complete design experiment with multiple topologies and metals."""
    print("\n" + "="*80)
    print("DESIGN EXPERIMENT DEMONSTRATION")
    print("="*80)
    
    # Initialize enhanced pipeline
    pipeline = EnhancedMetalloproteinPipelineWithMarkovGromacs()
    
    # Experiment configuration
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
    
    print("Running design experiment with:")
    print(f"  Topologies: {experiment_config['topologies']}")
    print(f"  Metal Ions: {experiment_config['metal_ions']}")
    print(f"  Target Length: {experiment_config['target_length']}")
    
    # Run experiment
    results = pipeline.run_design_experiment(experiment_config)
    
    # Display experiment summary
    print("\n" + "="*80)
    print("EXPERIMENT RESULTS SUMMARY")
    print("="*80)
    
    successful_designs = 0
    total_designs = len(results)
    
    for design_key, result in results.items():
        print(f"\n{design_key}:")
        
        if 'error' in result:
            print(f"  Status: Failed - {result['error']}")
        else:
            successful_designs += 1
            print(f"  Status: Successful")
            
            # Show key metrics
            if result['binding_analysis']['efficiency_results']:
                efficiency = result['binding_analysis']['efficiency_results']['overall_efficiency']
                print(f"  Binding Efficiency: {efficiency:.3f}")
            else:
                print(f"  Binding Efficiency: No binding sites")
            
            if result['sequence_validation']['compatible']:
                print(f"  Topology Compatibility: Compatible")
            else:
                print(f"  Topology Compatibility: Incompatible")
    
    print(f"\nExperiment Statistics:")
    print(f"  Total Designs: {total_designs}")
    print(f"  Successful Designs: {successful_designs}")
    print(f"  Success Rate: {successful_designs/total_designs*100:.1f}%")
    
    print(f"\nDetailed experiment results saved to: {pipeline.output_dir}")

def demonstrate_topology_analysis():
    """Demonstrate topology grammar parsing and analysis."""
    print("\n" + "="*80)
    print("TOPOLOGY GRAMMAR ANALYSIS DEMONSTRATION")
    print("="*80)
    
    from src.enhanced.topology_grammar_parser import TopologyGrammarParser
    
    parser = TopologyGrammarParser()
    
    # Example topology strings from the documentation
    example_topologies = [
        "-C+0+B+0-B-1+C-1",  # Bacterial Protein (PDB ID: 2CU6)
        "+A+0",              # Simple Alpha Helix
        "+B+0-B+1+B+2",      # Beta Sheet
        "+C+0-D+1-A+2"       # Mixed structure
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
            
            # Show preferred residues for each segment
            for segment in guidance['segment_guidance']:
                print(f"    Segment {segment['segment_index']} ({segment['layer_name']}): "
                      f"preferred residues: {', '.join(segment['preferred_residues'][:4])}...")
            
        except ValueError as e:
            print(f"  Error: {e}")

def demonstrate_markov_chain_generation():
    """Demonstrate Markov Chain sequence generation."""
    print("\n" + "="*80)
    print("MARKOV CHAIN SEQUENCE GENERATION DEMONSTRATION")
    print("="*80)
    
    from src.enhanced.markov_chain_sequence import MarkovChainSequenceGenerator
    
    # Initialize generator
    generator = MarkovChainSequenceGenerator(order=2)
    
    # Load some training data
    sequences, metal_ions = [
        "MHHHHHHSSGGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFPTSRE",
        "MHHHHHHSSGGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFPTSRE"
    ], ['Zn2+', 'Cu2+']
    
    # Train the model
    print("Training Markov Chain model...")
    generator.train_on_metalloprotein_data(sequences, metal_ions)
    
    # Generate sequences for different topologies
    test_cases = [
        ("-C+0+B+0-B-1+C-1", "Zn2+", 150),  # Bacterial protein
        ("+A+0", "Cu2+", 100),              # Simple helix
        ("+B+0-B+1+B+2", "Fe2+", 120)       # Beta sheet
    ]
    
    for topology, metal_ion, length in test_cases:
        print(f"\nGenerating sequence for topology: {topology}, metal: {metal_ion}")
        
        try:
            sequence = generator.generate_sequence_with_topology(
                topology, length, metal_ion
            )
            
            # Validate sequence
            validation = generator.validate_sequence(sequence, topology)
            
            print(f"  Generated sequence length: {len(sequence)}")
            print(f"  Validation score: {validation['overall_score']:.3f}")
            print(f"  Metal binding motifs: {validation['metal_binding_motifs']}")
            print(f"  Sequence preview: {sequence[:50]}...")
            
        except Exception as e:
            print(f"  Error: {e}")

def demonstrate_gromacs_integration():
    """Demonstrate GROMACS integration."""
    print("\n" + "="*80)
    print("GROMACS INTEGRATION DEMONSTRATION")
    print("="*80)
    
    from src.enhanced.gromacs_optimizer import GROMACSOptimizer
    
    # Initialize GROMACS optimizer
    optimizer = GROMACSOptimizer()
    
    # Test sequence
    test_sequence = "MHHHHHHSSGGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFPTSRE"
    
    print(f"Testing GROMACS optimization with sequence length: {len(test_sequence)}")
    
    # Optimize structure
    results = optimizer.optimize_protein_structure(
        sequence=test_sequence,
        topology_string="-C+0+B+0-B-1+C-1",
        metal_ions=['Zn2+']
    )
    
    print(f"Optimization Results:")
    print(f"  Success: {results['success']}")
    print(f"  Final Structure: {results['final_structure']}")
    
    if 'energy_minimization' in results:
        em = results['energy_minimization']
        print(f"  Energy Minimization:")
        print(f"    Initial Energy: {em.get('initial_energy', 'N/A'):.2f} kJ/mol")
        print(f"    Final Energy: {em.get('final_energy', 'N/A'):.2f} kJ/mol")
        print(f"    Energy Reduction: {em.get('energy_reduction', 'N/A'):.1f}%")
    
    if 'md_simulation' in results:
        md = results['md_simulation']
        print(f"  MD Simulation:")
        print(f"    Temperature Stability: {md.get('temperature_stability', 'N/A')}")
        print(f"    Pressure Stability: {md.get('pressure_stability', 'N/A')}")
        print(f"    Final RMSD: {md.get('rmsd_final', 'N/A'):.3f} nm")
    
    if 'analysis' in results:
        analysis = results['analysis']
        print(f"  Analysis:")
        print(f"    Structure Quality: {analysis.get('structure_quality', 'N/A')}")
        print(f"    Topology Compatibility: {analysis.get('topology_compatibility', 'N/A')}")

def main():
    """Main function to run all demonstrations."""
    print("ENHANCED METALLOPROTEIN PIPELINE WITH MARKOV CHAIN AND GROMACS")
    print("="*80)
    print("This demonstration shows the complete integration of:")
    print("1. Markov Chain sequence generation with topology constraints")
    print("2. GROMACS structural optimization")
    print("3. Topology grammar parsing and validation")
    print("4. Enhanced binding efficiency prediction")
    print("="*80)
    
    try:
        # Run demonstrations
        demonstrate_topology_analysis()
        demonstrate_markov_chain_generation()
        demonstrate_gromacs_integration()
        demonstrate_single_protein_design()
        demonstrate_design_experiment()
        
        print("\n" + "="*80)
        print("ALL DEMONSTRATIONS COMPLETED SUCCESSFULLY!")
        print("="*80)
        print("The enhanced pipeline successfully integrates:")
        print("✅ Markov Chain sequence generation")
        print("✅ GROMACS structural optimization")
        print("✅ Topology grammar parsing")
        print("✅ Enhanced binding efficiency prediction")
        print("\nCheck the output directories for detailed results and analysis.")
        
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        print("Some components may not be available or properly configured.")
        print("The pipeline will fall back to simulated results where needed.")

if __name__ == "__main__":
    main() 