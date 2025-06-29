#!/usr/bin/env python3
"""
Enhanced Example Usage for Metalloprotein Binding Efficiency Prediction Pipeline

This script demonstrates the enhanced pipeline with:
- Multi-algorithm binding site identification (MetalNet, Metal3D, bindEmbed21, AlphaFill)
- Environmental parameter coupling (temperature, pH, pressure, redox potential)
- Spatial discretization (1000-cube model)
- MESPEUS database integration
- CHED network analysis
- Advanced visualization (PyMOL, RF Diffusion)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import enhanced modules
from src.enhanced.enhanced_main import EnhancedMetalloproteinPipeline
from src.enhanced.enhanced_binding_kinetics import EnhancedBindingKinetics
from src.enhanced.enhanced_binding_site_identification import EnhancedBindingSiteIdentifier
from src.pdb_processor import PDBProcessor

def create_sample_pdb_with_metal_binding_sites():
    """Create a sample PDB file with metal-binding residues."""
    pdb_content = """ATOM      1  N   CYS A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  CYS A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   CYS A   1      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  O   CYS A   1      12.000  11.500  10.000  1.00 20.00           O
ATOM      5  CB  CYS A   1      11.000   9.000  10.000  1.00 20.00           C
ATOM      6  SG  CYS A   1      12.000   8.000  10.000  1.00 20.00           S
ATOM      7  N   HIS A   2      11.000  12.000  10.000  1.00 20.00           N
ATOM      8  CA  HIS A   2      11.500  13.000  10.000  1.00 20.00           C
ATOM      9  C   HIS A   2      12.000  14.000  10.000  1.00 20.00           C
ATOM     10  O   HIS A   2      12.500  14.500  10.000  1.00 20.00           O
ATOM     11  CB  HIS A   2      12.000  13.000  11.000  1.00 20.00           C
ATOM     12  CG  HIS A   2      13.000  13.000  11.000  1.00 20.00           C
ATOM     13  ND1 HIS A   2      13.500  13.000  12.000  1.00 20.00           N
ATOM     14  CD2 HIS A   2      14.000  13.000  11.000  1.00 20.00           C
ATOM     15  CE1 HIS A   2      14.500  13.000  12.000  1.00 20.00           C
ATOM     16  NE2 HIS A   2      14.500  13.000  12.000  1.00 20.00           N
ATOM     17  N   GLU A   3      12.000  15.000  10.000  1.00 20.00           N
ATOM     18  CA  GLU A   3      12.500  16.000  10.000  1.00 20.00           C
ATOM     19  C   GLU A   3      13.000  17.000  10.000  1.00 20.00           C
ATOM     20  O   GLU A   3      13.500  17.500  10.000  1.00 20.00           O
ATOM     21  CB  GLU A   3      13.000  16.000  11.000  1.00 20.00           C
ATOM     22  CG  GLU A   3      14.000  16.000  11.000  1.00 20.00           C
ATOM     23  CD  GLU A   3      14.500  16.000  12.000  1.00 20.00           C
ATOM     24  OE1 GLU A   3      15.000  16.000  12.500  1.00 20.00           O
ATOM     25  OE2 GLU A   3      14.000  16.000  13.000  1.00 20.00           O
ATOM     26  N   ASP A   4      13.000  18.000  10.000  1.00 20.00           N
ATOM     27  CA  ASP A   4      13.500  19.000  10.000  1.00 20.00           C
ATOM     28  C   ASP A   4      14.000  20.000  10.000  1.00 20.00           C
ATOM     29  O   ASP A   4      14.500  20.500  10.000  1.00 20.00           O
ATOM     30  CB  ASP A   4      14.000  19.000  11.000  1.00 20.00           C
ATOM     31  CG  ASP A   4      15.000  19.000  11.000  1.00 20.00           C
ATOM     32  OD1 ASP A   4      15.500  19.000  12.000  1.00 20.00           O
ATOM     33  OD2 ASP A   4      15.500  19.000  10.000  1.00 20.00           O
ATOM     34  N   CYS A   5      14.000  21.000  10.000  1.00 20.00           N
ATOM     35  CA  CYS A   5      14.500  22.000  10.000  1.00 20.00           C
ATOM     36  C   CYS A   5      15.000  23.000  10.000  1.00 20.00           C
ATOM     37  O   CYS A   5      15.500  23.500  10.000  1.00 20.00           O
ATOM     38  CB  CYS A   5      15.000  22.000  11.000  1.00 20.00           C
ATOM     39  SG  CYS A   5      16.000  22.000  11.000  1.00 20.00           S
ATOM     40  N   HIS A   6      15.000  24.000  10.000  1.00 20.00           N
ATOM     41  CA  HIS A   6      15.500  25.000  10.000  1.00 20.00           C
ATOM     42  C   HIS A   6      16.000  26.000  10.000  1.00 20.00           C
ATOM     43  O   HIS A   6      16.500  26.500  10.000  1.00 20.00           O
ATOM     44  CB  HIS A   6      16.000  25.000  11.000  1.00 20.00           C
ATOM     45  CG  HIS A   6      17.000  25.000  11.000  1.00 20.00           C
ATOM     46  ND1 HIS A   6      17.500  25.000  12.000  1.00 20.00           N
ATOM     47  CD2 HIS A   6      18.000  25.000  11.000  1.00 20.00           C
ATOM     48  CE1 HIS A   6      18.500  25.000  12.000  1.00 20.00           C
ATOM     49  NE2 HIS A   6      18.500  25.000  12.000  1.00 20.00           N
TER
END
"""
    
    with open("sample_protein.pdb", "w") as f:
        f.write(pdb_content)
    
    return "sample_protein.pdb"

def demonstrate_enhanced_pipeline():
    """Demonstrate the enhanced pipeline with all features."""
    print("="*80)
    print("ENHANCED METALLOPROTEIN BINDING EFFICIENCY PREDICTION")
    print("DEMONSTRATION")
    print("="*80)
    
    # Create sample PDB file
    print("\n1. Creating sample protein structure with metal-binding sites...")
    pdb_file = create_sample_pdb_with_metal_binding_sites()
    print(f"   Sample PDB file created: {pdb_file}")
    
    # Initialize enhanced pipeline
    print("\n2. Initializing enhanced pipeline...")
    pipeline = EnhancedMetalloproteinPipeline()
    print("   Enhanced pipeline initialized successfully")
    
    # Define protein sequence
    protein_sequence = "CHEDCH"  # Cysteine, Histidine, Glutamate, Aspartate, Cysteine, Histidine
    
    # Define metal ions and concentrations
    metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
    initial_concentrations = {
        'Zn2+': 1e-6,  # 1 μM
        'Cu2+': 1e-6,  # 1 μM
        'Fe2+': 1e-6   # 1 μM
    }
    
    print(f"\n3. Protein sequence: {protein_sequence}")
    print(f"   Metal ions: {metal_ions}")
    print(f"   Initial concentrations: {initial_concentrations}")
    
    # Run enhanced analysis
    print("\n4. Running enhanced analysis pipeline...")
    results = pipeline.run_enhanced_analysis(
        pdb_file=pdb_file,
        protein_sequence=protein_sequence,
        metal_ions=metal_ions,
        initial_concentrations=initial_concentrations,
        time_span=(0, 1000),
        save_results=True
    )
    
    # Display results
    print("\n5. Analysis Results:")
    print("-" * 50)
    
    # Binding site results
    binding_sites = results['binding_sites']
    print(f"   Binding Sites Identified: {binding_sites['total_sites']}")
    print(f"   Average Consensus Score: {binding_sites['average_consensus_score']:.3f}")
    
    print("\n   Algorithm Performance:")
    for algo, score in binding_sites['algorithm_scores'].items():
        print(f"     {algo}: {score:.3f}")
    
    # Display individual binding sites
    if binding_sites['binding_sites']:
        print(f"\n   Individual Binding Sites:")
        for i, site in enumerate(binding_sites['binding_sites']):
            print(f"     Site {i+1}:")
            print(f"       Center: [{site['center'][0]:.2f}, {site['center'][1]:.2f}, {site['center'][2]:.2f}]")
            print(f"       Radius: {site['radius']:.2f} Å")
            print(f"       Consensus Score: {site['consensus_score']:.3f}")
            print(f"       Algorithms Agreeing: {site['algorithm_count']}")
            print(f"       Residues: {site['residues']}")
    
    # Efficiency results
    if results['efficiency_results'] is not None:
        efficiency = results['efficiency_results']
        print(f"\n   Binding Efficiency: {efficiency['overall_efficiency']:.3f}")
        
        print(f"\n   Environmental Analysis:")
        env = efficiency['environmental_analysis']
        print(f"     Temperature: {env['temperature']['mean']-273.15:.1f}°C ± {env['temperature']['std']:.1f}°C")
        print(f"     pH: {env['pH']['mean']:.2f} ± {env['pH']['std']:.2f}")
        print(f"     Pressure: {env['pressure']['mean']:.2f} ± {env['pressure']['std']:.2f} atm")
        print(f"     Redox Potential: {env['redox_potential']['mean']:.3f} ± {env['redox_potential']['std']:.3f} V")
        
        print(f"\n   Final Concentrations:")
        for ion, conc in efficiency['final_free_concentrations'].items():
            print(f"     Free {ion}: {conc:.2e} M")
        for ion, conc in efficiency['final_bound_concentrations'].items():
            print(f"     Bound {ion}: {conc:.2e} M")
    
    # Brownian motion results
    if results['brownian_results'] is not None:
        brownian = results['brownian_results']
        print(f"\n   Brownian Motion Analysis:")
        for ion, diff_coeff in brownian['diffusion_coefficients'].items():
            print(f"     {ion} Diffusion Coefficient: {diff_coeff:.2e} Å²/s")
    
    print("\n" + "="*80)
    print("ENHANCED ANALYSIS COMPLETED SUCCESSFULLY!")
    print("="*80)
    print(f"Results saved to: {pipeline.output_dir}")
    print("\nGenerated files:")
    print("  - Plots: binding site identification, kinetics, environmental analysis")
    print("  - Data: JSON files with detailed results")
    print("  - Structures: PyMOL visualization scripts")
    print("  - Animations: RF Diffusion visualization scripts")
    print("  - Report: Comprehensive analysis report")
    
    return results

def demonstrate_environmental_parameter_effects():
    """Demonstrate the effects of environmental parameters on binding efficiency."""
    print("\n" + "="*80)
    print("ENVIRONMENTAL PARAMETER EFFECTS DEMONSTRATION")
    print("="*80)
    
    # Initialize components
    pipeline = EnhancedMetalloproteinPipeline()
    pdb_file = "sample_protein.pdb"
    protein_sequence = "CHEDCH"
    metal_ions = ['Zn2+']
    initial_concentrations = {'Zn2+': 1e-6}
    
    # Test different environmental conditions
    environmental_conditions = [
        {
            'name': 'Standard Conditions',
            'temperature': 298.15,  # 25°C
            'pH': 7.0,
            'pressure': 1.0,  # 1 atm
            'redox_potential': 0.0  # V
        },
        {
            'name': 'High Temperature',
            'temperature': 323.15,  # 50°C
            'pH': 7.0,
            'pressure': 1.0,
            'redox_potential': 0.0
        },
        {
            'name': 'Low pH',
            'temperature': 298.15,
            'pH': 5.0,
            'pressure': 1.0,
            'redox_potential': 0.0
        },
        {
            'name': 'High Pressure',
            'temperature': 298.15,
            'pH': 7.0,
            'pressure': 10.0,  # 10 atm
            'redox_potential': 0.0
        },
        {
            'name': 'Oxidizing Conditions',
            'temperature': 298.15,
            'pH': 7.0,
            'pressure': 1.0,
            'redox_potential': 0.3  # V
        }
    ]
    
    results_comparison = []
    
    for condition in environmental_conditions:
        print(f"\nTesting {condition['name']}...")
        
        # Update configuration for this condition
        pipeline.config['environmental_conditions']['temperature']['initial'] = condition['temperature']
        pipeline.config['environmental_conditions']['pH']['initial'] = condition['pH']
        pipeline.config['environmental_conditions']['pressure']['initial'] = condition['pressure']
        pipeline.config['environmental_conditions']['redox_potential']['initial'] = condition['redox_potential']
        
        # Reinitialize kinetics with new config
        pipeline.enhanced_kinetics = EnhancedBindingKinetics(pipeline.config)
        
        # Run analysis
        results = pipeline.run_enhanced_analysis(
            pdb_file=pdb_file,
            protein_sequence=protein_sequence,
            metal_ions=metal_ions,
            initial_concentrations=initial_concentrations,
            time_span=(0, 500),
            save_results=False
        )
        
        if results['efficiency_results'] is not None:
            efficiency = results['efficiency_results']['overall_efficiency']
            env_analysis = results['efficiency_results']['environmental_analysis']
            
            results_comparison.append({
                'condition': condition['name'],
                'efficiency': efficiency,
                'final_temp': env_analysis['temperature']['mean'] - 273.15,
                'final_ph': env_analysis['pH']['mean'],
                'final_pressure': env_analysis['pressure']['mean'],
                'final_redox': env_analysis['redox_potential']['mean']
            })
            
            print(f"   Binding Efficiency: {efficiency:.3f}")
    
    # Plot environmental effects
    if results_comparison:
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        conditions = [r['condition'] for r in results_comparison]
        efficiencies = [r['efficiency'] for r in results_comparison]
        temps = [r['final_temp'] for r in results_comparison]
        phs = [r['final_ph'] for r in results_comparison]
        pressures = [r['final_pressure'] for r in results_comparison]
        redox = [r['final_redox'] for r in results_comparison]
        
        # Binding efficiency
        ax1 = axes[0, 0]
        bars1 = ax1.bar(conditions, efficiencies, color='skyblue', alpha=0.7)
        ax1.set_ylabel('Binding Efficiency')
        ax1.set_title('Binding Efficiency Under Different Conditions')
        ax1.tick_params(axis='x', rotation=45)
        
        for bar, eff in zip(bars1, efficiencies):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{eff:.3f}', ha='center', va='bottom')
        
        # Final temperature
        ax2 = axes[0, 1]
        bars2 = ax2.bar(conditions, temps, color='red', alpha=0.7)
        ax2.set_ylabel('Final Temperature (°C)')
        ax2.set_title('Temperature Evolution')
        ax2.tick_params(axis='x', rotation=45)
        
        # Final pH
        ax3 = axes[1, 0]
        bars3 = ax3.bar(conditions, phs, color='blue', alpha=0.7)
        ax3.set_ylabel('Final pH')
        ax3.set_title('pH Evolution')
        ax3.tick_params(axis='x', rotation=45)
        
        # Final pressure
        ax4 = axes[1, 1]
        bars4 = ax4.bar(conditions, pressures, color='green', alpha=0.7)
        ax4.set_ylabel('Final Pressure (atm)')
        ax4.set_title('Pressure Evolution')
        ax4.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(pipeline.output_dir / "plots" / "environmental_effects.png", 
                   dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"\nEnvironmental effects plot saved to: {pipeline.output_dir / 'plots' / 'environmental_effects.png'}")

def demonstrate_spatial_discretization():
    """Demonstrate the spatial discretization features."""
    print("\n" + "="*80)
    print("SPATIAL DISCRETIZATION DEMONSTRATION")
    print("="*80)
    
    # Initialize enhanced kinetics
    pipeline = EnhancedMetalloproteinPipeline()
    kinetics = pipeline.enhanced_kinetics
    
    print(f"Spatial Grid: {kinetics.nx} × {kinetics.ny} × {kinetics.nz} = {kinetics.n_cubes} cubes")
    print(f"Cube Volume: {kinetics.cube_volume:.2e} m³")
    print(f"Chamber Dimensions: {pipeline.config['spatial_discretization']['chamber']['dimensions']}")
    
    # Show cube centers
    print(f"\nCube Centers (first 10):")
    for i in range(min(10, kinetics.n_cubes)):
        center = kinetics.cube_centers[i]
        print(f"  Cube {i+1}: [{center[0]:.2e}, {center[1]:.2e}, {center[2]:.2e}]")
    
    # Demonstrate environmental parameter distribution
    print(f"\nEnvironmental Parameter Distribution:")
    print(f"  Temperature: {np.mean(kinetics.T)-273.15:.1f}°C ± {np.std(kinetics.T):.1f}°C")
    print(f"  pH: {np.mean(kinetics.pH):.2f} ± {np.std(kinetics.pH):.2f}")
    print(f"  Pressure: {np.mean(kinetics.P):.2f} ± {np.std(kinetics.P):.2f} atm")
    print(f"  Redox Potential: {np.mean(kinetics.Eh):.3f} ± {np.std(kinetics.Eh):.3f} V")

def main():
    """Main function to run all demonstrations."""
    try:
        # Run main enhanced pipeline demonstration
        results = demonstrate_enhanced_pipeline()
        
        # Demonstrate environmental parameter effects
        demonstrate_environmental_parameter_effects()
        
        # Demonstrate spatial discretization
        demonstrate_spatial_discretization()
        
        print("\n" + "="*80)
        print("ALL DEMONSTRATIONS COMPLETED SUCCESSFULLY!")
        print("="*80)
        print("\nThe enhanced pipeline demonstrates:")
        print("✓ Multi-algorithm binding site identification")
        print("✓ Environmental parameter coupling")
        print("✓ Spatial discretization (1000-cube model)")
        print("✓ MESPEUS database integration")
        print("✓ CHED network analysis")
        print("✓ Advanced visualization capabilities")
        print("✓ Comprehensive reporting")
        
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 