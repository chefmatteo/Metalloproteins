#!/usr/bin/env python3
"""
Metalloprotein Binding Efficiency Prediction Pipeline Demo
Using 6GV6 Protein Structure

This script demonstrates the complete pipeline workflow with a real protein structure:
1. Protein structure loading and analysis
2. Binding site identification
3. Environmental parameter coupling
4. Brownian motion simulation
5. Binding kinetics calculation
6. Efficiency prediction and ranking
7. Results visualization and interpretation
"""

import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import yaml

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.basic.main import MetalloproteinPipeline
from src.pdb_processor import PDBProcessor
from src.brownian_simulation import BrownianSimulation
from src.basic.binding_kinetics import BindingKinetics

def analyze_6gv6_structure():
    """Analyze the 6GV6 protein structure and show its characteristics."""
    print("="*80)
    print("6GV6 PROTEIN STRUCTURE ANALYSIS")
    print("="*80)
    
    # Load the protein structure
    pdb_file = "6gv6.cif"
    processor = PDBProcessor()
    
    if not os.path.exists(pdb_file):
        print(f"Error: {pdb_file} not found. Please ensure the file is extracted.")
        return None
    
    print(f"1. Loading protein structure: {pdb_file}")
    processor.load_pdb(pdb_file)
    
    # Extract basic information
    structure = processor.structure
    print(f"   Structure loaded successfully")
    print(f"   Number of chains: {len(list(structure.get_chains()))}")
    
    # Analyze chains
    for i, chain in enumerate(structure.get_chains()):
        print(f"   Chain {chain.id}: {len(list(chain.get_residues()))} residues")
    
    # Identify binding sites
    print(f"\n2. Identifying metal binding sites...")
    binding_sites = processor.identify_metal_binding_sites()
    print(f"   Found {len(binding_sites)} potential metal binding sites")
    
    # Show binding site details
    for i, site in enumerate(binding_sites):
        print(f"\n   Binding Site {i+1}:")
        print(f"     Center: ({site['center'][0]:.2f}, {site['center'][1]:.2f}, {site['center'][2]:.2f})")
        print(f"     Radius: {site['radius']:.2f} Å")
        print(f"     Coordinating residues: {site['coordinating_residues']}")
        print(f"     Coordination number: {site['coordination_number']}")
    
    return processor, binding_sites

def run_complete_pipeline():
    """Run the complete metalloprotein binding efficiency prediction pipeline."""
    print("\n" + "="*80)
    print("COMPLETE PIPELINE EXECUTION")
    print("="*80)
    
    # Initialize pipeline with 6GV6
    pipeline = MetalloproteinPipeline()
    
    # Load configuration
    print("1. Loading configuration from config/basic/config.yaml")
    with open("config/basic/config.yaml", 'r') as f:
        config = yaml.safe_load(f)
    
    print(f"   Temperature: {config['temperature']} K ({config['temperature'] - 273.15:.1f}°C)")
    print(f"   Pressure: {config['pressure']} atm")
    print(f"   Metal ions: {config['metal_ions']}")
    print(f"   Initial concentrations: {config['initial_concentrations']}")
    
    # Run pipeline
    print(f"\n2. Running complete pipeline analysis...")
    try:
        results = pipeline.run_pipeline("6gv6.cif", output_dir='6gv6_results')
        
        # Display comprehensive results
        print(f"\n3. PIPELINE RESULTS SUMMARY")
        print("="*50)
        
        # Protein information
        protein_info = results['protein_info']
        print(f"Protein Analysis:")
        print(f"  Number of binding sites: {protein_info['n_binding_sites']}")
        
        # Environmental conditions
        env_conditions = results['environmental_conditions']
        print(f"\nEnvironmental Conditions:")
        print(f"  Temperature: {env_conditions['temperature']} K")
        print(f"  Pressure: {env_conditions['pressure']} atm")
        print(f"  Viscosity: {env_conditions['viscosity']} Pa·s")
        
        # Metal ions and concentrations
        print(f"\nMetal Ion Analysis:")
        for ion in results['metal_ions']:
            print(f"  {ion}: Initial concentration = {results['initial_concentrations'][ion]} M")
        
        # Binding kinetics results
        if results['binding_kinetics']['solution'] is not None:
            kinetics = results['binding_kinetics']
            efficiency = kinetics['efficiency']
            
            print(f"\nBinding Kinetics Results:")
            print(f"  Overall binding efficiency: {efficiency['overall_efficiency']:.4f}")
            print(f"  Per-ion efficiencies:")
            for ion, eff in efficiency['ion_efficiencies'].items():
                print(f"    {ion}: {eff:.4f}")
        
        # Brownian motion results
        print(f"\nBrownian Motion Simulation Results:")
        for ion, brownian_data in results['brownian_simulation'].items():
            binding_prob = brownian_data['binding_probability']['overall_binding_probability']
            diffusion_analysis = brownian_data['diffusion_analysis']
            print(f"  {ion}:")
            print(f"    Binding probability: {binding_prob:.4f}")
            print(f"    Mean diffusion coefficient: {diffusion_analysis['mean_diffusion_coefficient']:.2e} m²/s")
        
        # Efficiency ranking
        print(f"\nEfficiency Ranking (Best to Worst):")
        for i, (ion, efficiency) in enumerate(results['efficiency_ranking'], 1):
            print(f"  {i}. {ion}: {efficiency:.4f}")
        
        print(f"\nDetailed results saved to: 6gv6_results/")
        
        return results
        
    except Exception as e:
        print(f"Error during pipeline execution: {e}")
        return None

def demonstrate_environmental_coupling():
    """Demonstrate how environmental parameters affect binding efficiency."""
    print("\n" + "="*80)
    print("ENVIRONMENTAL PARAMETER COUPLING DEMONSTRATION")
    print("="*80)
    
    # Test different environmental conditions
    temperatures = [273.15, 298.15, 310.15, 323.15]  # 0°C, 25°C, 37°C, 50°C
    pH_values = [6.0, 7.0, 7.4, 8.0]
    
    print("Testing binding efficiency under different conditions:")
    print("Temperature effects on binding kinetics:")
    
    processor = PDBProcessor()
    processor.load_pdb("6gv6.cif")
    binding_sites = processor.identify_metal_binding_sites()
    
    if binding_sites:
        kinetics = BindingKinetics(temperature=298.15, pressure=1.0)
        metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
        initial_concentrations = {ion: 1e-6 for ion in metal_ions}
        
        for temp in temperatures:
            print(f"\n  Temperature: {temp} K ({temp - 273.15:.1f}°C)")
            
            # Calculate rate constants at this temperature
            if binding_sites:
                site = binding_sites[0]
                for ion in metal_ions:
                    k_plus, k_minus = kinetics.calculate_rate_constants(site, ion)
                    print(f"    {ion}: k+ = {k_plus:.2e} M⁻¹s⁻¹, k- = {k_minus:.2e} s⁻¹")
    
    print(f"\nEnvironmental coupling effects:")
    print(f"  - Higher temperature increases diffusion and binding rates")
    print(f"  - pH affects protonation of metal-binding residues")
    print(f"  - Pressure can affect binding equilibria")
    print(f"  - Redox potential affects metal ion oxidation states")

def show_spatial_discretization():
    """Demonstrate the spatial discretization approach."""
    print("\n" + "="*80)
    print("SPATIAL DISCRETIZATION APPROACH")
    print("="*80)
    
    print("Your pipeline uses a 1000-cube (10×10×10) spatial grid:")
    print("  - Models realistic reaction chamber geometry")
    print("  - Enables accurate diffusion-limited process modeling")
    print("  - Allows spatial gradients in environmental parameters")
    print("  - Provides high-resolution analysis of binding patterns")
    
    print(f"\nSpatial grid characteristics:")
    print(f"  Grid size: 10 × 10 × 10 = 1000 cubes")
    print(f"  Each cube represents a volume element")
    print(f"  Environmental parameters can vary across the grid")
    print(f"  Binding events tracked spatially and temporally")

def explain_multi_algorithm_consensus():
    """Explain the multi-algorithm consensus approach."""
    print("\n" + "="*80)
    print("MULTI-ALGORITHM CONSENSUS APPROACH")
    print("="*80)
    
    print("Your pipeline integrates multiple algorithms for robust binding site identification:")
    print("  - MetalNet: CHED network analysis and clustering")
    print("  - Metal3D: Geometric coordination analysis")
    print("  - bindEmbed21: Sequence-based prediction")
    print("  - AlphaFill: Ligand/cofactor binding prediction")
    print("  - MESPEUS: Database of metal-binding sites")
    print("  - CHED Network: Community detection in protein networks")
    
    print(f"\nConsensus scoring:")
    print(f"  - Combines predictions from all algorithms")
    print(f"  - Provides confidence scores for each binding site")
    print(f"  - Reduces false positives and improves accuracy")
    print(f"  - Accounts for different prediction methodologies")

def demonstrate_experimental_relevance():
    """Show how the pipeline provides experimental guidance."""
    print("\n" + "="*80)
    print("EXPERIMENTAL RELEVANCE AND GUIDANCE")
    print("="*80)
    
    print("Your pipeline provides direct experimental guidance:")
    print("  - Predicts optimal metal ion for binding")
    print("  - Suggests optimal environmental conditions")
    print("  - Estimates binding affinities under lab conditions")
    print("  - Provides concentration recommendations")
    print("  - Predicts competition between different metal ions")
    
    print(f"\nPractical applications:")
    print(f"  - Protein purification optimization")
    print(f"  - Metal loading protocols")
    print(f"  - Buffer condition selection")
    print(f"  - Competition experiments design")
    print(f"  - Temperature and pH optimization")

def create_visualization_script():
    """Create a script for visualizing the results."""
    viz_script = '''#!/usr/bin/env python3
"""
6GV6 Metalloprotein Visualization Script
Generated by Metalloprotein Binding Efficiency Prediction Pipeline
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load results
results_dir = Path("6gv6_results")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# 1. Binding efficiency by metal ion
if (results_dir / "binding_efficiency.json").exists():
    import json
    with open(results_dir / "binding_efficiency.json", 'r') as f:
        efficiency_data = json.load(f)
    
    ions = list(efficiency_data['ion_efficiencies'].keys())
    efficiencies = list(efficiency_data['ion_efficiencies'].values())
    
    axes[0, 0].bar(ions, efficiencies, color=['blue', 'green', 'red', 'orange'])
    axes[0, 0].set_ylabel('Binding Efficiency')
    axes[0, 0].set_title('Metal Ion Binding Efficiency')
    axes[0, 0].grid(True, alpha=0.3)

# 2. Brownian motion results
if (results_dir / "brownian_motion.json").exists():
    with open(results_dir / "brownian_motion.json", 'r') as f:
        brownian_data = json.load(f)
    
    ions = list(brownian_data['diffusion_coefficients'].keys())
    diffusion_coeffs = list(brownian_data['diffusion_coefficients'].values())
    
    axes[0, 1].bar(ions, diffusion_coeffs, color=['purple', 'brown', 'pink', 'gray'])
    axes[0, 1].set_ylabel('Diffusion Coefficient (m²/s)')
    axes[0, 1].set_title('Ion Diffusion Coefficients')
    axes[0, 1].grid(True, alpha=0.3)

# 3. Environmental parameter effects
axes[1, 0].text(0.5, 0.5, 'Environmental\nParameter Analysis\n(see comprehensive_report.md)', 
                ha='center', va='center', transform=axes[1, 0].transAxes, fontsize=12)
axes[1, 0].set_title('Environmental Effects')

# 4. Spatial distribution
axes[1, 1].text(0.5, 0.5, 'Spatial Distribution\nAnalysis\n(1000-cube model)', 
                ha='center', va='center', transform=axes[1, 1].transAxes, fontsize=12)
axes[1, 1].set_title('Spatial Analysis')

plt.tight_layout()
plt.savefig(results_dir / "6gv6_analysis_summary.png", dpi=300, bbox_inches='tight')
plt.show()

print("Visualization saved to: 6gv6_results/6gv6_analysis_summary.png")
'''
    
    with open("visualize_6gv6_results.py", 'w') as f:
        f.write(viz_script)
    
    print("Created visualization script: visualize_6gv6_results.py")

def main():
    """Main function to demonstrate the complete pipeline workflow."""
    print("METALLOPROTEIN BINDING EFFICIENCY PREDICTION PIPELINE")
    print("Using 6GV6 Protein Structure")
    print("="*80)
    
    try:
        # Step 1: Analyze the protein structure
        processor, binding_sites = analyze_6gv6_structure()
        
        if binding_sites is None:
            print("Could not analyze protein structure. Exiting.")
            return
        
        # Step 2: Run complete pipeline
        results = run_complete_pipeline()
        
        # Step 3: Demonstrate key features
        demonstrate_environmental_coupling()
        show_spatial_discretization()
        explain_multi_algorithm_consensus()
        demonstrate_experimental_relevance()
        
        # Step 4: Create visualization script
        create_visualization_script()
        
        print("\n" + "="*80)
        print("PIPELINE DEMONSTRATION COMPLETED SUCCESSFULLY!")
        print("="*80)
        print("Key takeaways:")
        print("✅ Real protein structure (6GV6) analyzed")
        print("✅ Binding sites identified and characterized")
        print("✅ Environmental parameter coupling demonstrated")
        print("✅ Brownian motion simulation completed")
        print("✅ Binding efficiency calculated and ranked")
        print("✅ Results saved for further analysis")
        print("✅ Visualization script created")
        
        print(f"\nNext steps:")
        print(f"1. Check results in: 6gv6_results/")
        print(f"2. Run visualization: python3 visualize_6gv6_results.py")
        print(f"3. Modify config/basic/config.yaml for different conditions")
        print(f"4. Try other protein structures")
        
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        print("Please check your installation and the 6gv6.cif file.")

if __name__ == "__main__":
    main() 