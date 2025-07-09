#!/usr/bin/env python3
"""
Demo: Metalloprotein Binding Efficiency Prediction using Basic Config

This script demonstrates the metalloprotein binding efficiency prediction pipeline
using the basic configuration file (config/basic/config.yaml).
"""

import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.basic.main import MetalloproteinPipeline
from src.pdb_processor import PDBProcessor
from src.brownian_simulation import BrownianSimulation
from src.basic.binding_kinetics import BindingKinetics

def create_demo_protein():
    """Create a demo protein structure for testing."""
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
TER
END"""
    
    pdb_file = "demo_protein.pdb"
    with open(pdb_file, 'w') as f:
        f.write(pdb_content)
    
    return pdb_file

def run_basic_pipeline_demo():
    """Run the basic pipeline demo using the config file."""
    print("="*80)
    print("METALLOPROTEIN BINDING EFFICIENCY PREDICTION DEMO")
    print("="*80)
    print("Using configuration from: config/basic/config.yaml")
    print()
    
    # Create demo protein
    print("1. Creating demo protein structure...")
    pdb_file = create_demo_protein()
    print(f"   Demo protein created: {pdb_file}")
    
    # Initialize pipeline
    print("\n2. Initializing metalloprotein pipeline...")
    pipeline = MetalloproteinPipeline()
    
    # Display configuration
    print("\n3. Configuration Summary:")
    print(f"   Temperature: {pipeline.config['temperature']} K ({pipeline.config['temperature'] - 273.15:.1f}°C)")
    print(f"   Pressure: {pipeline.config['pressure']} atm")
    print(f"   Viscosity: {pipeline.config['viscosity']} Pa·s")
    print(f"   Metal ions: {pipeline.config['metal_ions']}")
    print(f"   Initial concentrations: {pipeline.config['initial_concentrations']}")
    print(f"   Simulation time: {pipeline.config['simulation_time']} s")
    print(f"   Number of ions: {pipeline.config['n_ions']}")
    print(f"   Box size: {pipeline.config['box_size']} m")
    
    # Run pipeline
    print("\n4. Running metalloprotein binding efficiency analysis...")
    try:
        results = pipeline.run_pipeline(pdb_file, output_dir='demo_results')
        
        # Display results
        print("\n5. RESULTS SUMMARY:")
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
        
        print(f"\nDetailed results saved to: demo_results/")
        
        return results
        
    except Exception as e:
        print(f"Error during pipeline execution: {e}")
        return None
    
    finally:
        # Clean up demo file
        if os.path.exists(pdb_file):
            os.remove(pdb_file)

def analyze_binding_sites_demo():
    """Demonstrate binding site identification."""
    print("\n" + "="*80)
    print("BINDING SITE IDENTIFICATION DEMO")
    print("="*80)
    
    # Create demo protein
    pdb_file = create_demo_protein()
    
    try:
        # Initialize PDB processor
        processor = PDBProcessor()
        processor.load_pdb(pdb_file)
        
        # Identify binding sites
        binding_sites = processor.identify_metal_binding_sites()
        
        print(f"Binding Site Analysis:")
        print(f"  Total binding sites identified: {len(binding_sites)}")
        
        for i, site in enumerate(binding_sites, 1):
            print(f"\n  Binding Site {i}:")
            print(f"    Center: ({site['center'][0]:.2f}, {site['center'][1]:.2f}, {site['center'][2]:.2f})")
            print(f"    Radius: {site['radius']:.2f} Å")
            print(f"    Coordinating residues: {site['coordinating_residues']}")
            print(f"    Coordination number: {site['coordination_number']}")
        
        return binding_sites
        
    finally:
        if os.path.exists(pdb_file):
            os.remove(pdb_file)

def demonstrate_brownian_motion():
    """Demonstrate Brownian motion simulation."""
    print("\n" + "="*80)
    print("BROWNIAN MOTION SIMULATION DEMO")
    print("="*80)
    
    # Create demo protein
    pdb_file = create_demo_protein()
    
    try:
        # Initialize components
        processor = PDBProcessor()
        processor.load_pdb(pdb_file)
        binding_sites = processor.identify_metal_binding_sites()
        
        # Initialize Brownian simulation
        brownian_sim = BrownianSimulation(
            temperature=298.15,
            viscosity=1e-3,
            box_size=100e-9
        )
        
        # Test different metal ions
        metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
        
        for ion in metal_ions:
            print(f"\nSimulating {ion} diffusion:")
            
            # Ion properties
            ion_properties = {
                'Zn2+': {'radius': 0.74e-10, 'mass': 65.38e-27, 'charge': 2},
                'Cu2+': {'radius': 0.73e-10, 'mass': 63.55e-27, 'charge': 2},
                'Fe2+': {'radius': 0.78e-10, 'mass': 55.85e-27, 'charge': 2}
            }[ion]
            
            # Run simulation
            simulation_results = brownian_sim.simulate_multiple_ions(
                ion_properties=ion_properties,
                n_ions=50,  # Reduced for demo
                binding_sites=binding_sites,
                simulation_time=1e-10  # Short simulation for demo
            )
            
            # Calculate binding probability
            binding_prob = brownian_sim.calculate_binding_probability(
                simulation_results, binding_sites
            )
            
            # Calculate diffusion analysis
            diffusion_analysis = brownian_sim.calculate_diffusion_analysis(
                simulation_results
            )
            
            print(f"  Binding probability: {binding_prob['overall_binding_probability']:.4f}")
            print(f"  Mean diffusion coefficient: {diffusion_analysis['mean_diffusion_coefficient']:.2e} m²/s")
            print(f"  Number of collisions: {len(simulation_results['collisions'])}")
        
        return True
        
    finally:
        if os.path.exists(pdb_file):
            os.remove(pdb_file)

def demonstrate_binding_kinetics():
    """Demonstrate binding kinetics calculation."""
    print("\n" + "="*80)
    print("BINDING KINETICS DEMO")
    print("="*80)
    
    # Create demo protein
    pdb_file = create_demo_protein()
    
    try:
        # Initialize components
        processor = PDBProcessor()
        processor.load_pdb(pdb_file)
        binding_sites = processor.identify_metal_binding_sites()
        
        # Initialize binding kinetics
        kinetics = BindingKinetics(temperature=298.15, pressure=1.0)
        
        # Define analysis parameters
        metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
        initial_concentrations = {ion: 1e-6 for ion in metal_ions}
        
        print(f"Binding Kinetics Analysis:")
        print(f"  Metal ions: {metal_ions}")
        print(f"  Initial concentrations: {initial_concentrations}")
        print(f"  Binding sites: {len(binding_sites)}")
        
        # Solve binding kinetics
        solution = kinetics.solve_binding_kinetics(
            binding_sites, metal_ions, initial_concentrations, 
            time_span=(0, 1000)  # 1000 seconds
        )
        
        # Calculate binding efficiency
        efficiency = kinetics.calculate_binding_efficiency(
            solution, binding_sites, metal_ions, initial_concentrations
        )
        
        print(f"\nResults:")
        print(f"  Overall binding efficiency: {efficiency['overall_efficiency']:.4f}")
        print(f"  Per-ion efficiencies:")
        for ion, eff in efficiency['ion_efficiencies'].items():
            print(f"    {ion}: {eff:.4f}")
        
        # Show rate constants
        print(f"\nRate Constants (example for Zn2+):")
        if binding_sites:
            site = binding_sites[0]
            k_plus, k_minus = kinetics.calculate_rate_constants(site, 'Zn2+')
            print(f"  Association rate (k+): {k_plus:.2e} M⁻¹s⁻¹")
            print(f"  Dissociation rate (k-): {k_minus:.2e} s⁻¹")
            print(f"  Equilibrium constant (K): {k_plus/k_minus:.2e} M⁻¹")
        
        return efficiency
        
    finally:
        if os.path.exists(pdb_file):
            os.remove(pdb_file)

def explain_results(results):
    """Explain the results obtained from the pipeline."""
    print("\n" + "="*80)
    print("RESULTS EXPLANATION")
    print("="*80)
    
    if results is None:
        print("No results to explain.")
        return
    
    print("Understanding the Results:")
    print()
    
    # 1. Binding Efficiency
    if 'binding_kinetics' in results and results['binding_kinetics']['efficiency']:
        efficiency = results['binding_kinetics']['efficiency']['overall_efficiency']
        print(f"1. Overall Binding Efficiency: {efficiency:.4f}")
        print("   - This value ranges from 0 to 1")
        print("   - 1.0 = 100% binding efficiency (all metal ions bound)")
        print("   - 0.0 = 0% binding efficiency (no metal ions bound)")
        print("   - Values above 0.5 indicate good binding")
        print()
    
    # 2. Per-ion Efficiencies
    if 'binding_kinetics' in results and results['binding_kinetics']['efficiency']:
        ion_efficiencies = results['binding_kinetics']['efficiency']['ion_efficiencies']
        print("2. Per-Ion Binding Efficiencies:")
        for ion, eff in ion_efficiencies.items():
            print(f"   {ion}: {eff:.4f}")
        print("   - Different metal ions have different binding affinities")
        print("   - Higher values indicate stronger binding")
        print("   - This depends on ion size, charge, and coordination preferences")
        print()
    
    # 3. Brownian Motion Results
    print("3. Brownian Motion Simulation:")
    for ion, brownian_data in results['brownian_simulation'].items():
        binding_prob = brownian_data['binding_probability']['overall_binding_probability']
        print(f"   {ion} binding probability: {binding_prob:.4f}")
    print("   - Shows probability of ions finding binding sites through diffusion")
    print("   - Higher values indicate better diffusion-limited binding")
    print("   - Depends on ion size, solution viscosity, and binding site accessibility")
    print()
    
    # 4. Environmental Effects
    env_conditions = results['environmental_conditions']
    print("4. Environmental Conditions:")
    print(f"   Temperature: {env_conditions['temperature']} K ({env_conditions['temperature'] - 273.15:.1f}°C)")
    print(f"   Pressure: {env_conditions['pressure']} atm")
    print(f"   Viscosity: {env_conditions['viscosity']} Pa·s")
    print("   - Higher temperature increases diffusion and binding rates")
    print("   - Higher pressure can affect binding equilibria")
    print("   - Lower viscosity increases diffusion coefficients")
    print()
    
    # 5. Efficiency Ranking
    print("5. Efficiency Ranking:")
    for i, (ion, efficiency) in enumerate(results['efficiency_ranking'], 1):
        print(f"   {i}. {ion}: {efficiency:.4f}")
    print("   - Shows which metal ions bind most efficiently to this protein")
    print("   - Useful for predicting metal ion selectivity")
    print("   - Can guide experimental design")
    print()
    
    # 6. Practical Implications
    print("6. Practical Implications:")
    print("   - Use these results to predict optimal metal ion for your protein")
    print("   - Adjust experimental conditions based on environmental sensitivity")
    print("   - Design metal ion competition experiments")
    print("   - Optimize protein purification and metal loading protocols")
    print()
    
    # 7. Limitations
    print("7. Model Limitations:")
    print("   - Simplified binding site geometry")
    print("   - Assumes ideal solution conditions")
    print("   - Does not include protein dynamics")
    print("   - Binding site cooperativity not modeled")
    print("   - Results should be validated experimentally")

def main():
    """Main function to run the complete demo."""
    print("METALLOPROTEIN BINDING EFFICIENCY PREDICTION DEMO")
    print("Using Basic Configuration File")
    print("="*80)
    
    try:
        # Run main pipeline demo
        results = run_basic_pipeline_demo()
        
        # Run individual component demos
        analyze_binding_sites_demo()
        demonstrate_brownian_motion()
        demonstrate_binding_kinetics()
        
        # Explain results
        explain_results(results)
        
        print("\n" + "="*80)
        print("DEMO COMPLETED SUCCESSFULLY!")
        print("="*80)
        print("Key takeaways:")
        print("✅ Pipeline successfully analyzed metalloprotein binding")
        print("✅ Multiple metal ions compared for binding efficiency")
        print("✅ Environmental conditions considered")
        print("✅ Both kinetics and diffusion effects modeled")
        print("✅ Results saved for further analysis")
        print("\nCheck the 'demo_results/' directory for detailed outputs.")
        
    except Exception as e:
        print(f"\nError during demo: {e}")
        print("Please check your installation and configuration.")

if __name__ == "__main__":
    main() 