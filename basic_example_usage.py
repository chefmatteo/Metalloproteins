#!/usr/bin/env python3
"""
Basic Example Usage for Metalloprotein Binding Efficiency Prediction Pipeline

This script demonstrates the basic pipeline with fundamental ODE solving and
binding site identification capabilities.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

# Import basic modules
from src.basic.main import MetalloproteinPipeline
from src.basic.binding_kinetics import BindingKinetics
from src.pdb_processor import PDBProcessor
from src.brownian_simulation import BrownianSimulation

def create_example_pdb():
    """Create an example PDB file with multiple metal binding sites."""
    pdb_content = """ATOM      1  N   HIS A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  HIS A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   HIS A   1      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  O   HIS A   1      12.000  11.500  10.000  1.00 20.00           O
ATOM      5  CB  HIS A   1      11.500  10.000  11.000  1.00 20.00           C
ATOM      6  CG  HIS A   1      12.000  10.000  12.000  1.00 20.00           C
ATOM      7  ND1 HIS A   1      12.500  10.000  12.500  1.00 20.00           N
ATOM      8  CD2 HIS A   1      12.000  10.000  13.000  1.00 20.00           C
ATOM      9  CE1 HIS A   1      13.000  10.000  13.000  1.00 20.00           C
ATOM     10  NE2 HIS A   1      12.500  10.000  13.500  1.00 20.00           N
ATOM     11  N   CYS A   2      12.500  10.500  12.500  1.00 20.00           N
ATOM     12  CA  CYS A   2      13.000  10.500  12.500  1.00 20.00           C
ATOM     13  C   CYS A   2      13.500  11.500  12.500  1.00 20.00           C
ATOM     14  O   CYS A   2      14.000  12.000  12.500  1.00 20.00           O
ATOM     15  CB  CYS A   2      13.500  10.500  13.500  1.00 20.00           C
ATOM     16  SG  CYS A   2      14.000  10.500  14.500  1.00 20.00           S
ATOM     17  N   ASP A   3      15.000  15.000  15.000  1.00 20.00           N
ATOM     18  CA  ASP A   3      16.000  15.000  15.000  1.00 20.00           C
ATOM     19  C   ASP A   3      16.500  16.000  15.000  1.00 20.00           C
ATOM     20  O   ASP A   3      17.000  16.500  15.000  1.00 20.00           O
ATOM     21  CB  ASP A   3      16.500  15.000  16.000  1.00 20.00           C
ATOM     22  CG  ASP A   3      17.000  15.000  17.000  1.00 20.00           C
ATOM     23  OD1 ASP A   3      17.500  15.000  17.500  1.00 20.00           O
ATOM     24  OD2 ASP A   3      17.000  15.000  18.000  1.00 20.00           O
ATOM     25  N   GLU A   4      17.500  15.500  17.500  1.00 20.00           N
ATOM     26  CA  GLU A   4      18.000  15.500  17.500  1.00 20.00           C
ATOM     27  C   GLU A   4      18.500  16.500  17.500  1.00 20.00           C
ATOM     28  O   GLU A   4      19.000  17.000  17.500  1.00 20.00           O
ATOM     29  CB  GLU A   4      18.500  15.500  18.500  1.00 20.00           C
ATOM     30  CG  GLU A   4      19.000  15.500  19.500  1.00 20.00           C
ATOM     31  CD  GLU A   4      19.500  15.500  20.500  1.00 20.00           C
ATOM     32  OE1 GLU A   4      20.000  15.500  21.000  1.00 20.00           O
ATOM     33  OE2 GLU A   4      19.500  15.500  21.500  1.00 20.00           O
ATOM     34  N   HIS A   5      20.000  20.000  20.000  1.00 20.00           N
ATOM     35  CA  HIS A   5      21.000  20.000  20.000  1.00 20.00           C
ATOM     36  C   HIS A   5      21.500  21.000  20.000  1.00 20.00           C
ATOM     37  O   HIS A   5      22.000  21.500  20.000  1.00 20.00           O
ATOM     38  CB  HIS A   5      21.500  20.000  21.000  1.00 20.00           C
ATOM     39  CG  HIS A   5      22.000  20.000  22.000  1.00 20.00           C
ATOM     40  ND1 HIS A   5      22.500  20.000  22.500  1.00 20.00           N
ATOM     41  CD2 HIS A   5      22.000  20.000  23.000  1.00 20.00           C
ATOM     42  CE1 HIS A   5      23.000  20.000  23.000  1.00 20.00           C
ATOM     43  NE2 HIS A   5      22.500  20.000  23.500  1.00 20.00           N
ATOM     44  N   CYS A   6      22.500  20.500  23.000  1.00 20.00           N
ATOM     45  CA  CYS A   6      23.000  20.500  23.000  1.00 20.00           C
ATOM     46  C   CYS A   6      23.500  21.500  23.000  1.00 20.00           C
ATOM     47  O   CYS A   6      24.000  22.000  23.000  1.00 20.00           O
ATOM     48  CB  CYS A   6      23.500  20.500  24.000  1.00 20.00           C
ATOM     49  SG  CYS A   6      24.000  20.500  25.000  1.00 20.00           S
TER
END"""
    
    with open('example_protein.pdb', 'w') as f:
        f.write(pdb_content)
    
    return 'example_protein.pdb'

def example_1_basic_analysis():
    """Example 1: Basic binding site analysis."""
    print("Example 1: Basic Binding Site Analysis")
    print("-" * 40)
    
    # Create example PDB
    pdb_file = create_example_pdb()
    
    try:
        # Initialize PDB processor
        processor = PDBProcessor()
        
        # Load and analyze structure
        processor.load_pdb(pdb_file)
        binding_sites = processor.identify_metal_binding_sites()
        
        print(f"Found {len(binding_sites)} potential metal binding sites:")
        for i, site in enumerate(binding_sites):
            residue = site['residue']
            print(f"  Site {i+1}: {residue.get_resname()} {residue.get_id()[1]}")
            print(f"    Coordinating atoms: {[atom.get_name() for atom in site['coordinating_atoms']]}")
            print(f"    Center: {site['center']}")
            print(f"    Coordination number: {site['coordination_number']}")
        
        # Get binding site characteristics
        characteristics = processor.get_binding_site_characteristics()
        print(f"\nBinding site characteristics:")
        for i, char in enumerate(characteristics):
            print(f"  Site {i+1}:")
            print(f"    Volume: {char['volume']:.2e} m³")
            print(f"    Electrostatic potential: {char['electrostatic_potential']:.2e} V")
            print(f"    Accessible: {char['accessibility']}")
            print(f"    Geometry: {char['coordination_geometry']}")
        
    finally:
        # Clean up
        if os.path.exists(pdb_file):
            os.unlink(pdb_file)

def example_2_binding_kinetics():
    """Example 2: Binding kinetics analysis."""
    print("\nExample 2: Binding Kinetics Analysis")
    print("-" * 40)
    
    # Create example PDB
    pdb_file = create_example_pdb()
    
    try:
        # Initialize components
        processor = PDBProcessor()
        processor.load_pdb(pdb_file)
        binding_sites = processor.identify_metal_binding_sites()
        
        kinetics = BindingKinetics(temperature=298.15, pressure=1.0)
        
        # Define metal ions and concentrations
        metal_ions = ['Zn2+', 'Cu2+', 'Fe2+', 'Mg2+']
        initial_concentrations = {
            'Zn2+': 1e-6,
            'Cu2+': 1e-6,
            'Fe2+': 1e-6,
            'Mg2+': 1e-6
        }
        
        # Calculate rate constants for each ion-site combination
        print("Rate constants for different ion-site combinations:")
        for i, site in enumerate(binding_sites):
            print(f"\nBinding Site {i+1}:")
            for ion in metal_ions:
                k_plus, k_minus = kinetics.calculate_rate_constants(site, ion)
                print(f"  {ion}: k+ = {k_plus:.2e} M⁻¹s⁻¹, k- = {k_minus:.2e} s⁻¹")
        
        # Solve ODE system
        print(f"\nSolving binding kinetics ODE system...")
        solution = kinetics.solve_binding_kinetics(
            binding_sites=binding_sites,
            metal_ions=metal_ions,
            initial_concentrations=initial_concentrations,
            time_span=(0, 500)
        )
        
        # Calculate binding efficiency
        efficiency = kinetics.calculate_binding_efficiency(
            solution, binding_sites, metal_ions, initial_concentrations
        )
        
        print(f"\nBinding Efficiency Results:")
        print(f"Overall efficiency: {efficiency['overall_efficiency']:.4f}")
        print(f"\nPer-ion efficiencies:")
        for ion, eff in efficiency['ion_efficiencies'].items():
            print(f"  {ion}: {eff:.4f}")
        
        print(f"\nPer-site occupancies:")
        for site_name, site_data in efficiency['site_efficiencies'].items():
            print(f"  {site_name}: {site_data['total_occupancy']:.4f}")
        
        # Plot results
        kinetics.plot_binding_kinetics(solution, binding_sites, metal_ions)
        
    finally:
        if os.path.exists(pdb_file):
            os.unlink(pdb_file)

def example_3_brownian_simulation():
    """Example 3: Brownian motion simulation."""
    print("\nExample 3: Brownian Motion Simulation")
    print("-" * 40)
    
    # Create example PDB
    pdb_file = create_example_pdb()
    
    try:
        # Initialize components
        processor = PDBProcessor()
        processor.load_pdb(pdb_file)
        binding_sites = processor.identify_metal_binding_sites()
        
        brownian_sim = BrownianSimulation(
            temperature=298.15,
            viscosity=1e-3,
            box_size=100e-9
        )
        
        # Simulate different metal ions
        metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
        
        for ion in metal_ions:
            print(f"\nSimulating {ion}...")
            
            # Get ion properties
            ion_properties = {
                'Zn2+': {'radius': 0.74e-10, 'mass': 65.38e-27, 'charge': 2},
                'Cu2+': {'radius': 0.73e-10, 'mass': 63.55e-27, 'charge': 2},
                'Fe2+': {'radius': 0.78e-10, 'mass': 55.85e-27, 'charge': 2}
            }[ion]
            
            # Run simulation
            simulation_results = brownian_sim.simulate_multiple_ions(
                ion_properties=ion_properties,
                n_ions=50,
                binding_sites=binding_sites,
                simulation_time=1e-9
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
            print(f"  Total collisions: {binding_prob['total_collisions']}")
        
        # Plot trajectories for Zn2+
        print(f"\nGenerating trajectory plots for Zn2+...")
        ion_properties = {'radius': 0.74e-10, 'mass': 65.38e-27, 'charge': 2}
        simulation_results = brownian_sim.simulate_multiple_ions(
            ion_properties, 20, binding_sites, simulation_time=1e-9
        )
        brownian_sim.plot_trajectories(simulation_results, binding_sites)
        
    finally:
        if os.path.exists(pdb_file):
            os.unlink(pdb_file)

def example_4_temperature_dependence():
    """Example 4: Temperature dependence analysis."""
    print("\nExample 4: Temperature Dependence Analysis")
    print("-" * 40)
    
    # Create example PDB
    pdb_file = create_example_pdb()
    
    try:
        # Initialize components
        processor = PDBProcessor()
        processor.load_pdb(pdb_file)
        binding_sites = processor.identify_metal_binding_sites()
        
        # Test different temperatures
        temperatures = [273.15, 298.15, 323.15, 348.15]  # 0°C, 25°C, 50°C, 75°C
        metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
        initial_concentrations = {ion: 1e-6 for ion in metal_ions}
        
        temperature_results = {}
        
        for temp in temperatures:
            print(f"\nAnalyzing at {temp-273.15:.1f}°C ({temp:.1f}K)...")
            
            # Initialize kinetics solver for this temperature
            kinetics = BindingKinetics(temperature=temp, pressure=1.0)
            
            # Solve ODE system
            solution = kinetics.solve_binding_kinetics(
                binding_sites=binding_sites,
                metal_ions=metal_ions,
                initial_concentrations=initial_concentrations,
                time_span=(0, 200)
            )
            
            # Calculate efficiency
            efficiency = kinetics.calculate_binding_efficiency(
                solution, binding_sites, metal_ions, initial_concentrations
            )
            
            temperature_results[temp] = efficiency['ion_efficiencies']
        
        # Plot temperature dependence
        plt.figure(figsize=(12, 8))
        
        for ion in metal_ions:
            temps = list(temperature_results.keys())
            efficiencies = [temperature_results[temp][ion] for temp in temps]
            temps_celsius = [temp - 273.15 for temp in temps]
            
            plt.plot(temps_celsius, efficiencies, 'o-', linewidth=2, markersize=8, label=ion)
        
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Binding Efficiency')
        plt.title('Temperature Dependence of Binding Efficiency')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
        
        # Print summary
        print(f"\nTemperature Dependence Summary:")
        for temp in temperatures:
            print(f"\n{temp-273.15:.1f}°C:")
            for ion in metal_ions:
                eff = temperature_results[temp][ion]
                print(f"  {ion}: {eff:.4f}")
        
    finally:
        if os.path.exists(pdb_file):
            os.unlink(pdb_file)

def example_5_comprehensive_analysis():
    """Example 5: Comprehensive analysis with all components."""
    print("\nExample 5: Comprehensive Analysis")
    print("-" * 40)
    
    # Create example PDB
    pdb_file = create_example_pdb()
    
    try:
        # Initialize all components
        processor = PDBProcessor()
        processor.load_pdb(pdb_file)
        binding_sites = processor.identify_metal_binding_sites()
        
        kinetics = BindingKinetics(temperature=298.15, pressure=1.0)
        brownian_sim = BrownianSimulation(temperature=298.15, viscosity=1e-3)
        
        # Define analysis parameters
        metal_ions = ['Zn2+', 'Cu2+', 'Fe2+', 'Mg2+']
        initial_concentrations = {ion: 1e-6 for ion in metal_ions}
        
        print("Running comprehensive analysis...")
        
        # 1. Binding kinetics analysis
        print("1. Calculating binding kinetics...")
        kinetics_solution = kinetics.solve_binding_kinetics(
            binding_sites, metal_ions, initial_concentrations, time_span=(0, 300)
        )
        kinetics_efficiency = kinetics.calculate_binding_efficiency(
            kinetics_solution, binding_sites, metal_ions, initial_concentrations
        )
        
        # 2. Brownian motion analysis
        print("2. Running Brownian motion simulations...")
        brownian_results = {}
        for ion in metal_ions:
            ion_properties = kinetics.metal_properties.get(ion, {
                'radius': 0.75e-10, 'mass': 60e-27, 'charge': 2
            })
            
            simulation_results = brownian_sim.simulate_multiple_ions(
                ion_properties, 30, binding_sites, simulation_time=1e-9
            )
            
            binding_prob = brownian_sim.calculate_binding_probability(
                simulation_results, binding_sites
            )
            
            brownian_results[ion] = binding_prob
        
        # 3. Compile and compare results
        print("3. Compiling results...")
        print(f"\nComprehensive Binding Analysis Results:")
        print(f"{'Ion':<8} {'Kinetics Eff.':<15} {'Brownian Prob.':<15} {'Combined Score':<15}")
        print("-" * 60)
        
        for ion in metal_ions:
            kinetics_eff = kinetics_efficiency['ion_efficiencies'][ion]
            brownian_prob = brownian_results[ion]['overall_binding_probability']
            
            # Combined score (weighted average)
            combined_score = 0.7 * kinetics_eff + 0.3 * brownian_prob
            
            print(f"{ion:<8} {kinetics_eff:<15.4f} {brownian_prob:<15.4f} {combined_score:<15.4f}")
        
        # 4. Generate comprehensive plots
        print("4. Generating comprehensive plots...")
        
        # Kinetics plot
        kinetics.plot_binding_kinetics(kinetics_solution, binding_sites, metal_ions)
        
        # Brownian trajectories plot (for Zn2+)
        ion_properties = kinetics.metal_properties['Zn2+']
        simulation_results = brownian_sim.simulate_multiple_ions(
            ion_properties, 20, binding_sites, simulation_time=1e-9
        )
        brownian_sim.plot_trajectories(simulation_results, binding_sites)
        
        print("\nComprehensive analysis completed!")
        
    finally:
        if os.path.exists(pdb_file):
            os.unlink(pdb_file)

def main():
    """Run all examples."""
    print("Metalloprotein Binding Efficiency Prediction - Examples")
    print("=" * 60)
    
    try:
        # Run all examples
        example_1_basic_analysis()
        example_2_binding_kinetics()
        example_3_brownian_simulation()
        example_4_temperature_dependence()
        example_5_comprehensive_analysis()
        
        print("\n" + "=" * 60)
        print("All examples completed successfully!")
        print("\nTo use the pipeline with your own PDB file:")
        print("python src/main.py your_protein.pdb --output-dir results")
        
    except Exception as e:
        print(f"\nExample failed with error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 