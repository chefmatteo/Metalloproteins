#!/usr/bin/env python3
"""
Test script for the Metalloprotein Binding Efficiency Prediction Pipeline

This script demonstrates the pipeline functionality using a sample PDB structure.
"""

import os
import sys
import tempfile
import numpy as np

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.pdb_processor import PDBProcessor
from src.binding_kinetics import BindingKinetics
from src.brownian_simulation import BrownianSimulation

def create_sample_pdb():
    """Create a sample PDB file for testing."""
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
ATOM     11  N   CYS A   2      15.000  15.000  15.000  1.00 20.00           N
ATOM     12  CA  CYS A   2      16.000  15.000  15.000  1.00 20.00           C
ATOM     13  C   CYS A   2      16.500  16.000  15.000  1.00 20.00           C
ATOM     14  O   CYS A   2      17.000  16.500  15.000  1.00 20.00           O
ATOM     15  CB  CYS A   2      16.500  15.000  16.000  1.00 20.00           C
ATOM     16  SG  CYS A   2      17.000  15.000  17.000  1.00 20.00           S
ATOM     17  N   ASP A   3      20.000  20.000  20.000  1.00 20.00           N
ATOM     18  CA  ASP A   3      21.000  20.000  20.000  1.00 20.00           C
ATOM     19  C   ASP A   3      21.500  21.000  20.000  1.00 20.00           C
ATOM     20  O   ASP A   3      22.000  21.500  20.000  1.00 20.00           O
ATOM     21  CB  ASP A   3      21.500  20.000  21.000  1.00 20.00           C
ATOM     22  CG  ASP A   3      22.000  20.000  22.000  1.00 20.00           C
ATOM     23  OD1 ASP A   3      22.500  20.000  22.500  1.00 20.00           O
ATOM     24  OD2 ASP A   3      22.000  20.000  23.000  1.00 20.00           O
ATOM     25  N   GLU A   4      25.000  25.000  25.000  1.00 20.00           N
ATOM     26  CA  GLU A   4      26.000  25.000  25.000  1.00 20.00           C
ATOM     27  C   GLU A   4      26.500  26.000  25.000  1.00 20.00           C
ATOM     28  O   GLU A   4      27.000  26.500  25.000  1.00 20.00           O
ATOM     29  CB  GLU A   4      26.500  25.000  26.000  1.00 20.00           C
ATOM     30  CG  GLU A   4      27.000  25.000  27.000  1.00 20.00           C
ATOM     31  CD  GLU A   4      27.500  25.000  28.000  1.00 20.00           C
ATOM     32  OE1 GLU A   4      28.000  25.000  28.500  1.00 20.00           O
ATOM     33  OE2 GLU A   4      27.500  25.000  29.000  1.00 20.00           O
TER
END"""
    
    # Create temporary file
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
    temp_file.write(pdb_content)
    temp_file.close()
    
    return temp_file.name

def test_pdb_processor():
    """Test PDB processor functionality."""
    print("Testing PDB Processor...")
    
    # Create sample PDB
    pdb_file = create_sample_pdb()
    
    try:
        # Initialize processor
        processor = PDBProcessor()
        
        # Load PDB
        success = processor.load_pdb(pdb_file)
        assert success, "Failed to load PDB file"
        
        # Identify binding sites
        binding_sites = processor.identify_metal_binding_sites()
        print(f"   Found {len(binding_sites)} binding sites")
        
        # Get characteristics
        characteristics = processor.get_binding_site_characteristics()
        print(f"   Calculated characteristics for {len(characteristics)} sites")
        
        # Export results
        processor.export_binding_sites("test_binding_sites.txt")
        
        print("   PDB Processor test passed!")
        return binding_sites
        
    finally:
        # Clean up
        os.unlink(pdb_file)
        if os.path.exists("test_binding_sites.txt"):
            os.unlink("test_binding_sites.txt")

def test_binding_kinetics(binding_sites):
    """Test binding kinetics functionality."""
    print("\nTesting Binding Kinetics...")
    
    # Initialize kinetics solver
    kinetics = BindingKinetics(temperature=298.15, pressure=1.0)
    
    # Test metal ions
    metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
    initial_concentrations = {ion: 1e-6 for ion in metal_ions}
    
    # Calculate rate constants
    for site in binding_sites[:2]:  # Test first 2 sites
        for ion in metal_ions:
            k_plus, k_minus = kinetics.calculate_rate_constants(site, ion)
            print(f"   {ion} at site: k+ = {k_plus:.2e}, k- = {k_minus:.2e}")
    
    # Solve ODE system
    solution = kinetics.solve_binding_kinetics(
        binding_sites=binding_sites,
        metal_ions=metal_ions,
        initial_concentrations=initial_concentrations,
        time_span=(0, 100)
    )
    
    # Calculate efficiency
    efficiency = kinetics.calculate_binding_efficiency(
        solution, binding_sites, metal_ions, initial_concentrations
    )
    
    print(f"   Overall binding efficiency: {efficiency['overall_efficiency']:.4f}")
    print("   Binding Kinetics test passed!")

def test_brownian_simulation(binding_sites):
    """Test Brownian motion simulation."""
    print("\nTesting Brownian Motion Simulation...")
    
    # Initialize simulation
    brownian_sim = BrownianSimulation(
        temperature=298.15,
        viscosity=1e-3,
        box_size=50e-9
    )
    
    # Ion properties
    ion_properties = {
        'radius': 0.74e-10,
        'mass': 65.38e-27,
        'charge': 2,
        'diffusion_coeff': 7.0e-10
    }
    
    # Simulate single ion
    initial_position = np.array([25e-9, 25e-9, 25e-9])
    trajectory = brownian_sim.simulate_ion_trajectory(
        ion_properties=ion_properties,
        initial_position=initial_position,
        binding_sites=binding_sites,
        simulation_time=1e-10  # Short simulation for testing
    )
    
    print(f"   Simulated trajectory with {len(trajectory['collisions'])} collisions")
    
    # Simulate multiple ions
    simulation_results = brownian_sim.simulate_multiple_ions(
        ion_properties=ion_properties,
        n_ions=10,  # Small number for testing
        binding_sites=binding_sites,
        simulation_time=1e-10
    )
    
    # Calculate binding probability
    binding_prob = brownian_sim.calculate_binding_probability(
        simulation_results, binding_sites
    )
    
    print(f"   Overall binding probability: {binding_prob['overall_binding_probability']:.4f}")
    
    # Diffusion analysis
    diffusion_analysis = brownian_sim.calculate_diffusion_analysis(simulation_results)
    print(f"   Mean diffusion coefficient: {diffusion_analysis['mean_diffusion_coefficient']:.2e} m²/s")
    
    print("   Brownian Motion Simulation test passed!")

def test_integration():
    """Test integration of all components."""
    print("\nTesting Full Integration...")
    
    # Create sample PDB and get binding sites
    pdb_file = create_sample_pdb()
    
    try:
        # Initialize all components
        processor = PDBProcessor()
        processor.load_pdb(pdb_file)
        binding_sites = processor.identify_metal_binding_sites()
        
        kinetics = BindingKinetics(temperature=298.15, pressure=1.0)
        brownian_sim = BrownianSimulation(temperature=298.15)
        
        # Test parameters
        metal_ions = ['Zn2+', 'Cu2+']
        initial_concentrations = {ion: 1e-6 for ion in metal_ions}
        
        # Run kinetics
        solution = kinetics.solve_binding_kinetics(
            binding_sites, metal_ions, initial_concentrations, time_span=(0, 50)
        )
        
        # Run Brownian simulation
        ion_properties = kinetics.metal_properties['Zn2+']
        simulation_results = brownian_sim.simulate_multiple_ions(
            ion_properties, 5, binding_sites, simulation_time=1e-10
        )
        
        # Calculate efficiencies
        kinetics_efficiency = kinetics.calculate_binding_efficiency(
            solution, binding_sites, metal_ions, initial_concentrations
        )
        
        brownian_prob = brownian_sim.calculate_binding_probability(
            simulation_results, binding_sites
        )
        
        print(f"   Kinetics efficiency: {kinetics_efficiency['overall_efficiency']:.4f}")
        print(f"   Brownian probability: {brownian_prob['overall_binding_probability']:.4f}")
        
        print("   Integration test passed!")
        
    finally:
        os.unlink(pdb_file)

def main():
    """Run all tests."""
    print("Metalloprotein Binding Efficiency Prediction - Test Suite")
    print("=" * 60)
    
    try:
        # Test individual components
        binding_sites = test_pdb_processor()
        test_binding_kinetics(binding_sites)
        test_brownian_simulation(binding_sites)
        test_integration()
        
        print("\n" + "=" * 60)
        print("All tests passed successfully!")
        print("The pipeline is ready for use.")
        
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 