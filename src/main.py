#!/usr/bin/env python3
"""
Main Pipeline for Metalloprotein Binding Efficiency Prediction

This script integrates all components of the binding efficiency prediction model:
1. PDB structure processing
2. Metal binding site identification
3. Brownian motion simulation
4. Binding kinetics calculation
5. Efficiency prediction
"""

import argparse
import os
import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src.pdb_processor import PDBProcessor
from src.binding_kinetics import BindingKinetics
from src.brownian_simulation import BrownianSimulation

class MetalloproteinPipeline:
    """Main pipeline for metalloprotein binding efficiency prediction."""
    
    def __init__(self, config_file=None):
        """
        Initialize the pipeline.
        
        Parameters:
        -----------
        config_file : str, optional
            Path to configuration file
        """
        self.config = self._load_config(config_file)
        self.pdb_processor = PDBProcessor()
        self.binding_kinetics = BindingKinetics(
            temperature=self.config.get('temperature', 298.15),
            pressure=self.config.get('pressure', 1.0)
        )
        self.brownian_sim = BrownianSimulation(
            temperature=self.config.get('temperature', 298.15),
            viscosity=self.config.get('viscosity', 1e-3),
            box_size=self.config.get('box_size', 100e-9)
        )
        
    def _load_config(self, config_file):
        """Load configuration from file or use defaults."""
        default_config = {
            'temperature': 298.15,  # K
            'pressure': 1.0,        # atm
            'viscosity': 1e-3,      # Pa·s
            'box_size': 100e-9,     # m
            'simulation_time': 1e-9, # s
            'n_ions': 100,
            'metal_ions': ['Zn2+', 'Cu2+', 'Fe2+', 'Mg2+'],
            'initial_concentrations': {
                'Zn2+': 1e-6,
                'Cu2+': 1e-6,
                'Fe2+': 1e-6,
                'Mg2+': 1e-6
            }
        }
        
        if config_file and os.path.exists(config_file):
            with open(config_file, 'r') as f:
                user_config = yaml.safe_load(f)
                default_config.update(user_config)
        
        return default_config
    
    def run_pipeline(self, pdb_file, output_dir='results'):
        """
        Run the complete binding efficiency prediction pipeline.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        output_dir : str
            Output directory for results
        
        Returns:
        --------
        dict : Complete analysis results
        """
        print("Starting Metalloprotein Binding Efficiency Prediction Pipeline")
        print("=" * 60)
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Step 1: Process PDB structure
        print("\n1. Processing PDB structure...")
        if not self.pdb_processor.load_pdb(pdb_file):
            raise ValueError(f"Failed to load PDB file: {pdb_file}")
        
        # Step 2: Identify binding sites
        print("2. Identifying metal binding sites...")
        binding_sites = self.pdb_processor.identify_metal_binding_sites()
        print(f"   Found {len(binding_sites)} potential binding sites")
        
        # Export binding site information
        binding_sites_file = os.path.join(output_dir, 'binding_sites.txt')
        self.pdb_processor.export_binding_sites(binding_sites_file)
        print(f"   Binding sites exported to: {binding_sites_file}")
        
        # Step 3: Calculate binding site characteristics
        print("3. Calculating binding site characteristics...")
        characteristics = self.pdb_processor.get_binding_site_characteristics()
        
        # Step 4: Brownian motion simulation
        print("4. Running Brownian motion simulation...")
        brownian_results = self._run_brownian_simulation(binding_sites)
        
        # Step 5: Binding kinetics calculation
        print("5. Calculating binding kinetics...")
        kinetics_results = self._run_binding_kinetics(binding_sites)
        
        # Step 6: Generate comprehensive results
        print("6. Generating comprehensive results...")
        results = self._compile_results(
            binding_sites, characteristics, brownian_results, kinetics_results
        )
        
        # Step 7: Generate plots and reports
        print("7. Generating plots and reports...")
        self._generate_outputs(results, output_dir)
        
        print("\nPipeline completed successfully!")
        print(f"Results saved to: {output_dir}")
        
        return results
    
    def _run_brownian_simulation(self, binding_sites):
        """Run Brownian motion simulation for all metal ions."""
        brownian_results = {}
        
        for metal_ion in self.config['metal_ions']:
            print(f"   Simulating {metal_ion}...")
            
            # Get ion properties
            ion_properties = self.binding_kinetics.metal_properties.get(metal_ion, {})
            if not ion_properties:
                print(f"   Warning: No properties found for {metal_ion}, using defaults")
                ion_properties = {
                    'radius': 0.75e-10,
                    'mass': 60e-27,
                    'charge': 2,
                    'diffusion_coeff': 7.0e-10
                }
            
            # Run simulation
            simulation_results = self.brownian_sim.simulate_multiple_ions(
                ion_properties=ion_properties,
                n_ions=self.config['n_ions'],
                binding_sites=binding_sites,
                simulation_time=self.config['simulation_time']
            )
            
            # Calculate binding probabilities
            binding_prob = self.brownian_sim.calculate_binding_probability(
                simulation_results, binding_sites
            )
            
            # Calculate diffusion analysis
            diffusion_analysis = self.brownian_sim.calculate_diffusion_analysis(
                simulation_results
            )
            
            brownian_results[metal_ion] = {
                'simulation_results': simulation_results,
                'binding_probability': binding_prob,
                'diffusion_analysis': diffusion_analysis
            }
        
        return brownian_results
    
    def _run_binding_kinetics(self, binding_sites):
        """Run binding kinetics calculation."""
        print("   Solving ODE system for binding kinetics...")
        
        # Solve ODE system
        solution = self.binding_kinetics.solve_binding_kinetics(
            binding_sites=binding_sites,
            metal_ions=self.config['metal_ions'],
            initial_concentrations=self.config['initial_concentrations'],
            time_span=(0, 1000)  # 1000 seconds
        )
        
        # Calculate binding efficiency
        efficiency = self.binding_kinetics.calculate_binding_efficiency(
            solution, binding_sites, self.config['metal_ions'], 
            self.config['initial_concentrations']
        )
        
        return {
            'solution': solution,
            'efficiency': efficiency
        }
    
    def _compile_results(self, binding_sites, characteristics, brownian_results, kinetics_results):
        """Compile all results into a comprehensive summary."""
        results = {
            'protein_info': {
                'n_binding_sites': len(binding_sites),
                'binding_site_characteristics': characteristics
            },
            'environmental_conditions': {
                'temperature': self.config['temperature'],
                'pressure': self.config['pressure'],
                'viscosity': self.config['viscosity']
            },
            'brownian_simulation': brownian_results,
            'binding_kinetics': kinetics_results,
            'metal_ions': self.config['metal_ions'],
            'initial_concentrations': self.config['initial_concentrations']
        }
        
        # Calculate overall binding efficiency ranking
        ion_efficiencies = kinetics_results['efficiency']['ion_efficiencies']
        efficiency_ranking = sorted(
            ion_efficiencies.items(), 
            key=lambda x: x[1], 
            reverse=True
        )
        
        results['efficiency_ranking'] = efficiency_ranking
        
        return results
    
    def _generate_outputs(self, results, output_dir):
        """Generate plots and reports."""
        # Plot binding kinetics
        kinetics_plot = self.binding_kinetics.plot_binding_kinetics(
            results['binding_kinetics']['solution'],
            results['protein_info']['n_binding_sites'],
            results['metal_ions'],
            save_path=os.path.join(output_dir, 'binding_kinetics.png')
        )
        
        # Plot Brownian motion trajectories (for first metal ion)
        first_ion = results['metal_ions'][0]
        brownian_plot = self.brownian_sim.plot_trajectories(
            results['brownian_simulation'][first_ion]['simulation_results'],
            results['protein_info']['n_binding_sites'],
            save_path=os.path.join(output_dir, 'brownian_trajectories.png')
        )
        
        # Generate summary report
        self._generate_summary_report(results, output_dir)
        
        # Save results as JSON
        import json
        results_file = os.path.join(output_dir, 'results.json')
        
        # Convert numpy arrays to lists for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, dict):
                return {key: convert_numpy(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(item) for item in obj]
            else:
                return obj
        
        with open(results_file, 'w') as f:
            json.dump(convert_numpy(results), f, indent=2)
        
        print(f"   Results saved to: {results_file}")
    
    def _generate_summary_report(self, results, output_dir):
        """Generate a human-readable summary report."""
        report_file = os.path.join(output_dir, 'summary_report.txt')
        
        with open(report_file, 'w') as f:
            f.write("METALLOPROTEIN BINDING EFFICIENCY PREDICTION REPORT\n")
            f.write("=" * 60 + "\n\n")
            
            # Protein information
            f.write("PROTEIN ANALYSIS:\n")
            f.write(f"Number of binding sites identified: {results['protein_info']['n_binding_sites']}\n")
            f.write(f"Temperature: {results['environmental_conditions']['temperature']:.2f} K\n")
            f.write(f"Pressure: {results['environmental_conditions']['pressure']:.2f} atm\n\n")
            
            # Binding efficiency ranking
            f.write("BINDING EFFICIENCY RANKING:\n")
            for i, (ion, efficiency) in enumerate(results['efficiency_ranking'], 1):
                f.write(f"{i}. {ion}: {efficiency:.4f}\n")
            f.write("\n")
            
            # Brownian motion results
            f.write("BROWNIAN MOTION SIMULATION RESULTS:\n")
            for ion in results['metal_ions']:
                brownian_data = results['brownian_simulation'][ion]
                binding_prob = brownian_data['binding_probability']
                diffusion_data = brownian_data['diffusion_analysis']
                
                f.write(f"\n{ion}:\n")
                f.write(f"  Overall binding probability: {binding_prob['overall_binding_probability']:.4f}\n")
                f.write(f"  Mean diffusion coefficient: {diffusion_data['mean_diffusion_coefficient']:.2e} m²/s\n")
                f.write(f"  Total collisions: {binding_prob['total_collisions']}\n")
            
            f.write("\n" + "=" * 60 + "\n")
            f.write("Report generated by Metalloprotein Binding Efficiency Prediction Pipeline\n")
        
        print(f"   Summary report saved to: {report_file}")

def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Metalloprotein Binding Efficiency Prediction Pipeline"
    )
    parser.add_argument(
        'pdb_file', 
        help='Path to PDB file'
    )
    parser.add_argument(
        '--output-dir', 
        default='results',
        help='Output directory for results (default: results)'
    )
    parser.add_argument(
        '--config', 
        help='Path to configuration file (YAML format)'
    )
    parser.add_argument(
        '--temperature', 
        type=float, 
        default=298.15,
        help='Temperature in Kelvin (default: 298.15)'
    )
    parser.add_argument(
        '--pressure', 
        type=float, 
        default=1.0,
        help='Pressure in atm (default: 1.0)'
    )
    
    args = parser.parse_args()
    
    # Check if PDB file exists
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file not found: {args.pdb_file}")
        sys.exit(1)
    
    try:
        # Initialize pipeline
        pipeline = MetalloproteinPipeline(config_file=args.config)
        
        # Override config with command-line arguments
        if args.temperature != 298.15:
            pipeline.config['temperature'] = args.temperature
        if args.pressure != 1.0:
            pipeline.config['pressure'] = args.pressure
        
        # Run pipeline
        results = pipeline.run_pipeline(args.pdb_file, args.output_dir)
        
        # Print summary
        print("\nSUMMARY:")
        print(f"Overall binding efficiency ranking:")
        for i, (ion, efficiency) in enumerate(results['efficiency_ranking'], 1):
            print(f"  {i}. {ion}: {efficiency:.4f}")
        
    except Exception as e:
        print(f"Error running pipeline: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 