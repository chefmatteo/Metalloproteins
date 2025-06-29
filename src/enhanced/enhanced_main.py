"""
Enhanced Main Pipeline for Metalloprotein Binding Efficiency Prediction

This module integrates all enhanced components:
- Enhanced binding site identification with multiple algorithms
- Environmental parameter coupling (temperature, pH, pressure, redox)
- Spatial discretization (1000-cube model)
- Advanced visualization with PyMOL and RF Diffusion
- MESPEUS database integration
- CHED network analysis
"""

import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import enhanced modules
from .enhanced_binding_kinetics import EnhancedBindingKinetics
from .enhanced_binding_site_identification import EnhancedBindingSiteIdentifier
from ..pdb_processor import PDBProcessor
from ..brownian_simulation import BrownianSimulation

class EnhancedMetalloproteinPipeline:
    """Enhanced pipeline for metalloprotein binding efficiency prediction."""
    
    def __init__(self, config_path="config/enhanced_config.yaml"):
        """
        Initialize enhanced pipeline.
        
        Parameters:
        -----------
        config_path : str
            Path to configuration file
        """
        self.config = self._load_config(config_path)
        self.enhanced_kinetics = EnhancedBindingKinetics(self.config)
        self.binding_site_identifier = EnhancedBindingSiteIdentifier(self.config)
        self.pdb_processor = PDBProcessor()
        self.brownian_sim = BrownianSimulation()
        
        # Create output directories
        self.output_dir = Path("output/enhanced")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        (self.output_dir / "plots").mkdir(exist_ok=True)
        (self.output_dir / "structures").mkdir(exist_ok=True)
        (self.output_dir / "animations").mkdir(exist_ok=True)
        (self.output_dir / "data").mkdir(exist_ok=True)
        
    def _load_config(self, config_path):
        """Load configuration from YAML file."""
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    
    def run_enhanced_analysis(self, pdb_file, protein_sequence=None, 
                            metal_ions=None, initial_concentrations=None,
                            time_span=(0, 1000), save_results=True):
        """
        Run complete enhanced analysis pipeline.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        protein_sequence : str, optional
            Protein sequence string
        metal_ions : list, optional
            List of metal ion types to analyze
        initial_concentrations : dict, optional
            Initial concentrations of metal ions
        time_span : tuple, optional
            Time span for simulation
        save_results : bool, optional
            Whether to save results to files
        
        Returns:
        --------
        dict : Complete analysis results
        """
        print("="*80)
        print("ENHANCED METALLOPROTEIN BINDING EFFICIENCY PREDICTION PIPELINE")
        print("="*80)
        
        # Step 1: Load and process PDB structure
        print("\n1. Loading and processing PDB structure...")
        protein_structure = self.pdb_processor.load_pdb_structure(pdb_file)
        
        if protein_sequence is None:
            protein_sequence = self.pdb_processor.extract_sequence(protein_structure)
        
        print(f"   Protein sequence length: {len(protein_sequence)}")
        print(f"   Number of chains: {len(list(protein_structure.get_chains()))}")
        
        # Step 2: Enhanced binding site identification
        print("\n2. Enhanced binding site identification...")
        binding_site_results = self.binding_site_identifier.identify_binding_sites(
            protein_structure, protein_sequence
        )
        
        print(f"   Total binding sites identified: {binding_site_results['total_sites']}")
        print(f"   Average consensus score: {binding_site_results['average_consensus_score']:.3f}")
        
        # Step 3: Set up metal ions and concentrations
        if metal_ions is None:
            metal_ions = list(self.config['metal_ions'].keys())
        
        if initial_concentrations is None:
            initial_concentrations = {ion: 1e-6 for ion in metal_ions}
        
        print(f"\n3. Analyzing metal ions: {metal_ions}")
        print(f"   Initial concentrations: {initial_concentrations}")
        
        # Step 4: Enhanced binding kinetics simulation
        print("\n4. Running enhanced binding kinetics simulation...")
        if binding_site_results['binding_sites']:
            binding_sites = binding_site_results['binding_sites']
            kinetics_solution = self.enhanced_kinetics.solve_enhanced_kinetics(
                binding_sites, metal_ions, initial_concentrations, time_span
            )
            
            # Step 5: Calculate enhanced binding efficiency
            print("\n5. Calculating enhanced binding efficiency...")
            efficiency_results = self.enhanced_kinetics.calculate_enhanced_binding_efficiency(
                kinetics_solution, binding_sites, metal_ions, initial_concentrations
            )
            
            # Step 6: Brownian motion simulation
            print("\n6. Running Brownian motion simulation...")
            brownian_results = self.brownian_sim.simulate_brownian_motion(
                binding_sites, metal_ions, initial_concentrations,
                time_steps=1000, diffusion_coeffs=None
            )
            
            # Step 7: Generate visualizations
            print("\n7. Generating enhanced visualizations...")
            self._generate_enhanced_visualizations(
                binding_site_results, kinetics_solution, efficiency_results,
                brownian_results, protein_structure, protein_sequence
            )
            
            # Step 8: Save results
            if save_results:
                print("\n8. Saving results...")
                self._save_enhanced_results(
                    binding_site_results, kinetics_solution, efficiency_results,
                    brownian_results, protein_sequence
                )
            
            # Step 9: Generate comprehensive report
            print("\n9. Generating comprehensive report...")
            report = self._generate_comprehensive_report(
                binding_site_results, efficiency_results, brownian_results,
                protein_sequence, metal_ions, initial_concentrations
            )
            
            return {
                'binding_sites': binding_site_results,
                'kinetics_solution': kinetics_solution,
                'efficiency_results': efficiency_results,
                'brownian_results': brownian_results,
                'report': report
            }
        
        else:
            print("   No binding sites identified. Skipping kinetics simulation.")
            return {
                'binding_sites': binding_site_results,
                'kinetics_solution': None,
                'efficiency_results': None,
                'brownian_results': None,
                'report': None
            }
    
    def _generate_enhanced_visualizations(self, binding_site_results, kinetics_solution,
                                        efficiency_results, brownian_results,
                                        protein_structure, protein_sequence):
        """Generate enhanced visualizations."""
        
        # 1. Binding site identification plots
        print("   - Generating binding site identification plots...")
        binding_site_fig = self.binding_site_identifier.plot_binding_site_results(
            binding_site_results, 
            self.binding_site_identifier._extract_protein_info(protein_structure, protein_sequence),
            save_path=self.output_dir / "plots" / "binding_site_identification.png"
        )
        
        # 2. Enhanced kinetics plots
        if kinetics_solution is not None:
            print("   - Generating enhanced kinetics plots...")
            kinetics_fig = self.enhanced_kinetics.plot_enhanced_results(
                kinetics_solution,
                binding_site_results['binding_sites'],
                list(self.config['metal_ions'].keys()),
                save_path=self.output_dir / "plots" / "enhanced_kinetics.png"
            )
        
        # 3. Environmental parameter analysis
        if efficiency_results is not None:
            print("   - Generating environmental parameter analysis...")
            self._plot_environmental_analysis(efficiency_results)
        
        # 4. Spatial distribution plots
        if efficiency_results is not None:
            print("   - Generating spatial distribution plots...")
            self._plot_spatial_distributions(efficiency_results)
        
        # 5. Brownian motion visualization
        if brownian_results is not None:
            print("   - Generating Brownian motion visualization...")
            self._plot_brownian_motion(brownian_results)
        
        # 6. PyMOL integration (simulated)
        if self.config['visualization']['pymol']['enabled']:
            print("   - Generating PyMOL visualization scripts...")
            self._generate_pymol_scripts(binding_site_results, efficiency_results)
        
        # 7. RF Diffusion integration (simulated)
        if self.config['visualization']['rf_diffusion']['enabled']:
            print("   - Generating RF Diffusion visualization...")
            self._generate_rf_diffusion_visualization(efficiency_results)
    
    def _plot_environmental_analysis(self, efficiency_results):
        """Plot environmental parameter analysis."""
        env_data = efficiency_results['environmental_analysis']
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Temperature analysis
        ax1 = axes[0, 0]
        temp_data = env_data['temperature']
        ax1.bar(['Mean', 'Std', 'Min', 'Max'], 
               [temp_data['mean']-273.15, temp_data['std'], 
                temp_data['range'][0]-273.15, temp_data['range'][1]-273.15],
               color='red', alpha=0.7)
        ax1.set_ylabel('Temperature (°C)')
        ax1.set_title('Temperature Statistics')
        ax1.grid(True, alpha=0.3)
        
        # pH analysis
        ax2 = axes[0, 1]
        ph_data = env_data['pH']
        ax2.bar(['Mean', 'Std', 'Min', 'Max'], 
               [ph_data['mean'], ph_data['std'], 
                ph_data['range'][0], ph_data['range'][1]],
               color='blue', alpha=0.7)
        ax2.set_ylabel('pH')
        ax2.set_title('pH Statistics')
        ax2.grid(True, alpha=0.3)
        
        # Pressure analysis
        ax3 = axes[1, 0]
        pressure_data = env_data['pressure']
        ax3.bar(['Mean', 'Std', 'Min', 'Max'], 
               [pressure_data['mean'], pressure_data['std'], 
                pressure_data['range'][0], pressure_data['range'][1]],
               color='green', alpha=0.7)
        ax3.set_ylabel('Pressure (atm)')
        ax3.set_title('Pressure Statistics')
        ax3.grid(True, alpha=0.3)
        
        # Redox potential analysis
        ax4 = axes[1, 1]
        redox_data = env_data['redox_potential']
        ax4.bar(['Mean', 'Std', 'Min', 'Max'], 
               [redox_data['mean'], redox_data['std'], 
                redox_data['range'][0], redox_data['range'][1]],
               color='purple', alpha=0.7)
        ax4.set_ylabel('Redox Potential (V)')
        ax4.set_title('Redox Potential Statistics')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / "plots" / "environmental_analysis.png", 
                   dpi=300, bbox_inches='tight')
        plt.show()
    
    def _plot_spatial_distributions(self, efficiency_results):
        """Plot spatial distributions of environmental parameters."""
        spatial_data = efficiency_results['spatial_data']
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Temperature distribution
        ax1 = axes[0, 0]
        temp_dist = spatial_data['temperature'].reshape(10, 10, 10)
        temp_slice = temp_dist[:, :, 5]  # Middle slice
        im1 = ax1.imshow(temp_slice - 273.15, cmap='hot', aspect='auto')
        ax1.set_title('Temperature Distribution (°C)')
        plt.colorbar(im1, ax=ax1)
        
        # pH distribution
        ax2 = axes[0, 1]
        ph_dist = spatial_data['pH'].reshape(10, 10, 10)
        ph_slice = ph_dist[:, :, 5]  # Middle slice
        im2 = ax2.imshow(ph_slice, cmap='RdYlBu_r', aspect='auto')
        ax2.set_title('pH Distribution')
        plt.colorbar(im2, ax=ax2)
        
        # Pressure distribution
        ax3 = axes[1, 0]
        pressure_dist = spatial_data['pressure'].reshape(10, 10, 10)
        pressure_slice = pressure_dist[:, :, 5]  # Middle slice
        im3 = ax3.imshow(pressure_slice, cmap='Blues', aspect='auto')
        ax3.set_title('Pressure Distribution (atm)')
        plt.colorbar(im3, ax=ax3)
        
        # Redox potential distribution
        ax4 = axes[1, 1]
        redox_dist = spatial_data['redox_potential'].reshape(10, 10, 10)
        redox_slice = redox_dist[:, :, 5]  # Middle slice
        im4 = ax4.imshow(redox_slice, cmap='RdBu', aspect='auto')
        ax4.set_title('Redox Potential Distribution (V)')
        plt.colorbar(im4, ax=ax4)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / "plots" / "spatial_distributions.png", 
                   dpi=300, bbox_inches='tight')
        plt.show()
    
    def _plot_brownian_motion(self, brownian_results):
        """Plot Brownian motion results."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Ion trajectories
        ax1 = axes[0, 0]
        for ion_type, trajectories in brownian_results['trajectories'].items():
            for traj in trajectories[:5]:  # Plot first 5 trajectories
                ax1.plot(traj[:, 0], traj[:, 1], alpha=0.6, label=ion_type)
        ax1.set_xlabel('X Position (Å)')
        ax1.set_ylabel('Y Position (Å)')
        ax1.set_title('Ion Trajectories (Top View)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Mean squared displacement
        ax2 = axes[0, 1]
        for ion_type, msd in brownian_results['msd'].items():
            ax2.plot(brownian_results['time_points'], msd, label=ion_type, linewidth=2)
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('MSD (Å²)')
        ax2.set_title('Mean Squared Displacement')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Binding probability
        ax3 = axes[1, 0]
        for ion_type, prob in brownian_results['binding_probability'].items():
            ax3.plot(brownian_results['time_points'], prob, label=ion_type, linewidth=2)
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Binding Probability')
        ax3.set_title('Binding Probability vs Time')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Diffusion coefficient analysis
        ax4 = axes[1, 1]
        ion_types = list(brownian_results['diffusion_coefficients'].keys())
        diffusion_values = list(brownian_results['diffusion_coefficients'].values())
        ax4.bar(ion_types, diffusion_values, color=['blue', 'green', 'red', 'orange', 'purple'])
        ax4.set_ylabel('Diffusion Coefficient (Å²/s)')
        ax4.set_title('Diffusion Coefficients')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / "plots" / "brownian_motion.png", 
                   dpi=300, bbox_inches='tight')
        plt.show()
    
    def _generate_pymol_scripts(self, binding_site_results, efficiency_results):
        """Generate PyMOL visualization scripts."""
        pymol_script = f"""
# PyMOL Script for Enhanced Metalloprotein Visualization
# Generated by Enhanced Metalloprotein Pipeline

# Load protein structure
load protein.pdb

# Show protein as cartoon
show cartoon, all
color gray, all

# Color binding sites
"""
        
        for i, site in enumerate(binding_site_results['binding_sites']):
            confidence = site['consensus_score']
            color = self._get_confidence_color(confidence)
            pymol_script += f"""
# Binding site {i+1} (confidence: {confidence:.3f})
select binding_site_{i+1}, resi {','.join(map(str, site['residues']))}
show spheres, binding_site_{i+1}
color {color}, binding_site_{i+1}
"""
        
        # Add environmental parameter visualization
        if efficiency_results is not None:
            pymol_script += """
# Environmental parameter mapping
# Temperature gradient
color_by_temperature = True
# pH gradient  
color_by_ph = True
# Redox potential gradient
color_by_redox = True
"""
        
        pymol_script += """
# Set view
set_view (\
     0.123456789,    0.987654321,    0.000000000,\
    -0.987654321,    0.123456789,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000000000,    0.000000000, -100.000000000,\
    50.000000000,   50.000000000,   50.000000000,\
    50.000000000,  150.000000000,    0.000000000 )

# Save image
png enhanced_metalloprotein_visualization.png, 1920, 1080, dpi=300
"""
        
        # Save PyMOL script
        with open(self.output_dir / "structures" / "pymol_script.pml", 'w') as f:
            f.write(pymol_script)
        
        print(f"   PyMOL script saved to: {self.output_dir / 'structures' / 'pymol_script.pml'}")
    
    def _generate_rf_diffusion_visualization(self, efficiency_results):
        """Generate RF Diffusion visualization (simulated)."""
        # This is a simplified simulation of RF Diffusion visualization
        # In practice, you would use the actual RF Diffusion model
        
        diffusion_script = f"""
# RF Diffusion Visualization Script
# Generated by Enhanced Metalloprotein Pipeline

# Parameters
diffusion_steps = {self.config['visualization']['rf_diffusion']['diffusion_steps']}
noise_schedule = "{self.config['visualization']['rf_diffusion']['noise_schedule']}"
guidance_scale = {self.config['visualization']['rf_diffusion']['guidance_scale']}

# Generate diffusion process visualization
# This would integrate with the actual RF Diffusion model
# to show ion diffusion paths and binding probability fields

# Output: diffusion_animation.mp4
"""
        
        # Save RF Diffusion script
        with open(self.output_dir / "animations" / "rf_diffusion_script.py", 'w') as f:
            f.write(diffusion_script)
        
        print(f"   RF Diffusion script saved to: {self.output_dir / 'animations' / 'rf_diffusion_script.py'}")
    
    def _get_confidence_color(self, confidence):
        """Get color based on confidence score."""
        if confidence >= 0.8:
            return "red"
        elif confidence >= 0.6:
            return "orange"
        elif confidence >= 0.4:
            return "yellow"
        else:
            return "green"
    
    def _save_enhanced_results(self, binding_site_results, kinetics_solution,
                              efficiency_results, brownian_results, protein_sequence):
        """Save enhanced results to files."""
        import json
        
        # Save binding site results
        binding_site_data = {
            'total_sites': binding_site_results['total_sites'],
            'average_consensus_score': binding_site_results['average_consensus_score'],
            'algorithm_scores': binding_site_results['algorithm_scores'],
            'binding_sites': [
                {
                    'center': site['center'].tolist(),
                    'radius': site['radius'],
                    'residues': site['residues'],
                    'consensus_score': site['consensus_score'],
                    'algorithm_count': site['algorithm_count']
                }
                for site in binding_site_results['binding_sites']
            ]
        }
        
        with open(self.output_dir / "data" / "binding_site_results.json", 'w') as f:
            json.dump(binding_site_data, f, indent=2)
        
        # Save efficiency results
        if efficiency_results is not None:
            efficiency_data = {
                'overall_efficiency': efficiency_results['overall_efficiency'],
                'environmental_analysis': efficiency_results['environmental_analysis'],
                'final_free_concentrations': efficiency_results['final_free_concentrations'],
                'final_bound_concentrations': efficiency_results['final_bound_concentrations']
            }
            
            with open(self.output_dir / "data" / "efficiency_results.json", 'w') as f:
                json.dump(efficiency_data, f, indent=2)
        
        # Save Brownian motion results
        if brownian_results is not None:
            brownian_data = {
                'diffusion_coefficients': brownian_results['diffusion_coefficients'],
                'binding_probability': {
                    ion: prob.tolist() for ion, prob in brownian_results['binding_probability'].items()
                },
                'msd': {
                    ion: msd.tolist() for ion, msd in brownian_results['msd'].items()
                }
            }
            
            with open(self.output_dir / "data" / "brownian_results.json", 'w') as f:
                json.dump(brownian_data, f, indent=2)
        
        print(f"   Results saved to: {self.output_dir / 'data'}")
    
    def _generate_comprehensive_report(self, binding_site_results, efficiency_results,
                                     brownian_results, protein_sequence, metal_ions, 
                                     initial_concentrations):
        """Generate comprehensive analysis report."""
        report = f"""
# ENHANCED METALLOPROTEIN BINDING EFFICIENCY ANALYSIS REPORT

## Executive Summary
This report presents the results of an enhanced analysis of metalloprotein binding efficiency
using a multi-algorithm approach with environmental parameter coupling.

## Protein Information
- Sequence Length: {len(protein_sequence)}
- Metal Ions Analyzed: {', '.join(metal_ions)}
- Initial Concentrations: {initial_concentrations}

## Binding Site Identification Results
- Total Binding Sites Identified: {binding_site_results['total_sites']}
- Average Consensus Score: {binding_site_results['average_consensus_score']:.3f}

### Algorithm Performance
"""
        
        for algo, score in binding_site_results['algorithm_scores'].items():
            report += f"- {algo}: {score:.3f}\n"
        
        if efficiency_results is not None:
            report += f"""
## Binding Efficiency Results
- Overall Binding Efficiency: {efficiency_results['overall_efficiency']:.3f}

### Environmental Parameter Analysis
- Temperature Range: {efficiency_results['environmental_analysis']['temperature']['range'][0]-273.15:.1f}°C to {efficiency_results['environmental_analysis']['temperature']['range'][1]-273.15:.1f}°C
- pH Range: {efficiency_results['environmental_analysis']['pH']['range'][0]:.2f} to {efficiency_results['environmental_analysis']['pH']['range'][1]:.2f}
- Pressure Range: {efficiency_results['environmental_analysis']['pressure']['range'][0]:.2f} to {efficiency_results['environmental_analysis']['pressure']['range'][1]:.2f} atm
- Redox Potential Range: {efficiency_results['environmental_analysis']['redox_potential']['range'][0]:.3f} to {efficiency_results['environmental_analysis']['redox_potential']['range'][1]:.3f} V

### Final Concentrations
"""
            for ion, conc in efficiency_results['final_free_concentrations'].items():
                report += f"- Free {ion}: {conc:.2e} M\n"
            
            for ion, conc in efficiency_results['final_bound_concentrations'].items():
                report += f"- Bound {ion}: {conc:.2e} M\n"
        
        if brownian_results is not None:
            report += f"""
## Brownian Motion Analysis
### Diffusion Coefficients
"""
            for ion, diff_coeff in brownian_results['diffusion_coefficients'].items():
                report += f"- {ion}: {diff_coeff:.2e} Å²/s\n"
        
        report += f"""
## Technical Details
- Spatial Discretization: 1000 cubes (10×10×10)
- Environmental Parameters: Temperature, pH, Pressure, Redox Potential
- Algorithms Used: MetalNet, Metal3D, bindEmbed21, AlphaFill, MESPEUS, CHED Network
- Visualization: PyMOL, RF Diffusion, Matplotlib

## Conclusions
The enhanced analysis provides a comprehensive view of metalloprotein binding efficiency
under realistic environmental conditions, incorporating multiple algorithms and
environmental parameter coupling for improved accuracy.

---
Report generated by Enhanced Metalloprotein Pipeline
"""
        
        # Save report
        with open(self.output_dir / "comprehensive_report.md", 'w') as f:
            f.write(report)
        
        print(f"   Comprehensive report saved to: {self.output_dir / 'comprehensive_report.md'}")
        
        return report

def main():
    """Main function to run enhanced pipeline."""
    # Initialize enhanced pipeline
    pipeline = EnhancedMetalloproteinPipeline()
    
    # Example usage
    pdb_file = "sample_protein.pdb"  # Replace with your PDB file
    
    # Run enhanced analysis
    results = pipeline.run_enhanced_analysis(
        pdb_file=pdb_file,
        metal_ions=['Zn2+', 'Cu2+', 'Fe2+'],
        initial_concentrations={'Zn2+': 1e-6, 'Cu2+': 1e-6, 'Fe2+': 1e-6},
        time_span=(0, 1000),
        save_results=True
    )
    
    print("\n" + "="*80)
    print("ENHANCED ANALYSIS COMPLETED SUCCESSFULLY!")
    print("="*80)
    print(f"Results saved to: {pipeline.output_dir}")
    print("Check the output directory for plots, data files, and comprehensive report.")

if __name__ == "__main__":
    main() 