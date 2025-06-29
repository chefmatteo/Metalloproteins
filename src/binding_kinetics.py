"""
Binding Kinetics Module for Metalloprotein Analysis

This module implements the ODE system for predicting metal ion binding efficiency
based on the mathematical framework described in the documentation.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.constants import k as k_B, R
import matplotlib.pyplot as plt

class BindingKinetics:
    """Implements the binding kinetics ODE system."""
    
    def __init__(self, temperature=298.15, pressure=1.0):
        """
        Initialize binding kinetics solver.
        
        Parameters:
        -----------
        temperature : float
            Temperature in Kelvin
        pressure : float
            Pressure in atm
        """
        self.T = temperature
        self.P = pressure
        self.k_B = k_B  # Boltzmann constant
        self.R = R      # Gas constant
        
        # Metal ion properties (example values)
        self.metal_properties = {
            'Zn2+': {
                'radius': 0.74e-10,  # meters
                'mass': 65.38e-27,   # kg
                'charge': 2,
                'diffusion_coeff': 7.0e-10  # m²/s at 298K
            },
            'Cu2+': {
                'radius': 0.73e-10,
                'mass': 63.55e-27,
                'charge': 2,
                'diffusion_coeff': 7.1e-10
            },
            'Fe2+': {
                'radius': 0.78e-10,
                'mass': 55.85e-27,
                'charge': 2,
                'diffusion_coeff': 7.2e-10
            },
            'Mg2+': {
                'radius': 0.72e-10,
                'mass': 24.31e-27,
                'charge': 2,
                'diffusion_coeff': 7.0e-10
            }
        }
    
    def calculate_rate_constants(self, binding_site, metal_ion, 
                               activation_energy_assoc=50e3, 
                               activation_energy_dissoc=80e3):
        """
        Calculate temperature-dependent rate constants using Arrhenius equation.
        
        Parameters:
        -----------
        binding_site : dict
            Binding site information from PDBProcessor
        metal_ion : str
            Metal ion type (e.g., 'Zn2+', 'Cu2+')
        activation_energy_assoc : float
            Association activation energy (J/mol)
        activation_energy_dissoc : float
            Dissociation activation energy (J/mol)
        
        Returns:
        --------
        tuple : (k_plus, k_minus) rate constants
        """
        # Pre-exponential factors (collision frequency)
        A_plus = 1e9  # M⁻¹s⁻¹ (typical for diffusion-limited reactions)
        A_minus = 1e3  # s⁻¹ (typical for dissociation)
        
        # Arrhenius equation
        k_plus = A_plus * np.exp(-activation_energy_assoc / (self.R * self.T))
        k_minus = A_minus * np.exp(-activation_energy_dissoc / (self.R * self.T))
        
        # Pressure correction
        k_plus = self._apply_pressure_correction(k_plus, binding_site, 'association')
        k_minus = self._apply_pressure_correction(k_minus, binding_site, 'dissociation')
        
        return k_plus, k_minus
    
    def _apply_pressure_correction(self, rate_constant, binding_site, reaction_type):
        """Apply pressure correction to rate constants."""
        # Simplified pressure effect through activation volume
        if reaction_type == 'association':
            delta_V = -1e-6  # m³/mol (typical for association)
        else:
            delta_V = 1e-6   # m³/mol (typical for dissociation)
        
        # Pressure correction factor
        P_correction = np.exp(-delta_V * self.P * 101325 / (self.R * self.T))  # Convert atm to Pa
        return rate_constant * P_correction
    
    def solve_binding_kinetics(self, binding_sites, metal_ions, initial_concentrations, 
                             time_span=(0, 1000), method='RK45'):
        """
        Solve the coupled ODE system for binding kinetics.
        
        Parameters:
        -----------
        binding_sites : list
            List of binding site dictionaries
        metal_ions : list
            List of metal ion types
        initial_concentrations : dict
            Initial concentrations of metal ions
        time_span : tuple
            Time span for integration (start, end) in seconds
        method : str
            ODE solver method
        
        Returns:
        --------
        scipy.integrate.OdeSolution
            Solution object containing time points and concentrations
        """
        # Set up the ODE system
        def binding_ode_system(t, y):
            return self._binding_ode_rhs(t, y, binding_sites, metal_ions, initial_concentrations)
        
        # Initial conditions
        y0 = self._setup_initial_conditions(binding_sites, metal_ions, initial_concentrations)
        
        # Solve ODE system
        solution = solve_ivp(
            binding_ode_system,
            time_span,
            y0,
            method=method,
            t_eval=np.linspace(time_span[0], time_span[1], 1000)
        )
        
        return solution
    
    def _binding_ode_rhs(self, t, y, binding_sites, metal_ions, initial_concentrations):
        """
        Right-hand side of the binding kinetics ODE system.
        
        The ODE system includes:
        1. Free ion concentrations
        2. Bound ion concentrations
        3. Binding site occupancies
        """
        n_sites = len(binding_sites)
        n_ions = len(metal_ions)
        
        # Extract variables from y vector
        free_ions = y[:n_ions]  # Free ion concentrations
        bound_ions = y[n_ions:n_ions + n_sites * n_ions]  # Bound ion concentrations
        occupancies = y[n_ions + n_sites * n_ions:]  # Site occupancies
        
        # Reshape bound_ions and occupancies
        bound_ions = bound_ions.reshape(n_sites, n_ions)
        occupancies = occupancies.reshape(n_sites, n_ions)
        
        # Initialize derivatives
        d_free_ions = np.zeros(n_ions)
        d_bound_ions = np.zeros((n_sites, n_ions))
        d_occupancies = np.zeros((n_sites, n_ions))
        
        # Calculate derivatives for each binding site and metal ion
        for i, site in enumerate(binding_sites):
            for j, ion in enumerate(metal_ions):
                # Get rate constants
                k_plus, k_minus = self.calculate_rate_constants(site, ion)
                
                # Association term
                association_term = k_plus * free_ions[j] * (1 - occupancies[i, j])
                
                # Dissociation term
                dissociation_term = k_minus * occupancies[i, j]
                
                # Update derivatives
                d_occupancies[i, j] = association_term - dissociation_term
                d_bound_ions[i, j] = association_term - dissociation_term
                d_free_ions[j] -= association_term - dissociation_term
        
        # Return flattened derivatives
        return np.concatenate([d_free_ions, d_bound_ions.flatten(), d_occupancies.flatten()])
    
    def _setup_initial_conditions(self, binding_sites, metal_ions, initial_concentrations):
        """Set up initial conditions for the ODE system."""
        n_sites = len(binding_sites)
        n_ions = len(metal_ions)
        
        # Initial free ion concentrations
        free_ions = [initial_concentrations.get(ion, 1e-6) for ion in metal_ions]
        
        # Initial bound ion concentrations (all zero)
        bound_ions = np.zeros(n_sites * n_ions)
        
        # Initial occupancies (all zero)
        occupancies = np.zeros(n_sites * n_ions)
        
        return np.concatenate([free_ions, bound_ions, occupancies])
    
    def calculate_binding_efficiency(self, solution, binding_sites, metal_ions, initial_concentrations):
        """
        Calculate overall binding efficiency from ODE solution.
        
        Parameters:
        -----------
        solution : scipy.integrate.OdeSolution
            Solution from ODE solver
        binding_sites : list
            List of binding site dictionaries
        metal_ions : list
            List of metal ion types
        initial_concentrations : dict
            Initial concentrations of metal ions
        
        Returns:
        --------
        dict : Binding efficiency metrics
        """
        n_sites = len(binding_sites)
        n_ions = len(metal_ions)
        
        # Extract final values
        final_y = solution.y[:, -1]
        free_ions = final_y[:n_ions]
        bound_ions = final_y[n_ions:n_ions + n_sites * n_ions].reshape(n_sites, n_ions)
        occupancies = final_y[n_ions + n_sites * n_ions:].reshape(n_sites, n_ions)
        
        # Calculate total initial and final bound ions
        total_initial = sum(initial_concentrations.get(ion, 1e-6) for ion in metal_ions)
        total_bound = np.sum(bound_ions)
        
        # Overall binding efficiency
        efficiency = total_bound / total_initial if total_initial > 0 else 0
        
        # Per-site and per-ion efficiencies
        site_efficiencies = {}
        ion_efficiencies = {}
        
        for i, site in enumerate(binding_sites):
            site_efficiencies[f"site_{i}"] = {
                'total_occupancy': np.sum(occupancies[i, :]),
                'max_occupancy': min(1.0, len(site['coordinating_atoms']))
            }
        
        for j, ion in enumerate(metal_ions):
            initial_conc = initial_concentrations.get(ion, 1e-6)
            final_bound = np.sum(bound_ions[:, j])
            ion_efficiencies[ion] = final_bound / initial_conc if initial_conc > 0 else 0
        
        return {
            'overall_efficiency': efficiency,
            'site_efficiencies': site_efficiencies,
            'ion_efficiencies': ion_efficiencies,
            'final_free_concentrations': dict(zip(metal_ions, free_ions)),
            'final_bound_concentrations': dict(zip(metal_ions, np.sum(bound_ions, axis=0)))
        }
    
    def plot_binding_kinetics(self, solution, binding_sites, metal_ions, save_path=None):
        """Plot binding kinetics results."""
        n_sites = len(binding_sites)
        n_ions = len(metal_ions)
        
        # Handle case with no binding sites
        if n_sites == 0:
            print("Warning: No binding sites found. Skipping kinetics plots.")
            return None
        
        # Extract solution components
        times = solution.t
        free_ions = solution.y[:n_ions, :]
        bound_ions = solution.y[n_ions:n_ions + n_sites * n_ions, :].reshape(n_sites, n_ions, -1)
        occupancies = solution.y[n_ions + n_sites * n_ions:, :].reshape(n_sites, n_ions, -1)
        
        # Create subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Plot 1: Free ion concentrations
        ax1 = axes[0, 0]
        for i, ion in enumerate(metal_ions):
            ax1.plot(times, free_ions[i, :], label=ion, linewidth=2)
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Free Ion Concentration (M)')
        ax1.set_title('Free Ion Concentration vs Time')
        ax1.legend()
        ax1.grid(True)
        
        # Plot 2: Total bound ions
        ax2 = axes[0, 1]
        total_bound = np.sum(bound_ions, axis=0)  # Sum over sites
        for i, ion in enumerate(metal_ions):
            ax2.plot(times, total_bound[i, :], label=ion, linewidth=2)
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Total Bound Ions (M)')
        ax2.set_title('Total Bound Ions vs Time')
        ax2.legend()
        ax2.grid(True)
        
        # Plot 3: Site occupancies
        ax3 = axes[1, 0]
        for i in range(n_sites):
            for j, ion in enumerate(metal_ions):
                ax3.plot(times, occupancies[i, j, :], 
                        label=f'Site {i+1} - {ion}', linewidth=2, alpha=0.7)
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Site Occupancy')
        ax3.set_title('Binding Site Occupancy vs Time')
        ax3.legend()
        ax3.grid(True)
        
        # Plot 4: Binding efficiency
        ax4 = axes[1, 1]
        efficiency = total_bound / (total_bound + free_ions)
        for i, ion in enumerate(metal_ions):
            ax4.plot(times, efficiency[i, :], label=ion, linewidth=2)
        ax4.set_xlabel('Time (s)')
        ax4.set_ylabel('Binding Efficiency')
        ax4.set_title('Binding Efficiency vs Time')
        ax4.legend()
        ax4.grid(True)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
        
        return fig 