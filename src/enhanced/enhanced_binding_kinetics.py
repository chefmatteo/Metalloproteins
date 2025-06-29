"""
Enhanced Binding Kinetics Module for Metalloprotein Analysis

This module implements the enhanced ODE/PDE system for predicting metal ion binding efficiency
with environmental parameter coupling (temperature, pH, pressure, redox potential).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.constants import k as k_B, R, F
from scipy.sparse import csr_matrix, lil_matrix
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import multiprocessing as mp
from functools import partial
import warnings
warnings.filterwarnings('ignore')

class EnhancedBindingKinetics:
    """Implements the enhanced binding kinetics with environmental parameter coupling."""
    
    def __init__(self, config):
        """
        Initialize enhanced binding kinetics solver.
        
        Parameters:
        -----------
        config : dict
            Configuration dictionary containing all parameters
        """
        self.config = config
        self.k_B = k_B  # Boltzmann constant
        self.R = R      # Gas constant
        self.F = F      # Faraday constant
        
        # Extract configuration parameters
        self.env_config = config['environmental_conditions']
        self.spatial_config = config['spatial_discretization']
        self.metal_config = config['metal_ions']
        
        # Initialize spatial grid
        self._initialize_spatial_grid()
        
        # Initialize environmental parameters
        self._initialize_environmental_parameters()
        
    def _initialize_spatial_grid(self):
        """Initialize the 1000-cube spatial discretization."""
        grid_size = self.spatial_config['chamber']['grid_size']
        self.nx, self.ny, self.nz = grid_size
        self.n_cubes = self.nx * self.ny * self.nz
        
        # Create coordinate arrays
        x = np.linspace(0, self.spatial_config['chamber']['dimensions'][0], self.nx)
        y = np.linspace(0, self.spatial_config['chamber']['dimensions'][1], self.ny)
        z = np.linspace(0, self.spatial_config['chamber']['dimensions'][2], self.nz)
        
        self.X, self.Y, self.Z = np.meshgrid(x, y, z, indexing='ij')
        self.cube_centers = np.column_stack([
            self.X.flatten(), self.Y.flatten(), self.Z.flatten()
        ])
        
        # Calculate cube volume
        self.cube_volume = self.spatial_config['cube']['volume']
        
    def _initialize_environmental_parameters(self):
        """Initialize environmental parameters for each cube."""
        # Temperature
        self.T = np.full(self.n_cubes, self.env_config['temperature']['initial'])
        
        # pH
        self.pH = np.full(self.n_cubes, self.env_config['pH']['initial'])
        self.H_plus = 10**(-self.pH)  # Proton concentration
        
        # Pressure
        self.P = np.full(self.n_cubes, self.env_config['pressure']['initial'])
        
        # Redox potential
        self.Eh = np.full(self.n_cubes, self.env_config['redox_potential']['initial'])
        
    def calculate_temperature_dependent_diffusion(self, metal_ion, T, P):
        """
        Calculate temperature and pressure dependent diffusion coefficient.
        
        Parameters:
        -----------
        metal_ion : str
            Metal ion type
        T : np.array
            Temperature array for each cube
        P : np.array
            Pressure array for each cube
        
        Returns:
        --------
        np.array : Diffusion coefficients for each cube
        """
        metal_props = self.metal_config[metal_ion]
        D0 = metal_props['diffusion_coeff_0']
        Ea_D = metal_props['activation_energy_diffusion']
        Delta_V_D = metal_props['activation_volume_diffusion']
        T0 = 298.15  # Reference temperature
        
        # Temperature dependence (Arrhenius)
        D_T = D0 * (T / T0) * np.exp(-Ea_D / self.R * (1/T - 1/T0))
        
        # Pressure dependence
        D_P = D_T * np.exp(-P * Delta_V_D / (self.R * T))
        
        return D_P
    
    def calculate_enhanced_rate_constants(self, binding_site, metal_ion, T, P, pH, Eh):
        """
        Calculate enhanced rate constants with environmental parameter dependence.
        
        Parameters:
        -----------
        binding_site : dict
            Binding site information
        metal_ion : str
            Metal ion type
        T : np.array
            Temperature array
        P : np.array
            Pressure array
        pH : np.array
            pH array
        Eh : np.array
            Redox potential array
        
        Returns:
        --------
        tuple : (k_plus, k_minus) rate constant arrays
        """
        metal_props = self.metal_config[metal_ion]
        
        # Base rate constants (Arrhenius)
        A_plus = 1e9  # M⁻¹s⁻¹
        A_minus = 1e3  # s⁻¹
        Ea_plus = 50e3  # J/mol (association)
        Ea_minus = 80e3  # J/mol (dissociation)
        
        # Temperature dependence
        k_plus_T = A_plus * np.exp(-Ea_plus / (self.R * T))
        k_minus_T = A_minus * np.exp(-Ea_minus / (self.R * T))
        
        # Pressure dependence
        Delta_V_plus = metal_props.get('activation_volume_association', -1e-6)
        Delta_V_minus = metal_props.get('activation_volume_dissociation', 1e-6)
        
        k_plus_P = k_plus_T * np.exp(-P * Delta_V_plus / (self.R * T))
        k_minus_P = k_minus_T * np.exp(-P * Delta_V_minus / (self.R * T))
        
        # pH dependence
        pKa = 9.0  # Typical pKa for metal-binding residues
        f_pH = 1 / (1 + 10**(pH - pKa))
        pKa_effect = metal_props.get('pKa_effect', 0.5)
        
        k_plus_pH = k_plus_P * (1 - pKa_effect * (1 - f_pH))
        k_minus_pH = k_minus_P * (1 + pKa_effect * (1 - f_pH))
        
        # Redox dependence
        Eh0 = 0.0  # Standard redox potential
        n = 1  # Number of electrons
        redox_sensitivity = metal_props.get('redox_sensitivity', 0.1)
        
        f_Eh = np.exp(-n * self.F * (Eh - Eh0) / (self.R * T))
        
        k_plus_Eh = k_plus_pH * (1 + redox_sensitivity * (f_Eh - 1))
        k_minus_Eh = k_minus_pH * (1 - redox_sensitivity * (f_Eh - 1))
        
        return k_plus_Eh, k_minus_Eh
    
    def solve_enhanced_kinetics(self, binding_sites, metal_ions, initial_concentrations, 
                              time_span=(0, 1000), method='RK45'):
        """
        Solve the enhanced coupled ODE/PDE system.
        
        Parameters:
        -----------
        binding_sites : list
            List of binding site dictionaries
        metal_ions : list
            List of metal ion types
        initial_concentrations : dict
            Initial concentrations of metal ions
        time_span : tuple
            Time span for integration
        method : str
            ODE solver method
        
        Returns:
        --------
        dict : Solution and analysis results
        """
        print("Solving enhanced coupled ODE/PDE system...")
        
        # Set up the enhanced ODE system
        def enhanced_ode_system(t, y):
            return self._enhanced_ode_rhs(t, y, binding_sites, metal_ions, initial_concentrations)
        
        # Initial conditions
        y0 = self._setup_enhanced_initial_conditions(binding_sites, metal_ions, initial_concentrations)
        
        # Solve ODE system
        solution = solve_ivp(
            enhanced_ode_system,
            time_span,
            y0,
            method=method,
            t_eval=np.linspace(time_span[0], time_span[1], 1000),
            rtol=1e-6,
            atol=1e-9
        )
        
        return solution
    
    def _enhanced_ode_rhs(self, t, y, binding_sites, metal_ions, initial_concentrations):
        """
        Right-hand side of the enhanced ODE system with environmental coupling.
        
        The enhanced system includes:
        1. Free ion concentrations in each cube
        2. Bound ion concentrations
        3. Binding site occupancies
        4. Environmental parameters (T, pH, P, Eh)
        """
        n_sites = len(binding_sites)
        n_ions = len(metal_ions)
        
        # Extract variables from y vector
        # Structure: [free_ions (n_cubes * n_ions), bound_ions (n_sites * n_ions), 
        #            occupancies (n_sites * n_ions), T (n_cubes), pH (n_cubes), 
        #            P (n_cubes), Eh (n_cubes)]
        
        start_idx = 0
        free_ions = y[start_idx:start_idx + self.n_cubes * n_ions].reshape(self.n_cubes, n_ions)
        start_idx += self.n_cubes * n_ions
        
        bound_ions = y[start_idx:start_idx + n_sites * n_ions].reshape(n_sites, n_ions)
        start_idx += n_sites * n_ions
        
        occupancies = y[start_idx:start_idx + n_sites * n_ions].reshape(n_sites, n_ions)
        start_idx += n_sites * n_ions
        
        T = y[start_idx:start_idx + self.n_cubes]
        start_idx += self.n_cubes
        
        pH = y[start_idx:start_idx + self.n_cubes]
        start_idx += self.n_cubes
        
        P = y[start_idx:start_idx + self.n_cubes]
        start_idx += self.n_cubes
        
        Eh = y[start_idx:start_idx + self.n_cubes]
        
        # Calculate H+ concentration from pH
        H_plus = 10**(-pH)
        
        # Initialize derivatives
        d_free_ions = np.zeros((self.n_cubes, n_ions))
        d_bound_ions = np.zeros((n_sites, n_ions))
        d_occupancies = np.zeros((n_sites, n_ions))
        d_T = np.zeros(self.n_cubes)
        d_pH = np.zeros(self.n_cubes)
        d_P = np.zeros(self.n_cubes)
        d_Eh = np.zeros(self.n_cubes)
        
        # Calculate diffusion and reaction terms for each cube
        for cube_idx in range(self.n_cubes):
            # Diffusion terms
            d_free_ions[cube_idx, :] = self._calculate_diffusion_terms(
                cube_idx, free_ions, metal_ions, T[cube_idx], P[cube_idx]
            )
            
            # Reaction terms
            for i, site in enumerate(binding_sites):
                for j, ion in enumerate(metal_ions):
                    # Get rate constants for this cube
                    k_plus, k_minus = self.calculate_enhanced_rate_constants(
                        site, ion, T[cube_idx], P[cube_idx], pH[cube_idx], Eh[cube_idx]
                    )
                    
                    # Binding reaction rate
                    binding_rate = k_plus * free_ions[cube_idx, j] * (1 - occupancies[i, j])
                    unbinding_rate = k_minus * occupancies[i, j]
                    
                    # Update derivatives
                    d_free_ions[cube_idx, j] -= binding_rate - unbinding_rate
                    d_bound_ions[i, j] += binding_rate - unbinding_rate
                    d_occupancies[i, j] = binding_rate - unbinding_rate
                    
                    # Heat generation from binding
                    metal_props = self.metal_config[ion]
                    Delta_H = metal_props.get('binding_enthalpy', -50e3)
                    d_T[cube_idx] += -Delta_H * (binding_rate - unbinding_rate) / (
                        self.env_config['temperature']['heat_transfer']['density'] *
                        self.env_config['temperature']['heat_transfer']['specific_heat']
                    )
                    
                    # Proton release from binding
                    d_pH[cube_idx] += 0.1 * (binding_rate - unbinding_rate)  # Simplified
        
        # Environmental parameter evolution
        d_T += self._calculate_heat_diffusion(T)
        d_pH += self._calculate_proton_diffusion(H_plus)
        d_P += self._calculate_pressure_evolution(P, T, free_ions)
        d_Eh += self._calculate_redox_evolution(Eh, bound_ions)
        
        # Return flattened derivatives
        return np.concatenate([
            d_free_ions.flatten(), d_bound_ions.flatten(), d_occupancies.flatten(),
            d_T, d_pH, d_P, d_Eh
        ])
    
    def _calculate_diffusion_terms(self, cube_idx, free_ions, metal_ions, T, P):
        """Calculate diffusion terms for a specific cube."""
        diffusion_terms = np.zeros(len(metal_ions))
        
        for j, ion in enumerate(metal_ions):
            # Get diffusion coefficient
            D = self.calculate_temperature_dependent_diffusion(ion, T, P)
            
            # Simple finite difference approximation
            # This is a simplified version - in practice, you'd use proper finite differences
            diffusion_terms[j] = D * 1e6  # Simplified diffusion term
        
        return diffusion_terms
    
    def _calculate_heat_diffusion(self, T):
        """Calculate heat diffusion terms."""
        # Simplified heat diffusion
        kappa = self.env_config['temperature']['heat_transfer']['thermal_conductivity']
        return kappa * 1e6 * np.gradient(T)  # Simplified
    
    def _calculate_proton_diffusion(self, H_plus):
        """Calculate proton diffusion terms."""
        D_H = self.env_config['pH']['proton_diffusion_coeff']
        return D_H * 1e6 * np.gradient(H_plus)  # Simplified
    
    def _calculate_pressure_evolution(self, P, T, free_ions):
        """Calculate pressure evolution terms."""
        beta_T = self.env_config['pressure']['compressibility']['thermal']
        beta_C = self.env_config['pressure']['compressibility']['chemical']
        
        dP = beta_T * np.gradient(T) + beta_C * np.sum(np.gradient(free_ions, axis=0), axis=1)
        return dP
    
    def _calculate_redox_evolution(self, Eh, bound_ions):
        """Calculate redox potential evolution terms."""
        D_ox = self.env_config['redox_potential']['redox_diffusion_coeff']
        return D_ox * 1e6 * np.gradient(Eh)  # Simplified
    
    def _setup_enhanced_initial_conditions(self, binding_sites, metal_ions, initial_concentrations):
        """Set up initial conditions for the enhanced ODE system."""
        n_sites = len(binding_sites)
        n_ions = len(metal_ions)
        
        # Initial free ion concentrations (uniform across all cubes)
        free_ions = np.array([initial_concentrations.get(ion, 1e-6) for ion in metal_ions])
        free_ions_all_cubes = np.tile(free_ions, (self.n_cubes, 1))
        
        # Initial bound ion concentrations (all zero)
        bound_ions = np.zeros(n_sites * n_ions)
        
        # Initial occupancies (all zero)
        occupancies = np.zeros(n_sites * n_ions)
        
        # Initial environmental parameters
        T_init = np.full(self.n_cubes, self.env_config['temperature']['initial'])
        pH_init = np.full(self.n_cubes, self.env_config['pH']['initial'])
        P_init = np.full(self.n_cubes, self.env_config['pressure']['initial'])
        Eh_init = np.full(self.n_cubes, self.env_config['redox_potential']['initial'])
        
        return np.concatenate([
            free_ions_all_cubes.flatten(), bound_ions, occupancies,
            T_init, pH_init, P_init, Eh_init
        ])
    
    def calculate_enhanced_binding_efficiency(self, solution, binding_sites, metal_ions, initial_concentrations):
        """
        Calculate enhanced binding efficiency with environmental parameter analysis.
        
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
        dict : Enhanced binding efficiency metrics
        """
        n_sites = len(binding_sites)
        n_ions = len(metal_ions)
        
        # Extract final values
        final_y = solution.y[:, -1]
        
        # Parse solution vector
        start_idx = 0
        free_ions = final_y[start_idx:start_idx + self.n_cubes * n_ions].reshape(self.n_cubes, n_ions)
        start_idx += self.n_cubes * n_ions
        
        bound_ions = final_y[start_idx:start_idx + n_sites * n_ions].reshape(n_sites, n_ions)
        start_idx += n_sites * n_ions
        
        occupancies = final_y[start_idx:start_idx + n_sites * n_ions].reshape(n_sites, n_ions)
        start_idx += n_sites * n_ions
        
        T_final = final_y[start_idx:start_idx + self.n_cubes]
        start_idx += self.n_cubes
        
        pH_final = final_y[start_idx:start_idx + self.n_cubes]
        start_idx += self.n_cubes
        
        P_final = final_y[start_idx:start_idx + self.n_cubes]
        start_idx += self.n_cubes
        
        Eh_final = final_y[start_idx:start_idx + self.n_cubes]
        
        # Calculate efficiency metrics
        total_initial = sum(initial_concentrations.get(ion, 1e-6) for ion in metal_ions) * self.n_cubes
        total_bound = np.sum(bound_ions)
        overall_efficiency = total_bound / total_initial if total_initial > 0 else 0
        
        # Per-cube analysis
        cube_efficiencies = []
        for cube_idx in range(self.n_cubes):
            cube_bound = np.sum(bound_ions[:, :]) / self.n_cubes  # Average bound ions per cube
            cube_initial = sum(initial_concentrations.get(ion, 1e-6) for ion in metal_ions)
            cube_efficiency = cube_bound / cube_initial if cube_initial > 0 else 0
            cube_efficiencies.append(cube_efficiency)
        
        # Environmental parameter analysis
        env_analysis = {
            'temperature': {
                'mean': np.mean(T_final),
                'std': np.std(T_final),
                'range': [np.min(T_final), np.max(T_final)]
            },
            'pH': {
                'mean': np.mean(pH_final),
                'std': np.std(pH_final),
                'range': [np.min(pH_final), np.max(pH_final)]
            },
            'pressure': {
                'mean': np.mean(P_final),
                'std': np.std(P_final),
                'range': [np.min(P_final), np.max(P_final)]
            },
            'redox_potential': {
                'mean': np.mean(Eh_final),
                'std': np.std(Eh_final),
                'range': [np.min(Eh_final), np.max(Eh_final)]
            }
        }
        
        return {
            'overall_efficiency': overall_efficiency,
            'cube_efficiencies': cube_efficiencies,
            'environmental_analysis': env_analysis,
            'final_free_concentrations': dict(zip(metal_ions, np.mean(free_ions, axis=0))),
            'final_bound_concentrations': dict(zip(metal_ions, np.sum(bound_ions, axis=0))),
            'spatial_data': {
                'free_ions': free_ions,
                'bound_ions': bound_ions,
                'temperature': T_final,
                'pH': pH_final,
                'pressure': P_final,
                'redox_potential': Eh_final
            }
        }
    
    def plot_enhanced_results(self, solution, binding_sites, metal_ions, save_path=None):
        """Plot enhanced binding kinetics results with environmental parameters."""
        n_sites = len(binding_sites)
        n_ions = len(metal_ions)
        
        # Extract solution components
        times = solution.t
        
        # Parse solution for plotting
        start_idx = 0
        free_ions = solution.y[start_idx:start_idx + self.n_cubes * n_ions, :].reshape(self.n_cubes, n_ions, -1)
        start_idx += self.n_cubes * n_ions
        
        bound_ions = solution.y[start_idx:start_idx + n_sites * n_ions, :].reshape(n_sites, n_ions, -1)
        start_idx += n_sites * n_ions
        
        occupancies = solution.y[start_idx:start_idx + n_sites * n_ions, :].reshape(n_sites, n_ions, -1)
        start_idx += n_sites * n_ions
        
        T = solution.y[start_idx:start_idx + self.n_cubes, :]
        start_idx += self.n_cubes
        
        pH = solution.y[start_idx:start_idx + self.n_cubes, :]
        start_idx += self.n_cubes
        
        P = solution.y[start_idx:start_idx + self.n_cubes, :]
        start_idx += self.n_cubes
        
        Eh = solution.y[start_idx:start_idx + self.n_cubes, :]
        
        # Create subplots
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        
        # Plot 1: Average free ion concentrations
        ax1 = axes[0, 0]
        for i, ion in enumerate(metal_ions):
            avg_concentration = np.mean(free_ions[:, i, :], axis=0)
            ax1.plot(times, avg_concentration, label=ion, linewidth=2)
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Average Free Ion Concentration (M)')
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
        
        # Plot 3: Temperature evolution
        ax3 = axes[0, 2]
        avg_temp = np.mean(T, axis=0)
        ax3.plot(times, avg_temp - 273.15, linewidth=2, color='red')
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Average Temperature (°C)')
        ax3.set_title('Temperature Evolution')
        ax3.grid(True)
        
        # Plot 4: pH evolution
        ax4 = axes[1, 0]
        avg_ph = np.mean(pH, axis=0)
        ax4.plot(times, avg_ph, linewidth=2, color='blue')
        ax4.set_xlabel('Time (s)')
        ax4.set_ylabel('Average pH')
        ax4.set_title('pH Evolution')
        ax4.grid(True)
        
        # Plot 5: Pressure evolution
        ax5 = axes[1, 1]
        avg_pressure = np.mean(P, axis=0)
        ax5.plot(times, avg_pressure, linewidth=2, color='green')
        ax5.set_xlabel('Time (s)')
        ax5.set_ylabel('Average Pressure (atm)')
        ax5.set_title('Pressure Evolution')
        ax5.grid(True)
        
        # Plot 6: Redox potential evolution
        ax6 = axes[1, 2]
        avg_redox = np.mean(Eh, axis=0)
        ax6.plot(times, avg_redox, linewidth=2, color='purple')
        ax6.set_xlabel('Time (s)')
        ax6.set_ylabel('Average Redox Potential (V)')
        ax6.set_title('Redox Potential Evolution')
        ax6.grid(True)
        
        # Plot 7: Binding efficiency
        ax7 = axes[2, 0]
        efficiency = total_bound / (total_bound + np.mean(free_ions, axis=0))
        for i, ion in enumerate(metal_ions):
            ax7.plot(times, efficiency[i, :], label=ion, linewidth=2)
        ax7.set_xlabel('Time (s)')
        ax7.set_ylabel('Binding Efficiency')
        ax7.set_title('Binding Efficiency vs Time')
        ax7.legend()
        ax7.grid(True)
        
        # Plot 8: Spatial temperature distribution (final time)
        ax8 = axes[2, 1]
        final_temp = T[:, -1].reshape(self.nx, self.ny, self.nz)
        temp_slice = final_temp[:, :, self.nz//2]  # Middle slice
        im8 = ax8.imshow(temp_slice - 273.15, cmap='hot', aspect='auto')
        ax8.set_title('Final Temperature Distribution (°C)')
        plt.colorbar(im8, ax=ax8)
        
        # Plot 9: Spatial pH distribution (final time)
        ax9 = axes[2, 2]
        final_ph = pH[:, -1].reshape(self.nx, self.ny, self.nz)
        ph_slice = final_ph[:, :, self.nz//2]  # Middle slice
        im9 = ax9.imshow(ph_slice, cmap='RdYlBu_r', aspect='auto')
        ax9.set_title('Final pH Distribution')
        plt.colorbar(im9, ax=ax9)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
        
        return fig 