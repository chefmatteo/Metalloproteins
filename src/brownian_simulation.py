"""
Brownian Motion Simulation for Metal Ion Diffusion

This module implements stochastic Brownian motion simulation for metal ions
in solution, including collision detection with protein binding sites.
"""

import numpy as np
from scipy.constants import k as k_B
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class BrownianSimulation:
    """Simulates Brownian motion of metal ions in solution."""
    
    def __init__(self, temperature=298.15, viscosity=1e-3, box_size=100e-9):
        """
        Initialize Brownian motion simulation.
        
        Parameters:
        -----------
        temperature : float
            Temperature in Kelvin
        viscosity : float
            Solvent viscosity in Pa·s
        box_size : float
            Simulation box size in meters
        """
        self.T = temperature
        self.eta = viscosity
        self.box_size = box_size
        self.k_B = k_B
        
        # Time step for simulation
        self.dt = 1e-12  # 1 ps
        
    def simulate_ion_trajectory(self, ion_properties, initial_position, 
                              binding_sites, simulation_time=1e-9):
        """
        Simulate Brownian motion trajectory of a single ion.
        
        Parameters:
        -----------
        ion_properties : dict
            Ion properties (mass, radius, charge)
        initial_position : np.array
            Initial position [x, y, z]
        binding_sites : list
            List of binding site dictionaries
        simulation_time : float
            Total simulation time in seconds
        
        Returns:
        --------
        dict : Trajectory and collision data
        """
        # Extract ion properties
        mass = ion_properties['mass']
        radius = ion_properties['radius']
        
        # Calculate friction coefficient (Stokes law)
        gamma = 6 * np.pi * self.eta * radius
        
        # Calculate diffusion coefficient
        D = self.k_B * self.T / gamma
        
        # Number of time steps
        n_steps = int(simulation_time / self.dt)
        
        # Initialize trajectory
        positions = np.zeros((n_steps, 3))
        velocities = np.zeros((n_steps, 3))
        positions[0] = initial_position
        
        # Track collisions
        collisions = []
        
        # Run simulation
        for step in range(1, n_steps):
            # Current position and velocity
            r_old = positions[step-1]
            v_old = velocities[step-1]
            
            # Random force (Gaussian white noise)
            random_force = np.sqrt(2 * gamma * self.k_B * self.T / self.dt) * np.random.randn(3)
            
            # External forces (from binding sites)
            external_force = self._calculate_binding_forces(r_old, binding_sites)
            
            # Total force
            total_force = external_force + random_force
            
            # Update velocity (Langevin equation)
            v_new = v_old + (total_force - gamma * v_old) * self.dt / mass
            
            # Update position
            r_new = r_old + v_new * self.dt
            
            # Apply periodic boundary conditions
            r_new = self._apply_boundary_conditions(r_new)
            
            # Store new position and velocity
            positions[step] = r_new
            velocities[step] = v_new
            
            # Check for collisions with binding sites
            collision = self._check_collision(r_new, binding_sites, radius)
            if collision:
                collisions.append({
                    'time': step * self.dt,
                    'position': r_new,
                    'binding_site': collision['site'],
                    'distance': collision['distance']
                })
        
        return {
            'positions': positions,
            'velocities': velocities,
            'times': np.arange(n_steps) * self.dt,
            'collisions': collisions,
            'diffusion_coefficient': D
        }
    
    def _calculate_binding_forces(self, position, binding_sites):
        """Calculate forces from binding sites."""
        total_force = np.zeros(3)
        
        for site in binding_sites:
            site_center = np.array(site['center'])
            distance = np.linalg.norm(position - site_center)
            
            if distance > 0:
                # Simplified force model (Lennard-Jones like)
                r0 = 3e-10  # Equilibrium distance
                epsilon = 1e-20  # Interaction strength
                
                # Force magnitude
                force_mag = 12 * epsilon * ((r0/distance)**13 - (r0/distance)**7)
                
                # Force direction
                force_dir = (site_center - position) / distance
                
                total_force += force_mag * force_dir
        
        return total_force
    
    def _check_collision(self, position, binding_sites, ion_radius):
        """Check if ion collides with any binding site."""
        for site in binding_sites:
            site_center = np.array(site['center'])
            distance = np.linalg.norm(position - site_center)
            
            # Collision threshold (ion radius + binding site radius)
            collision_threshold = ion_radius + 2e-10  # 2 Å binding site radius
            
            if distance < collision_threshold:
                return {
                    'site': site,
                    'distance': distance
                }
        
        return None
    
    def _apply_boundary_conditions(self, position):
        """Apply periodic boundary conditions."""
        return position % self.box_size
    
    def simulate_multiple_ions(self, ion_properties, n_ions, binding_sites, 
                             simulation_time=1e-9, initial_positions=None):
        """
        Simulate multiple ions simultaneously.
        
        Parameters:
        -----------
        ion_properties : dict
            Ion properties
        n_ions : int
            Number of ions to simulate
        binding_sites : list
            List of binding site dictionaries
        simulation_time : float
            Total simulation time
        initial_positions : np.array, optional
            Initial positions for all ions
        
        Returns:
        --------
        dict : Multi-ion simulation results
        """
        if initial_positions is None:
            # Random initial positions
            initial_positions = np.random.rand(n_ions, 3) * self.box_size
        
        # Simulate each ion
        trajectories = []
        all_collisions = []
        
        for i in range(n_ions):
            trajectory = self.simulate_ion_trajectory(
                ion_properties, 
                initial_positions[i], 
                binding_sites, 
                simulation_time
            )
            trajectories.append(trajectory)
            all_collisions.extend(trajectory['collisions'])
        
        return {
            'trajectories': trajectories,
            'all_collisions': all_collisions,
            'n_ions': n_ions,
            'simulation_time': simulation_time
        }
    
    def calculate_binding_probability(self, simulation_results, binding_sites):
        """
        Calculate binding probability from simulation results.
        
        Parameters:
        -----------
        simulation_results : dict
            Results from multi-ion simulation
        binding_sites : list
            List of binding site dictionaries
        
        Returns:
        --------
        dict : Binding probability statistics
        """
        n_ions = simulation_results['n_ions']
        collisions = simulation_results['all_collisions']
        
        # Count collisions per binding site
        site_collisions = {i: 0 for i in range(len(binding_sites))}
        
        for collision in collisions:
            site_index = binding_sites.index(collision['binding_site'])
            site_collisions[site_index] += 1
        
        # Calculate probabilities
        total_collisions = len(collisions)
        binding_probabilities = {
            f'site_{i}': count / n_ions if n_ions > 0 else 0
            for i, count in site_collisions.items()
        }
        
        overall_probability = total_collisions / n_ions if n_ions > 0 else 0
        
        return {
            'overall_binding_probability': overall_probability,
            'site_binding_probabilities': binding_probabilities,
            'total_collisions': total_collisions,
            'collisions_per_site': site_collisions
        }
    
    def plot_trajectories(self, simulation_results, binding_sites, save_path=None):
        """Plot ion trajectories and binding sites."""
        fig = plt.figure(figsize=(15, 10))
        
        # 3D trajectory plot
        ax1 = fig.add_subplot(221, projection='3d')
        
        for i, trajectory in enumerate(simulation_results['trajectories']):
            positions = trajectory['positions']
            ax1.plot(positions[:, 0], positions[:, 1], positions[:, 2], 
                    alpha=0.7, linewidth=1, label=f'Ion {i+1}')
        
        # Plot binding sites
        for i, site in enumerate(binding_sites):
            center = np.array(site['center'])
            ax1.scatter(center[0], center[1], center[2], 
                       s=100, c='red', marker='o', label=f'Binding Site {i+1}')
        
        ax1.set_xlabel('X (m)')
        ax1.set_ylabel('Y (m)')
        ax1.set_zlabel('Z (m)')
        ax1.set_title('Ion Trajectories and Binding Sites')
        ax1.legend()
        
        # 2D projections
        ax2 = fig.add_subplot(222)
        for trajectory in simulation_results['trajectories']:
            positions = trajectory['positions']
            ax2.plot(positions[:, 0], positions[:, 1], alpha=0.7, linewidth=1)
        
        for site in binding_sites:
            center = np.array(site['center'])
            ax2.scatter(center[0], center[1], s=100, c='red', marker='o')
        
        ax2.set_xlabel('X (m)')
        ax2.set_ylabel('Y (m)')
        ax2.set_title('XY Projection')
        ax2.grid(True)
        
        # Collision time distribution
        ax3 = fig.add_subplot(223)
        collision_times = [collision['time'] for collision in simulation_results['all_collisions']]
        if collision_times:
            ax3.hist(collision_times, bins=20, alpha=0.7, edgecolor='black')
            ax3.set_xlabel('Collision Time (s)')
            ax3.set_ylabel('Frequency')
            ax3.set_title('Collision Time Distribution')
            ax3.grid(True)
        
        # Mean squared displacement
        ax4 = fig.add_subplot(224)
        for i, trajectory in enumerate(simulation_results['trajectories']):
            positions = trajectory['positions']
            times = trajectory['times']
            
            # Calculate MSD
            msd = []
            for t_idx in range(len(times)):
                if t_idx == 0:
                    msd.append(0)
                else:
                    displacement = positions[t_idx] - positions[0]
                    msd.append(np.sum(displacement**2))
            
            ax4.plot(times, msd, alpha=0.7, linewidth=2, label=f'Ion {i+1}')
        
        ax4.set_xlabel('Time (s)')
        ax4.set_ylabel('Mean Squared Displacement (m²)')
        ax4.set_title('Mean Squared Displacement')
        ax4.legend()
        ax4.grid(True)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
        
        return fig
    
    def calculate_diffusion_analysis(self, simulation_results):
        """
        Analyze diffusion behavior from simulation results.
        
        Parameters:
        -----------
        simulation_results : dict
            Results from multi-ion simulation
        
        Returns:
        --------
        dict : Diffusion analysis results
        """
        diffusion_coefficients = []
        msd_slopes = []
        
        for trajectory in simulation_results['trajectories']:
            positions = trajectory['positions']
            times = trajectory['times']
            
            # Calculate MSD
            msd = []
            for t_idx in range(len(times)):
                if t_idx == 0:
                    msd.append(0)
                else:
                    displacement = positions[t_idx] - positions[0]
                    msd.append(np.sum(displacement**2))
            
            # Fit MSD to linear function (MSD = 6Dt)
            if len(times) > 1:
                slope = np.polyfit(times, msd, 1)[0]
                msd_slopes.append(slope)
                diffusion_coefficients.append(slope / 6)  # D = slope / 6
        
        return {
            'mean_diffusion_coefficient': np.mean(diffusion_coefficients),
            'std_diffusion_coefficient': np.std(diffusion_coefficients),
            'msd_slopes': msd_slopes,
            'diffusion_coefficients': diffusion_coefficients
        } 