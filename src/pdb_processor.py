"""
PDB Structure Processor for Metalloprotein Analysis

This module handles PDB file parsing, structure analysis, and metal binding site identification.
"""

import numpy as np
from Bio import PDB
from Bio.PDB import *
import warnings
warnings.filterwarnings('ignore')

class PDBProcessor:
    """Process PDB files and identify metal binding sites."""
    
    def __init__(self):
        self.parser = PDB.PDBParser(QUIET=True)
        self.structure = None
        self.binding_sites = []
        
    def load_pdb(self, pdb_file):
        """Load PDB file and extract structure information."""
        try:
            self.structure = self.parser.get_structure('protein', pdb_file)
            return True
        except Exception as e:
            print(f"Error loading PDB file: {e}")
            return False
    
    def identify_metal_binding_sites(self, distance_threshold=3.0):
        """
        Identify potential metal binding sites based on:
        1. Histidine, Cysteine, Aspartate, Glutamate residues
        2. Proximity to other potential ligands
        3. Geometric arrangement
        """
        if not self.structure:
            raise ValueError("No structure loaded. Call load_pdb() first.")
        
        metal_ligands = ['HIS', 'CYS', 'ASP', 'GLU', 'MET', 'TYR']
        binding_sites = []
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in metal_ligands:
                        # Get potential coordinating atoms
                        coord_atoms = self._get_coordinating_atoms(residue)
                        
                        # Find nearby potential ligands
                        nearby_ligands = self._find_nearby_ligands(
                            residue, coord_atoms, distance_threshold
                        )
                        
                        if len(nearby_ligands) >= 2:  # At least 2 coordinating groups
                            binding_site = {
                                'residue': residue,
                                'coordinating_atoms': coord_atoms,
                                'nearby_ligands': nearby_ligands,
                                'center': self._calculate_binding_center(coord_atoms),
                                'coordination_number': len(nearby_ligands)
                            }
                            binding_sites.append(binding_site)
        
        self.binding_sites = binding_sites
        return binding_sites
    
    def _get_coordinating_atoms(self, residue):
        """Get atoms that can coordinate metal ions."""
        coord_atoms = []
        
        if residue.get_resname() == 'HIS':
            # Histidine can coordinate through Nδ1 or Nε2
            for atom in residue:
                if atom.get_name() in ['ND1', 'NE2']:
                    coord_atoms.append(atom)
                    
        elif residue.get_resname() == 'CYS':
            # Cysteine coordinates through Sγ
            for atom in residue:
                if atom.get_name() == 'SG':
                    coord_atoms.append(atom)
                    
        elif residue.get_resname() in ['ASP', 'GLU']:
            # Aspartate/Glutamate coordinate through carboxylate oxygens
            for atom in residue:
                if atom.get_name() in ['OD1', 'OD2', 'OE1', 'OE2']:
                    coord_atoms.append(atom)
        
        return coord_atoms
    
    def _find_nearby_ligands(self, residue, coord_atoms, distance_threshold):
        """Find nearby residues that could coordinate the same metal ion."""
        nearby = []
        
        for model in self.structure:
            for chain in model:
                for other_residue in chain:
                    if other_residue != residue:
                        other_coord_atoms = self._get_coordinating_atoms(other_residue)
                        
                        for atom1 in coord_atoms:
                            for atom2 in other_coord_atoms:
                                distance = atom1 - atom2
                                if distance < distance_threshold:
                                    nearby.append({
                                        'residue': other_residue,
                                        'atom': atom2,
                                        'distance': distance
                                    })
        
        return nearby
    
    def _calculate_binding_center(self, coord_atoms):
        """Calculate the geometric center of coordinating atoms."""
        if not coord_atoms:
            return None
        
        coords = np.array([atom.get_coord() for atom in coord_atoms])
        return np.mean(coords, axis=0)
    
    def get_binding_site_characteristics(self):
        """Calculate characteristics of identified binding sites."""
        characteristics = []
        
        for site in self.binding_sites:
            # Calculate binding pocket volume
            volume = self._estimate_binding_volume(site)
            
            # Calculate electrostatic potential (simplified)
            electrostatic = self._calculate_electrostatic_potential(site)
            
            # Calculate accessibility
            accessibility = self._calculate_accessibility(site)
            
            char = {
                'volume': volume,
                'electrostatic_potential': electrostatic,
                'accessibility': accessibility,
                'coordination_geometry': self._determine_coordination_geometry(site)
            }
            characteristics.append(char)
        
        return characteristics
    
    def _estimate_binding_volume(self, binding_site):
        """Estimate the volume of the binding pocket."""
        # Simplified volume estimation based on coordinating atoms
        coords = np.array([atom.get_coord() for atom in binding_site['coordinating_atoms']])
        
        if len(coords) >= 3:
            # Use convex hull to estimate volume
            from scipy.spatial import ConvexHull
            try:
                hull = ConvexHull(coords)
                return hull.volume
            except:
                # Fallback to simple sphere approximation
                center = np.mean(coords, axis=0)
                max_dist = np.max([np.linalg.norm(coord - center) for coord in coords])
                return (4/3) * np.pi * max_dist**3
        else:
            return 0.0
    
    def _calculate_electrostatic_potential(self, binding_site):
        """Calculate simplified electrostatic potential at binding site."""
        # This is a simplified calculation
        # In practice, you would use Poisson-Boltzmann solvers
        potential = 0.0
        center = binding_site['center']
        
        for atom in binding_site['coordinating_atoms']:
            charge = self._get_atom_charge(atom)
            distance = np.linalg.norm(atom.get_coord() - center)
            if distance > 0:
                potential += charge / distance
        
        return potential
    
    def _get_atom_charge(self, atom):
        """Get approximate atomic charge."""
        # Simplified charge assignment
        charge_dict = {
            'ND1': -0.5, 'NE2': -0.5,  # Histidine
            'SG': -0.5,   # Cysteine
            'OD1': -0.5, 'OD2': -0.5,  # Aspartate
            'OE1': -0.5, 'OE2': -0.5   # Glutamate
        }
        return charge_dict.get(atom.get_name(), 0.0)
    
    def _calculate_accessibility(self, binding_site):
        """Calculate accessibility of binding site."""
        # Simplified accessibility calculation
        center = binding_site['center']
        accessible = True
        
        # Check if binding site is buried or exposed
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom not in binding_site['coordinating_atoms']:
                            distance = np.linalg.norm(atom.get_coord() - center)
                            if distance < 2.0:  # Very close atoms might block access
                                accessible = False
                                break
        
        return accessible
    
    def _determine_coordination_geometry(self, binding_site):
        """Determine the coordination geometry."""
        n_ligands = len(binding_site['coordinating_atoms'])
        
        if n_ligands == 2:
            return "linear"
        elif n_ligands == 3:
            return "trigonal_planar"
        elif n_ligands == 4:
            return "tetrahedral"
        elif n_ligands == 5:
            return "trigonal_bipyramidal"
        elif n_ligands == 6:
            return "octahedral"
        else:
            return "unknown"
    
    def export_binding_sites(self, output_file):
        """Export binding site information to file."""
        with open(output_file, 'w') as f:
            f.write("Binding Site Analysis Results\n")
            f.write("=" * 50 + "\n\n")
            
            for i, site in enumerate(self.binding_sites):
                f.write(f"Binding Site {i+1}:\n")
                f.write(f"  Residue: {site['residue'].get_resname()} {site['residue'].get_id()[1]}\n")
                f.write(f"  Coordinating atoms: {[atom.get_name() for atom in site['coordinating_atoms']]}\n")
                f.write(f"  Coordination number: {site['coordination_number']}\n")
                f.write(f"  Center coordinates: {site['center']}\n")
                f.write(f"  Nearby ligands: {len(site['nearby_ligands'])}\n\n") 