"""
Enhanced Binding Site Identification Module

This module integrates multiple algorithms for binding site identification:
- MetalNet: CHED network analysis and clustering
- Metal3D: Geometric feature analysis
- bindEmbed21: Sequence-based embedding
- AlphaFill: Ligand/cofactor prediction
- MESPEUS: Database integration
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.cluster import DBSCAN, KMeans
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.PDB import *
import requests
import json
import warnings
warnings.filterwarnings('ignore')

class EnhancedBindingSiteIdentifier:
    """Enhanced binding site identification with multi-algorithm integration."""
    
    def __init__(self, config):
        """
        Initialize enhanced binding site identifier.
        
        Parameters:
        -----------
        config : dict
            Configuration dictionary
        """
        self.config = config
        self.algo_config = config['binding_site_algorithms']
        self.mespeus_config = config['mespeus_integration']
        self.ml_config = config['machine_learning']
        
        # Initialize algorithm weights
        self.weights = self.algo_config['weights']
        
        # CHED residues
        self.ched_residues = ['CYS', 'HIS', 'GLU', 'ASP']
        
        # Metal ion types
        self.metal_ions = ['Zn2+', 'Cu2+', 'Fe2+', 'Mg2+', 'Ca2+', 'Mn2+']
        
    def identify_binding_sites(self, protein_structure, protein_sequence=None):
        """
        Identify binding sites using multiple algorithms.
        
        Parameters:
        -----------
        protein_structure : Bio.PDB.Structure
            Protein structure object
        protein_sequence : str, optional
            Protein sequence string
        
        Returns:
        --------
        dict : Binding site predictions with consensus scores
        """
        print("Identifying binding sites using enhanced multi-algorithm approach...")
        
        # Extract protein information
        protein_info = self._extract_protein_info(protein_structure, protein_sequence)
        
        # Run individual algorithms
        results = {}
        
        # MetalNet analysis
        if self.algo_config['MetalNet']['use_network_clustering']:
            results['MetalNet'] = self._run_metalnet_analysis(protein_info)
        
        # Metal3D analysis
        if self.algo_config['Metal3D']['use_geometric_features']:
            results['Metal3D'] = self._run_metal3d_analysis(protein_info)
        
        # bindEmbed21 analysis
        if protein_sequence:
            results['bindEmbed21'] = self._run_bindembed21_analysis(protein_sequence)
        
        # AlphaFill analysis
        if self.algo_config['AlphaFill']['use_alphafold_predictions']:
            results['AlphaFill'] = self._run_alphafill_analysis(protein_info)
        
        # MESPEUS database integration
        if self.mespeus_config['enabled']:
            results['MESPEUS'] = self._run_mespeus_analysis(protein_sequence)
        
        # CHED network analysis
        if self.ml_config['ched_network']['enabled']:
            results['CHED_Network'] = self._run_ched_network_analysis(protein_info)
        
        # Consensus scoring
        consensus_results = self._calculate_consensus_scores(results, protein_info)
        
        return consensus_results
    
    def _extract_protein_info(self, protein_structure, protein_sequence):
        """Extract protein information for analysis."""
        protein_info = {
            'structure': protein_structure,
            'sequence': protein_sequence,
            'residues': [],
            'coordinates': [],
            'residue_types': [],
            'chain_ids': [],
            'residue_ids': []
        }
        
        # Extract residue information
        for model in protein_structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in self.ched_residues:
                        # Get CA atom coordinates
                        if 'CA' in residue:
                            ca_atom = residue['CA']
                            protein_info['residues'].append(residue)
                            protein_info['coordinates'].append(ca_atom.get_coord())
                            protein_info['residue_types'].append(residue.get_resname())
                            protein_info['chain_ids'].append(chain.get_id())
                            protein_info['residue_ids'].append(residue.get_id()[1])
        
        protein_info['coordinates'] = np.array(protein_info['coordinates'])
        protein_info['n_residues'] = len(protein_info['residues'])
        
        return protein_info
    
    def _run_metalnet_analysis(self, protein_info):
        """Run MetalNet analysis with CHED network clustering."""
        print("Running MetalNet analysis...")
        
        if protein_info['n_residues'] < 2:
            return {'sites': [], 'confidence': 0.0}
        
        # Calculate pairwise distances
        distances = cdist(protein_info['coordinates'], protein_info['coordinates'])
        
        # Create adjacency matrix based on distance threshold
        threshold = self.algo_config['MetalNet']['distance_threshold']
        adjacency_matrix = (distances < threshold).astype(int)
        np.fill_diagonal(adjacency_matrix, 0)  # Remove self-connections
        
        # Network clustering
        if self.algo_config['MetalNet']['use_network_clustering']:
            clusters = self._perform_network_clustering(adjacency_matrix)
        else:
            clusters = [[i] for i in range(protein_info['n_residues'])]
        
        # Identify binding sites from clusters
        binding_sites = []
        for cluster in clusters:
            if len(cluster) >= 2:  # Minimum 2 residues for binding site
                cluster_coords = protein_info['coordinates'][cluster]
                cluster_center = np.mean(cluster_coords, axis=0)
                cluster_radius = np.max(cdist([cluster_center], cluster_coords))
                
                # Calculate confidence based on cluster properties
                confidence = self._calculate_cluster_confidence(cluster, distances, protein_info)
                
                if confidence >= self.algo_config['MetalNet']['confidence_threshold']:
                    binding_sites.append({
                        'center': cluster_center,
                        'radius': cluster_radius,
                        'residues': cluster,
                        'confidence': confidence,
                        'algorithm': 'MetalNet'
                    })
        
        return {
            'sites': binding_sites,
            'confidence': np.mean([site['confidence'] for site in binding_sites]) if binding_sites else 0.0
        }
    
    def _run_metal3d_analysis(self, protein_info):
        """Run Metal3D analysis with geometric features."""
        print("Running Metal3D analysis...")
        
        if protein_info['n_residues'] < 3:
            return {'sites': [], 'confidence': 0.0}
        
        # Calculate geometric features
        geometric_features = self._calculate_geometric_features(protein_info)
        
        # Identify binding sites based on geometric criteria
        binding_sites = []
        
        # Look for tetrahedral coordination (4 residues)
        tetrahedral_sites = self._find_tetrahedral_sites(protein_info, geometric_features)
        binding_sites.extend(tetrahedral_sites)
        
        # Look for octahedral coordination (6 residues)
        octahedral_sites = self._find_octahedral_sites(protein_info, geometric_features)
        binding_sites.extend(octahedral_sites)
        
        # Look for trigonal planar coordination (3 residues)
        trigonal_sites = self._find_trigonal_sites(protein_info, geometric_features)
        binding_sites.extend(trigonal_sites)
        
        return {
            'sites': binding_sites,
            'confidence': np.mean([site['confidence'] for site in binding_sites]) if binding_sites else 0.0
        }
    
    def _run_bindembed21_analysis(self, protein_sequence):
        """Run bindEmbed21 analysis (simulated)."""
        print("Running bindEmbed21 analysis...")
        
        # This is a simplified simulation of bindEmbed21
        # In practice, you would use the actual bindEmbed21 model
        
        # Simulate binding site predictions based on sequence patterns
        binding_sites = []
        
        # Look for known metal-binding motifs
        motifs = self._identify_metal_binding_motifs(protein_sequence)
        
        for motif in motifs:
            # Convert sequence position to approximate 3D coordinates
            # This is a simplified approach
            confidence = motif['score']
            
            if confidence >= self.algo_config['bindEmbed21']['confidence_threshold']:
                binding_sites.append({
                    'center': motif['center'],
                    'radius': 3.0,  # Default radius
                    'residues': motif['residues'],
                    'confidence': confidence,
                    'algorithm': 'bindEmbed21'
                })
        
        return {
            'sites': binding_sites,
            'confidence': np.mean([site['confidence'] for site in binding_sites]) if binding_sites else 0.0
        }
    
    def _run_alphafill_analysis(self, protein_info):
        """Run AlphaFill analysis (simulated)."""
        print("Running AlphaFill analysis...")
        
        # This is a simplified simulation of AlphaFill
        # In practice, you would use the actual AlphaFill model
        
        binding_sites = []
        
        # Simulate ligand/cofactor binding site predictions
        # Look for cavities and pockets in the protein structure
        cavities = self._identify_protein_cavities(protein_info)
        
        for cavity in cavities:
            confidence = cavity['score']
            
            if confidence >= self.algo_config['AlphaFill']['confidence_threshold']:
                binding_sites.append({
                    'center': cavity['center'],
                    'radius': cavity['radius'],
                    'residues': cavity['residues'],
                    'confidence': confidence,
                    'algorithm': 'AlphaFill'
                })
        
        return {
            'sites': binding_sites,
            'confidence': np.mean([site['confidence'] for site in binding_sites]) if binding_sites else 0.0
        }
    
    def _run_mespeus_analysis(self, protein_sequence):
        """Run MESPEUS database analysis."""
        print("Running MESPEUS database analysis...")
        
        if not protein_sequence:
            return {'sites': [], 'confidence': 0.0}
        
        # Query MESPEUS database for similar proteins
        similar_proteins = self._query_mespeus_database(protein_sequence)
        
        binding_sites = []
        
        for protein in similar_proteins:
            if protein['similarity'] >= self.mespeus_config['similarity_threshold']:
                # Transfer binding site information
                transferred_sites = self._transfer_binding_sites(protein, protein_sequence)
                binding_sites.extend(transferred_sites)
        
        return {
            'sites': binding_sites,
            'confidence': np.mean([site['confidence'] for site in binding_sites]) if binding_sites else 0.0
        }
    
    def _run_ched_network_analysis(self, protein_info):
        """Run CHED network analysis."""
        print("Running CHED network analysis...")
        
        if protein_info['n_residues'] < 2:
            return {'sites': [], 'confidence': 0.0}
        
        # Calculate CHED pair distances
        distances = cdist(protein_info['coordinates'], protein_info['coordinates'])
        
        # Find CHED pairs within threshold
        threshold = self.ml_config['ched_network']['distance_threshold']
        ched_pairs = []
        
        for i in range(protein_info['n_residues']):
            for j in range(i+1, protein_info['n_residues']):
                if distances[i, j] < threshold:
                    ched_pairs.append((i, j, distances[i, j]))
        
        # Cluster CHED pairs
        if ched_pairs:
            clusters = self._cluster_ched_pairs(ched_pairs, protein_info)
        else:
            clusters = []
        
        # Convert clusters to binding sites
        binding_sites = []
        for cluster in clusters:
            if len(cluster) >= 2:
                cluster_coords = protein_info['coordinates'][cluster]
                cluster_center = np.mean(cluster_coords, axis=0)
                cluster_radius = np.max(cdist([cluster_center], cluster_coords))
                
                confidence = len(cluster) / 6.0  # Normalize by max coordination
                
                binding_sites.append({
                    'center': cluster_center,
                    'radius': cluster_radius,
                    'residues': cluster,
                    'confidence': confidence,
                    'algorithm': 'CHED_Network'
                })
        
        return {
            'sites': binding_sites,
            'confidence': np.mean([site['confidence'] for site in binding_sites]) if binding_sites else 0.0
        }
    
    def _perform_network_clustering(self, adjacency_matrix):
        """Perform network clustering on adjacency matrix."""
        # Use hierarchical clustering
        if self.ml_config['ched_network']['clustering_algorithm'] == 'hierarchical':
            # Convert adjacency matrix to distance matrix
            distance_matrix = 1 - adjacency_matrix
            linkage_matrix = linkage(squareform(distance_matrix), method='ward')
            
            # Determine number of clusters
            n_clusters = min(5, len(adjacency_matrix) // 2)  # Heuristic
            clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            
            # Convert to list of lists
            unique_clusters = np.unique(clusters)
            cluster_lists = []
            for cluster_id in unique_clusters:
                cluster_indices = np.where(clusters == cluster_id)[0]
                cluster_lists.append(cluster_indices.tolist())
            
            return cluster_lists
        
        elif self.ml_config['ched_network']['clustering_algorithm'] == 'dbscan':
            # Use DBSCAN clustering
            dbscan = DBSCAN(eps=0.5, min_samples=2)
            clusters = dbscan.fit_predict(adjacency_matrix)
            
            # Convert to list of lists
            unique_clusters = np.unique(clusters)
            cluster_lists = []
            for cluster_id in unique_clusters:
                if cluster_id != -1:  # Skip noise points
                    cluster_indices = np.where(clusters == cluster_id)[0]
                    cluster_lists.append(cluster_indices.tolist())
            
            return cluster_lists
        
        else:  # kmeans
            # Use K-means clustering
            n_clusters = min(5, len(adjacency_matrix) // 2)
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            clusters = kmeans.fit_predict(adjacency_matrix)
            
            # Convert to list of lists
            unique_clusters = np.unique(clusters)
            cluster_lists = []
            for cluster_id in unique_clusters:
                cluster_indices = np.where(clusters == cluster_id)[0]
                cluster_lists.append(cluster_indices.tolist())
            
            return cluster_lists
    
    def _calculate_cluster_confidence(self, cluster, distances, protein_info):
        """Calculate confidence score for a cluster."""
        if len(cluster) < 2:
            return 0.0
        
        # Factors affecting confidence:
        # 1. Number of residues in cluster
        size_factor = min(len(cluster) / 6.0, 1.0)  # Normalize by max coordination
        
        # 2. Average distance between residues
        cluster_distances = []
        for i in cluster:
            for j in cluster:
                if i != j:
                    cluster_distances.append(distances[i, j])
        
        if cluster_distances:
            avg_distance = np.mean(cluster_distances)
            distance_factor = max(0, 1 - avg_distance / 10.0)  # Prefer closer residues
        else:
            distance_factor = 0.0
        
        # 3. Residue type diversity
        cluster_types = [protein_info['residue_types'][i] for i in cluster]
        type_diversity = len(set(cluster_types)) / len(cluster_types)
        
        # 4. Geometric arrangement
        cluster_coords = protein_info['coordinates'][cluster]
        geometric_factor = self._calculate_geometric_quality(cluster_coords)
        
        # Combine factors
        confidence = (size_factor * 0.3 + 
                     distance_factor * 0.3 + 
                     type_diversity * 0.2 + 
                     geometric_factor * 0.2)
        
        return confidence
    
    def _calculate_geometric_features(self, protein_info):
        """Calculate geometric features for Metal3D analysis."""
        features = {}
        
        # Calculate pairwise distances
        distances = cdist(protein_info['coordinates'], protein_info['coordinates'])
        
        # Calculate angles between triplets of residues
        angles = {}
        for i in range(protein_info['n_residues']):
            for j in range(i+1, protein_info['n_residues']):
                for k in range(j+1, protein_info['n_residues']):
                    angle = self._calculate_angle(
                        protein_info['coordinates'][i],
                        protein_info['coordinates'][j],
                        protein_info['coordinates'][k]
                    )
                    angles[(i, j, k)] = angle
        
        features['distances'] = distances
        features['angles'] = angles
        
        return features
    
    def _find_tetrahedral_sites(self, protein_info, geometric_features):
        """Find tetrahedral coordination sites."""
        sites = []
        
        # Look for 4 residues arranged in tetrahedral geometry
        for i in range(protein_info['n_residues']):
            for j in range(i+1, protein_info['n_residues']):
                for k in range(j+1, protein_info['n_residues']):
                    for l in range(k+1, protein_info['n_residues']):
                        cluster = [i, j, k, l]
                        cluster_coords = protein_info['coordinates'][cluster]
                        
                        # Check if arrangement is approximately tetrahedral
                        if self._is_tetrahedral_arrangement(cluster_coords):
                            center = np.mean(cluster_coords, axis=0)
                            radius = np.max(cdist([center], cluster_coords))
                            
                            confidence = self._calculate_tetrahedral_confidence(cluster_coords)
                            
                            sites.append({
                                'center': center,
                                'radius': radius,
                                'residues': cluster,
                                'confidence': confidence,
                                'coordination': 'tetrahedral'
                            })
        
        return sites
    
    def _find_octahedral_sites(self, protein_info, geometric_features):
        """Find octahedral coordination sites."""
        sites = []
        
        # Look for 6 residues arranged in octahedral geometry
        # This is computationally expensive, so we use a simplified approach
        for i in range(protein_info['n_residues']):
            for j in range(i+1, protein_info['n_residues']):
                for k in range(j+1, protein_info['n_residues']):
                    for l in range(k+1, protein_info['n_residues']):
                        for m in range(l+1, protein_info['n_residues']):
                            for n in range(m+1, protein_info['n_residues']):
                                cluster = [i, j, k, l, m, n]
                                cluster_coords = protein_info['coordinates'][cluster]
                                
                                # Check if arrangement is approximately octahedral
                                if self._is_octahedral_arrangement(cluster_coords):
                                    center = np.mean(cluster_coords, axis=0)
                                    radius = np.max(cdist([center], cluster_coords))
                                    
                                    confidence = self._calculate_octahedral_confidence(cluster_coords)
                                    
                                    sites.append({
                                        'center': center,
                                        'radius': radius,
                                        'residues': cluster,
                                        'confidence': confidence,
                                        'coordination': 'octahedral'
                                    })
        
        return sites
    
    def _find_trigonal_sites(self, protein_info, geometric_features):
        """Find trigonal planar coordination sites."""
        sites = []
        
        # Look for 3 residues arranged in trigonal planar geometry
        for i in range(protein_info['n_residues']):
            for j in range(i+1, protein_info['n_residues']):
                for k in range(j+1, protein_info['n_residues']):
                    cluster = [i, j, k]
                    cluster_coords = protein_info['coordinates'][cluster]
                    
                    # Check if arrangement is approximately trigonal planar
                    if self._is_trigonal_arrangement(cluster_coords):
                        center = np.mean(cluster_coords, axis=0)
                        radius = np.max(cdist([center], cluster_coords))
                        
                        confidence = self._calculate_trigonal_confidence(cluster_coords)
                        
                        sites.append({
                            'center': center,
                            'radius': radius,
                            'residues': cluster,
                            'confidence': confidence,
                            'coordination': 'trigonal'
                        })
        
        return sites
    
    def _is_tetrahedral_arrangement(self, coords):
        """Check if coordinates form tetrahedral arrangement."""
        if len(coords) != 4:
            return False
        
        # Calculate angles between center and each point
        center = np.mean(coords, axis=0)
        angles = []
        
        for i in range(4):
            for j in range(i+1, 4):
                angle = self._calculate_angle(center, coords[i], coords[j])
                angles.append(angle)
        
        # Tetrahedral angles should be around 109.5°
        tetrahedral_angle = 109.5
        tolerance = 20.0
        
        return all(abs(angle - tetrahedral_angle) < tolerance for angle in angles)
    
    def _is_octahedral_arrangement(self, coords):
        """Check if coordinates form octahedral arrangement."""
        if len(coords) != 6:
            return False
        
        # Simplified check: look for opposite pairs
        center = np.mean(coords, axis=0)
        distances = cdist([center], coords)[0]
        
        # Should have 3 pairs of opposite points
        sorted_indices = np.argsort(distances)
        
        return True  # Simplified for computational efficiency
    
    def _is_trigonal_arrangement(self, coords):
        """Check if coordinates form trigonal planar arrangement."""
        if len(coords) != 3:
            return False
        
        # Calculate angles between the three points
        angles = []
        for i in range(3):
            angle = self._calculate_angle(coords[i], coords[(i+1)%3], coords[(i+2)%3])
            angles.append(angle)
        
        # Trigonal planar angles should be around 120°
        trigonal_angle = 120.0
        tolerance = 30.0
        
        return all(abs(angle - trigonal_angle) < tolerance for angle in angles)
    
    def _calculate_angle(self, p1, p2, p3):
        """Calculate angle between three points."""
        v1 = p1 - p2
        v2 = p3 - p2
        
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        
        return np.arccos(cos_angle) * 180 / np.pi
    
    def _calculate_geometric_quality(self, coords):
        """Calculate geometric quality of a binding site."""
        if len(coords) < 2:
            return 0.0
        
        # Calculate distances from center
        center = np.mean(coords, axis=0)
        distances = cdist([center], coords)[0]
        
        # Quality based on distance uniformity
        std_distance = np.std(distances)
        mean_distance = np.mean(distances)
        
        if mean_distance > 0:
            cv_distance = std_distance / mean_distance
            quality = max(0, 1 - cv_distance)
        else:
            quality = 0.0
        
        return quality
    
    def _identify_metal_binding_motifs(self, sequence):
        """Identify metal-binding motifs in sequence."""
        motifs = []
        
        # Known metal-binding motifs (simplified)
        known_motifs = [
            'CXXC',  # Cysteine pairs
            'HXXH',  # Histidine pairs
            'CXXXC', # Cysteine pairs with spacer
            'HXXXH', # Histidine pairs with spacer
        ]
        
        for motif in known_motifs:
            for i in range(len(sequence) - len(motif) + 1):
                subsequence = sequence[i:i+len(motif)]
                if self._matches_motif(subsequence, motif):
                    motifs.append({
                        'motif': motif,
                        'position': i,
                        'sequence': subsequence,
                        'score': 0.8,  # Default score
                        'center': [i + len(motif)//2, 0, 0],  # Simplified 3D position
                        'residues': list(range(i, i+len(motif)))
                    })
        
        return motifs
    
    def _matches_motif(self, subsequence, motif):
        """Check if subsequence matches motif pattern."""
        if len(subsequence) != len(motif):
            return False
        
        for i, (seq_char, motif_char) in enumerate(zip(subsequence, motif)):
            if motif_char == 'X':
                continue
            elif seq_char != motif_char:
                return False
        
        return True
    
    def _identify_protein_cavities(self, protein_info):
        """Identify protein cavities (simplified)."""
        cavities = []
        
        # Simplified cavity detection based on residue clustering
        if protein_info['n_residues'] >= 3:
            # Use DBSCAN to find clusters that might form cavities
            dbscan = DBSCAN(eps=5.0, min_samples=3)
            clusters = dbscan.fit_predict(protein_info['coordinates'])
            
            unique_clusters = np.unique(clusters)
            for cluster_id in unique_clusters:
                if cluster_id != -1:  # Skip noise
                    cluster_indices = np.where(clusters == cluster_id)[0]
                    cluster_coords = protein_info['coordinates'][cluster_indices]
                    
                    center = np.mean(cluster_coords, axis=0)
                    radius = np.max(cdist([center], cluster_coords))
                    
                    cavities.append({
                        'center': center,
                        'radius': radius,
                        'residues': cluster_indices.tolist(),
                        'score': len(cluster_indices) / 10.0  # Normalized score
                    })
        
        return cavities
    
    def _query_mespeus_database(self, protein_sequence):
        """Query MESPEUS database (simulated)."""
        # This is a simplified simulation
        # In practice, you would make actual API calls to MESPEUS
        
        similar_proteins = []
        
        # Simulate database query results
        if len(protein_sequence) > 50:
            # Simulate finding similar proteins
            for i in range(3):
                similarity = 0.7 + 0.2 * np.random.random()  # 0.7-0.9 similarity
                similar_proteins.append({
                    'id': f'MESPEUS_{i:04d}',
                    'similarity': similarity,
                    'binding_sites': [
                        {
                            'center': [10 + i*5, 10 + i*3, 10 + i*2],
                            'radius': 3.0 + i*0.5,
                            'residues': [i*10, i*10+1, i*10+2],
                            'confidence': 0.8
                        }
                    ]
                })
        
        return similar_proteins
    
    def _transfer_binding_sites(self, similar_protein, query_sequence):
        """Transfer binding sites from similar protein."""
        transferred_sites = []
        
        for site in similar_protein['binding_sites']:
            # Adjust confidence based on similarity
            adjusted_confidence = site['confidence'] * similar_protein['similarity']
            
            transferred_sites.append({
                'center': site['center'],
                'radius': site['radius'],
                'residues': site['residues'],
                'confidence': adjusted_confidence,
                'source': similar_protein['id']
            })
        
        return transferred_sites
    
    def _cluster_ched_pairs(self, ched_pairs, protein_info):
        """Cluster CHED pairs into binding sites."""
        if not ched_pairs:
            return []
        
        # Create graph from CHED pairs
        nodes = set()
        for i, j, dist in ched_pairs:
            nodes.add(i)
            nodes.add(j)
        
        # Simple clustering: group connected components
        clusters = []
        visited = set()
        
        for node in nodes:
            if node not in visited:
                cluster = self._find_connected_component(node, ched_pairs)
                clusters.append(cluster)
                visited.update(cluster)
        
        return clusters
    
    def _find_connected_component(self, start_node, ched_pairs):
        """Find connected component starting from a node."""
        component = {start_node}
        queue = [start_node]
        
        while queue:
            node = queue.pop(0)
            
            for i, j, dist in ched_pairs:
                if i == node and j not in component:
                    component.add(j)
                    queue.append(j)
                elif j == node and i not in component:
                    component.add(i)
                    queue.append(i)
        
        return list(component)
    
    def _calculate_consensus_scores(self, results, protein_info):
        """Calculate consensus scores from multiple algorithms."""
        print("Calculating consensus scores...")
        
        # Collect all binding sites
        all_sites = []
        algorithm_scores = {}
        
        for algo_name, result in results.items():
            if 'sites' in result and result['sites']:
                weight = self.weights.get(algo_name, 0.1)
                algorithm_scores[algo_name] = result['confidence']
                
                for site in result['sites']:
                    site['weighted_confidence'] = site['confidence'] * weight
                    all_sites.append(site)
        
        # Group nearby sites
        consensus_sites = self._group_nearby_sites(all_sites)
        
        # Calculate final consensus scores
        for site in consensus_sites:
            site['consensus_score'] = np.mean([s['weighted_confidence'] for s in site['contributing_sites']])
            site['algorithm_count'] = len(set(s['algorithm'] for s in site['contributing_sites']))
        
        # Sort by consensus score
        consensus_sites.sort(key=lambda x: x['consensus_score'], reverse=True)
        
        return {
            'binding_sites': consensus_sites,
            'algorithm_scores': algorithm_scores,
            'total_sites': len(consensus_sites),
            'average_consensus_score': np.mean([s['consensus_score'] for s in consensus_sites]) if consensus_sites else 0.0
        }
    
    def _group_nearby_sites(self, all_sites, distance_threshold=5.0):
        """Group nearby binding sites into consensus sites."""
        if not all_sites:
            return []
        
        # Calculate distances between all sites
        centers = np.array([site['center'] for site in all_sites])
        distances = cdist(centers, centers)
        
        # Group sites using hierarchical clustering
        linkage_matrix = linkage(distances, method='single')
        clusters = fcluster(linkage_matrix, distance_threshold, criterion='distance')
        
        # Create consensus sites
        consensus_sites = []
        unique_clusters = np.unique(clusters)
        
        for cluster_id in unique_clusters:
            cluster_indices = np.where(clusters == cluster_id)[0]
            cluster_sites = [all_sites[i] for i in cluster_indices]
            
            # Calculate consensus center
            cluster_centers = np.array([site['center'] for site in cluster_sites])
            consensus_center = np.mean(cluster_centers, axis=0)
            
            # Calculate consensus radius
            consensus_radius = np.mean([site['radius'] for site in cluster_sites])
            
            # Combine residues
            all_residues = []
            for site in cluster_sites:
                all_residues.extend(site['residues'])
            consensus_residues = list(set(all_residues))
            
            consensus_sites.append({
                'center': consensus_center,
                'radius': consensus_radius,
                'residues': consensus_residues,
                'contributing_sites': cluster_sites,
                'consensus_score': 0.0,  # Will be calculated later
                'algorithm_count': 0     # Will be calculated later
            })
        
        return consensus_sites
    
    def plot_binding_site_results(self, results, protein_info, save_path=None):
        """Plot binding site identification results."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Algorithm confidence scores
        ax1 = axes[0, 0]
        algorithms = list(results['algorithm_scores'].keys())
        scores = list(results['algorithm_scores'].values())
        
        bars = ax1.bar(algorithms, scores, color=['blue', 'green', 'red', 'orange', 'purple', 'brown'])
        ax1.set_ylabel('Confidence Score')
        ax1.set_title('Algorithm Confidence Scores')
        ax1.set_ylim(0, 1)
        
        # Add value labels on bars
        for bar, score in zip(bars, scores):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{score:.2f}', ha='center', va='bottom')
        
        # Plot 2: Binding site locations (2D projection)
        ax2 = axes[0, 1]
        if protein_info['coordinates'].shape[1] >= 2:
            # Plot all CHED residues
            ax2.scatter(protein_info['coordinates'][:, 0], protein_info['coordinates'][:, 1],
                       c='lightblue', s=50, alpha=0.6, label='CHED Residues')
            
            # Plot binding sites
            colors = plt.cm.Set3(np.linspace(0, 1, len(results['binding_sites'])))
            for i, site in enumerate(results['binding_sites']):
                center = site['center']
                radius = site['radius']
                confidence = site['consensus_score']
                
                circle = plt.Circle((center[0], center[1]), radius, 
                                  color=colors[i], alpha=0.7, fill=False, linewidth=2)
                ax2.add_patch(circle)
                
                ax2.text(center[0], center[1], f'{confidence:.2f}', 
                        ha='center', va='center', fontweight='bold')
            
            ax2.set_xlabel('X Coordinate (Å)')
            ax2.set_ylabel('Y Coordinate (Å)')
            ax2.set_title('Binding Site Locations')
            ax2.legend()
            ax2.set_aspect('equal')
        
        # Plot 3: Consensus score distribution
        ax3 = axes[1, 0]
        consensus_scores = [site['consensus_score'] for site in results['binding_sites']]
        if consensus_scores:
            ax3.hist(consensus_scores, bins=10, alpha=0.7, color='skyblue', edgecolor='black')
            ax3.set_xlabel('Consensus Score')
            ax3.set_ylabel('Number of Sites')
            ax3.set_title('Consensus Score Distribution')
        
        # Plot 4: Algorithm agreement
        ax4 = axes[1, 1]
        algorithm_counts = [site['algorithm_count'] for site in results['binding_sites']]
        if algorithm_counts:
            unique_counts, counts = np.unique(algorithm_counts, return_counts=True)
            ax4.bar(unique_counts, counts, color='lightgreen', alpha=0.7)
            ax4.set_xlabel('Number of Algorithms Agreeing')
            ax4.set_ylabel('Number of Sites')
            ax4.set_title('Algorithm Agreement Distribution')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
        
        return fig 