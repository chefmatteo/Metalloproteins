# Enhanced Metalloprotein Binding Efficiency Prediction Pipeline

A comprehensive computational framework for predicting metal ion binding efficiency in metalloproteins under realistic environmental conditions.

## 🚀 Key Features

### 🌡️ Environmental Parameter Coupling
- **Temperature effects**: Arrhenius behavior and heat generation from binding reactions
- **pH dynamics**: Protonation effects on metal-binding residues
- **Pressure effects**: Activation volume effects on binding kinetics
- **Redox coupling**: Oxidation state effects on binding affinity

### 🔬 Multi-Algorithm Integration
- **MetalNet**: CHED network analysis and clustering
- **Metal3D**: Geometric coordination analysis
- **bindEmbed21**: Sequence-based binding site prediction
- **AlphaFill**: Ligand/cofactor binding site prediction
- **MESPEUS**: Database integration for experimental validation
- **CHED Network**: Machine learning-based clustering

### 📊 Spatial Discretization
- **1000-cube model**: 10×10×10 spatial grid for high-resolution analysis
- **Finite difference methods**: Proper discretization of coupled PDEs
- **Parallel computing**: Multi-processor support for efficient computation

### 🎨 Advanced Visualization
- **PyMOL integration**: 3D structure visualization with binding sites
- **RF Diffusion**: Ion diffusion and binding process animation
- **Environmental mapping**: Color-coded parameter visualization
- **Spatial plots**: 2D slices of 3D environmental distributions

## 📋 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Enhanced Features](#enhanced-features)
- [Configuration](#configuration)
- [Usage Examples](#usage-examples)
- [Output and Results](#output-and-results)
- [Validation](#validation)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## 🛠️ Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Install Dependencies
```bash
# Clone the repository
git clone https://github.com/yourusername/metalloproteins.git
cd metalloproteins

# Install dependencies
pip install -r requirements.txt
```

### Optional Dependencies
For advanced features, you can install optional packages:

```bash
# Deep learning support (for future enhancements)
pip install torch torch-geometric transformers

# Advanced visualization
pip install pymol-open-source vtk

# High-performance computing
pip install mpi4py cupy-cuda11x
```

## 🚀 Quick Start

### Basic Usage
```python
from src.enhanced_main import EnhancedMetalloproteinPipeline

# Initialize the enhanced pipeline
pipeline = EnhancedMetalloproteinPipeline()

# Run enhanced analysis
results = pipeline.run_enhanced_analysis(
    pdb_file="your_protein.pdb",
    protein_sequence="CHEDCH",  # Cysteine, Histidine, Glutamate, Aspartate, Cysteine, Histidine
    metal_ions=['Zn2+', 'Cu2+', 'Fe2+'],
    initial_concentrations={'Zn2+': 1e-6, 'Cu2+': 1e-6, 'Fe2+': 1e-6},
    time_span=(0, 1000),
    save_results=True
)

# Access results
print(f"Binding efficiency: {results['efficiency_results']['overall_efficiency']:.3f}")
print(f"Binding sites identified: {results['binding_sites']['total_sites']}")
```

### Run Enhanced Example
```bash
python enhanced_example_usage.py
```

## 🔬 Enhanced Features

### 1. Environmental Parameter Analysis

The pipeline incorporates realistic environmental conditions:

```python
# Test different environmental conditions
conditions = [
    {'name': 'Standard', 'T': 298.15, 'pH': 7.0, 'P': 1.0, 'Eh': 0.0},
    {'name': 'High Temperature', 'T': 323.15, 'pH': 7.0, 'P': 1.0, 'Eh': 0.0},
    {'name': 'Low pH', 'T': 298.15, 'pH': 5.0, 'P': 1.0, 'Eh': 0.0},
    {'name': 'High Pressure', 'T': 298.15, 'pH': 7.0, 'P': 10.0, 'Eh': 0.0},
    {'name': 'Oxidizing', 'T': 298.15, 'pH': 7.0, 'P': 1.0, 'Eh': 0.3}
]

for condition in conditions:
    # Update configuration
    pipeline.config['environmental_conditions']['temperature']['initial'] = condition['T']
    pipeline.config['environmental_conditions']['pH']['initial'] = condition['pH']
    pipeline.config['environmental_conditions']['pressure']['initial'] = condition['P']
    pipeline.config['environmental_conditions']['redox_potential']['initial'] = condition['Eh']
    
    # Run analysis
    results = pipeline.run_enhanced_analysis(...)
```

### 2. Multi-Algorithm Binding Site Identification

Consensus scoring from multiple algorithms:

```python
# Algorithm weights can be customized
pipeline.config['binding_site_algorithms']['weights'] = {
    'MetalNet': 0.3,
    'Metal3D': 0.25,
    'bindEmbed21': 0.25,
    'AlphaFill': 0.2
}

# Results include consensus scores
binding_sites = results['binding_sites']['binding_sites']
for site in binding_sites:
    print(f"Consensus score: {site['consensus_score']:.3f}")
    print(f"Algorithms agreeing: {site['algorithm_count']}")
```

### 3. Spatial Discretization

1000-cube spatial model for high-resolution analysis:

```python
# Access spatial data
spatial_data = results['efficiency_results']['spatial_data']
temperature_dist = spatial_data['temperature'].reshape(10, 10, 10)
ph_dist = spatial_data['pH'].reshape(10, 10, 10)

# Analyze spatial gradients
print(f"Temperature gradient: {np.gradient(temperature_dist)}")
print(f"pH gradient: {np.gradient(ph_dist)}")
```

## ⚙️ Configuration

### Enhanced Configuration File
```yaml
# config/enhanced_config.yaml
environmental_conditions:
  temperature:
    initial: 298.15  # Kelvin (25°C)
    range: [273.15, 373.15]  # 0°C to 100°C
    heat_transfer:
      thermal_conductivity: 0.6  # W/m·K
      specific_heat: 4186  # J/kg·K
      density: 997  # kg/m³
  
  pH:
    initial: 7.0
    range: [4.0, 10.0]
    buffer_system: "phosphate"
    buffer_concentration: 0.1  # M
  
  pressure:
    initial: 1.0  # atm
    range: [0.1, 100.0]  # atm
    compressibility:
      thermal: 2.1e-4  # 1/K
      chemical: 1e-6  # m³/mol
  
  redox_potential:
    initial: 0.0  # V vs SHE
    range: [-0.5, 0.5]  # V
    buffer_system: "glutathione"
    buffer_concentration: 1e-3  # M

spatial_discretization:
  chamber:
    dimensions: [10e-6, 10e-6, 10e-6]  # 10×10×10 μm³
    grid_size: [10, 10, 10]  # 1000 cubes total
    boundary_conditions: "periodic"

binding_site_algorithms:
  weights:
    MetalNet: 0.3
    Metal3D: 0.25
    bindEmbed21: 0.25
    AlphaFill: 0.2
```

## 📊 Usage Examples

### Example 1: Basic Analysis
```python
# Simple metalloprotein analysis
results = pipeline.run_enhanced_analysis(
    pdb_file="metallothionein.pdb",
    metal_ions=['Zn2+'],
    initial_concentrations={'Zn2+': 1e-6},
    time_span=(0, 500)
)
```

### Example 2: Environmental Effects Study
```python
# Study temperature effects on binding
temperatures = [273.15, 298.15, 323.15, 348.15]  # 0°C, 25°C, 50°C, 75°C
efficiencies = []

for T in temperatures:
    pipeline.config['environmental_conditions']['temperature']['initial'] = T
    pipeline.enhanced_kinetics = EnhancedBindingKinetics(pipeline.config)
    
    results = pipeline.run_enhanced_analysis(...)
    efficiencies.append(results['efficiency_results']['overall_efficiency'])

# Plot temperature dependence
plt.plot([T-273.15 for T in temperatures], efficiencies, 'o-')
plt.xlabel('Temperature (°C)')
plt.ylabel('Binding Efficiency')
plt.title('Temperature Dependence of Binding Efficiency')
plt.show()
```

### Example 3: Multi-Metal Competition
```python
# Study competition between different metal ions
results = pipeline.run_enhanced_analysis(
    pdb_file="protein.pdb",
    metal_ions=['Zn2+', 'Cu2+', 'Fe2+'],
    initial_concentrations={
        'Zn2+': 1e-6,
        'Cu2+': 1e-6,
        'Fe2+': 1e-6
    },
    time_span=(0, 1000)
)

# Analyze competition
final_bound = results['efficiency_results']['final_bound_concentrations']
for metal, concentration in final_bound.items():
    print(f"{metal}: {concentration:.2e} M")
```

## 📈 Output and Results

### Generated Files
The pipeline generates comprehensive output:

```
output/enhanced/
├── plots/
│   ├── binding_site_identification.png
│   ├── enhanced_kinetics.png
│   ├── environmental_analysis.png
│   ├── spatial_distributions.png
│   ├── brownian_motion.png
│   └── environmental_effects.png
├── data/
│   ├── binding_site_results.json
│   ├── efficiency_results.json
│   └── brownian_results.json
├── structures/
│   └── pymol_script.pml
├── animations/
│   └── rf_diffusion_script.py
└── comprehensive_report.md
```

### Key Metrics
- **Binding efficiency**: Overall and per-cube analysis
- **Environmental analysis**: Temperature, pH, pressure, redox statistics
- **Algorithm performance**: Individual and consensus scores
- **Spatial distribution**: Parameter gradients and distributions

### Results Structure
```python
results = {
    'binding_sites': {
        'total_sites': 3,
        'average_consensus_score': 0.75,
        'algorithm_scores': {'MetalNet': 0.8, 'Metal3D': 0.7, ...},
        'binding_sites': [
            {
                'center': [x, y, z],
                'radius': 3.2,
                'residues': [1, 2, 3],
                'consensus_score': 0.85,
                'algorithm_count': 4
            }
        ]
    },
    'efficiency_results': {
        'overall_efficiency': 0.65,
        'environmental_analysis': {
            'temperature': {'mean': 298.5, 'std': 0.1, 'range': [298.0, 299.0]},
            'pH': {'mean': 7.1, 'std': 0.05, 'range': [7.0, 7.2]},
            'pressure': {'mean': 1.0, 'std': 0.01, 'range': [0.99, 1.01]},
            'redox_potential': {'mean': 0.0, 'std': 0.01, 'range': [-0.01, 0.01]}
        },
        'final_free_concentrations': {'Zn2+': 3.5e-7},
        'final_bound_concentrations': {'Zn2+': 6.5e-7}
    },
    'brownian_results': {
        'diffusion_coefficients': {'Zn2+': 7.0e-10},
        'binding_probability': {'Zn2+': [0.1, 0.2, ..., 0.65]},
        'msd': {'Zn2+': [0.0, 1.0, ..., 100.0]}
    }
}
```

## ✅ Validation

### Cross-Validation
```python
# Perform cross-validation on binding site prediction
from sklearn.model_selection import KFold

kf = KFold(n_splits=5, shuffle=True, random_state=42)
scores = []

for train_idx, test_idx in kf.split(protein_dataset):
    # Train and test
    score = pipeline.cross_validate_binding_sites(protein_dataset[train_idx])
    scores.append(score)

print(f"Cross-validation score: {np.mean(scores):.3f} ± {np.std(scores):.3f}")
```

### Experimental Validation
```python
# Compare with experimental data
experimental_data = [
    {'pdb_file': 'protein1.pdb', 'metal_ions': ['Zn2+'], 'binding_efficiency': 0.7},
    {'pdb_file': 'protein2.pdb', 'metal_ions': ['Cu2+'], 'binding_efficiency': 0.8}
]

validation_results = pipeline.validate_with_experimental_data(experimental_data)
for result in validation_results.values():
    print(f"Error: {result['error']:.3f}")
```

## 📚 Documentation

### Comprehensive Documentation
- **[Enhanced Mathematical Framework](docs/enhanced_mathematical_framework.md)**: Detailed mathematical formulation
- **[Architecture Summary](docs/enhanced_architecture_summary.md)**: System architecture overview
- **[Demo1](docs/Demo1.md)**: Original environmental parameter equations

### API Documentation
```python
# Enhanced Binding Kinetics
from src.enhanced_binding_kinetics import EnhancedBindingKinetics
kinetics = EnhancedBindingKinetics(config)

# Enhanced Binding Site Identification
from src.enhanced_binding_site_identification import EnhancedBindingSiteIdentifier
identifier = EnhancedBindingSiteIdentifier(config)

# Main Pipeline
from src.enhanced_main import EnhancedMetalloproteinPipeline
pipeline = EnhancedMetalloproteinPipeline()
```

### Complete File Strcture Created: 
```text 
Metalloproteins/
├── config/
│   └── enhanced_config.yaml          # Comprehensive configuration
├── docs/
│   ├── enhanced_mathematical_framework.md    # Detailed equations
│   ├── enhanced_architecture_summary.md      # Architecture overview
│   └── Demo1.md                              # Original equations
├── src/
│   ├── enhanced_binding_kinetics.py          # Enhanced ODE/PDE solver
│   ├── enhanced_binding_site_identification.py # Multi-algorithm integration
│   ├── enhanced_main.py                       # Main pipeline
│   ├── pdb_processor.py                       # PDB processing
│   └── brownian_simulation.py                 # Brownian motion
├── enhanced_example_usage.py                  # Comprehensive demonstration
├── requirements.txt                           # Updated dependencies
└── README.md                                  # Enhanced documentation
```
## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup
```bash
# Clone repository
git clone https://github.com/yourusername/metalloproteins.git
cd metalloproteins

# Install development dependencies
pip install -r requirements.txt
pip install -e .

# Run tests
pytest tests/

# Run linting
black src/
flake8 src/
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **MESPEUS Database**: For experimental metalloprotein data
- **MetalNet, Metal3D, bindEmbed21, AlphaFill**: For binding site prediction algorithms
- **PyMOL**: For 3D structure visualization
- **RF Diffusion**: For diffusion process visualization

## 📞 Support

For questions and support:
- Create an issue on GitHub
- Check the documentation in the `docs/` folder
- Run the example scripts for usage demonstrations

## Differences between CHARMM-GUI and this repository: 
- CHARMM-GUI:
    - GUI-driven, multi-platform input generation 
    - Advance system preparation (mutation, glycosylation, ligand parameterization)
    - Visualization using JSmol
    - High-throughput, automation, and enhanced sampling support

---

**The Enhanced Metalloprotein Binding Efficiency Prediction Pipeline** - Making metalloprotein analysis more accurate and comprehensive than ever before! 🧬⚗️🔬 