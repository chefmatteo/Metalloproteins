# Enhanced Metalloprotein Pipeline: Installation and Setup Guide

## Abstract

This guide provides comprehensive instructions for installing and setting up the enhanced metalloprotein binding efficiency prediction pipeline using UV package manager. It includes dependency installation, environment configuration, parameter modification, and troubleshooting steps.

## 1. Prerequisites

### 1.1 System Requirements

- **Operating System**: Linux, macOS, or Windows (WSL recommended for Windows)
- **Python**: 3.8 or higher
- **Memory**: Minimum 4GB RAM, 8GB recommended for large simulations
- **Storage**: 2GB free space for installation and data
- **Internet**: Required for package downloads and database access

### 1.2 Required Software

- **UV Package Manager**: Modern Python package manager
- **Git**: For cloning the repository
- **C Compiler**: For some scientific packages (gcc on Linux/macOS, Visual Studio on Windows)

## 2. UV Installation

### 2.1 Install UV Package Manager

**Linux/macOS:**
```bash
# Install UV using curl
curl -LsSf https://astral.sh/uv/install.sh | sh

# Add UV to PATH (add to ~/.bashrc or ~/.zshrc)
export PATH="$HOME/.cargo/bin:$PATH"

# Verify installation
uv --version
```

**Windows:**
```powershell
# Install UV using PowerShell
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

# Verify installation
uv --version
```

**Alternative Installation Methods:**
```bash
# Using pip (if you have Python already)
pip install uv

# Using Homebrew (macOS)
brew install uv

# Using Cargo (if you have Rust)
cargo install uv
```

### 2.2 Verify UV Installation

```bash
# Check UV version
uv --version

# Check Python version
uv python --version

# List available Python versions
uv python list
```

## 3. Project Setup

### 3.1 Clone Repository

```bash
# Clone the repository
git clone https://github.com/your-username/metalloproteins.git
cd metalloproteins

# Verify project structure
ls -la
```

### 3.2 Create Virtual Environment

```bash
# Create virtual environment with UV
uv venv

# Activate virtual environment
# On Linux/macOS:
source .venv/bin/activate

# On Windows:
.venv\Scripts\activate

# Verify activation
which python  # Should point to .venv/bin/python
```

### 3.3 Install Dependencies

```bash
# Install all dependencies from requirements.txt
uv pip install -r requirements.txt

# Or install dependencies individually for better control
uv pip install numpy scipy matplotlib seaborn
uv pip install biopython scikit-learn
uv pip install requests pandas
uv pip install pyyaml
uv pip install multiprocessing-logging

# Verify installations
python -c "import numpy, scipy, matplotlib, seaborn, Bio, sklearn; print('All packages installed successfully!')"
```

### 3.4 Install Development Dependencies (Optional)

```bash
# Install development tools
uv pip install pytest pytest-cov
uv pip install black flake8 mypy
uv pip install jupyter notebook
uv pip install ipython

# Install documentation tools
uv pip install sphinx sphinx-rtd-theme
uv pip install myst-parser
```

## 4. Configuration Setup

### 4.1 Initial Configuration

```bash
# Create configuration directory if it doesn't exist
mkdir -p config/enhanced

# Copy default configuration
cp config/enhanced/enhanced_config.yaml config/enhanced/my_config.yaml
```

### 4.2 Environment Variables

Create a `.env` file in the project root:

```bash
# Create environment file
cat > .env << EOF
# Project settings
PROJECT_ROOT=$(pwd)
DATA_DIR=\${PROJECT_ROOT}/data
RESULTS_DIR=\${PROJECT_ROOT}/results
LOGS_DIR=\${PROJECT_ROOT}/logs

# Database settings
MESPEUS_URL=http://mespeus.nchu.edu.tw/
PDB_URL=https://files.rcsb.org/download/

# Computational settings
N_PROCESSES=4
MEMORY_LIMIT=8GB
CACHE_DIR=\${PROJECT_ROOT}/cache

# Logging
LOG_LEVEL=INFO
LOG_FILE=\${LOGS_DIR}/metalloprotein.log
EOF

# Load environment variables
source .env
```

### 4.3 Directory Structure Setup

```bash
# Create necessary directories
mkdir -p data/{pdb,sequences,experimental}
mkdir -p results/{plots,data,reports}
mkdir -p logs
mkdir -p cache
mkdir -p temp

# Set permissions
chmod 755 data results logs cache temp
```

## 5. Parameter Configuration

### 5.1 Basic Configuration File

Create a custom configuration file `config/enhanced/my_config.yaml`:

```yaml
# Enhanced Configuration for Metalloprotein Binding Efficiency Prediction

# =============================================================================
# ENVIRONMENTAL PARAMETERS
# =============================================================================

environmental_conditions:
  # Temperature settings
  temperature:
    initial: 298.15  # Kelvin (25°C)
    range: [273.15, 373.15]  # 0°C to 100°C
    gradient: false  # Enable temperature gradients
  
  # pH settings
  pH:
    initial: 7.0
    range: [4.0, 10.0]
    buffer_system: "phosphate"
    buffer_concentration: 0.1  # M
  
  # Pressure settings
  pressure:
    initial: 1.0  # atm
    range: [0.1, 100.0]  # atm
  
  # Redox potential settings
  redox_potential:
    initial: 0.0  # V vs SHE
    range: [-0.5, 0.5]  # V

# =============================================================================
# SPATIAL DISCRETIZATION
# =============================================================================

spatial_discretization:
  # Reaction chamber dimensions
  chamber:
    dimensions: [10e-6, 10e-6, 10e-6]  # 10×10×10 μm³
    grid_size: [10, 10, 10]  # 1000 cubes total
    boundary_conditions: "periodic"  # Options: periodic, reflective, absorbing
  
  # Cube properties
  cube:
    volume: 1e-18  # m³ (1 μm³)
    diffusion_time_step: 1e-12  # s (1 ps)
    reaction_time_step: 1e-9  # s (1 ns)

# =============================================================================
# METAL ION PROPERTIES
# =============================================================================

metal_ions:
  Zn2+:
    radius: 0.74e-10  # m
    mass: 65.38e-27  # kg
    charge: 2
    diffusion_coeff_0: 7.0e-10  # m²/s at 298K
    activation_energy_diffusion: 15e3  # J/mol
    binding_enthalpy: -50e3  # J/mol
    pKa_effect: 0.5
    redox_sensitivity: 0.1
    
  Cu2+:
    radius: 0.73e-10
    mass: 63.55e-27
    charge: 2
    diffusion_coeff_0: 7.1e-10
    activation_energy_diffusion: 16e3
    binding_enthalpy: -60e3
    pKa_effect: 0.8
    redox_sensitivity: 0.2
    
  Fe2+:
    radius: 0.78e-10
    mass: 55.85e-27
    charge: 2
    diffusion_coeff_0: 7.2e-10
    activation_energy_diffusion: 14e3
    binding_enthalpy: -45e3
    pKa_effect: 0.3
    redox_sensitivity: 0.15

# =============================================================================
# BINDING SITE IDENTIFICATION ALGORITHMS
# =============================================================================

binding_site_algorithms:
  # Algorithm weights for consensus scoring
  weights:
    MetalNet: 0.30
    Metal3D: 0.25
    bindEmbed21: 0.25
    AlphaFill: 0.20
  
  # MetalNet settings
  MetalNet:
    distance_threshold: 3.0  # Å
    confidence_threshold: 0.7
    use_network_clustering: true
    
  # Metal3D settings
  Metal3D:
    resolution: 1.0  # Å
    confidence_threshold: 0.6
    use_geometric_features: true
    
  # bindEmbed21 settings
  bindEmbed21:
    embedding_dimension: 21
    sequence_window: 15
    confidence_threshold: 0.5
    
  # AlphaFill settings
  AlphaFill:
    fill_ligands: true
    confidence_threshold: 0.8
    use_alphafold_predictions: true

# =============================================================================
# NUMERICAL SETTINGS
# =============================================================================

numerical_settings:
  # ODE solver configuration
  ode_solver:
    method: "RK45"  # Options: RK45, RK23, BDF, LSODA
    rtol: 1e-6      # Relative tolerance
    atol: 1e-8      # Absolute tolerance
    max_step: 1e-3  # Maximum time step
  
  # Parallel processing
  parallel_processing: true
  n_processes: 4
  
  # Memory management
  use_sparse_matrices: true
  compress_results: false
  compression_ratio: 0.1
  
  # Performance optimization
  cache_intermediate_results: true
  max_memory_usage: "8GB"

# =============================================================================
# OUTPUT SETTINGS
# =============================================================================

output_settings:
  # File formats
  save_plots: true
  plot_format: "png"  # Options: png, pdf, svg
  save_data: true
  data_format: "json"  # Options: json, csv, hdf5
  
  # Plot settings
  plot_style: "seaborn-v0_8"
  figure_dpi: 300
  figure_size: [10, 6]
  
  # Output directories
  plots_dir: "results/plots"
  data_dir: "results/data"
  reports_dir: "results/reports"
  
  # Logging
  log_level: "INFO"
  log_file: "logs/metalloprotein.log"
```

### 5.2 Parameter Modification Guide

#### 5.2.1 Environmental Parameters

**Temperature Settings:**
```yaml
environmental_conditions:
  temperature:
    initial: 310.15  # Change to 37°C for physiological conditions
    range: [298.15, 323.15]  # Narrow range for specific studies
    gradient: true  # Enable temperature gradients
```

**pH Settings:**
```yaml
environmental_conditions:
  pH:
    initial: 6.5  # Change for acidic conditions
    range: [5.0, 8.0]  # Adjust range for specific pH studies
    buffer_system: "tris"  # Change buffer system
    buffer_concentration: 0.05  # Adjust buffer concentration
```

**Pressure Settings:**
```yaml
environmental_conditions:
  pressure:
    initial: 2.0  # Change for high-pressure studies
    range: [1.0, 10.0]  # Adjust range for pressure studies
```

**Redox Settings:**
```yaml
environmental_conditions:
  redox_potential:
    initial: 0.1  # Change for oxidizing conditions
    range: [-0.2, 0.3]  # Adjust range for redox studies
    buffer_system: "glutathione"  # Change redox buffer
```

#### 5.2.2 Spatial Discretization

**Grid Resolution:**
```yaml
spatial_discretization:
  chamber:
    grid_size: [20, 20, 20]  # Increase for higher resolution (8000 cubes)
    # or
    grid_size: [5, 5, 5]     # Decrease for faster computation (125 cubes)
```

**Chamber Size:**
```yaml
spatial_discretization:
  chamber:
    dimensions: [20e-6, 20e-6, 20e-6]  # Larger reaction chamber
    # or
    dimensions: [5e-6, 5e-6, 5e-6]     # Smaller reaction chamber
```

**Boundary Conditions:**
```yaml
spatial_discretization:
  chamber:
    boundary_conditions: "reflective"  # Change from "periodic"
    # Options: "periodic", "reflective", "absorbing"
```

#### 5.2.3 Metal Ion Properties

**Add New Metal Ions:**
```yaml
metal_ions:
  # Existing ions...
  
  Mg2+:
    radius: 0.72e-10
    mass: 24.31e-27
    charge: 2
    diffusion_coeff_0: 7.0e-10
    activation_energy_diffusion: 13e3
    binding_enthalpy: -30e3
    pKa_effect: 0.1
    redox_sensitivity: 0.05
    
  Ca2+:
    radius: 1.00e-10
    mass: 40.08e-27
    charge: 2
    diffusion_coeff_0: 7.9e-10
    activation_energy_diffusion: 12e3
    binding_enthalpy: -25e3
    pKa_effect: 0.05
    redox_sensitivity: 0.02
```

**Modify Existing Properties:**
```yaml
metal_ions:
  Zn2+:
    # ... existing properties ...
    binding_enthalpy: -60e3  # Change binding strength
    pKa_effect: 0.8          # Increase pH sensitivity
    redox_sensitivity: 0.2   # Increase redox sensitivity
```

#### 5.2.4 Algorithm Weights

**Adjust Algorithm Importance:**
```yaml
binding_site_algorithms:
  weights:
    MetalNet: 0.40      # Increase MetalNet weight
    Metal3D: 0.30       # Increase Metal3D weight
    bindEmbed21: 0.20   # Decrease bindEmbed21 weight
    AlphaFill: 0.10     # Decrease AlphaFill weight
```

**Algorithm-Specific Parameters:**
```yaml
binding_site_algorithms:
  MetalNet:
    distance_threshold: 3.5    # Increase distance threshold
    confidence_threshold: 0.6  # Lower confidence threshold
    
  Metal3D:
    resolution: 0.5            # Increase geometric precision
    confidence_threshold: 0.5  # Lower confidence threshold
    
  bindEmbed21:
    sequence_window: 20        # Increase sequence window
    confidence_threshold: 0.4  # Lower confidence threshold
    
  AlphaFill:
    confidence_threshold: 0.7  # Lower confidence threshold
    template_similarity_threshold: 0.6  # Adjust similarity threshold
```

#### 5.2.5 Numerical Settings

**Solver Configuration:**
```yaml
numerical_settings:
  ode_solver:
    method: "BDF"        # Change to stiff solver
    rtol: 1e-8          # Increase precision
    atol: 1e-10         # Increase precision
    max_step: 1e-4      # Decrease maximum step size
```

**Performance Settings:**
```yaml
numerical_settings:
  parallel_processing: true
  n_processes: 8        # Increase for more cores
  
  # Memory management
  use_sparse_matrices: true
  compress_results: true
  compression_ratio: 0.05  # Higher compression
  
  max_memory_usage: "16GB"  # Increase memory limit
```

## 6. Running the Pipeline

### 6.1 Basic Usage

```bash
# Activate virtual environment
source .venv/bin/activate

# Run basic analysis
python enhanced_example_usage.py
```

### 6.2 Custom Configuration

```bash
# Run with custom configuration
python -c "
from src.enhanced.enhanced_main import EnhancedMetalloproteinPipeline
import yaml

# Load custom configuration
with open('config/enhanced/my_config.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Initialize pipeline with custom config
pipeline = EnhancedMetalloproteinPipeline(config=config)

# Run analysis
results = pipeline.run_enhanced_analysis(
    pdb_file='data/pdb/protein.pdb',
    protein_sequence='CHEDCH',
    metal_ions=['Zn2+', 'Cu2+'],
    initial_concentrations={'Zn2+': 1e-6, 'Cu2+': 1e-6},
    time_span=(0, 1000),
    save_results=True
)

print('Analysis completed successfully!')
"
```

### 6.3 Batch Processing

Create a batch script `run_batch_analysis.py`:

```python
#!/usr/bin/env python3
"""
Batch analysis script for multiple proteins
"""

import os
import yaml
from src.enhanced.enhanced_main import EnhancedMetalloproteinPipeline

def run_batch_analysis():
    # Load configuration
    with open('config/enhanced/my_config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    
    # Initialize pipeline
    pipeline = EnhancedMetalloproteinPipeline(config=config)
    
    # Define protein files
    protein_files = [
        'data/pdb/protein1.pdb',
        'data/pdb/protein2.pdb',
        'data/pdb/protein3.pdb'
    ]
    
    # Define metal ions
    metal_ions = ['Zn2+', 'Cu2+', 'Fe2+']
    initial_concentrations = {metal: 1e-6 for metal in metal_ions}
    
    # Run analysis for each protein
    for pdb_file in protein_files:
        if os.path.exists(pdb_file):
            print(f"Analyzing {pdb_file}...")
            
            try:
                results = pipeline.run_enhanced_analysis(
                    pdb_file=pdb_file,
                    metal_ions=metal_ions,
                    initial_concentrations=initial_concentrations,
                    time_span=(0, 1000),
                    save_results=True
                )
                
                print(f"Completed analysis of {pdb_file}")
                print(f"Binding efficiency: {results['binding_efficiency']['overall']:.3f}")
                
            except Exception as e:
                print(f"Error analyzing {pdb_file}: {e}")
        else:
            print(f"File not found: {pdb_file}")

if __name__ == "__main__":
    run_batch_analysis()
```

Run the batch script:
```bash
python run_batch_analysis.py
```

## 7. Troubleshooting

### 7.1 Common Installation Issues

**UV Installation Problems:**
```bash
# If UV installation fails, try alternative methods
pip install uv

# Or use system package manager
# Ubuntu/Debian:
sudo apt install uv

# macOS:
brew install uv
```

**Python Version Issues:**
```bash
# Check Python version
python --version

# If Python version is too old, install newer version with UV
uv python install 3.11
uv python use 3.11
```

**Package Installation Issues:**
```bash
# Clear cache and reinstall
uv cache clean
uv pip install --force-reinstall -r requirements.txt

# Install packages individually to identify problematic packages
uv pip install numpy
uv pip install scipy
# ... continue for each package
```

### 7.2 Runtime Issues

**Memory Issues:**
```yaml
# Reduce memory usage in config
numerical_settings:
  n_processes: 2        # Reduce number of processes
  max_memory_usage: "4GB"  # Reduce memory limit
  compress_results: true
  compression_ratio: 0.1

spatial_discretization:
  chamber:
    grid_size: [5, 5, 5]  # Reduce grid size
```

**Convergence Issues:**
```yaml
# Improve numerical stability
numerical_settings:
  ode_solver:
    method: "BDF"        # Use stiff solver
    rtol: 1e-4          # Relax tolerance
    atol: 1e-6          # Relax tolerance
    max_step: 1e-2      # Increase step size
```

**Performance Issues:**
```yaml
# Optimize performance
numerical_settings:
  parallel_processing: true
  n_processes: 8        # Increase processes
  cache_intermediate_results: true
  use_sparse_matrices: true
```

### 7.3 Debugging

**Enable Debug Logging:**
```yaml
output_settings:
  log_level: "DEBUG"
  log_file: "logs/debug.log"
```

**Run with Verbose Output:**
```bash
python -u enhanced_example_usage.py 2>&1 | tee debug_output.log
```

**Check System Resources:**
```bash
# Monitor memory usage
htop

# Monitor disk usage
df -h

# Check Python process
ps aux | grep python
```

## 8. Best Practices

### 8.1 Configuration Management

1. **Version Control**: Keep configuration files in version control
2. **Backup Configurations**: Save working configurations with descriptive names
3. **Documentation**: Document parameter changes and their rationale
4. **Testing**: Test configurations on small datasets before large runs

### 8.2 Performance Optimization

1. **Grid Resolution**: Start with coarse grids and refine as needed
2. **Time Spans**: Use shorter time spans for initial testing
3. **Parallel Processing**: Adjust number of processes based on system capabilities
4. **Memory Management**: Monitor memory usage and adjust settings accordingly

### 8.3 Data Management

1. **Organize Data**: Use consistent directory structure
2. **Backup Results**: Regularly backup important results
3. **Cleanup**: Remove temporary files and old results periodically
4. **Documentation**: Keep detailed records of analysis parameters and results

## 9. Conclusion

This installation and setup guide provides comprehensive instructions for:

- Installing the enhanced metalloprotein pipeline using UV
- Configuring the environment and dependencies
- Modifying calculation parameters for different research scenarios
- Running analyses with custom configurations
- Troubleshooting common issues
- Following best practices for optimal performance

By following this guide, users can successfully set up and configure the pipeline for their specific research needs, ensuring accurate and efficient metalloprotein binding efficiency predictions.

## References

[1] UV Package Manager Documentation: https://docs.astral.sh/uv/
[2] Python Virtual Environments: https://docs.python.org/3/tutorial/venv.html
[3] YAML Configuration: https://yaml.org/spec/
[4] Scientific Python Stack: https://scipy.org/ 