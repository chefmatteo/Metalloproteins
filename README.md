# Metalloprotein Binding Efficiency Prediction Model

## Project Overview
This project implements a predictive model for estimating the binding efficiency of metalloproteins under various environmental conditions including temperature, pressure, and ion concentrations.

## Key Components

### 1. Structure Analysis Pipeline
- **PDB Structure Processing**: Parse and analyze protein structures from PDB files
- **Metal Binding Site Identification**: Use multiple algorithms to predict binding sites
- **Geometric Analysis**: Calculate binding pocket characteristics

### 2. Binding Efficiency Prediction
- **Brownian Motion Simulation**: Model ion diffusion in solution
- **Binding Probability Calculation**: Estimate likelihood of successful binding
- **Environmental Parameter Integration**: Temperature, pressure, concentration effects

### 3. Mathematical Framework
- **ODE System**: Coupled differential equations for binding kinetics
- **Stochastic Processes**: Brownian motion and collision theory
- **Thermodynamic Integration**: Free energy calculations

## Dependencies
- Python 3.8+
- NumPy, SciPy, Matplotlib
- BioPython (for PDB processing)
- PyTorch (for deep learning models)
- MDAnalysis (for molecular dynamics)

## Usage
```bash
# Install dependencies
pip install -r requirements.txt

# Run the main pipeline
python src/main.py --pdb_file protein.pdb --temperature 298 --pressure 1.0
```

## Mathematical Framework
See `docs/mathematical_framework.md` for detailed equations and derivations. 