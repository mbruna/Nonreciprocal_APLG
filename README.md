# Active-Passive Lattice Gas

This repository contains Julia code for the active-passive lattice gas model described in [Dynamical patterns and nonreciprocal effective interactions in an active-passive mixture through exact hydrodynamic analysis](https://doi.org/10.48550/arXiv.2408.03932).

## 1. System Requirements

### Software Dependencies
- **Operating System**: Compatible with Windows, macOS, and Linux
- **Julia**: Version 1.11.2 or later
- **Required Packages**: see `Project.toml`

### Tested Versions
- Successfully tested on Julia 1.11.2
- Verified on macOS 14.6, and Ubuntu 22.04

### Hardware Requirements
- Standard desktop computer (4+ CPU cores recommended for faster simulations)
- Minimum 8GB RAM (16GB recommended for large system sizes)
- No specialized hardware required

## 2. Installation Guide

### Instructions

1. Clone this repository:
   ```bash
   git clone https://github.com/mbruna/Nonreciprocal_APLG
   cd Nonreciprocal_APLG
   ```

2. Activate the project in Julia:
   ```julia
   using DrWatson
   @quickactivate "APLG"
   ```

3. Install dependencies:
   ```julia
   using Pkg
   Pkg.instantiate()
   ```

### Installation Time
- Typical installation time: [MEASURE AND INSERT ACTUAL TIME] on a standard desktop computer
- Package precompilation may take an additional on first use

## 3. Demo

In the `notebooks/` directory, you can find the notebook `demo.ipynb` that demonstrates the key features of the code.

### Running the Demo
Using IJulia, Jupyter Notebook or JupyterLab, open the notebook `demo.ipynb` in the `notebooks/` directory.

### Expected Output
- Particle simulation: Density profiles showing phase separation
- PDE simulation: Density profiles showing phase separation
- Travelling solution: Converged solution with wave speed c
- Outer region travelling solutions: three different types depending on the number of interfaces

### Expected Run Time
- Phase diagram: ~6 minutes
- Particle simulation: 2-3 minutes
- PDE simulation: 2-3 seconds
- Travelling solution (full problem): 2-3 minutes
- Outer region travelling solutions: 1-3 seconds

## 4. Instructions for Use

### Parameter Configuration

The `new_param` function generates a dictionary of parameters compatible with all simulation and analysis functions:

```julia
param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; 
                  T=T, name=param_name, save_interval=save_interval, 
                  save_on=true, pert=pert)
```

Key parameters:
- `DT`: Translational diffusion coefficient
- `v0`: Self-propulsion velocity (Péclet number)
- `ϕa`, `ϕp`: Active and passive particle densities
- `Lx`, `Ly`: System dimensions

### Particle Simulations

Run a particle-based simulation with the Gillespie algorithm:

```julia
run_new_sim(param)
```

### Time-Dependent Solutions

Solve the time-dependent hydrodynamic PDE (Eq. 1 in the paper):

```julia
run_new_pde(param)

# Load and analyze results
loaded, f, t = load_last_pde(param)
plot_density_profile(f)
```

### Travelling Wave Solutions

For finite domains, use the Newton-Raphson method:

```julia
f, u, c = solve_full(Lx, Nx, ϕa, ϕp, v0, u0; tol=1e-8, maxiters=100)
```

For infinite domains, use one of the following:

```julia
# No interfaces
f, u, c = solve_outer0(Lx, Nx, ϕa, ϕp, v0, u0)

# One interface
f, u, c = solve_outer1(Lx, Nx, ϕa, ϕp, v0, ϕ, γ, u0)

# Two interfaces
f, u, c = solve_outer2(Lx, Nx, ϕa, ϕp, v0, ind, γ, u0)
```

### Phase Diagram Analysis

Analyze stability and find coexisting phases:

```julia
# Check stability
max_eigenvalue = is_stable_value(ϕa, ϕp, v0)

# Find coexisting phases
find_sol, lower_limits, upper_limits = colapse_sol_interval(Pe=v0, γ=γ)
```

### Reproduction Instructions

The notebooks in the `notebooks/` directory reproduce the figures from our paper, taking the data from the `data/` directory (it can be found in the FigShare repository [here](https://figshare.com/articles/dataset/Data_for_Dynamical_patterns_and_nonreciprocal_effective_interactions_in_an_active-passive_mixture_through_exact_hydrodynamic_analysis_/28844123)).


## Contributing

This project was developed by James Mason, Robert L. Jack and [Maria Bruna](https://www.maths.ox.ac.uk/people/maria.bruna).

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use this code in your research, please cite:
```
@article{mason2024dynamical,
  title={Dynamical patterns and nonreciprocal effective interactions in an active-passive mixture through exact hydrodynamic analysis},
  author={Mason, James and Jack, Robert L. and Bruna, Maria},
  journal={arXiv preprint arXiv:2408.03932},
  year={2024},
  doi={10.48550/arXiv.2408.03932}
}
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE.md) file for details.
