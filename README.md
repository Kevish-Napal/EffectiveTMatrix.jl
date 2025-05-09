# EffectiveTMatrix.jl


This package implements the Effective Waves Method for calculating the effective T-matrix of a random particulate system, as detailed in [[1]]((https://royalsocietypublishing.org/doi/full/10.1098/rspa.2023.0660)). It provides tools for simulating wave scattering and validating results using Monte Carlo simulations.


---
**Reference:** 

[[1]](https://royalsocietypublishing.org/doi/full/10.1098/rspa.2023.0660)  K. K. Napal, P. S. Piva, and A. L. Gower. Effective T-matrix of a cylinder filled with a random two-dimensional particulate. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, 480(2292):20230660, 2024.


[![][doi-img]][doi-url] 

## Installation

[Install Julia v1.6.1 or later](https://julialang.org/downloads/) then run in the Julia REPL:

```julia
using Pkg
Pkg.add https://github.com/Kevish-Napal/EffectiveTMatrix.jl.git
```

---
## Purpose of this package

<ol>
  <li> Compute the effective T-matrix by using the Effective Waves Method as described in <a href="http://arxiv.org/abs/2308.13338">[1]</a> 
  <li> Validate the result with Monte Carlo simulations</li>
</ol>

**__Long term goal:__** The package currently focuses on the acoustic scattering of particles in a homogeneous host medium, confined in a cylindrical area. The code is setup to easily integrate higher dimensions and other physics such as elastodynamics or electromagnetics in the future. The package will also take into account more complex microstructures. For example, the case of particles contained in a cylinder having different properties from the host medium, known as the two medium problem.


## Learn through examples

**Example 1:** compute the coefficients of the effective T-matrix and the resulting averaged scattered field from a given incident field.

*--> examples/pressure_field* 

**Example 2:** covers the main functionality of the package. That is validation of the Effective Waves Method (fast computation of the effective T-matrix) against the Monte Carlo method (slow computation of the effective T-matrix based on averaging of simulated scattering data).

*--> examples/monte_carlo_validation*

**Example 3:** attempts to clarify the theory by exhibiting the "mode to mode scattering" of the averaged cylinder.

*--> examples/modal_scattering*

