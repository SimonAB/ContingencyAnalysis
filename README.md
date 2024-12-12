# ContingencyAnalysis.jl

A Julia package for analysing contingency tables using exact and Monte Carlo versions of Fisher's exact test and Chi-square tests. The package automatically selects appropriate tests based on table characteristics and provides Monte Carlo alternatives when exact computations are infeasible.

## Installation

```julia
using Pkg
Pkg.add("ContingencyAnalysis")
```

## Usage

```julia
using ContingencyAnalysis

# Example: Treatment effectiveness
treatment_table = [
    10  2;  # Treatment group: success, failure
    1   8   # Control group: success, failure
]

# Automatic analysis - selects appropriate test
analyse_contingency(treatment_table)

# Manual test selection
p_value = monte_carlo_fisher(treatment_table, n_simulations=10000)
chi_stat = calculate_chi_statistic(treatment_table)
```

## Features

- Automatic test selection based on:
  - Table dimensions
  - Sample size
  - Expected cell frequencies
  - Computational feasibility
- Monte Carlo simulation for larger tables
- One-sided and two-sided tests for 2×2 tables
- Handling of tables with zero marginals
- Chi-square statistics and expected frequencies

## Test Selection Logic

The package chooses between:

- Fisher's exact test (for small samples or sparse tables)
- Chi-square test (for larger samples with sufficient expected frequencies)
- Monte Carlo versions when exact computation is impractical

## Examples

See the `examples` directory for:

- `iris_example.jl`: Species classification
- `Sim_examples.jl`: Simulation demonstrations

## Dependencies

This package builds upon several core Julia packages:

- [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) - Provides the underlying statistical tests
- [Statistics.jl](https://github.com/JuliaLang/Statistics.jl) - Core statistical functionality
- [Random.jl](https://github.com/JuliaLang/Random.jl) - Random number generation for Monte Carlo simulations

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## License

This package is licensed under the MIT License.
