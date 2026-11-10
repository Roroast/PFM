# Counting The Dimension of Fréchet Mean Sets induced by Polytope Norms

Code for the exact computation of the **Fréchet mean set** of a finite dataset when the distance is induced by the $\ell_1$, $\ell_{\infty}$ or tropical norms.
This code was used to generate the numerical computation results included in the paper ["On the Uniqueness of Fréchet Means for Polytope Norms"](https://arxiv.org/abs/2603.15785).

These computations were run on the Hamilton HPC Service of Durham University, using SLURM array-job scripts.

## Citation

To cite this work, please use the following bibtex:
```
@misc{talbut2026uniquenessfrechetmeanspolytope,
      title={On the Uniqueness of Fr\'echet Means for Polytope Norms}, 
      author={Roan Talbut and Andrew McCormack and Anthea Monod},
      year={2026},
      eprint={2603.15785},
      archivePrefix={arXiv},
      primaryClass={math.PR},
      url={https://arxiv.org/abs/2603.15785}, 
}
```

## How it works

We formulate the Fr\'echet mean set via the constrained quadratic optimisation problem:

$$
  \text{minimise } \Vert\vec d\Vert^2_2, \quad \text{given } \forall i \leq n \text{ and }\forall a \in A : \, a^{\top}( x_i -  \theta) \leq d_i,
$$

where $A$ is the set of normal vectors of the norm's unit ball and $\{ x_1, \dots, x_n\}$ is the dataset.

We then implement a variant of [Wolfe's combinatorial algorithm](https://link.springer.com/article/10.1007/BF01580381) to solve this problem exactly. For full technical details on our implementation, see ["On the Uniqueness of Fréchet Means for Polytope Norms"](https://arxiv.org/abs/2603.15785).

The pipeline in `exact_frechet_computation.jl` is:

1. **`lifted_Q`** — builds the feasible set from the data and the normal vectors, over `CDDLib` with exact rational arithmetic.
2. **`lifted_wolfe`** — a **Wolfe-style minimum-norm-point algorithm** on the lifted polytope, using exact linear solvers. It returns the optimal distances and a minimiser $x^*$.
3. **`Wolfe_FM_dimension`** — computes the Fr\'echet mean set from the optimal distances, then computes the dimension of the resulting Fréchet mean set.
4. **`solution_checker`** — asserts that $x^*$ lies in the mean set, that the distances match, and that $\mathbf 0$ lies in the convex hull of the active subgradients (the first-order optimality condition).

## Repository layout

```
exact_frechet_computation.jl   Core: alpha sets, lifted LP, Wolfe algorithm, dimension, checks
L1_Wolfe.jl                    Julia script for the ell_1 experiments
Linf_Wolfe.jl                  Julia script for the ell_inf experiments
Trop_Wolfe.jl                  Julia script for the tropical experiments
L1_Wolfe.sh, Linf_Wolfe.sh,    SLURM array-job scripts (1000 repetitions each)
Trop_Wolfe.sh
file_reading.jl                Parses job output, writes LaTeX tables, draws figures
Project.toml, Manifest.toml    Julia environment (Julia >= 1.9)
L_1_13_13_1000/                Raw SLURM .out/.err files for ell_1 norm
L_{inf}_13_13_1000/            Raw SLURM .out/.err files for ell_inf norm
Trop_13_13_1000/               Raw SLURM .out/.err files for tropical norm
LatexTables/                   Generated tables: average times, expected dimension, uniqueness probability
Figures/                       Generated PDFs of the two plots per problem
```

## Requirements

- **Julia 1.9 or newer** with the packages pinned in `Project.toml`: `JuMP`, `Polyhedra`, `CDDLib`, `Distributions`, `JLD`, `SlurmClusterManager`, plus `Plots`, `LaTeXTabulars` and `LaTeXStrings` for the post-processing step.
- A **SLURM cluster** for the bundled scripts. The drivers use `SlurmClusterManager` to add workers, so they expect to be launched as batch jobs; on a single machine you can replace `addprocs(SlurmManager())` with `addprocs(n)`.

## Running the experiments

Each Julia script takes six integers — `k_min k_max n_min n_max rep_min rep_max` — and runs every \((k, n, \text{rep})\) combination in parallel over the worker pool. For each run it draws a \(k \times n\) matrix of standard normal samples (converted to exact rationals) and prints one line:

```
k = 8, n = 6, rep = 137, dim = 0, time = 0.42 seconds
```

The bundled scripts sweep \(k = 2{:}13\), \(n = 2{:}13\), 1000 repetitions, with one array task per repetition:

```bash
sbatch L1_Wolfe.sh      # ell_1
sbatch Linf_Wolfe.sh    # ell_inf
sbatch Trop_Wolfe.sh    # tropical
```

To run a single configuration locally:

```bash
julia --project=. L1_Wolfe.jl 2 13 2 13 1 1
```

## Post-processing

`file_reading.jl` reads the `.out` files from a results folder, flags any non-empty `.err` files, and produces:

- `LatexTables/<problem>_average_times.tex` — mean runtime in minutes per \((n, k)\);
- `LatexTables/<problem>_ex_dim.tex` — average mean-set dimension divided by \(k-1\);
- `LatexTables/<problem>_uniqueness_prob.tex` — empirical probability that the Fréchet mean is unique;
- `Figures/<problem>_expected_dimension.pdf` and `Figures/<problem>_uniqueness_probability.pdf`.

The committed tables and figures are the output of the full 1000-repetition sweep.
