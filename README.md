# gaSampling

Genetic algorithm allocation for stratified two-stage cluster sampling in R.

gaSampling searches for the number of clusters to sample in each stratum
and the number of units to sample within each selected cluster. The aim is
to minimise the variance of the overall mean while satisfying sample-size,
allocation and budget constraints.

The project is under development and currently runs by sourcing R scripts.

## Overview

In a stratified two-stage design, the population is divided into strata.
Clusters are selected within each stratum, followed by units within each
selected cluster.

The joint allocation procedure searches over:

- `c_h`: the number of sampled clusters in stratum h.
- `k_h`: the number of sampled units per cluster in stratum h.

The resulting stratum sample size is `n_h = c_h * k_h`.

The project provides functions to prepare stratum summaries, construct
starting allocations, repair candidate allocations, run the genetic
algorithm and check the resulting constraints.

## Setup

Install the GA package if it is not already available:

```r
install.packages("GA")
```

From the project directory, load the functions needed for the example:

```r
source("variance_twostage.R")
source("cost_twostage.R")
source("repair_joint.R")
source("check_constraints.R")
source("fitness_ch_kh.R")
source("allocate_twostage_joint.R")
```

## Quick example

This example uses illustrative stratum summaries and costs. It does not
require an external dataset.

```r
strata <- data.frame(
  stratum = c("A", "B"),
  C_h = c(20, 25),
  W_h = c(0.4, 0.6),
  sigma2_bh = c(2, 3),
  sigma2_wh = c(6, 9),
  cost_cluster = c(100, 120),
  cost_unit = c(5, 6)
)

result <- allocate_twostage_ga_joint(
  strata = strata,
  total_n = 150,
  max_cost = 3900,
  popSize = 50,
  maxiter = 100,
  run = 20,
  seed = 123,
  suggestion = c(10, 15, 5, 8)
)
```

The starting suggestion contains cluster counts for all strata followed
by unit counts per cluster: `c(c_h, k_h)`.

Check feasibility before interpreting the allocation:

```r
result$feasible
result$constraint_check
result$allocation
```

Inspect the objective and achieved totals:

```r
result$variance
result$total_n
result$total_cost
result$fallback_used
```

The settings above provide a small demonstration. Assess stability across
multiple random seeds before choosing settings for a substantive analysis.

## Input data

The joint allocation function expects one row per stratum with these columns:

| Column | Meaning |
|--------|---------|
| `stratum` | Stratum identifier |
| `C_h` | Number of available clusters |
| `W_h` | Stratum population weight |
| `sigma2_bh` | Between-cluster variance component |
| `sigma2_wh` | Within-cluster variance component |
| `cost_cluster` | Cost per sampled cluster |
| `cost_unit` | Cost per sampled secondary unit |

Population weights should sum to one. Variance components and costs must
refer to the same outcome, design and cost assumptions used in the analysis.
Use a consistent currency for all costs.

For workflows beginning with unit-level data, `read_survey_data()` and
`prepare_sampling_frame()` provide data preparation helpers.

## Objective and constraints

The variance objective is:

```r
sum(W_h^2 * (
  sigma2_bh / c_h +
  sigma2_wh / (c_h * k_h)
))
```

Total cost is:

```r
sum(cost_cluster * c_h + cost_unit * c_h * k_h)
```

The current joint allocator applies these constraints:

- Cluster and unit counts are integers.
- Each stratum has at least two sampled clusters.
- Cluster counts do not exceed `C_h`.
- Units per cluster range from one to `total_n`.
- The achieved total sample size meets or exceeds `total_n`.
- When a budget is supplied, spending falls between 95% and 102%
  of `max_cost`, using the current helper defaults.

`total_n` is a minimum requirement, so the returned allocation may exceed it.
Similarly, `max_cost` is a reference budget with an upper tolerance, rather
than a strict spending ceiling. Supply `max_cost = NULL` to disable budget
constraints.

## Understanding the results

`allocate_twostage_ga_joint()` returns a list containing:

| Component | Meaning |
|-----------|---------|
| `feasible` | Whether the returned allocation passes the constraints |
| `allocation` | Stratum allocations and variance contributions |
| `variance` | Total variance for the returned allocation |
| `between_variance` | Between-cluster contribution to total variance |
| `within_variance` | Within-cluster contribution to total variance |
| `total_n` | Achieved total sample size |
| `total_cost` | Achieved total cost |
| `constraint_check` | Separate indicators for each constraint |
| `repair` | Allocation repair results and diagnostics |
| `fallback_used` | Whether a repaired starting suggestion replaced an infeasible GA result |
| `ga_object` | The fitted genetic algorithm object |

The final allocation can differ from the solution stored in `ga_object`
because of repair or fallback. Use `result$allocation` when reporting
the returned design.

## Main source files

| File | Purpose |
|------|---------|
| `read_survey_data.R` | Standardise unit-level survey data |
| `sampling_frame.R` | Construct stratum-level sampling summaries |
| `feasible_joint_start.R` | Construct and assess starting allocations |
| `allocate_twostage_joint.R` | Run the joint genetic algorithm |
| `fitness_ch_kh.R` | Evaluate candidate allocation fitness |
| `repair_joint.R` | Attempt to repair cluster and unit allocations |
| `check_constraints.R` | Check allocation feasibility |
| `variance_twostage.R` | Calculate the variance objective |
| `cost_twostage.R` | Calculate total sampling cost |
| `frame_c.R` | Calculate an available and affordable total cluster count |
| `sample_size_feasibility.R` | Provide sample-size feasibility and greatest-common-divisor helpers |

## Assumptions and limitations

The genetic algorithm and repair procedure are heuristic. They do not
guarantee a global optimum or a feasible result. Always inspect `feasible`
and the constraint diagnostics.

The variance calculation does not include finite population corrections.

Within each stratum, the allocation uses a common number of sampled units
per selected cluster. The current joint allocator does not enforce observed
unit capacity for each individual cluster.

Results depend on variance estimates, cost assumptions, starting values
and genetic algorithm settings. Record these inputs and the random seed
when reporting an analysis.

The procedure returns allocation counts, it does not select the actual
clusters or units.

## Vignettes
 For a complete walkthrough of the workflow, see the [gaSampling vignette (PDF)](https://github.com/roxycampher/gaSampling/blob/main/vignettes/gaSampling.pdf)
