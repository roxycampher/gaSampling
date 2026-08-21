# gaSampling

Multistage stratified cluster sampling package for two-stage allocation using
genetic algorithms.

## Main Diagnostic Script

Use `final_diagnostic_checks_present.R` when you want to test the full joint
GA workflow on the ECDC data and produce stability summaries, plots and CSV
outputs.

Before running it, make sure these files are in the same working directory:

- `ecdc-2021-public-v2.csv`
- `read_survey_data.R`
- `sampling_frame.R`
- `feasible_joint_start.R`
- `repair_joint.R`
- `check_constraints.R`
- `variance_twostage.R`
- `cost_twostage.R`
- `fitness_ch_kh.R`
- `allocate_twostage_joint.R`

Run the script line by line in RStudio if you want to inspect each object as it
is created. The multiple-run section can take longer, so reduce
`joint_ga_seed_count` if you only want a quick check.

## Final Diagnostics Flow

The final diagnostics script follows this flow:

1. Load packages and source helper functions.
2. Read `ecdc-2021-public-v2.csv`.
3. Clean the survey data with `read_survey_data()`.
4. Build the stratum-level sampling frame with `prepare_sampling_frame()`.
5. Set the optimization target:
   - `total_n <- 5000`
   - `max_cost <- 1400000`
   - GA settings such as population size, iterations, and seeds.
6. Find a feasible starting allocation with `find_feasible_joint_start()`.
7. Run one joint GA from the repaired Neyman starting allocation.
8. Run repeated joint GAs across starting suggestions and random seeds.
9. Summarise feasibility, variance, total cost, and allocation stability.
10. Save plots and CSV summaries to `joint_ga_stability_plots/`.

## Function Flow

`read_survey_data.R`

Standardises the raw survey data. It creates the common columns used by the rest
of the workflow, especially `strata`, `cluster` and the analysis variable.

`sampling_frame.R`

Uses the cleaned data to create a stratum-level frame. The important output is:

```r
frame$strata
```

This contains one row per stratum with population size, cluster counts, stratum
weights, variance components and cost variables.

`feasible_joint_start.R`

Finds a feasible starting allocation for the joint GA. The key function is:

```r
find_feasible_joint_start()
```

It tries several candidate total cluster counts, allocates clusters across
strata using different starting rules, repairs each allocation and checks
whether the repaired allocation satisfies the sample-size and cost constraints.

The main returned values are:

```r
joint_start$total_c
joint_start$c_h
joint_start$k_h
joint_start$suggestion
joint_start$all_starts
joint_start$diagnostics
```

Use `joint_start$suggestion` as the starting vector for
`allocate_twostage_ga_joint()`.

`repair_joint.R`

Contains `repair_ckh()`, which adjusts proposed `c_h` and `k_h` values so that
the allocation is closer to the requested sample size and budget constraints.

`check_constraints.R`

Checks whether an allocation is feasible. It verifies:

- integer `c_h` and `k_h`,
- lower and upper bounds,
- total sample size,
- total cost,
- minimum budget use and upper budget tolerance.

`fitness_ch_kh.R`

Defines the constrained fitness function used by the joint GA. Infeasible
solutions are penalised; feasible solutions are evaluated by variance.

`allocate_twostage_joint.R`

Runs the joint GA over both allocation vectors:

```r
c_h = sampled clusters per stratum
k_h = sampled units per selected cluster
```

It then repairs the final GA solution and returns the final allocation,
variance, total sample size, total cost, feasibility status and repair details.

## Interpreting Key Objects

`joint_starting_suggestions`

Shows the repaired starting allocations for each starting rule. Use this to see
whether the starting allocations are feasible before the GA runs.

`starting_suggestion_diagnostics`

Summarises each starting rule at the total-allocation level, including repaired
sample size, repaired total cost, feasibility and repair iterations.

`ga_result`

The single joint GA run from the repaired Neyman proportional start.

`run_summary`

One row per repeated GA run. This is the main table for checking stability
across seeds and starting suggestions.

`suggestion_run_summary`

Summarises variance, feasibility rate, cost and repair behaviour by starting
suggestion.

`allocation_summary`

Shows how stable each stratum's `c_h`, `k_h` and `n_h` allocations are across
repeated runs.

## Important Notes

The feasible-start helper finds a feasible starting allocation, not the final
optimal allocation. The joint GA can still move both `c_h` and `k_h`.

If cost parameters change, update them consistently in `prepare_sampling_frame()`
and in any affordability calculation. With the current cost model:

```r
total_cost = sum(c_h * cost_cluster) + sum(c_h * k_h * cost_unit)
```

A rough affordable cluster cap is:

```r
affordable_c <- floor((max_cost - cost_unit * total_n) / cost_cluster)
```

For repeated runs, start with a small `joint_ga_seed_count` while developing.
Increase it only when you need stronger stability evidence.