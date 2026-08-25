After the multiple test run, everything will be saved in the folder:

"joint_ga_stability_plots"

It will have these CSV result files:

- starting_suggestion_diagnostics.csv 
  Initial/repaired totals for each of the four GA starting suggestions.

- run_summary.csv 
  One row per GA run: seed, starting suggestion, variance, total sample size, total cost, feasibility, repair info.

- allocation_summary.csv  
  Average c_h, k_h, n_h, shares and sampling fractions by stratum and starting suggestion.

- suggestion_run_summary.csv  
  Summary by GA starting suggestion: min/median/mean/max variance, SD, CV, feasibility rate, cost, sample-size gap.

- basic_allocation_summary.csv 
  Non-GA benchmark summaries for equal, proportional and Neyman allocations.

- basic_allocation_allocations.csv  
  The actual benchmark c_h, k_h and n_h values by stratum.

- basic_benchmark_variances.csv  
  Clean benchmark variance table for equal, proportional and Neyman.

- ga_vs_basic_variance_comparison.csv  
  GA results compared with the basic benchmarks, including improvement vs equal, proportional and Neyman.

- suggestion_metric_comparison.csv  
  Long-format metrics table for plotting/inspection.

- stability_summary.csv  
  Final numerical stability summary by GA starting suggestion.

It will also have these PNG plots:

- 01_variance_by_seed.png
- 02_variance_distribution.png
- 03_variance_boxplot.png
- 04_starting_suggestion_metric_comparison.png
- 04b_ga_vs_basic_variance_comparison.png
- 04c_basic_benchmark_variance_histogram.png
- 05_ch_stability_by_stratum.png
- 06_kh_stability_by_stratum.png
- 07_nh_stability_by_stratum.png
- 08_combined_allocation_stability.png
- 10_repair_diagnostics.png