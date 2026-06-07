# Numerical Experiments for Higher-Order Stochastic Dominance Thesis

This repository contains the datasets, Julia code, and output files used in the numerical experiments of the thesis:

**Higher-Order Stochastic Dominance and Application in Portfolio Optimization**

The thesis is a rewrite and numerical replication of the paper on higher-order stochastic dominance constraints in optimization. The experiments here use Fama--French industry portfolio returns instead of the original dataset.

## Data

The experiments use five Fama--French industry portfolios:

* NoDur
* Manuf
* BusEq
* Shops
* Hlth

The benchmark is constructed as the equally weighted return of these five industries.

The main data files are:

```text
ff12_selected5_all_months.csv
ff12_selected5_last22.csv
ff12_selected5_last36.csv
```

## Julia scripts

The main experiment scripts are:

```text
step8_ff12_run_22m_order_inf.jl
step7_ff12_run_36m_order_inf.jl
step9_ff12_risk_36m_beta_scan.jl
step13_compare_tsd_36m_v2_with_counts.jl
```

The plotting and table-generation scripts are:

```text
plot_figure6_style.jl
make_final_tables.jl
```

## Experiments

### Experiment 1: 22-month higher-order stochastic dominance experiment

Script:

```text
step8_ff12_run_22m_order_inf.jl
```

This experiment solves the expected-return maximization problem under stochastic dominance constraints for the orders

```text
2, 3, 5, 10, 15, 20, infinity.
```

### Experiment 2: 36-month robustness experiment

Script:

```text
step7_ff12_run_36m_order_inf.jl
```

This experiment repeats the first experiment using a 36-month window.

### Experiment 3: Risk-function experiment

Script:

```text
step9_ff12_risk_36m_beta_scan.jl
```

This experiment fixes the stochastic order at the second order and compares different values of the risk level:

```text
beta = 0.1, 0.5, 0.8.
```

### Experiment 4: Comparison with a Post--Kopa type formulation

Script:

```text
step13_compare_tsd_36m_v2_with_counts.jl
```

This experiment compares the third-order stochastic dominance solution with a Post--Kopa type baseline using 25 equally spaced threshold levels.

## Reproducing the results

To reproduce the experiments, activate the Julia project environment:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Then run the scripts in the following order:

```bash
julia --project=. step8_ff12_run_22m_order_inf.jl
julia --project=. step7_ff12_run_36m_order_inf.jl
julia --project=. step9_ff12_risk_36m_beta_scan.jl
julia --project=. step13_compare_tsd_36m_v2_with_counts.jl
julia --project=. plot_figure6_style.jl
julia --project=. make_final_tables.jl
```

## Output files

The main output files are:

```text
step8_summary_with_inf.csv
step8_weights_with_inf.csv
step7_summary_with_inf.csv
step7_weights_with_inf.csv
step9_summary.csv
step9_table1_like.csv
compare_tsd_36m_v2_with_counts_summary.csv
compare_tsd_36m_v2_with_counts_weights.csv
Table1_final.csv
Table3_final.csv
fig6_style_22m_final.png
fig6_style_36m_final.png
```

These files correspond to the numerical results, tables, and figures reported in the thesis.
# Higher-Order-Stochastic-Dominance-and-Application-in-Portfolio-Optimization-thesis-experiment
