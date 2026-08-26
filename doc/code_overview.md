# Code overview

This project models life-course trajectories with semi-Markov processes and uses permutation-based tree methods to identify covariates associated with trajectory differences.

## 1. Main data model

The core dataset used throughout the project is a long-format trajectory table with one row per visited state episode.

Required columns:

- `id`: individual identifier
- `state`: current state index
- `time`: duration spent in that state

The last row for each individual is treated as right-censored, meaning the duration is known only to be at least `time`.

This representation is defined and used in:

- [src/semi_markov/mle_estimation.R](../src/semi_markov/mle_estimation.R)
- [src/two_samples_test.R](../src/two_samples_test.R)
- [src/random_forest/tree_construction.R](../src/random_forest/tree_construction.R)

## 2. Core statistical engine

### Semi-Markov estimation

The file [src/semi_markov/mle_estimation.R](../src/semi_markov/mle_estimation.R) contains the core likelihood-based estimation routines.

Key functions:

- `last_row_per_id()`: identifies the final censored observation for each individual
- `log_likelihood_alpha()`: contribution of the initial-state distribution
- `log_likelihood_P()`: contribution of transition probabilities
- `log_likelihood_omega()`: contribution of sojourn-time densities
- `mle_alpha()`: MLE for the initial distribution `alpha`
- `mle_P()`: MLE for the transition matrix `P`
- `mle_omega_gamma()`, `mle_omega_weibull()`, `mle_omega_exponential()`: MLEs for state-dependent sojourn times
- `mle_fit()`: combines the above into a complete semi-Markov estimator

The function `mle_fit()` returns:

- `estimator`: a named list with `alpha`, `P`, `omega`
- `log_likelihood`: total log-likelihood

This is the foundational object used throughout the rest of the project.

## 3. Two-sample testing

The file [src/two_samples_test.R](../src/two_samples_test.R) implements tests for comparing two trajectory populations.

Functions:

- `log_likelihood_ratio()`: compares global vs group-specific fits
- `lambda_statistic()`: $-2 \times$ log-likelihood ratio
- `likelihood_ratio_test()`: asymptotic chi-square approximation
- `parametric_bootstrap()`: bootstrap approximation under the null
- `permutation_test()`: permutation-based comparison

These functions are used as split criteria in the tree construction step.

## 4. Decision trees

The file [src/random_forest/tree_construction.R](../src/random_forest/tree_construction.R) builds a regression tree whose split criterion is based on the two-sample test statistic.

Main ideas:

- candidate splits are tested over each covariate
- numeric covariates are split by threshold
- categorical covariates are split by category subsets
- a split is retained only if it is statistically significant at level `alpha`
- recursive splitting continues until a stopping condition is reached

Important functions:

- `generate_bipartitions()`: enumerates candidate categorical partitions
- `best_split_categorical()`: scores categorical splits
- `best_split_numeric()`: scores numeric splits
- `find_best_split()`: picks the best overall split at a node
- `build_tree()`: recursively builds the tree
- `get_leaf()`: maps an observation to its terminal leaf

The result is a tree where each leaf stores an estimated semi-Markov model for the subset of trajectories in that node.

## 5. Random forest and variable importance

The file [src/random_forest/random_forest.R](../src/random_forest/random_forest.R) creates a bootstrap forest of trees.

Key parts:

- `random_forest()`: builds many trees from bootstrap samples
- each tree includes out-of-bag IDs (`oob_ids`) for later importance calculations

The file [src/random_forest/variable_importance.R](../src/random_forest/variable_importance.R) computes variable importance using two measures:

- `MDI_tree()`, `MDI()`, `MDI_all()`: mean decrease in impurity
- `oob_score()`, `MDA()`, `MDA_parallelized()`, `MDA_all()`: mean decrease in accuracy on out-of-bag samples

## 6. Data preparation and analysis scripts

The project also contains analysis scripts:

- [src/data/extract_data.R](../src/data/extract_data.R): data loading and preparation
- [src/data/data_analysis.R](../src/data/data_analysis.R): exploratory analysis of the raw dataset
- [src/data/covariates_analysis.R](../src/data/covariates_analysis.R): covariate screening and preprocessing
- [src/visualization.R](../src/visualization.R): plotting utilities

These scripts are not the core statistical engine, but they support data inspection, visualization, and exploratory analysis.

## 7. Entry points and workflow

The main workflow is driven by:

- [src/main.R](../src/main.R): script that loads data, builds a tree, and runs forest-based variable importance computations

Typical execution flow:

1. Load trajectory data and covariates
2. Fit a semi-Markov model using `mle_fit()`
3. Use `permutation_test()` or related tests to compare groups
4. Split the data with `build_tree()`
5. Build a forest with `random_forest()`
6. Rank covariates with `MDI_all()` or `MDA_all()`

## 8. Simulation experiments

The folder [simulations](../simulations) contains experiment scripts used to benchmark statistical properties and computation time.

Examples:

- [simulations/nominal_level_test.R](../simulations/nominal_level_test.R)
- [simulations/power_test.R](../simulations/power_test.R)
- [simulations/time_test.R](../simulations/time_test.R)
- [simulations/tree_RF.R](../simulations/tree_RF.R)

These scripts are research/validation utilities rather than core library code.

## 9. Good reading order

If you are new to the project, read in this order:

1. [src/semi_markov/mle_estimation.R](../src/semi_markov/mle_estimation.R)
2. [src/two_samples_test.R](../src/two_samples_test.R)
3. [src/random_forest/tree_construction.R](../src/random_forest/tree_construction.R)
4. [src/random_forest/random_forest.R](../src/random_forest/random_forest.R)
5. [src/random_forest/variable_importance.R](../src/random_forest/variable_importance.R)
6. [src/main.R](../src/main.R)

