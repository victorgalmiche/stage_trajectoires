# Life course trajectory segmentation with semi-Markov processes

Report available at: https://www.overleaf.com/read/kspwhnrtcdtg#9d5b53

## Code documentation

This repository contains the implementation of a semi-Markov trajectory model with tree-based feature selection and random-forest variable importance.

For a structured overview of the project architecture, data model, and main functions, see [doc/code_overview.md](doc/code_overview.md).

## Repository structure

- [src/semi_markov/mle_estimation.R](src/semi_markov/mle_estimation.R): semi-Markov likelihood and parameter estimation
- [src/two_samples_test.R](src/two_samples_test.R): two-sample tests for trajectory comparison
- [src/random_forest/tree_construction.R](src/random_forest/tree_construction.R): regression tree construction
- [src/random_forest/random_forest.R](src/random_forest/random_forest.R): random forest generation
- [src/random_forest/variable_importance.R](src/random_forest/variable_importance.R): variable importance measures
- [simulations](simulations): experimental scripts for nominal-level, power, and timing studies

## Typical workflow

1. Prepare the trajectory dataframe with columns `id`, `state`, and `time`.
2. Estimate the semi-Markov parameters with `mle_fit()`.
3. Compare groups using the likelihood-ratio or permutation-based tests.
4. Build and analyze trees or random forests to identify important covariates.

See [doc/code_overview.md](doc/code_overview.md) for the detailed code map and reading order.

