# NMAPower: Power Evaluation for Bayesian Network Meta-Analysis

`NMAPower` provides a simulation-based framework to evaluate the empirical operating characteristics of Bayesian Network Meta-Analyses (NMAs). 
Whether you are designing a new set of clinical trials or evaluating the robustness of an existing evidence network, this package helps you estimate statistical power, treatment efficacy ranking, and estimation bias.

The package integrates seamlessly with the Bayesian modeling ecosystem (for example `gemtc` for CB-NMA) and utilizes parallel processing to handle intensive Markov chain Monte Carlo (MCMC) simulations efficiently.

📦 Installation
Since `NMAPower` relies on Bayesian MCMC sampling, you must install JAGS (Just Another Gibbs Sampler) on your system before installing the package in R.

* Download and install JAGS from mcmc-jags.sourceforge.io.

* Install `NMAPower` and its core dependencies in R:

```r
# Install dependencies
install.packages(c("gemtc", "foreach", "doParallel", "igraph"))

# Install NMAPower 
# devtools::install_github("yourusername/NMAPower")
library(NMAPower)
library(gemtc)
```



## Quick Start Guide

### Prospective Power Analysis (Arbitrary Networks)
If you are planning an NMA and want to know how much power a specific network design will yield, use `nma_power_sim()`. This function supports any arbitrary network topology.  

You simply define your network as an edge-list dataframe. If you leave an odds ratio (OR) as NA, the function will mathematically infer it assuming statistical consistency across the network.

```r
# Define a 4-node network: Standard of Care (A) vs Active Treatments (B, C, D)
design <- data.frame(
  t1 = c("A", "A", "A", "B"),
  t2 = c("B", "C", "D", "C"),
  k  = c(6,   6,   2,   3),         # Number of studies per direct comparison
  OR = c(1.2, 1.4, 1.1, NA)         # True Odds Ratio (B vs C inferred)
)

# Run the simulation evaluating the B vs C contrast
results_general <- nma_power_sim(
  S = 100,                          # Number of iterations
  target_contrast = c("B", "C"),    # Contrast of interest
  network_design = design,
  pi_base = 0.5,                    # Baseline event probability for Node A
  tau = 0.2                         # Between-study heterogeneity
)

print(results_general)
```


### Retrospective / Post-Hoc Analysis
If you already have a dataset and want to evaluate if your NMA was adequately powered, use `nma_power_posthoc()`.

First, let's simulate a raw dataset using the built-in `generate_NMA()` function, which creates trial data based on a specified number of studies, treatments, and underlying success probabilities.

```r
set.seed(123)

# Generate a mock dataset with 3 treatments and 10 studies
# Sample sizes are drawn from a default range of 50 to 200
mock_data <- generate_NMA(
  num_studies = 10,
  num_treatments = 3,
  treatment_effects = c(0.3, 0.5, 0.7)
)
```



Before evaluating power, you can visualize your dataset's network structure using `gemtc`'s built-in plotting capabilities:

```r
# Convert the mock data into a gemtc network format
# The dataset maps to expected columns: study.id, treatment.id, sample.size, response
dat_gemtc <- data.frame(
  study = as.character(mock_data$study.id),
  treatment = as.character(mock_data$treatment.id),
  sampleSize = as.numeric(mock_data$sample.size),
  responders = as.numeric(mock_data$response)
)

network <- mtc.network(data.ab = dat_gemtc)

# Plot the network topology
plot(network, 
     vertex.color = "lightblue", 
     edge.color = "darkgray",
     main = "NMA Network Topology")
```


Now, pass the dataset into `nma_power_posthoc()`. 
This function will automatically fit a Bayesian NMA to your data to extract the underlying parameters (effect sizes and heterogeneity), and then perform simulations that perfectly preserve your original network topology, arm configurations, and sample sizes.

```r
# Evaluate power for comparing Treatment 1 vs Treatment 2
posthoc_results <- nma_power_posthoc(
  data = mock_data, 
  target_contrast = c(1, 2), # Evaluates Treatment 2 vs Treatment 1[cite: 9]
  S = 50,
  verbose = FALSE # Suppresses JAGS console output[cite: 9]
)

print(posthoc_results)
```

