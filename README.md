# Optimal Learning in Fairness-Aware Contextual Bandits

Code for the paper by Gengyu Xue, Fan Wang, Feiyu Jiang, Zifeng Zhao and Yi Yu (2026).

## Files

| File | Purpose |
|---|---|
| `algorithms.R` | Fair greedy, polyhedral UCB, their one-dimensional implementations, Standard UCB, policy optimisation and training checkpoints. |
| `cost_of_fairness.R` | Oracle policies, population reward integration, cost of fairness and reward improvement. |
| `simulation_studies.R` | Baseline simulations, varying L and policy visualisation. |
| `mixed_design.R` | Mixed discrete–continuous simulations, cumulative and relative regret. |
| `real_data_analysis.R` | Adult and Jobs experiments and cumulative online accuracy. |

All computational helpers, oracle construction and dataset source information
are included in the five R files. Each experiment file loads its required functions automatically;
sourcing `algorithms.R` alone only loads functions.

Running the experiments creates `results/` for checkpoints, summaries and figures.
Real-data preparation also creates `data/` for downloaded and transformed records.
These directories are generated during execution.

## Installation

Use R 4.5 or newer. Install the R packages with:

```r
install.packages(c("Rcpp", "Rglpk", "slam", "ggplot2", "jsonlite",
                   "pbivnorm", "filelock", "scales"))
```

Rcpp requires a C++ compiler, and Rglpk requires GLPK. On macOS, install the Xcode
command-line tools and use compatible CRAN binaries where available. On
Debian/Ubuntu, the system packages are `build-essential` and `libglpk-dev`.
Windows requires the Rtools version matching R. The code has been checked on
macOS; Windows compilation has not been tested.

The R environment used was R 4.5.2 with Rcpp 1.1.1, Rglpk 0.6-5.1,
slam 0.1-55, ggplot2 4.0.3, jsonlite 2.0.0, pbivnorm 0.6.0,
filelock 1.0.3 and scales 1.4.0.

Python is needed for real-data preparation and the mixed numerical oracle.
These programs are included in the corresponding R files and invoked
automatically when needed. Use Python 3.12 with:

```sh
python3 -m venv .venv
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install "numpy>=1.24,<3" "pandas>=2.0,<3" "highspy>=1.15,<2"
```

Data preparation was checked with Python 3.9.6, NumPy 2.0.2 and pandas 2.3.3;
oracle construction used Python 3.12.14, NumPy 2.3.5 and HiGHS 1.15.1.

## Running experiments

From this directory:

```sh
Rscript cost_of_fairness.R
Rscript simulation_studies.R
Rscript mixed_design.R
Rscript real_data_analysis.R
```

Alternatively, in an R console, select an experiment file with:

```r
source(file.choose(), echo = FALSE)
```

To load functions without starting an experiment, set
`options(mab.release.autorun = FALSE)` before sourcing. Set it back to `TRUE`
to run experiments.

The default experiments are:

| Study | Settings |
|---|---|
| Cost | Uniform and truncated Gaussian covariates, homogeneous and heterogeneous rewards, d = 1, 2, 3, 4; 81 logarithmically spaced L values from 0.02 to 10. Deterministic population calculation. |
| Baseline | The same 16 settings; L = 0.5; T = 1024, 2048, 4096, 8192, 16384, 32768; 50 replications each. |
| Varying L | Uniform, homogeneous, d = 2; L = 0.25, 0.5, 1, 2, 4; the same six horizons and 50 replications. L = 0.5 reuses baseline runs. |
| Policy illustration | Arm 1 small, uniform, d = 2, L = 0.5; T = 8192 and 16384; replication 15, comparing greedy and Standard UCB with the oracle. |
| Mixed | Six categories and three continuous covariates, d = 4, L = 0.5; six horizons and 50 replications, scheduled in stages of 5, 10, 20 and 50 across horizons. |
| Real data | Adult and Jobs; six horizons, five L values and 50 replications per dataset and (T,L) pair. |

There are 6,004 simulation learning runs, 300 mixed runs and 3,000 real-data
runs. High-dimensional experiments can require substantial computation.
The paper figures use greedy learning. To run mixed polyhedral UCB instead:

```sh
Rscript mixed_design.R --algorithm=polyhedral_ucb
```

UCB results use a separate output directory. Its exploration rule is specified
separately from greedy's.

## Resume and parallel execution

Completed runs and intermediate training and reward-evaluation states are saved.
Training checkpoints are written every 128 rounds and before and after expensive
fair-policy updates. Interrupt a process and run the same command again to resume
from its last checkpoint. Use a new output directory if scientific settings change.

To partition a study between two terminals:

```sh
Rscript mixed_design.R --worker-id=1 --worker-count=2
Rscript mixed_design.R --worker-id=2 --worker-count=2
```

These worker arguments also apply to simulation and real-data learning. Use the
same output directory and worker count with distinct worker IDs. Seeds and
integration points are independent of worker assignment. Cost calculation uses
one process.

`--output=/path/to/results` selects an output directory. `--max-tasks=1` completes
at most one new task before returning. The corresponding R options are
`mab.release.output`, `mab.release.max_tasks`, `mab.release.worker_id` and
`mab.release.worker_count`.

## Figures

After experiments have saved results, redraw their figures with:

```sh
Rscript cost_of_fairness.R --plot-only
Rscript simulation_studies.R --plot-only
Rscript mixed_design.R --plot-only
Rscript real_data_analysis.R --plot-only
```

Figures and summaries are written to the corresponding `results/` subdirectory,
or the selected output directory. The commands use completed local runs or their
saved summaries. They do not start training; if no results are available, run the
corresponding experiment first.

| Study | PDF and PNG figures |
|---|---|
| Cost | `cost_of_fairness`, `cost_reduction` |
| Simulation | `ordinary_relative`, `ordinary_cumulative`, `vary_L_relative`, `vary_L_cumulative`, `policy` |
| Mixed | `mixed_d4_relative_regret`, `mixed_d4_cumulative_loglog` |
| Real data | `real_accuracy` |

Means and sample standard deviations are reported across replications. Shading
shows one sample standard deviation. The mixed and baseline figures use the same
d = 4 colour and formatting.

Each setting schedules 50 replications. Intermediate figures use the completed
replications, and summary tables record the actual count. Floating-point results
can vary with solver versions, platforms and numerical tolerances.

The policy illustration uses replication 15, selected among 20 pairs for its
large improvement, rather than a replication mean. Its two horizons are separate
experiments.

Relative regret is cumulative population reward loss divided by T times the
fair-oracle reward. Initial exploration is included. Baseline reward evaluation
uses analytical integration in d = 1, 2 and shared Halton points in d = 3, 4:
16,384 and 32,768 points respectively, skipping the first 128. Uniform cost
calculations use analytical formulas in all four dimensions. The mixed design
uses 32,768 shared Halton points and exact summation over six categories.

## Adult and Jobs: sources and preprocessing

Adult is downloaded from the [UCI Machine Learning Repository](https://archive.ics.uci.edu/dataset/2/adult).
Credit Barry Becker and Ronny Kohavi (1996), *Adult*,
[doi:10.24432/C5XW20](https://doi.org/10.24432/C5XW20), distributed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Jobs uses version 1 of Ayush Tankha's (2023)
[70k+ Job Applicants Data (Human Resource)](https://www.kaggle.com/datasets/ayushtankha/70k-job-applicants-data-human-resource).
If automatic download is unavailable, use the page's Download button, sign in
if requested, and extract `stackoverflow_full.csv` to `data/raw/`, then rerun
`Rscript real_data_analysis.R`. Kaggle's
[dataset metadata](https://www.kaggle.com/api/v1/datasets/view/ayushtankha/70k-job-applicants-data-human-resource)
lists its license as Unknown. Raw records are not redistributed here; the
providers' terms apply. Dataset URLs, input-file SHA-256 values and reference
processed-array hashes are included in `real_data_analysis.R`.

The real-data entry prepares missing data automatically. To prepare data alone,
or to choose different data directories, use:

```r
options(mab.release.autorun = FALSE)
source("real_data_analysis.R")
real_prepare_data(c("adult", "jobs"))
# Optional: raw_dir = "/path/to/raw", output = "/path/to/processed"
options(mab.release.autorun = TRUE)
```

Set `options(mab.release.data = "/path/to/processed")` before running experiments
if prepared files are stored elsewhere. `MAB_PYTHON` selects the Python executable.

**Sample selection.** For Adult, concatenate `adult.data` and `adult.test`,
treat `?` as missing and remove incomplete records. Retain United States, Mexico,
Philippines, Germany, Puerto Rico, Canada, El Salvador and India. The response
is one if income exceeds $50,000 and zero otherwise. The 12 covariates are age,
work class, `fnlwgt`, numeric education level, marital status, occupation,
relationship, race, sex, capital gain, capital loss and weekly working hours.
`fnlwgt` is a covariate, not an observation weight. Country and the text
education-category field are excluded.

For Jobs, remove `Unnamed: 0` and `HaveWorkedWith`, then discard incomplete records.
Use the binary `Employed` response and the 11 covariates `Age`, `Accessibility`,
`EdLevel`, `Employment`, `Gender`, `MentalHealth`, `MainBranch`, `YearsCode`,
`YearsCodePro`, `PreviousSalary` and `ComputerSkills`. Country is excluded and all
observed gender categories are retained. `Employment` and `Employed` are separate
supplied fields. The task classifies the supplied label.

**Transformation.** Select `floor(0.1*n)` preprocessing records without replacement
using NumPy `default_rng` with seeds 20260721 for Adult and 20260722 for Jobs.
The subset is fixed; its labels are not used to fit the transformation.

1. Identify numeric fields and categorical levels on the preprocessing subset.
   One-hot encode every observed level, without dropping a reference level.
   Unseen levels map to zeros in that field's indicator columns. Adult's
   categorical fields are work class, marital status, occupation, relationship,
   race and sex. Jobs uses age group, accessibility, education level, gender,
   mental health and main branch. Encoded dimensions are 46 and 21.
2. Subtract preprocessing means and divide by preprocessing sample standard
   deviations (`ddof=1`), using scale one for constant coordinates. Apply PCA
   to the covariance matrix and retain the three leading eigenvectors. Fix
   each sign so its largest absolute loading is positive.
3. Project every record onto these directions. For each component subtract its
   preprocessing minimum, divide by its preprocessing range, and clip to `[0,1]`.
   The encoding, means, scales, PCA directions and ranges are fitted only on the
   preprocessing subset. Fairness uses L1 distance in the resulting `[0,1]^3`.

Before the R splits, Adult rows are ordered as United States followed by the
other countries in the order listed above; Jobs uses non-US followed by US.
Original row order is retained within groups. This order determines seeded
sampling; the groups have no separate weighting or learning role.

The preparation outputs are `records.csv.gz`, `covariates.bin`,
`preprocessing_model.npz` and `manifest.json` for each dataset. Covariates are
stored as column-major little-endian 64-bit doubles. The manifest records the
source hashes, transformation and software versions. PCA scores can differ at
roundoff level across BLAS/LAPACK versions.

**Data splits and online sampling.** Exclude preprocessing records, sample a
fixed disjoint holdout of 512 records, and use the remainder as the online pool:

| Dataset | Eligible | Preprocessing | Held out | Online pool |
|---|---:|---:|---:|---:|
| Adult | 43,303 | 4,330 | 512 | 38,461 |
| Jobs | 73,462 | 7,346 | 512 | 65,604 |

For dataset index `ds` (Adult 1, Jobs 2), replication `r` and stream index `k`
(holdout 1, permutation 2, action 3), the R seed is
`(20260924 + ds*10000019 + r*100003 + k*7919) %% 2147483646 + 1`.
The holdout uses replication zero. R uses Mersenne-Twister, Inversion and
Rejection. Each replication permutes the online pool and uses the first T
records. Record order and action uniforms are shared across T and L within that
replication, while every (T,L) pair trains a separate learner.

Arm 1 predicts one and arm 2 predicts zero; rewards are the label and one minus
the label. Only the selected arm's reward updates the learner. No artificial
reward noise is added. Online accuracy averages realised rewards over all rounds,
including the initial `2*ceiling(sqrt(T*log(4*T)))` exploration rounds.
Set `options(mab.release.heldout=TRUE)` to additionally evaluate final-policy and
time-averaged expected accuracy on the holdout; it does not enter online updates.

## Mixed numerical oracle

The mixed experiment computes its numerical oracle before training. It uses the
six-category grid `{0,1/64,...,1}^3`, with one arm-selection probability in `[0,1]`
for each category and grid point. The linear programme maximises population reward
subject to adjacent continuous-grid differences bounded by `L/64` and adjacent
category differences bounded by `L/5`. The objective uses two-point tensor Gauss
quadrature, splitting at grid boundaries and discontinuities of the category
probabilities. Trilinear interpolation evaluates the policy between grid points
and preserves the L1-Lipschitz constraint.

The first mixed run solves this programme and saves its probabilities and
metadata under `results/oracle/`. Subsequent runs reuse the saved solution.
Parallel workers share a lock so only one computes the oracle. The programme has
1,647,750 variables and 6,240,325 edge constraints and can require substantial
memory and time. An interruption before the oracle has been saved requires
restarting that solve; completed oracle and learning results are reused.

To compute the oracle separately:

```r
options(mab.release.autorun = FALSE)
source("mixed_design.R")
mx_build_oracle(grid = 64L, threads = 1L)
options(mab.release.autorun = TRUE)
```

To use a computed oracle stored elsewhere:

```sh
Rscript mixed_design.R --oracle=/path/to/mixed_d4_grid64.bin
```

Its matching JSON metadata is required. The driver checks the environment,
grid, probability bounds, fairness and finite-grid objective gap, then records
the oracle checksum with the experiment settings. Numerical solutions may differ
across solvers or platforms, and an LP can have multiple optimal policies.

This grid policy approximates the population oracle; it is not an exact solution
to the infinite-dimensional problem. The mixed distribution does not satisfy
the continuous-density assumption, so comparison of its empirical slope with
2/3 is descriptive.

## Numerical settings and license

`MAB_NUMERICAL_TOLERANCE = 1e-6` bounds accepted probability, fairness,
boundary-gap and feasibility residuals. Optimisation uses the specified L.
Models with cache keys longer than 9,000 bytes are computed without caching.

The code is distributed under the MIT license in `LICENSE`. Dataset terms are
separate. Please cite *Optimal Learning in Fairness-Aware Contextual Bandits*,
Gengyu Xue, Fan Wang, Feiyu Jiang, Zifeng Zhao and Yi Yu (2026), when using the code.
