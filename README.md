# ts_clustering

> **Time-Series Clustering Framework for [GENeSYS-MOD](https://github.com/GENeSYS-MOD/GENeSYS_MOD.jl)**

Aggregates high-resolution temporal profiles (solar PV capacity factors, wind capacity factors, electricity demand) into a compact set of representative periods, making large-scale energy system optimisation tractable without sacrificing temporal diversity.

---

## Table of Contents

- [Background](#background)
- [Features](#features)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Methodology](#methodology)
- [Outputs](#outputs)
- [Dependencies](#dependencies)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)

---

## Background

Energy system models like GENeSYS-MOD require full-year hourly profiles (8,760 time steps) as inputs. Solving an optimisation at this resolution is computationally prohibitive. Time-series clustering reduces those 8,760 hours to a small set of *k* representative periods while preserving key statistical properties of the original profiles — seasonal patterns, peak demand, and the correlation between generation and consumption.

---

## Features

| Feature | Options |
|---|---|
| **Clustering algorithm** | Hierarchical – Ward's linkage, Complete linkage |
| **Distance metric** | Euclidean, Dynamic Time Warping (DTW) |
| **DTW implementation** | FastDTW (approximate), Dynamic Programming (exact) |
| **Normalisation** | z-score (optional) |
| **Pre-processing** | Savitzky–Golay smoothing (optional) |
| **Visualisation** | Interactive plots via PlotlyJS |
| **Input format** | XLSX |
| **Configuration** | YAML |

---

## Repository Structure

```
ts_clustering/
├── TSClustering.jl          # Module entry point
├── Project.toml             # Dependency declarations
├── Manifest.toml            # Fully-resolved lockfile
├── src/
│   ├── clustering.jl        # Core clustering logic
│   └── plot_results.jl      # Visualisation (PlotlyJS)
├── utils/
│   ├── load_data.jl         # XLSX data ingestion
│   ├── datastructs.jl       # Custom types and data structures
│   └── post_processing.jl   # Representative profile extraction & weighting
└── data/                    # Input time-series files (XLSX)
```

---

## Installation

**Prerequisites:** Julia 1.8+

```julia
# 1. Clone the repo
# git clone https://github.com/danareu/ts_clustering.git
# cd ts_clustering

# 2. Activate the environment and install dependencies
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

`Pkg.instantiate()` will download and precompile all required packages from the locked versions in `Manifest.toml`. This may take a few minutes on first run.

---

## Configuration

Runs are controlled by a YAML configuration file:

```yaml
n_clusters:    12          # Number of representative periods to generate
normalization: true        # Apply z-score normalisation before clustering
algorithm:     ward        # "ward" | "complete"
distance:      euclidean   # "euclidean" | "fastdtw" | "dp"
data_path:     data/profiles.xlsx
output_path:   results/
```

| Parameter | Description |
|---|---|
| `n_clusters` | Number of clusters *k* |
| `normalization` | Standardise profiles to zero mean / unit variance before computing distances |
| `algorithm` | Linkage method for hierarchical clustering |
| `distance` | Distance metric used to build the dissimilarity matrix |
| `data_path` | Path to the XLSX input file |
| `output_path` | Directory for output files and plots |

---

## Usage

```julia
using Pkg; Pkg.activate(".")

include("TSClustering.jl")
using .TSClustering

TSClustering.run_clustering("config.yaml")
```

### Choosing an algorithm

| Method | When to use |
|---|---|
| `ward` | Default. Minimises within-cluster variance; produces compact, evenly sized clusters. |
| `complete` | Maximises separation between clusters; useful when outlier robustness matters more than compactness. |

### Choosing a distance metric

| Metric | Speed | Best for |
|---|---|---|
| `euclidean` | Fast | Profiles that are already temporally aligned |
| `fastdtw` | Moderate | Large datasets with mild temporal shifts |
| `dp` | Slow | Small datasets where exact DTW accuracy is required |

### Normalisation

Enable `normalization: true` when your input contains profiles on different scales (e.g., demand in MW alongside 0–1 capacity factors). Disable it when all profiles are already on the same scale and you want to cluster on absolute values.

---

## Methodology

### Hierarchical Clustering

The framework uses **agglomerative** hierarchical clustering (via [Clustering.jl](https://github.com/JuliaStats/Clustering.jl)). Starting with every time period as its own cluster, it iteratively merges the two closest clusters until *k* remain. The full merge history is captured in a dendrogram, which can guide selection of *k* via the elbow method.

### Dynamic Time Warping

DTW allows non-linear alignment between time series, matching peaks that occur at slightly different times. This makes it more robust than Euclidean distance when comparing renewable generation profiles that are structurally similar but temporally shifted. FastDTW achieves linear time and space complexity through a multi-resolution approximation; the DP variant computes exact distances using standard dynamic programming.

### Savitzky–Golay Smoothing

An optional Savitzky–Golay filter (polynomial fitting over a sliding window) can be applied before clustering to remove high-frequency noise from wind and PV profiles while preserving the underlying shape.

---

## Outputs

After a run, the following are written to `output_path`:

- **Representative profiles** — one aggregated time-series per cluster (medoid or mean)
- **Weights** — number of original hours each representative period stands for (sum = 8,760 for annual data)
- **XLSX tables** — formatted for direct import into GENeSYS-MOD
- **Interactive plots** — HTML files (open in any browser) including:
  - Dendrogram
  - Cluster assignments over time
  - Original vs. representative profile overlay
  - Within-cluster variance vs. *k* (elbow plot)

---

## Dependencies

| Package | Role |
|---|---|
| [Clustering.jl](https://github.com/JuliaStats/Clustering.jl) | Hierarchical clustering algorithms |
| [Distances.jl](https://github.com/JuliaStats/Distances.jl) | Euclidean and other distance metrics |
| [DynamicAxisWarping.jl](https://github.com/baggepinnen/DynamicAxisWarping.jl) | FastDTW and DP-DTW |
| [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) | Tabular data manipulation |
| [XLSX.jl](https://github.com/felipenoris/XLSX.jl) | Excel I/O |
| [YAML.jl](https://github.com/JuliaData/YAML.jl) | Configuration parsing |
| [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl) | Multivariate statistics |
| [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) | Interactive HTML visualisations |
| [SavitzkyGolay.jl](https://github.com/lnacquaroli/SavitzkyGolay.jl) | Profile smoothing |
| [JuMP.jl](https://github.com/jump-dev/JuMP.jl) + [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) | Optimisation in post-processing |
| [Statistics.jl](https://docs.julialang.org/en/v1/stdlib/Statistics/) | Standard statistical functions |
| [DelimitedFiles.jl](https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/) | CSV-style file I/O |

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `Pkg.instantiate()` fails | Check internet access; run `Pkg.update()` to refresh the General registry and retry. |
| XLSX file not found | Verify `data_path` in your config matches the actual file location relative to the project root. |
| Clustering returns 1 cluster | Confirm `n_clusters > 1` and that the input contains more than one distinct time period. |
| Plots don't display | PlotlyJS generates HTML files — open them in a web browser. |
| DTW is very slow | Switch from `dp` to `fastdtw`, or aggregate to daily resolution before clustering. |
| Out-of-memory errors | Reduce the number of simultaneous profiles, or switch to `euclidean` distance. |

---

## Contributing

1. Fork the repository
2. Create a feature branch from `main`
3. Make your changes and add tests where applicable
4. Open a pull request with a description of what changed and why

Please follow [Julia style conventions](https://docs.julialang.org/en/v1/manual/style-guide/) and add any new dependencies to `Project.toml`.

---

## Citation

If you use this framework in your research, please cite:

```bibtex
@software{reulein2023tsclustering,
  author  = {Reulein, Dana},
  title   = {ts\_clustering: Time-Series Clustering Framework for GENeSYS-MOD},
  year    = {2023},
  url     = {https://github.com/danareu/ts_clustering}
}
```

---

*Developed by [Dana Reulein](https://github.com/danareu), 2023.*
