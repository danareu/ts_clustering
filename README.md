# ts_clustering

> **Time-Series Clustering Framework for [GENeSYS-MOD](https://github.com/GENeSYS-MOD/GENeSYS_MOD.jl)**

Aggregates hourly temporal profiles (e.g., solar PV capacity factors, wind capacity factors, electricity demand) into a set of representative days. Importantly, the input data must follow the given format as outlined in [GENeSYS-MOD.data]( https://github.com/GENeSYS-MOD/GENeSYS_MOD.data)**

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

---

## Background

Energy system models like GENeSYS-MOD require full-year hourly profiles (8,760 time steps) as inputs. Solving an optimisation at this resolution is computationally expensive. Time-series clustering reduces those 8,760 hours to a small set of *k* representative periods while preserving key statistical properties of the original profiles.

---

## Features

| Feature | Options |
|---|---|
| **Clustering algorithm** | Hierarchical – Ward's linkage, Complete linkage |
| **Distance metric** | Euclidean, Dynamic Time Warping (DTW) |
| **DTW implementation** | FastDTW (approximate), Dynamic Programming (exact) |
| **Normalisation** | z-score (optional) |
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
---

## Configuration

Runs are controlled by a YAML configuration file:

```yaml
Country_Data_Entries: # the attributes in the xlsx files for which representative periods should be found
  - TS_LOAD
  - TS_PV_AVG          
Clustering: # the attributes in the xlsx files that should be used for the clustering algorithm
  - TS_LOAD
  - TS_PV_AVG
countries: # the country profiles for each attribute that should be considered
  - DE
SCTOLERANCE: 10.0e-6 # the max. accuracy between the mean of the aggregated and the original time-series
Load: # the profiles that represent a load and are therefore normalized 
  - TS_LOAD
```

## Usage

```julia
using Pkg; Pkg.activate(".")

include("TSClustering.jl")
using .TSClustering
```


### Step-by-step pipeline

**1. Load data**


```julia
# Hourly profiles are read. A dictionary with the time-series attributes is generated using key_mappings which maps the sheet names to attributes 𝓣
hourly_data = XLSX.readxlsx(joinpath(inputdir, hourly_data_file * ".xlsx"))
CountryData = Dict(t => DataFrame(XLSX.gettable(hourly_data[keys_mapping[t]])) for t ∈ 𝓣)
```

**2. Normalise & build the clustering matrix**

```julia
data = TSClustering.normalize_data(config=config, CountryData=CountryData)

# With PCA (if pca_path is set):
pca_res = TSClustering.derive_principal_components(config=config, CountryData=CountryData, technology=pca_ts)
data_clustering = TSClustering.create_clustering_matrix(
    technology=["PCA"],
    CountryData=Dict("PCA" => DataFrame(Matrix(pca_res), :auto))
)

# Without PCA — seasonal profiles excluded from the clustering matrix:
data_clustering = TSClustering.create_clustering_matrix(technology=setdiff(𝓣, seasonal), CountryData=data)
```

**3. Compute the distance matrix**

```julia
D = TSClustering.define_distance(w=warping_window, data_clustering=data_clustering, fast_dtw=false)
```

**4. Cluster**

```julia
# Hierarchical Ward:
result = hclust(D, linkage=:ward)
cl = cutree(result, k=clusters)

# K-Means (when "Kmeans" ∈ resultdir):
Random.Seed(1)
R = kmeans(data_clustering, clusters; maxiter=200, display=:iter)
cl = assignments(R)
```

**5. Compute representative profiles**

Three methods are available, controlled by `hoffmann` and `resultdir`:

```julia
# Hoffmann — optimised hourly distribution per cluster:
sc1 = TSClustering.calculate_representative_value_distribution(
    data_org=filter(kv -> kv[1] ∉ seasonal, CountryData), cl=cl, config=config, K=clusters
)

# Medoid — most central day per cluster (default):
cluster_dict_org = TSClustering.calculate_medoid(
    data_org=CountryData, cl=cl, config=config, K=clusters, technology=𝓣
)
sc = TSClustering.scaling(
    data_org=CountryData, scaled_clusters=cluster_dict_org,
    k=clusters, weights=weights, config=config, technology=𝓣
)

# Centroid — mean profile per cluster:
cluster_dict_org = TSClustering.calculate_centroid(
    data_org=CountryData, cl=cl, config=config, K=clusters, technology=𝓣
)
```

> **Seasonal profiles:** Technologies listed in `seasonal` are always represented via medoid and are excluded from the DTW clustering step. They are handled separately via `TSClustering.scaling` after the main clustering is done.

> **Heat pumps:** `HLR_Heatpump_Aerial` and `HLR_Heatpump_Ground` are written to `TimeDepEfficiency` rather than `CapacityFactor`, reflecting their time-dependent COP.

**6. Write outputs**

Set `write_reduced_timeserie = 1` to write results to:

```
{inputdir}/input_reduced_timeseries_{k}_{resultdir_stem}.xlsx
```

Sheets written: `SpecifiedDemandProfile`, `CapacityFactor`, `TimeDepEfficiency`, `YearSplit`.

---



## Outputs

The pipeline returns the following JuMP `DenseAxisArray` objects:

| Return value | Axes | Description |
|---|---|---|
| `SpecifiedDemandProfile` | `[Region, Fuel, Timeslice, Year]` | Normalised hourly demand profiles; each cluster's 24 values sum to `weights[k] / 365` |
| `CapacityFactor` | `[Region, Technology, Timeslice, Year]` | Hourly capacity factors for renewable technologies |
| `TimeDepEfficiency` | `[Region, Technology, Timeslice, Year]` | Time-dependent COP for heat pumps |
| `YearSplit` | `[Timeslice, Year]` | Fraction of the year per timeslice (`weights[k] / 8760`, repeated 24× per cluster) |
| `cl` | `Vector{Int}` | Cluster assignment for each of the 365 input days |
| `weights` | `Dict{Int, Int}` | Number of days assigned to each cluster (sums to 365) |

---

## Methodology

### Hierarchical clustering (default)

Agglomerative clustering with Ward's minimum-variance linkage via [Clustering.jl](https://github.com/JuliaStats/Clustering.jl). A full 365×365 DTW distance matrix is computed once and the resulting tree is cut at depth *k*.

### Dynamic Time Warping

DTW allows non-linear temporal alignment, making it more robust than Euclidean distance for renewable profiles that are structurally similar but temporally shifted. `warping_window` constrains the maximum allowed shift, balancing accuracy against compute time. DTW is implemented via [DynamicAxisWarping.jl](https://github.com/baggepinnen/DynamicAxisWarping.jl).

### Representative value methods

- **Medoid** — selects the actual day closest to each cluster centre; always produces physically plausible profiles.
- **Centroid** — averages all days in a cluster; smoother but can produce unrealistic values.
- **Hoffmann** — selects the representative which finds the hourly distribution that best preserves the statistical properties of the original data. The method was published in https://www.sciencedirect.com/science/article/abs/pii/S0306261922004342

### Demand profile normalisation

`SpecifiedDemandProfile` values within each cluster are normalised so the 24 hourly values sum to `weights[k] / 365`, preserving each cluster's correct share of annual energy. Any `NaN` values from zero-sum clusters are set to zero with a warning.

---

## Dependencies

| Package | Role |
|---|---|
| [Clustering.jl](https://github.com/JuliaStats/Clustering.jl) | `hclust`, `cutree`, `kmeans` |
| [Distances.jl](https://github.com/JuliaStats/Distances.jl) | Distance metric primitives |
| [DynamicAxisWarping.jl](https://github.com/baggepinnen/DynamicAxisWarping.jl) | DTW distance computation |
| [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) | Tabular data manipulation |
| [XLSX.jl](https://github.com/felipenoris/XLSX.jl) | Input/output Excel files |
| [YAML.jl](https://github.com/JuliaData/YAML.jl) | Configuration parsing |
| [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl) | PCA |
| [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) | Interactive result visualisation |
| [Statistics.jl](https://docs.julialang.org/en/v1/stdlib/Statistics/) | Mean, std, etc. |
| [DelimitedFiles.jl](https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/) | Text file I/O |
| [Dates.jl](https://docs.julialang.org/en/v1/stdlib/Dates/) | Runtime benchmarking |

---
