using Pkg
include("/cluster/home/danare/git/Clustering/TSClustering.jl")
cd("/cluster/home/danare/git/Clustering")
Pkg.activate(".")
using .TSClustering
using Distances
using Clustering
using DelimitedFiles
using DataFrames
using DynamicAxisWarping
using Dates
using Statistics
using CSV
using PlotlyJS

# config
path = "/cluster/home/danare/git/GENeSYS_MOD.data/Output/output_excel/Timeseries.xlsx"
config = "/cluster/home/danare/git/Clustering/data/config3.yml"

config = TSClustering.read_yaml_file(file=config);
data_org = TSClustering.read_data(path=path, config=config)


p = make_subplots(rows=3, 
cols=2, 
)

for (w, i, j, k) in zip(["TS_WIND_OFFSHORE_DEEP" , "TS_WIND_ONSHORE_AVG", "TS_PV_AVG", "TS_HYDRO_ROR", "TS_LOAD"], [1,1,2,2,3], [1,2,1,2,1], [true,false,false, false, false])
    z = cor(Matrix(data_org[w]))

    add_trace!(p,
    heatmap(z=z,
    colorscale="PdBu",
    showscale=k,
    y=names(data_org[w]),
    x=names(data_org[w]),),
    row=i, col=j)
end


relayout!(p, 
    yaxis_title = "TS_WIND_OFFSHORE_DEEP", 
    yaxis2_title = "TS_WIND_ONSHORE_AVG",
    yaxis3_title = "TS_PV_AVG",
    yaxis4_title = "TS_HYDRO_ROR",
    yaxis5_title = "TS_LOAD",
    xaxis_title = "TS_WIND_OFFSHORE_DEEP", 
    xaxis2_title = "TS_WIND_ONSHORE_AVG",
    xaxis3_title = "TS_PV_AVG",
    xaxis4_title = "TS_HYDRO_ROR",
    xaxis5_title = "TS_LOAD",
    ) 



savefig(p, "spatial_correlation_all_countries.html")