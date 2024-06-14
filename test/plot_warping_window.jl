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



#### info: script to compare warping window size on chronology

# config
path = "/cluster/home/danare/git/GENeSYS_MOD.data/Output/output_excel/Timeseries.xlsx"
country = "DE"
list_window = [0,1,5,10] 
config = "/cluster/home/danare/git/Clustering/data/config5.yml"

# read in data
config = TSClustering.read_yaml_file(file=config);
data_org = TSClustering.read_data(path=path, config=config);
data_clustering_org = TSClustering.create_clustering_matrix(config=config, CountryData=data_org);
# apply normalization
data = TSClustering.normalize_data(CountryData=data_org, config=config);
data_clustering = TSClustering.create_clustering_matrix(config=config, CountryData=data);


for k in [6]
    list_results = []
    for w in list_window
        D = TSClustering.define_distance(w=w, data_clustering=data_clustering, fast_dtw=false)
        result = hclust(D, linkage=:complete)
        cl = cutree(result, k=k)
        push!(list_results, cl)
    end
    TSClustering.map_to_org(list_df=list_results,list_window=list_window, write_html=true, name="Chronology_$(Dates.now()).html")
end