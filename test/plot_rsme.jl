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


# config
path = "/cluster/home/danare/git/GENeSYS_MOD.data/Output/output_excel/Timeseries.xlsx"
config_list = ["/cluster/home/danare/git/Clustering/data/config$i.yml" for i in [3, 5]]
K = [1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,20,30, 40,50,100]
W = Dict(
    "Euclidean" => 0,
    "DTW (1)" => 1,
    "DTW (2)" => 2,
    "DTW (3)" => 3,
    "DTW (4)" => 4,
    "DTW (5)" => 5,
    "DTW (6)" => 6,
    "DTW (7)" => 7,
    "DTW (8)"  => 8,
    "DTW (9)" => 9,
    "DTW (10)"  => 10,
    "DTW (15)"  => 15,
    
);


for (config, name) in zip(config_list, ["spatial", "onenode"])
    # empty DataFrame
    df = DataFrame(Cluster = Float64[], Window = Float64[], Value = Float64[])

    # read in data
    config = TSClustering.read_yaml_file(file=config);

    data_org = TSClustering.read_data(path=path, config=config);
    data_clustering_org = TSClustering.create_clustering_matrix(config=config, CountryData=data_org);
    # apply normalization
    data = TSClustering.normalize_data(CountryData=data_org, config=config);
    data_clustering = TSClustering.create_clustering_matrix(config=config, CountryData=data);

    # calculate the distance matrix
    for w ∈ values(W)
        D = TSClustering.define_distance(w=w, data_clustering=data_clustering, fast_dtw=false)
        result = hclust(D, linkage=:ward)
        for k ∈ K 
            cl = cutree(result, k=k)
            weights = Dict()
            for i in cl
                weights[i] = get(weights, i, 0) + 1
            end
            
            # bring data in Jump format
            m_cluster_org = TSClustering.calculate_representative(representative=:medoid, data_clustering=data_clustering_org, cl=cl, weights=weights, k=k);
            cluster_dict_org = TSClustering.convert_data(k=k, config=config, M=m_cluster_org);
            data_org = TSClustering.read_data(path=path, config=config);
            sc = TSClustering.scaling(data_org=data_org, scaled_clusters=cluster_dict_org, k=k, weights=weights, config=config)
            # calculate rmse
            println(data_org)
            rmse = TSClustering.calculate_rmse(data_org=data_org, cluster_dict=sc, cl=cl)
            # insert new row dataframe 
            push!(df, [k w rmse])
        end 
    end 

    # write file as csv 
    CSV.write("results/rmse_$name.csv", df)
end

