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
using JuMP


#### info: script to compare warping window size on chronology
# config
path = "/cluster/home/danare/git/GENeSYS_MOD.data/Output/output_excel/Timeseries.xlsx"
config =  "/cluster/home/danare/git/Clustering/data/config5.yml"
k = 10
technology = ["TS_LOAD", "TS_WIND_ONSHORE_AVG", "TS_WIND_OFFSHORE_DEEP", "TS_PV_AVG"]
region = "DE"
warping_window = Dict(
    "Euclidean" => 0,
    #"DTW (1)" => 1,
    "DTW (5)"  => 5,
    #"DTW (7)"  => 7,
    #"DTW (9)"  => 9,
    "DTW (10)"=> 10,
    
);
total_ts = true

# read in data
config = TSClustering.read_yaml_file(file=config);
data_org = TSClustering.read_data(path=path, config=config);
data_clustering_org = TSClustering.create_clustering_matrix(config=config, CountryData=data_org);
# apply normalization
data = TSClustering.normalize_data(CountryData=data_org, config=config);
data_clustering = TSClustering.create_clustering_matrix(config=config, CountryData=data);

p = collect(keys(warping_window))
append!(p, ["Original_Data"])


# define empty array
data_array = JuMP.Containers.DenseAxisArray(zeros(length(p),length(technology),8760), p, technology, 1:8760) 

# add the remaining data
for (j, w) in enumerate(keys(warping_window))
    # define distance matrix
    D = TSClustering.define_distance(w=warping_window[w], data_clustering=data_clustering, fast_dtw=false)
    result = hclust(D, linkage=:complete)
    cl = cutree(result, k=k)
    weights = Dict{Int64, Int64}()
    for i in cl
        weights[i] = get(weights, i, 0) + 1
    end

    # bring data in Jump format
    m_cluster_org = TSClustering.calculate_representative(representative=:medoid, data_clustering=data_clustering_org, cl=cl, weights=weights, k=k);
    cluster_dict_org = TSClustering.convert_data(k=k, config=config, M=m_cluster_org);
    data_org = TSClustering.read_data(path=path, config=config);
    sc = TSClustering.scaling(data_org=data_org, scaled_clusters=cluster_dict_org, k=k, weights=weights, config=config)

    for t ∈ technology
        data_array[w,t,:] = sort(TSClustering.upsample_time_series(weight=weights, cluster_dict=sc, technology=t, region=region, config=config, tot_sum=total_ts), rev=true)
    end
end


data_org_n = TSClustering.read_data(path=path, config=config);

# add original data
for t ∈ technology
    if total_ts
        # sum per region
        ts = 0
        for r in config["countries"]
            if sum(data_org_n[t][:,r]) > 0.0000000000
                ts+=1
            end
        end
        data_array["Original_Data",t,:]  = sort(sum(data_org_n[t][:,r] for r in config["countries"])/ts,  rev=true)
    else
        data_array["Original_Data",t,:] = sort(data_org_n[t][:,region], rev=true)
    end
end


# plot the results
TSClustering.plot_duration_curve(list_data=data_array, technology=technology, write_html=true)