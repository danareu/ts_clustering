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
config = "/cluster/home/danare/git/Clustering/data/config.yml"
K = [6]
W = Dict(
    "Euclidean" => 0,
    "DTW (10)"  => 10,
);

# empty DataFrame
df = DataFrame(Cluster = Float64[], Window = Float64[], Technology =String[], Value = Vector{Vector{}}())

# read in data
config = TSClustering.read_yaml_file(file=config);
data_org = TSClustering.read_data(path=path, config=config);
data_clustering_org = TSClustering.create_clustering_matrix(config=config, CountryData=data_org);
# apply normalization
data = TSClustering.normalize_data(CountryData=data_org, config=config);
data_clustering = TSClustering.create_clustering_matrix(config=config, CountryData=data);

# calculate the distance matrix
for w ∈ values(W)
    D = TSClustering.define_distance(w=w, data_clustering=data_clustering)
    result = hclust(D, linkage=:ward)
    for k ∈ K 
        cl = cutree(result, k=k)
        weights = Dict()
        for i in cl
            weights[i] = get(weights, i, 0) + 1
        end
    
        # bring data in Jump format
        m_cluster_org, mapping_org_data = TSClustering.calculate_representative(representative=:medoid, data_clustering=data_clustering_org, cl=cl, weights=weights, k=k);
        cluster_dict_org = TSClustering.convert_data(k=k, config=config, M=m_cluster_org);
        data_org = TSClustering.read_data(path=path, config=config);
        sc = TSClustering.scaling(data_org=data_org, scaled_clusters=cluster_dict_org, k=k, weights=weights, config=config)

        # sum the region

        technology_map = Dict(
            "TS_WIND_ONSHORE_AVG" => ["TS_WIND_ONSHORE_OPT", "TS_WIND_ONSHORE_AVG", "TS_WIND_ONSHORE_INF", "TS_WIND_ONSHORE_OPT"],
            "TS_WIND_OFFSHORE" => ["TS_WIND_OFFSHORE", "TS_WIND_OFFSHORE_SHALLOW", "TS_WIND_OFFSHORE_DEEP"],
            "TS_PV_AVG" => ["TS_PV_AVG", "TS_PV_INF", "TS_PV_TRA", "TS_PV_OPT"]
        )

        z = zeros(length(axes(sc)[3])*length(axes(sc)[4]),1)
        for ts ∈ ["TS_WIND_OFFSHORE", "TS_WIND_ONSHORE_AVG", "TS_PV_AVG"] 
            for t ∈ technology_map[ts] 
                for r ∈ axes(sc)[1]
                    m = vcat(sc[r, t, :,:])
                    z += reshape(m, length(axes(sc)[3])*length(axes(sc)[4]))
                end
            end
            push!(df, [k, w, ts, vec(z)])
        end
    end 
end 

# add org data
technology_map = Dict(
    "TS_WIND_ONSHORE_AVG" => ["TS_WIND_ONSHORE_OPT", "TS_WIND_ONSHORE_AVG", "TS_WIND_ONSHORE_INF", "TS_WIND_ONSHORE_OPT"],
    "TS_WIND_OFFSHORE" => ["TS_WIND_OFFSHORE", "TS_WIND_OFFSHORE_SHALLOW", "TS_WIND_OFFSHORE_DEEP"],
    "TS_PV_AVG" => ["TS_PV_AVG", "TS_PV_INF", "TS_PV_TRA", "TS_PV_OPT"]
)


data_org = TSClustering.read_data(path=path, config=config);
for ts ∈ ["TS_WIND_OFFSHORE", "TS_WIND_ONSHORE_AVG", "TS_PV_AVG"] 
    z = zeros(8760,1)
    for t ∈ technology_map[ts] 
        for r ∈ eachcol(data_org[t])
            m = data_org[t][:,r]
            z += m
        end
    end
    push!(df, [6, 99999, ts, vec(z)])
end



TSClustering.plot_box_technology(df=df)