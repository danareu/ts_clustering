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
config_list = "/cluster/home/danare/git/Clustering/data/config5.yml"


    # empty DataFrame
    df = DataFrame(Cluster = Float64[], Window = Float64[], Value = Float64[])

    # read in data
    config = TSClustering.read_yaml_file(file=config);



    # calculate the distance matrix
    for w ∈ values(W)
        data_org = TSClustering.read_data(path=path, config=config);
        ClusteredData = JuMP.Containers.DenseAxisArray(zeros(length(config["countries"]), length(config["Country_Data_Entries"]), k,24), config["countries"], config["Country_Data_Entries"], 1:k, 1:24) 




            # calculate rmse
            rmse = TSClustering.calculate_rmse(data_org=data_org, cluster_dict=sc, cl=cl)
            # insert new row dataframe 
            push!(df, [k w rmse])
        end 
    end 

    # write file as csv 
    CSV.write("results/rmse_$name.csv", df)
end

