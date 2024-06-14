function cluster_kmeans(;k::Integer, x::Matrix{Float64})
    k = k
    R = kmeans(x, k)
    a = assignments(R) # get the assignments of points to clusters
    c = counts(R) # get the cluster sizes
    M = R.centers
    return ClusterResults(R,a,c,M)
end


function convert_data(; k::Integer, config::Dict, M)

    #clustered data

    ClusteredData = JuMP.Containers.DenseAxisArray(zeros(length(config["countries"]), length(config["Country_Data_Entries"]), k,24), config["countries"], config["Country_Data_Entries"], 1:k, 1:24) 
    for d in 1:k
        for (m,t) in enumerate(config["Country_Data_Entries"])
            for (b,c) in enumerate(config["countries"])
                a = max((m-1)*length(config["countries"])*48+(((b-1)*48)+1),1)
                ClusteredData[c,t,d,:] = M[a:a+23,d]
            end
        end
    end
    return ClusteredData
end




"""
Calculate representative data points based on the given method.

Parameters:
- `representative::Symbol`: A symbol indicating the method to calculate representatives. Possible values are :medoid.
- `data_clustering::Matrix{Float64}`: A matrix containing the clustered data points.
- `cl::Vector{Int64}`: A vector containing the cluster assignments for each data point.
- `weights::Dict`: A dictionary containing weights for each cluster.
- `k::Integer`: The number of clusters.

Returns:
- `m_cluster::Matrix{Float64}`: A matrix containing the representative data points.

"""
function calculate_representative(; representative::Symbol, data_clustering::Matrix{Float64}, cl::Vector{Int64}, weights::Dict, k::Integer)

    m_cluster = zeros(size(data_clustering)[1], k)

    #TODO facilitate
    for (j, i) in enumerate(cl)
        m_cluster[:,i] += data_clustering[:,j]
    end
    for i in 1:k
        m_cluster[:,i] /= weights[i]
    end

    if representative == :medoid
        df = DataFrame(Day=1:365, Cluster=cl, Distance=[euclidean(m_cluster[:,i], data_clustering[:,j]) for (j,i) in enumerate(cl)])
        #d = DTW(radius=5)
        #df = DataFrame(Day=1:365, Cluster=cl, Distance=[d(m_cluster[:,i], data_clustering[:,j]) for (j,i) in enumerate(cl)])

        for i in 1:k
            filtered_df = filter(row -> row[:Cluster] == i, df)
            # Sort values by column A in descending order
            min_day = sort(filtered_df, :Distance, rev=false)[1,:Day]
            m_cluster[:,i] = data_clustering[:,min_day]
        end
    end   
    return m_cluster
end



function define_distance(; w::Integer, data_clustering::Matrix{Float64}, fast_dtw=true)
    z = zeros(365,365)
    start = Dates.now()
    # iterate pairwise
    for i in 1:365
        for j in i+1:365  # Avoid redundant calculations (distance between i and j is the same as between j and i)
            if i == j
                z[i, j] = 0
                z[j, i] = 0
            else
                if w == 0
                    z[i,j] = sqeuclidean(data_clustering[:,i],data_clustering[:,j])
                else
                    if fast_dtw
                        z[i, j] = fastdtw(data_clustering[:,i],data_clustering[:,j],SqEuclidean(),w)[1]
                    else
                        d = DTW(radius=w)
                        z[i, j] = d(data_clustering[:,i],data_clustering[:,j])
                    end
                end
                z[j, i] = z[i, j]  # Symmetric matrix, so we fill both sides
            end
        end
    end
    ends = Dates.now()
    
    # compute and storage of runtime
    string = [Dict("Time"=>ends-start, "Size"=>size(data_clustering)[1]*size(data_clustering)[2], "warping path"=> w, "Method"=> fast_dtw)]
    open("/cluster/home/danare/git/Clustering/data/computation_time.txt", "a") do io
        writedlm(io, string)
        end
    return z
end








function denormalized_data!(col)
    # normalize the data
    mean = mean(col)
    std = std(col)
    return (col .* std) .+ mean
end
