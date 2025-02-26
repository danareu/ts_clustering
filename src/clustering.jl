function cluster_kmeans(;k::Integer, x::Matrix{Float64})
    k = k
    R = kmeans(x, k)
    a = assignments(R) # get the assignments of points to clusters
    c = counts(R) # get the cluster sizes
    M = R.centers
    return ClusterResults(R,a,c,M)
end


function convert_data(; 
    k::Integer, 
    config::Dict, M,  
    technology::Vector{String},)

    #clustered data

    ClusteredData = JuMP.Containers.DenseAxisArray(zeros(length(config["countries"]), length(config["Country_Data_Entries"]), k,24), config["countries"], config["Country_Data_Entries"], 1:k, 1:24) 
    for d in 1:k
        for (m,t) in enumerate(technology)
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
- `representative::Symbol`: A symbol indicating the method to calculate representatives. Possible values are :medoid. The default is the mean .
- `data_clustering::Matrix{Float64}`: A matrix containing the clustered data points.
- `cl::Vector{Int64}`: A vector containing the cluster assignments for each data point.
- `weights::Dict`: A dictionary containing weights for each cluster.
- `k::Integer`: The number of clusters.

Returns:
- `m_cluster::Matrix{Float64}`: A matrix containing the representative data points.

"""
function calculate_representative(; 
    representative::Symbol, 
    data_clustering::Matrix{Float64}, 
    cl::Vector{Int64}, 
    weights::Dict, 
    k::Integer)

    m_cluster = zeros(size(data_clustering)[1], k)

    #TODO facilitate
    for (j, i) in enumerate(cl)
        m_cluster[:,i] += data_clustering[:,j]
    end
    for i in 1:k
        m_cluster[:,i] /= weights[i]
    end

    if representative == :medoid
        df = DataFrame(Day=1:365, Cluster=cl, Distance=[sqeuclidean(m_cluster[:,i], data_clustering[:,j]) for (j,i) in enumerate(cl)])
        for i in 1:k
            filtered_df = filter(row -> row[:Cluster] == i, df)
            #Sort values by column A in descending order
            min_day = sort(filtered_df, :Distance, rev=false)[1,:Day]
            m_cluster[:,i] = data_clustering[:,min_day]
        end
    end   
    return m_cluster
end



function define_distance(; 
    w::Integer, 
    data_clustering::Matrix{Float64}, 
    fast_dtw=true,
    #kwargs...
    )

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
                        #z[i, j] = abs(cor(data_clustering[:,i],data_clustering[:,j]))
                        #z[i, j] = fastdtw(data_clustering[:,i],data_clustering[:,j],corr_dist(),w)[1]
                    # elseif :correlation in keys(kwargs)
                    #     corr = cor(data_clustering[:,i],data_clustering[:,j])
                    #     corr = (corr + corr') / 2
                    #     fill!(diagm(corr), 1.0)
                    #     z[i, j] =  1 - abs(corr)
                    else
                        d = DTW(radius=w)
                        z[i, j] = d(data_clustering[:,i],data_clustering[:,j])
                    end
                end
                z[j, i] = z[i, j] 
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



"""
Calculate representative value distribution for clustered data.

Args:
    data_clustering (Matrix{Float64}): Matrix of clustered data.
    cl (Vector{Int64}): Vector indicating cluster assignments.
    config (Dict): Configuration dictionary containing parameters.
    K (Integer): Number of clusters.

Returns:
ClusteredData. JuMP.Containers.DenseAxisArray

The function calculates representative value distributions for each cluster (`k`) 
according to 10.1016/j.apenergy.2022.119029

"""

function calculate_representative_value_distribution(;    
    data_org::Dict, 
    cl::Vector{Int64}, 
    config::Dict,
    K::Integer)

    ClusteredData = JuMP.Containers.DenseAxisArray(zeros(length(config["countries"]), length(keys(data_org)), K,24), config["countries"], keys(data_org), 1:K, 1:24) 

    for t ∈ keys(data_org), c ∈ config["countries"], k ∈ 1:K

        m = zeros(24, count(x -> x == k, cl))
        for (i,j) ∈ enumerate(findall(x -> x == k, cl))
            m[:,i] = data_org[t][(j-1)*24+1:j*24, c]
        end
        # first load duration curve & then mean value
        duration = sort(vec(m), rev=true)
        dur_avg = mean.(Iterators.partition(duration, size(m)[2]))
        # first mean of clusters and duration curve
        m_avg = vec(mean(m, dims=2))
        m_avg_index = sortperm(m_avg, rev=true)
        # calculate representative 
        z = hcat(m_avg_index, dur_avg)
        value = z[sortperm(z[:, 1]),:][:,2]
        ClusteredData[c,t,k,:] = value
    end
    return ClusteredData
end






function denormalized_data!(col)
    # normalize the data
    mean = mean(col)
    std = std(col)
    return (col .* std) .+ mean
end



function calculate_medoid(; 
    data_org::Dict,
    cl::Vector{Int64}, 
    config,
    K::Integer,
    technology::Vector{String})

    ClusteredData = JuMP.Containers.DenseAxisArray(zeros(length(config["countries"]), length(technology), K,24), config["countries"], technology, 1:K, 1:24) 

    ## alternative: calculate medoid for each feature

    for t ∈ technology, c ∈ config["countries"], k ∈ 1:K
        # Indices for the current cluster
        indices = findall(x -> x == k, cl)
        num_points = length(indices)
        m = zeros(24, num_points)

        for (i, idx) ∈ enumerate(indices)
            m[:, i] = data_org[t][(idx - 1) * 24 + 1 : idx * 24, c]
        end

        # Mean vector of the cluster
        m_avg = vec(mean(m, dims=2))

        # Compute squared Euclidean distances
        distances = [sqeuclidean(m[:, i], m_avg) for i in 1:num_points]

        # Find the medoid
        min_day_idx = argmin(distances)
        ClusteredData[c, t, k, :] = m[:, min_day_idx]
    end 
    return ClusteredData
end


