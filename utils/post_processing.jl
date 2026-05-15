"""
    scaling(; data_org, scaled_clusters, k, technology, weights, config) -> DenseAxisArray

Scale representative cluster profiles so that their weighted mean matches the
original annual mean for each `(country, technology)` pair.

For load-type technologies the profiles are normalised by their total sum.
For all others, a scaling factor `σ` is applied iteratively: values are first
scaled to preserve the annual sum, capped at 1, then rescaled again on the
remaining sub-unity values to recover the lost energy. A warning is printed
when the residual deviation exceeds `config["SCTOLERANCE"]`.

# Arguments
- `data_org::Dict`: Original hourly data keyed by technology.
- `scaled_clusters`: Representative profiles `[country, technology, cluster, hour]`, modified in-place.
- `k::Integer`: Number of clusters.
- `technology::Vector{String}`: Technologies to scale.
- `weights::Dict`: Number of days assigned to each cluster.
- `config::Dict`: Must contain `"Load"` and `"SCTOLERANCE"`.
"""
function scaling(; 
    data_org::Dict, 
    scaled_clusters, 
    k::Integer, 
    technology::Vector{String},
    weights::Dict, 
    config::Dict)
    # idea: weighted time-series has the same avergae as the original time-series

    list_diff = []

    for t in technology
        if t ∈ config["Load"]
            # normalize to max value
            for c in names(data_org[t])
                scaled_clusters[c,t,:,:] /=  sum(scaled_clusters[c,t,:,:])
            end
        else
            # make mean value
            for c in names(data_org[t])
                sum_org = sum(data_org[t][:,c])
                max_org = maximum(data_org[t][:,c])
                min_org = minimum(data_org[t][:,c])

                # determine sum clustering
                sum_cl = 0
                for d in axes(scaled_clusters)[3]
                    sum_cl += sum(scaled_clusters[c,t,d,:]) * weights[d]
                end
                
                # scale time-series such that it has the original average
                σ = (sum_org / sum_cl)
                scaled_clusters[c,t,:,:] *= σ

                # remove values above 1
                for d in axes(scaled_clusters)[3], h in axes(scaled_clusters)[4]
                    if scaled_clusters[c,t,d,h] > 1
                        scaled_clusters[c,t,d,h] = 1
                    end
                end   
                
                # scale remaining values
                sum_cl_new = 0
                for d in axes(scaled_clusters)[3]
                    sum_cl_new += sum(scaled_clusters[c,t,d,:]) * weights[d]
                end

                σ = (sum_org / sum_cl_new) 

                for d in axes(scaled_clusters)[3], h in axes(scaled_clusters)[4]
                    if scaled_clusters[c,t,d,h] < 1
                        scaled_clusters[c,t,d,h] *= σ
                    elseif isnan(scaled_clusters[c,t,d,h])
                        scaled_clusters[c,t,d,h] = 0
                    end
                end  

                # check tolerance level
                mean_cl=0
                for d in 1:k
                    mean_cl += sum(scaled_clusters[c, t,d,:]*(weights[d]/8760))
                end
                if abs(mean_cl - mean(data_org[t][:,c])) > config["SCTOLERANCE"]
                    println("Deviation abs.: $(abs(mean_cl - mean(data_org[t][:,c]))) for $c and $t")
                end
                push!(list_diff, abs(mean_cl/8760-sum(data_org[t][:,c])/8760))
            end
        end
    end
    return scaled_clusters
end



function make_vector(; array:: JuMP.Containers.DenseAxisArray)
    return vec(reshape(array', 1, :))
end


"""
    calculate_variance_per_cluster(; CountryData::Dict, k::Integer, country::String) -> Dict

Calculate the variance for different time series data per cluster and return a dictionary containing these variances.

# Arguments
- `CountryData::Dict`: A dictionary where the keys are strings representing different time series (e.g., "TS_LOAD", "TS_HEAT_LOW", "TS_HEAT_HIGH") and the values are 2D arrays. Each array's rows correspond to different time points, and the columns correspond to different countries.
- `k::Integer`: The number of clusters.
- `cl::Vector{Int64}`: Mapping which days belong to which clusters.
- `country::String`: The name of the country for which the variance is to be calculated.

# Returns
- `Dict`: A dictionary where each key is an integer representing a cluster and the value is another dictionary. This nested dictionary has the same keys as `CountryData`, and each key maps to a list of variances calculated for the specified country.

"""

function calculate_variance_per_cluster(; CountryData:: Dict, k::Integer, country::String, cl::Vector{Int64})

    v = Dict(t => [] for t in collect(keys(CountryData)))
    temp_dict = Dict(d => Dict(t => [] for t in collect(keys(CountryData))) for d in 1:k)
    
    # Calculate and append variance
    for (k, m) in enumerate(0:24:8735)
        for t in collect(keys(CountryData))
            y = var(CountryData[t][m+1:m+24, country])
            push!(temp_dict[cl[k]][t], y)
        end
    end
    
    return temp_dict
end



function upsample_time_series(; weight::Dict, cluster_dict:: JuMP.Containers.DenseAxisArray, technology:: String, region, tot_sum=false, config::Dict)
    # average time-series if not zeros
    lst_tmp = []       
    if tot_sum
        # count number of zeros
        t = 0
        for r in config["countries"]
            if sum(cluster_dict[r, technology,:,:]) > 0.00000000
                t+=1
            end
        end
        cluster_dict = sum(cluster_dict[r, technology,:,:] for r in config["countries"])/t

        for i in 1:length(keys(weight))
            append!(lst_tmp, repeat(Array(cluster_dict[i, :]), weight[i]))
        end
    else   
        for i in 1:length(keys(weight))
            append!(lst_tmp, repeat(Array(cluster_dict[region, technology,i ,:]), weight[i]))
        end

        if technology ∈ config["Load"]
            tmp_sum = sum(lst_tmp)
            lst_tmp = lst_tmp/tmp_sum
        end

    end
    return lst_tmp
end




"""
calculate_rmse(data_org, cluster_dict)

This function calculates the rmse according to 10.1016/j.energy.2016.06.081.

# Arguments
- `data_org::Dict`: The original time-series data.
- `cluster_dict::JuMP.Containers.DenseAxisArray`: The cluster representatives.
- `cl::Vector{Int64}`: The mapping of the original days to the clusters.

# Returns
- `result::Float64`: The rmse.
"""

function calculate_rmse(; data_org:: Dict, cluster_dict:: JuMP.Containers.DenseAxisArray, cl:: Vector{Int64})
    tmp=0.0
    for t ∈ keys(data_org), r ∈ names(data_org[t])
        z = reshape(data_org[t][:,r], (24,365))
        for i ∈ 1:365
            tmp += sqeuclidean(convert(Vector{Float64}, vcat(z[:,i])),vcat(cluster_dict[r,t,cl[i],:]))
            #tmp += sum((vcat(z[:,i])) - vcat(cluster_dict[r,t,cl[i],:])^2)
        end
    end
    return (1/(length(keys(data_org))*length(axes(cluster_dict)[1])))*sqrt((1/8760)*tmp)
end



"""
calculate_variance(data_org, cluster_dict)

This function calculates the variance according to 10.1016/j.energy.2011.08.021. 

# Arguments
- `cluster_dict::JuMP.Containers.DenseAxisArray`: The cluster representatives.
- `cl::Vector{Int64}`: The mapping of the original days to the clusters.
- `weight::Dict`: Weight of clusters.
- `config::Dict`: config as dictionary.

# Returns
- `result::Float64`: The variance.
"""


function calculate_variance(; 
    cluster_dict:: JuMP.Containers.DenseAxisArray, 
    cl:: Vector{Int64}, 
    weight:: Dict, 
    config:: Dict
    )

    tmp = 0.0
    counter = 0.0
    for t ∈ axes(cluster_dict)[2], r ∈ axes(cluster_dict)[1]
        if sum(vcat(cluster_dict[r,t,:,:])) > 0
            z = upsample_time_series(weight=weight, cluster_dict=cluster_dict, technology=t, region=r, tot_sum=false, config=config)
            tmp +=var(z)
            counter += 1
        end
    end
    return tmp/counter
end


"""
calculate_variance_daily(cluster_dict, weight)

This function calculates the daily variance. 

# Arguments
- `cluster_dict::JuMP.Containers.DenseAxisArray`: The cluster representatives.
- `weight::Dict`: Weight of clusters.

# Returns
- `result::Float64`: The variance.
"""


function calculate_variance_daily(; cluster_dict:: JuMP.Containers.DenseAxisArray, weight:: Dict)

    tmp = 0.0
    for k ∈ axes(cluster_dict)[3]
        tmp += (var(cluster_dict[:,:,k,:])weight[k])
    end
    return tmp/365
end



"""
calculate_correlation(data_org, cluster_dict)

This function calculates the variance according to 10.1016/j.energy.2011.08.021. 

# Arguments
- `cluster_dict::JuMP.Containers.DenseAxisArray`: The cluster representatives.
- `cl::Vector{Int64}`: The mapping of the original days to the clusters.
- `weight::Dict`: Weight of clusters.
- `config::Dict`: config as dictionary.

# Returns
- `result::Float64`: The variance.
"""


function calculate_correlation(; 
    cluster_dict:: JuMP.Containers.DenseAxisArray, 
    data_org:: Dict, 
    cl:: Vector{Int64}, 
    weight:: Dict, 
    config:: Dict)
    
    tmp = 0.0
    counter = 0.0
    for t ∈ axes(cluster_dict)[2], r ∈ axes(cluster_dict)[1]
        if sum(vcat(cluster_dict[r,t,:,:])) > 0
            z = upsample_time_series(weight=weight, cluster_dict=cluster_dict, technology=t, region=r, tot_sum=false, config=config)
            tmp += cor(z, data_org[t][:,r])
            counter += 1
        end
    end
    return tmp/counter
end




function calculate_tail_rmse(;
    cluster_dict::JuMP.Containers.DenseAxisArray,
    data_org::Dict,
    weight::Dict,
    techs::Vector{},
    config::Dict)

    mae_low_list  = Float64[]

    println("Tail MAE diagnostics (normalized by range):")
    for t ∈ techs, r ∈ axes(cluster_dict)[1]
        if sum(vcat(cluster_dict[r,t,:,:])) > 0

            z_hat = upsample_time_series(
                weight=weight,
                cluster_dict=cluster_dict,
                technology=t,
                region=r,
                tot_sum=false,
                config=config
            )
            z_org = data_org[t][:, r]
            n = length(z_org)

            # Sort both series in descending order for load duration curve comparison
            z_org = sort(z_org, rev=true)
            #z_org = z_org[n - 438 + 1:n]
            z_org = z_org[1:438]
            z_hat = sort(z_hat, rev=true)
            #z_hat = z_hat[n - 438 + 1:n]
            z_hat = z_hat[1:438]

            # compute euclidean Distance
            distance = mean(abs.(z_org - z_hat))
            push!(mae_low_list, distance)
        end
    end

    avg_mae_low  = mean(mae_low_list)

    return avg_mae_low
end