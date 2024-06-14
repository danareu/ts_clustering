struct ClusterResults
    R:: Clustering.KmeansResult{Matrix{Float64}, Float64, Int64}
    a::Vector{Int64}
    c::Vector{Int64}
    M::Matrix{Float64}
end