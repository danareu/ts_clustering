##### LOAD THE DATA


function read_data(;path::String, config::Dict)
# returns dictionary with technology time series
    CountryData = Dict()
    hourly_data = XLSX.readxlsx(path)

    for cde ∈ config["Country_Data_Entries"]
        CountryData[cde] = DataFrame(XLSX.gettable(hourly_data[cde]))
        select!(CountryData[cde], (config["countries"]))
    end

    # reformat the load as share of the total load
    for t ∈ intersect(config["Load"], config["Country_Data_Entries"])
        for c in config["countries"]
            CountryData[t][:,c]  = CountryData[t][:,c] / 8760
        end
    end
    return CountryData
end


function normalize_data(; config::Dict, CountryData::Dict)

    for cde ∈ config["Country_Data_Entries"]
        for col in names(CountryData[cde])
            try
                CountryData[cde][!, col] .= zscore_column!(CountryData[cde][!, col])
            catch
                println(col, cde)
            end
        end
    end
    return CountryData
end


function zscore_column!(col)
    # normalize the data
    mean_col = mean(col)
    std_col = std(col)
    return (col .- mean_col) / std_col
end




function read_yaml_file(; file)
    return YAML.load_file(file)
end


"""
    create_clustering_matrix(; technology::Vector{String}, CountryData::Dict)

Create a clustering matrix for k-means algorithms.

# Arguments
- `technology::Vector{String}`: A vector of strings representing different technologies.
- `CountryData::Dict`: A dictionary where the keys are technology names and the values are data matrices for different countries.

# Returns
- `x::Matrix{Float64}`: A matrix, where each technology-country combination's data is reshaped and concatenated, interleaved with zero matrices.

"""
function create_clustering_matrix(;technology::Vector{String}, CountryData::Dict)
    i = 1
    x = 0
    for cde in technology
        for n in names(CountryData[cde])
            if i == 1
                x = reshape(CountryData[cde][:,n], (24,365))
            else
                x = vcat(x, reshape(CountryData[cde][:,n], (24,365)))
            end
            # create zero ts for dynamic warping window
            x = vcat(x, zeros(Int, 24, 365))
            i += 1
        end
    end
    try
        x = convert(Matrix{Float64}, x)
        x[isnan.(x)] .= 0.0
    catch
        println(technology)
    end
    return x
end


