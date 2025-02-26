## plot the results


function plot_cluster_centers(;K::Integer, config::Dict, FullData, CountryData, country, a, weights, hoffmann)

    p = make_subplots(rows=length(collect(keys(CountryData))), 
    cols=K, 
    vertical_spacing=0.04, 

    column_titles=["Day $(i)" for i in 1:K], row_titles=[split(i, "TS_")[end] for i in collect(keys(CountryData))])

    # real da
    for (i,t) in enumerate(collect(keys(CountryData)))
        for (k,m) in enumerate(0:24:8735)
            if t ∈ ["TS_LOAD", "TS_HEAT_LOW", "TS_HEAT_HIGH", "TS_MOBILITY_PSNG"]
                y=CountryData[t][m+1:m+24,country]/sum(CountryData[t][m+1:m+24,country])
            else
                y=CountryData[t][m+1:m+24,country]
            end

            add_trace!(p, 
            scatter(y=y, 
            marker=attr(
                color="grey"
            ),
            showlegend=false,
            ), 
            row=i, 
            col=a[k+1])
        end
    end


    # clustered data
    for d in 1:K
        for (i,t) in enumerate(collect(keys(CountryData)))
            show_legend = (i == 1 ) && (d == 1)
            if t ∈ ["TS_LOAD", "TS_HEAT_LOW", "TS_HEAT_HIGH", "TS_MOBILITY_PSNG"]
                y=FullData[country,t,d,:]/sum(FullData[country,t,d,:])
                c = hoffmann[country,t,d,:]/sum(hoffmann[country,t,d,:]) 
            else
                y=FullData[country,t,d,:]
                c=hoffmann[country,t,d,:]
            end

            add_trace!(p, 
            scatter(y=y, 
            marker=attr(
                color="#bc203e"
            ),
            name="medoid",
            showlegend=show_legend,
            ), 
            row=i, 
            col=d)

            add_trace!(p, 
            scatter(y=c, 
            marker=attr(
                color="blue"
            ),
            name="hoffmann",
            showlegend=show_legend,
            ), 
            row=i, 
            col=d)
        end
    end


    relayout!(p)
    p.plot.layout.height = 200*length(collect(keys(CountryData)))  
    p.plot.layout.width = 200 * K
    return p
end


function plot_variance_cluster(;K::Integer, var_temp::Dict,)

    p = make_subplots(rows=K, 
    cols=1, 
    vertical_spacing=0.02, 
    #column_titles=["Day$(i)_w_$(round(weights[i]/365,digits=2))" for i in 1:K], row_titles=[split(i, "TS_")[end] for i in collect(keys(CountryData))]
    )

    color_dict = Dict(t => ColorSchemes.tab20.colors[i] for (i,t) in enumerate(keys(var_temp[K])))


    for d in keys(var_temp)
        for t in keys(var_temp[d])
            add_trace!(p, 
            box(y=var_temp[d][t],
            name=t, 
            showlegend=false,
            jitter=0.5,
            boxpoints="all",
            marker_color=color_dict[t],
            fillcolor=color_dict[t],
            whiskerwidth=0.2,
            ), 
            row=d, 
            col=1)
        end
    end

    p.plot.layout.height = 200*K  # Set the width of the plot (in pixels)
    return p
end






function plot_heatmaps(; full_data, clustered_data)
    p = make_subplots(rows=1, 
    cols=2, 
    vertical_spacing=0.02, 
    )

    for (i, df) in enumerate([full_data, clustered_data])
        add_trace!(p,
        heatmap(z=Matrix(reshape(df, 24, 365)'),
        y=1:365,
        x=1:24),
        #showscale=true if i==1 else false
        row=1, col=i)
    end

    relayout!(p)
    return p

end




function map_to_org(; list_df, list_window, write_html=false, name="Chronology")
    traces = GenericTrace[]

    #color_dict = create_palette(list_data=list_window)
    color_dict = Dict(0 => "rgb(181,48,57)", 5 => "rgb(247,144,111)", 10 => "rgb(1,151,239)", 1 => "rgb(212,144,202)", 3 => "rgb(1,151,239)")

    for (df, w) in zip(list_df, list_window)
        n = w == 999 ? "Original" : (w == 0 ? "Euclidean" : "DTW ($w)")
        push!(traces, scatter(
            x=1:365, 
            y=df, 
            mode="markers", 
            name=n,
            marker_size=10,
            marker_color=color_dict[w],
        ))
    end

    scatter_plot = plot(traces)

    relayout!(scatter_plot,
    xaxis_title="Number of Days",
    yaxis_title="Number of Representative Days")


    if write_html
        format_layout(p=scatter_plot, max_value=0)
        savefig(scatter_plot, name)
    end
    return scatter_plot
end




function plot_dtw_heatmaps(; matrix, x,y)
    p = make_subplots(rows=1, 
    cols=1, 
    vertical_spacing=0.02
    )

    add_trace!(p,
    heatmap(z=matrix,
    y=1:size(matrix)[1],
    x=1:size(matrix)[2]),
    row=1, col=1)


    add_trace!(p,
    scatter(
    y=y,
    marker_color="black",
    x=x),
    row=1, col=1) 
    #format_layout(p=p, max_value=0)
    relayout!(p)
    return p

end


function plot_basic_heatmap(; matrix)

    return plot(heatmap(z=matrix))

end


function plot_time_series_dtw(; ts1, ts2)
    p = make_subplots(rows=1, 
    cols=1, 
    vertical_spacing=0.02
    )

    add_trace!(p,
    scatter(
    marker_color="red",
    y=ts2),
    row=1, col=1) 


    add_trace!(p,
    scatter(
    marker_color="blue",
    y=ts1),
    row=1, col=1) 

    relayout!(p)
    return p
end





"""
plot_duration_curve(list_data, technology)

This function calculates the rmse according to 10.1016/j.energy.2016.06.081.

# Arguments
- `list_data::JuMP.Containers.DenseAxisArray`: Array with 1:8760 scenario data.
- `technology::Vector`: Technologies to be plotted.
- `write_html::Boolean`: Write html file 

# Returns
- `p::figure`: Subplots for each technology.
"""

function plot_duration_curve(; 
    list_data::JuMP.Containers.DenseAxisArray, 
    technology::Vector, 
    write_html::Bool=true,
    kwargs...)

    p = make_subplots(rows=length(axes(list_data)[2]), 
    cols=length(technology), 
    vertical_spacing=0.02,
    #specs=[Spec() Spec() Spec() Spec()],
    subplot_titles=["TS_LOAD" "TS_WIND_ONSHORE_AVG" "TS_WIND_OFFSHORE" "TS_PV_AVG"]
    
    )

    palette = ["rgb(1,151,239)", "rgb(247,144,111)", "rgb(82,169,105)", "rgb(212,144,202)", "rgb(231,202,145)", "rgb(181,48,57)", "rgb(252,184,32)", "grey"]
    color_dict = Dict(key => palette[i] for (i, key) in enumerate(axes(list_data)[1]))

    #color_dict = Dict("Original_Data"=> "grey", "Euclidean" => "rgb(181,48,57)", "DTW (5)" => "rgb(247,144,111)", "DTW (10)" => "rgb(1,151,239)", "DTW (1)" => "rgb(212,144,202)")

    for (g,r) ∈ enumerate(axes(list_data)[2])
        for (j,s) ∈ enumerate(axes(list_data)[1])
            for (i,t) ∈ enumerate(technology)
                add_trace!(p,
                scatter(x=1:8760,
                y=list_data[s,r,t,:],
                line=attr(
                        color=color_dict[s]),
                name="$(r)_$(s)",
                showlegend=(i == 1),
                legendgroup=s),
                row=g, 
                col=i) 
            end
        end
    end

    # plot genesys timeseries 
    # probably wont workk because of region row is missing
    if :genesys in keys(kwargs)
        ts_gerb = Dict(kwargs)[:genesys]
        for (i,t) ∈ enumerate(technology)
            z = sort(Vector(ts_gerb["DE", t, 1:8760]), rev=true)
            add_trace!(p,
            scatter(x=1:8760,
            y=z,
            line=attr(color="red"),
            name="old_method",
            legendgroup="GENESYS"),
            row=1, 
            col=i) 
        end
    end
        
    if write_html
        format_layout(p=p, max_value=0)
        savefig(p, "LoadDuration_$(Dates.now()).html")
    end
    p.plot.layout.height = 200*length(axes(list_data)[2])  # Set the width of the plot (in pixels)

    return p
end




"""
plot_rmse_clusters(df, K, W)

This function calculates the rmse according to 10.1016/j.energy.2016.06.081.

# Arguments
- `df::DataFrame`: Dataframe with Cluster, Weight and Values for RMSE.
- `name::String`: name of html file.

# Returns
- `scatter_plot::figure`: scatterplot.
"""

function plot_rmse_clusters(; df:: DataFrame, write_html=false, name:: String)

    traces = GenericTrace[]

    if "Window" ∈ names(df)

        W = unique(df.Window)

        palette = ["rgb(1,151,239)", "rgb(247,144,111)", "rgb(82,169,105)", "rgb(212,144,202)", "rgb(231,202,145)", "rgb(252,184,32)"]

        if length(unique(df.Window)) <= length(palette)
            color_dict = Dict((key == 0 ? key => "rgb(181,48,57)" : key => palette[i]) for (i, key) in enumerate(W))
        else
            color_dict = Dict((key == 0 ? key => "rgb(181,48,57)" : key => ColorSchemes.tab20.colors[i]) for (i, key) in enumerate(W))
        end

        for w in W
            n = w == 999 ? "Original" : (w == 0 ? "Euclidean" : "DTW ($w)")
            l = w == 4 ? write_html : 2
            dash_value = (w == 0) ? "dash" : "solid"
            push!(traces, scatter(
                x=filter(row -> row.Window == w, df)[:, :Cluster], 
                y=filter(row -> row.Window == w, df)[:, :Value], 
                mode="lines", 
                name=w,
                line=attr(width=l,
                    color=color_dict[w],
                    dash=dash_value),
            ))
        end
    else
        push!(traces, scatter(
                x=df[:, :TimeSlice]/24, 
                y=df[:, :Value], 
                mode="lines", 
                line=attr(
                    color="rgb(1,151,239)")))
    end


    scatter_plot = plot(traces)
    relayout!(scatter_plot,
    xaxis_title="Number of Clusters",
    yaxis_title="RMSE")

    # Show the scatter plot
    if write_html
        #format_layout(p=scatter_plot, max_value=maximum(df[!, :Value]))
        savefig(scatter_plot, name)
    end
    return scatter_plot
end


function format_layout(; p, max_value)
    return relayout!(
        p,
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin=attr(l=10, r=10, t=10, b=10),  # Adjust the margins for spacing
        font_family="Calibri",
        font=attr(size=25, color="white"),
        legend_traceorder="normal",
        font_color="white",
        
        xaxis=attr(
            showline=true,
            mirror="allticks",  # Mirrors the axis lines across all subplots
            showgrid=false,
            tickfont=attr(color="white"),
            linecolor="black",  # Border color
            linewidth=2,  # Border thickness
        ),
        
        yaxis=attr(
            showline=true,
            mirror="allticks",  # Mirrors the axis lines across all subplots
            showgrid=false,
            tickfont=attr(color="white"),
            linecolor="black",  # Border color
            linewidth=2,  # Border thickness
        ),
        
        # Apply to all subplots
        xaxis2=attr(
            showline=true,
            mirror="allticks",  
            showgrid=false,
            tickfont=attr(color="white"),
            linecolor="black",  
            linewidth=2,
        ),
        
        yaxis2=attr(
            showline=true,
            mirror="allticks",  
            showgrid=false,
            tickfont=attr(color="white"),
            linecolor="black",  
            linewidth=2,
        ),
        
        # Add additional x/y axes as needed for more subplots
    )
    
    
end



"""
plot_box_technology(df)

This function plots the box for clustered technology data.

# Arguments
- `df::DataFrame`: Dataframe with Cluster, Weight and Values for RMSE.

# Returns
- `box_plot::figure`: boxplot.
"""

function plot_box_technology(; df:: DataFrame)

    traces = GenericTrace[]

    W = unique(df.Window)
    T = unique(df.Technology)

    palette = ["rgb(1,151,239)", "rgb(247,144,111)", "rgb(82,169,105)", "rgb(212,144,202)", "rgb(231,202,145)", "rgb(181,48,57)", "rgb(252,184,32)"]
    color_dict = Dict(key => palette[i] for (i, key) in enumerate(W))

    for t ∈ T, w ∈ W
        push!(traces, box(
            y=filter(row -> (row.Window == w) & (row.Technology == t), df)[:, :Value][1],
            name="$(t)_$w",
            jitter=0.5,
            whiskerwidth=0.2,
            fillcolor=color_dict[w],
            marker_size=2,
            line_width=1
        ))
    end


    p = plot(traces)
    format_layout(p=p)

    # Show the scatter plot
    savefig(p, "results/box_$(Dates.now()).html")
    return p
end



"""
create_palette(list_data)

This function plots the box for clustered technology data.

# Arguments
- `list_data`: Iterable item to assign colors to.

# Returns
- `color_dict::Dict`: key,values with colors.
"""

function create_palette(; list_data)
    palette = ["rgb(1,151,239)", "rgb(247,144,111)", "rgb(82,169,105)", "rgb(212,144,202)", "rgb(231,202,145)", "rgb(181,48,57)", "rgb(252,184,32)", "grey"]
    color_dict = Dict((key == 0 ? key => "rgb(181,48,57)" : key => palette[i]) for (i, key) in enumerate(list_data))
    return color_dict
end




function plot_cluster_centers_converence(;K::Integer, config::Dict, FullData, CountryData, country, a, weights)

    p = make_subplots(rows=1, 
    cols=4, 
    vertical_spacing=0.02, 
    column_titles=["Wind Offshore", "Wind Onshore", "Load", "PV"]
    )

    # real da
    for (i,t) in enumerate(collect(keys(CountryData)))
        for (k,m) in enumerate(0:24:8735)
            if t ∈ ["TS_LOAD", "TS_HEAT_LOW", "TS_HEAT_HIGH"]
                y=CountryData[t][m+1:m+24,country]/sum(CountryData[t][m+1:m+24,country])
            else
                y=CountryData[t][m+1:m+24,country]
            end
            if a[k+1] == 3
                add_trace!(p, 
                scatter(y=y, 
                marker=attr(
                    color="grey"
                ),
                showlegend=false,
                ), 
                row=1, 
                col=i)
            end
        end
    end


    # clustered data
    for d in [3]
        for (i,t) in enumerate(collect(keys(CountryData)))
            if t ∈ ["TS_LOAD", "TS_HEAT_LOW", "TS_HEAT_HIGH"]
                y=FullData[country,t,d,:]/sum(FullData[country,t,d,:])
            else
                y=FullData[country,t,d,:]
            end

            add_trace!(p, 
            scatter(y=y, 
            marker=attr(
                color="#bc203e"
            ),
            name="medoid",
            showlegend=false,
            ), 
            row=1, 
            col=i)


        end
    end

    update_xaxes!(p, range=[0,24], tickvals=[0,12,24])
    update_yaxes!(p, range=[0,1], tickvals=[0,1])


    relayout!(p)
    format_layout(p=p, max_value=0)
    return p
end