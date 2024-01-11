####################################################################################################
# Convenience methods to plot data and results.
####################################################################################################





# Loading data.
####################################################################################################

export load_results

"""
	load_results(path::AbstractString)
Load results and data as saved by [`fit_condition`](@ref). Returns `(result, data,replicates)` as `AdaptiveResult`, `FittingData` and vector of `FittingData` objects.
"""
function load_results(path::AbstractString)
	# Full result data must be present at the path!
	result = deserialize(joinpath(path,"result.jld"))
	grid = deserialize(joinpath(path,"grid.jld"))
	optimizer = deserialize(joinpath(path,"optimizer.jld"))
	objective_value = deserialize(joinpath(path,"objective_value.jld"))
	time = deserialize(joinpath(path,"time.jld"))
	
	if isfile(joinpath(path,"data.jld"))
		data = deserialize(joinpath(path,"data.jld"))
	else
		data = nothing
	end

	if isfile(joinpath(path,"replicates.jld"))
		replicates = deserialize(joinpath(path,"replicates.jld"))
	else
		replicates = nothing
	end
	
	return AdaptiveResult(result,grid,optimizer,objective_value,time), data, replicates
end





















# Default options as functions.
####################################################################################################


export dr_base_plot, density_base_plot,data_options, replicate_options, fit_options, density_options, eu_options, du_options


"""
	function dr_base_plot(keywords...)
Function to return a modified dose-response base plot. Most `Plots.jl` keywords are available. By default, the following keyword arguments are set:

* `xaxis = :log`
* `xlabel = "dilution"`
* `ylabel = "response"`
* `legend = :topleft`
* `density = 300`
"""
function dr_base_plot(;xaxis = :log, xlabel = "dilution", ylabel = "response", legend = :topleft, density = 300, args...)
	return plot(;xaxis = xaxis, xlabel = xlabel, ylabel = ylabel, legend = legend, density = density, args...)
end


"""
	density_base_plot(keywords...) 
Function to return a modified K_τ density base plot. Most `Plots.jl` keywords are available. By default, the following keyword arguments are set:

* `xaxis = :log`
* `xlabel = "K_τ"`
* `ylabel = "density"`
* `legend = :topleft`
* `density = 300`
"""
function density_base_plot(;xaxis = :log, xlabel = "K_τ", ylabel = "density", legend = :topleft, density = 300, args...) 
	return plot(;xaxis = xaxis, xlabel = xlabel, ylabel = ylabel, legend = legend, density = density, args...)
end


"""
	data_options(keywords...) 
Function to return a modified tuple of plotting keyword arguments for the `FittingData` data series. Most `Plots.jl` keywords are available. By default, the following keyword arguments are set:

* `seriestype = :scatter`
* `color = 1`
* `label = "mean values"`
* `yerrors = nothing`
"""
function data_options(;seriestype = :scatter,color = 1, label = "mean values", yerrors = nothing, args...)
	return (seriestype = seriestype,color = color, label = label, yerrors = yerrors, args...)
end


"""
	replicate_options(keywords...) 
Function to return a modified tuple of plotting keyword arguments for the replicate data series. Most `Plots.jl` keywords are available. By default, the following keyword arguments are set:

* `seriestype = :scatter`
* `color = :black`
* `opacity = 0.2`
* `label = "replicates"`
* `yerrors = nothing`
"""
function replicate_options(;seriestype = :scatter ,color = :black, opacity = 0.2, label = "replicates", yerror = nothing, args...)
	return (seriestype = seriestype ,color = color, opacity = opacity, label = label, yerror = yerror, args...)
end


"""
	fit_options(keywords...) 
Function to return a modified tuple of plotting keyword arguments for the `DoseResponseResult` data series. Most `Plots.jl` keywords are available. By default, the following keyword arguments are set:

* `seriestype = :path`
* `color = 2`
* `label = "fit result"`
"""
function fit_options(;seriestype = :path, color = 2, label = "fit result", args...)
	return (seriestype = seriestype, color = color, label = label, args...)
end


"""
	density_options(keywords...) 
Function to return a modified tuple of plotting keyword arguments for the `DoseResponseResult` data series. Most `Plots.jl` keywords are available. By default, the following keyword arguments are set:

* `seriestype = :path`
* `color = 2`
* `fillalpha = 0.5`
* `label = "fitted density"`
"""
function density_options(;color = 2, fill = 0, fillalpha = 0.5, label = "fitted density", args...)
	return (color = color, fill = fill, fillalpha = fillalpha, label = label, args...)
end


"""
	eu_options(n::Integer, bins = nothing; keywords...) 
Function to return a modified tuple of plotting keyword arguments for the `EpitopeUncertainty` data series. 

`n` must be the number of levels (or larger) for the `EpitopeUncertainty` object that is plotted. If not `bins = nothing` the passed bins are marked in the plot with dashed lines (the color can be changed with the keyword `bin_color`).

Most `Plots.jl` keywords are available. By default, the following keyword arguments are set:

* `colors = colormap("RdBu",n)[end:-1:1]`
* `opacities = [1]`
"""
function eu_options(n::Integer, bins = nothing; colors =colormap("RdBu", n)[end:-1:1], opacities = [1], args...)
	return (colors = colors, opacities = opacities, bins = bins, args...)
end


"""
	du_options(n::Integer; keywords...) 
Function to return a modified tuple of plotting keyword arguments for the `EpitopeUncertainty` data series. 

`n` must be the number of levels (or larger) for the `EpitopeUncertainty` object that is plotted. 

Most `Plots.jl` keywords are available. By default, the following keyword arguments are set:

* `colors = colormap("RdBu",n)[end:-1:1]`
* `opacities = [1]`
"""
function du_options(n::Integer; colors =colormap("RdBu", n)[end:-1:1], opacities = [1], args...)
	return (colors = colors, opacities = opacities, args...)
end

























# Plotting
####################################################################################################

export bin_analysis_plot, peak_analysis_plot, uncertainty_plot


"""
	bin_analysis_plot(results::Union{AdaptiveResult,Nothing},data = nothing,replicates = nothing; keywords...)
Create and return basic `DoseResponseResult` and K_τ-density `gird` plots with the option to evaluate the epitope number in specified bins. 

If `results` is an AdaptiveResult, both the gird is plotted into the density plot and the corresponding resulting curve is plotted into the dose-response plot.

If `data` is a `FittingData` object the data points are plotted in the dose-response plot. Similarly, if `replicates` is an array of `FittingData` objects, the data points are plotted as replicates in the dose-response plot.

**Keywords**

* `dr_plot = dr_base_plot()`: The base plot onto which the `AdaptiveResult.result` and the `FittingData` objects (`data` and `replicates`) are plotted. It can be any `Plots.jl` plot (e.g. another dose-response plot). The `AdaptiveResult.result` and `FittingData` objects are then plotted on top.
* `density_plot = density_base_plot()`: The base plot onto which the `AdaptiveResult.grid` is plotted. It can be any `Plots.jl` plot (e.g. another `K_τ` density plot). The `AdaptiveResult.grid` is then plotted on top.
* `fit_arguments = fit_options()`: Keyword arguments for the `AdaptiveResult.result` data-series.
* `data_arguments = data_options()`: Keyword arguments for the `FittingData` data-series.
* `replicate_arguments = replicate_options()`: Keyword arguments for the replicate data-series.
* `density_arguments = density_options()`: Keyword arguments for the `AdaptiveResult.grid` data-series.
* `annotation_arguments = NamedTuple()`: Keyword arguments for bin-annotations in the density plot.

**Bin annotations**

The annotation bins are bins for the K_τ density. The total number of epitopes (in the unit of the density) within the respective bin is calculated and added as annotation to the plot.

The following keywords can be used for the `annotation_arguments`:

* `annotation_bins = []`: The bins.
* `annotation_size = 10`: Size of the annotation font.
* `annotation_offset = 0.05`: Relative vertical offset for the annotation from the top of the plot (can also be a vector of offsets that matches the number of bins).
* `annotation_color = :black`: Color of the annotations.
* `hover_points = false`: If true, adds scatter-points with tool-tips for the `Plotly.jl` backend.

"""
function bin_analysis_plot(results::Union{AdaptiveResult,Nothing},data::Union{FittingData,Nothing} = nothing,replicates::Union{Nothing,AbstractArray{FittingData}} = nothing;
	dr_plot = dr_base_plot(),
	density_plot = density_base_plot(),
	fit_arguments = fit_options(),
	data_arguments = data_options(),
	replicate_arguments = replicate_options(),
	density_arguments = density_options(),
	annotation_arguments = NamedTuple(),
	)
	
	dr_plot = plot(dr_plot)

	if !isnothing(replicates)
		concentrations = deepcopy(replicates[1].independent)
		responses = deepcopy(replicates[1].dependent)
		for i in 2:length(replicates)
			push!(concentrations, replicates[i].independent...)
			push!(responses, replicates[i].dependent...)
		end
		plot!(FittingData(concentrations,responses); replicate_arguments...)
	end
	

	if !isnothing(results)
		plot!(results.result; fit_arguments...)
	end

	if !isnothing(data)
		plot!(data; data_arguments...)
	end

	density_plot = plot(density_plot)

	if !isnothing(results)
		plot!(DensityPlot(results.grid); density_arguments...)
	end

	if !isempty(annotation_arguments)
		plot!(DensityAnnotations(results.grid); annotation_arguments...)
	end

	return dr_plot, density_plot
end




"""
	peak_analysis_plot(results::AdaptiveResult,data = nothing; keywords...)
Create and return plots to analyze the impact of peaks in the K_τ density on the corresponding dose-response curve.

Returns `(individual_dr_plot, cumulative_dr_plot, density_plot)`, where

* `density_plot` Contains a plot of the K_τ-density with different colors for the different peaks (specified by the bins).
* `individual_dr_plot`: Contains the individual dose-response curves (color matched) that originate from the different peaks alone.
* `cumulative_dr_plot`: Contains the cumulative dose-response curves, i.e. dose-responses include the response increases caused by earlier peaks. Again the curves are color matched with the peaks.

If `data` is a `FittingData` object the data points are plotted in the `cumulative_dr_plot`.

**Keywords**

* `individual_dr_plot = dr_base_plot()`: The base plot onto which the individual dose-response curves are plotted. It can be any `Plots.jl` plot (e.g. another dose-response plot). The individual dose-response curves are then plotted on top.
* `cumulative_dr_plot = dr_base_plot()`: The base plot onto which the cumulative dose-response curves are plotted. It can be any `Plots.jl` plot (e.g. another dose-response plot). The individual dose-response curves are then plotted on top.
* `density_plot = density_base_plot()`: The base plot onto which the K_τ-peaks are plotted. It can be any `Plots.jl` plot (e.g. another dose-response plot). The K_τ-peaks are then plotted on top.
* `bins = peak_detection(results.grid, fill = false)[2]`: The bins (domain ranges) that define the K_τ-peaks.
* `colors = collect(1:length(bins))`: The colors for the different bins.
* `join_bins = true`: If true, extends the bins if needed, s.t. there remains no gap between the bins.
"""
function peak_analysis_plot(results::AdaptiveResult,data::Union{FittingData, Nothing} = nothing;
	individual_dr_plot = dr_base_plot(),
	cumulative_dr_plot = dr_base_plot(),
	density_plot = density_base_plot(),
	bins = peak_detection(results.grid, fill = false)[2],
	colors = collect(1:length(bins)),
	join_bins = true
	)
	
	
	if !isnothing(data)
		cumulative_dr_plot = plot(cumulative_dr_plot)
		scatter!(data, color = :black, label = "data")
	end


	sorted_bins = sort(bins)

	if join_bins && length(sorted_bins) > 1
		sorted_bins[1] = [-Inf,(sorted_bins[1][2]+sorted_bins[2][1])/2]
		for i in 2:length(sorted_bins)-1
			sorted_bins[i] = [(sorted_bins[i-1][2]+sorted_bins[i][1])/2,(sorted_bins[i][2]+sorted_bins[i+1][1])/2 ]
		end
		sorted_bins[end] = [(sorted_bins[end-1][2]+sorted_bins[end][1])/2,Inf]
	elseif join_bins
		sorted_bins[1] = [-Inf,Inf]
	end

	for i in length(sorted_bins):-1:1
		single_peak_grid = restrict_domain!(deepcopy(results.grid),lower = sorted_bins[i][1], upper = sorted_bins[i][2], weight_distribution = :log)

		density_plot = plot(density_plot)
		plot!(DensityPlot(single_peak_grid), color = colors[i], fill = 0, fillalpha = 0.5, label = "peak $i")

		if length(results.optimizer) > length(results.grid)
			dose_response = DoseResponseResult(single_peak_grid,results.result.concentrations, offset = results.optimizer[end])
		else
			dose_response = DoseResponseResult(single_peak_grid,results.result.concentrations)
		end
		individual_dr_plot = plot(individual_dr_plot)
		plot!(dose_response, color = colors[i], label = "peak $i")

		cumulative_peak_grid = restrict_domain!(deepcopy(results.grid), upper = sorted_bins[i][2], weight_distribution = :log)
		if length(results.optimizer) > length(results.grid)
			dose_response = DoseResponseResult(cumulative_peak_grid,results.result.concentrations, offset = results.optimizer[end])
		else
			dose_response = DoseResponseResult(cumulative_peak_grid,results.result.concentrations)
		end
		cumulative_dr_plot = plot(cumulative_dr_plot)
		plot!(dose_response, color = colors[i], label = "peak $i")
	end
	return individual_dr_plot, cumulative_dr_plot, density_plot
end












"""
	uncertainty_plot(e_uncertainty::EpitopeUncertainty,d_uncertainty::DoseResponseUncertainty,grid::OneDimGrid; keywords...)
Create and return uncertainty visualizations `(dr_uncertainty_plot, density_uncertainty_plot)`.

The estimated bounds of the `DoseResponseUncertainty` and `EpitopeUncertainty` are plotted as color-matched ribbons.

**Keywords**

* `dr_plot = dr_base_plot()`: The base plot onto which the `AdaptiveResult.result` and the `FittingData` objects (`data` and `replicates`) are plotted. It can be any `Plots.jl` plot (e.g. another dose-response plot). The `AdaptiveResult.result` and `FittingData` objects are then plotted on top.
* `density_plot = density_base_plot()`: The base plot onto which the `AdaptiveResult.grid` is plotted. It can be any `Plots.jl` plot (e.g. another `K_τ` density plot). The `AdaptiveResult.grid` is then plotted on top.
* `eu_arguments = eu_options(length(e_uncertainty.levels))`: Keyword arguments for the `EpitopeUncertainty` data series.
* `du_arguments = du_options(length(d_uncertainty.levels))`: Keyword arguments for the `DoseResponseUncertainty` data series.

"""
function uncertainty_plot(e_uncertainty::EpitopeUncertainty,d_uncertainty::DoseResponseUncertainty,grid::AdaptiveDensityApproximation.OneDimGrid;
	dr_plot = dr_base_plot(),
	density_plot = density_base_plot(),
	eu_arguments = eu_options(length(e_uncertainty.levels)),
	du_arguments = du_options(length(d_uncertainty.levels))
	)
	
	dr_plot = plot(dr_plot)
	plot!(d_uncertainty; du_arguments...)

	density_plot = plot(density_plot)
	plot!(grid,e_uncertainty; eu_arguments...)

	return dr_plot, density_plot
end