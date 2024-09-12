####################################################################################################
# Convenience methods to plot data and results.
####################################################################################################





# Loading data.
####################################################################################################

export load_results

"""
	load_results(path::AbstractString)
Load results and data as saved by [`fit_condition`](@ref). Returns `(result, data,replicates)` as [`AdaptiveResult`](@ref), [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) and vector of [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects.
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
Function to return a modified dose-response base plot. Most [`Plots.jl`](https://docs.juliaplots.org/stable/) keywords are available. By default, the following keyword arguments are set:

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
Function to return a modified K_τ density base plot. Most [`Plots.jl`](https://docs.juliaplots.org/stable/) keywords are available. By default, the following keyword arguments are set:

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
Function to return a modified tuple of plotting keyword arguments for the [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object. Most [`Plots.jl`](https://docs.juliaplots.org/stable/) keywords are available. By default, the following keyword arguments are set:

* `seriestype = :scatter`
* `color = 1`
* `label = "mean values"`
* `yerrors = nothing`

In addition, the following keyword ia available:

* `filter_zeros = [true,false]`: Select to omit data points from the plot when the x-value is zero (if `true` for first element of `filter_zeros`) or the y-value is zero (if `true` for second element of `filter_zeros`). See [Measurement data - Plotting](@ref measurement_data_plotting) for further information.
"""
function data_options(;seriestype = :scatter,color = 1, label = "mean values", yerrors = nothing, args...)
	return (seriestype = seriestype,color = color, label = label, yerrors = yerrors, args...)
end


"""
	replicate_options(keywords...) 
Function to return a modified tuple of plotting keyword arguments for the replicate data. Most [`Plots.jl`](https://docs.juliaplots.org/stable/) keywords are available. By default, the following keyword arguments are set:

* `seriestype = :scatter`
* `color = :black`
* `opacity = 0.2`
* `label = "replicates"`
* `yerrors = nothing`

In addition, the following keyword ia available:

* `filter_zeros = [true,false]`: Select to omit data points from the plot when the x-value is zero (if `true` for first element of `filter_zeros`) or the y-value is zero (if `true` for second element of `filter_zeros`). See [Measurement data - Plotting](@ref measurement_data_plotting) for further information.
"""
function replicate_options(;seriestype = :scatter ,color = :black, opacity = 0.2, label = "replicates", yerror = nothing, args...)
	return (seriestype = seriestype ,color = color, opacity = opacity, label = label, yerror = yerror, args...)
end


"""
	fit_options(keywords...) 
Function to return a modified tuple of plotting keyword arguments for the [`DoseResponseResult`](@ref) response curve. Most [`Plots.jl`](https://docs.juliaplots.org/stable/) keywords are available. By default, the following keyword arguments are set:

* `seriestype = :path`
* `color = 2`
* `label = "fit result"`

In addition, the following keyword ia available:

* `filter_zeros = [true,false]`: Select to omit data points from the plot when the x-value is zero (if `true` for first element of `filter_zeros`) or the y-value is zero (if `true` for second element of `filter_zeros`). See [Measurement data - Plotting](@ref measurement_data_plotting) for further information.
"""
function fit_options(;seriestype = :path, color = 2, label = "fit result", args...)
	return (seriestype = seriestype, color = color, label = label, args...)
end


"""
	density_options(keywords...) 
Function to return a modified tuple of plotting keyword arguments for the [`DoseResponseResult`](@ref) density. Most [`Plots.jl`](https://docs.juliaplots.org/stable/) keywords are available. By default, the following keyword arguments are set:

* `seriestype = :path`
* `color = 2`
* `fillalpha = 0.5
* `label = "fitted density"`


In addition, the following keyword ia available:

* `volume_normalization = :log`: Normalizes the grid weights for plotting. If `:none` the weights are not normalized, if `:linear` the weights are divided by their corresponding interval length and if `:log` the weights are divided by the interval length as it appears in a logarithmic plot. See [Background: log-volume normalization](@ref log_volume_normalization) for further information.

"""
function density_options(;color = 2, fill = 0, fillalpha = 0.5, label = "fitted density", args...)
	return (color = color, fill = fill, fillalpha = fillalpha, label = label, args...)
end


"""
	eu_options(n::Integer, bins = nothing; keywords...) 
Function to return a modified tuple of plotting keyword arguments for the [`EpitopeUncertainty`](@ref) data series. 

`n` must be equal to the number of levels (or larger) of the `EpitopeUncertainty` object that is plotted. If not `bins = nothing` the passed bins are marked in the plot with dashed lines (the color can be changed with the keyword `bin_color`). The bins must specify the indices of the gird intervals, e.g. `[[1,2,3],[4,5], [9,10,11]]`. The function [`select_indices`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.select_indices) can be used to obtain grid indices from grid domain ranges.

Most [`Plots.jl`](https://docs.juliaplots.org/stable/) keywords are available. In addition, the following keyword arguments are available:

* `colors = colormap("RdBu",n)[end:-1:1]`: Array of colors (that [`Plots.jl`](https://docs.juliaplots.org/stable/) accepts for the `color` keyword) corresponding to the different uncertainty levels. If the array contains less colors than uncertainty levels, the last color is repeated for the remaining levels.
* `opacities = [1]`: Array of opacities (number between `0` and `1`) that correspond to the different uncertainty levels. Again, the last opacity is repeated if there are more uncertainty levels than opacities.
* `reverse = false`: If `true` the plotting order of the uncertainty levels is reversed. Since the uncertainty ranges are plotted on top of each other, this can become necessary when the [`EpitopeUncertainty`](@ref) constructor for samples is used, where larger levels correspond to larger uncertainty (as opposed to the bin-wise shifting constructor). 
* `volume_normalization = :log`: Normalizes the grid weights for plotting. If `:none` the weights are not normalized, if `:linear` the weights are divided by their corresponding interval length and if `:log` the weights are divided by the interval length as it appears in a logarithmic plot. See [Background: log-volume normalization](@ref log_volume_normalization) for further information.
* `hide_labels = true`: If `true` the labels are omitted. Can become necessary when a large number of uncertainty levels is used.
* `bins = nothing`: Specifies the positions for the bin-markers (dashed lines). The bins must specify the interval indices, e.g. `[[1,2,3], [5,6]]`. Ideally, the bins used for the [`EpitopeUncertainty`](@ref) construction should be used. If `bins = nothing`, bin markers are omitted.
* `bin_color = :gray`: Color of the bin markers.
"""
function eu_options(n::Integer, bins = nothing; colors =colormap("RdBu", n)[end:-1:1], opacities = [1], args...)
	return (colors = colors, opacities = opacities, bins = bins, args...)
end


"""
	du_options(n::Integer; keywords...) 
Function to return a modified tuple of plotting keyword arguments for the [`DoseResponseUncertainty`](@ref) data series. 

`n` must be equal to the number of levels (or larger) of the [`DoseResponseUncertainty`](@ref) object that is plotted. 

Most [`Plots.jl`](https://docs.juliaplots.org/stable/) keywords are available. In addition, the following keyword arguments are available:

* `colors = colormap("RdBu",n)[end:-1:1]`: Array of colors (that [`Plots.jl`](https://docs.juliaplots.org/stable/) accepts for the `color` keyword) corresponding to the different uncertainty levels. If the array contains less colors than uncertainty levels, the last color is repeated for the remaining levels.
* `opacities = [1]`: Array of opacities (number between `0` and `1`) that correspond to the different uncertainty levels. Again, the last opacity is repeated if there are more uncertainty levels than opacities.
* `reverse = false`: If `true` the plotting order of the uncertainty levels is reversed. Since the uncertainty ranges are plotted on top of each other, this can become necessary when the [`EpitopeUncertainty`](@ref) constructor for samples is used, where larger levels correspond to larger uncertainty (as opposed to the bin-wise shifting constructor). 
* `hide_labels = true`: If `true` the labels are omitted. Can become necessary when a large number of uncertainty levels is used.
* `filter_zeros = [true,false]`: Select to omit data points from the plot when the x-value is zero (if `true` for first element of `filter_zeros`) or the y-value is zero (if `true` for second element of `filter_zeros`). See [Measurement data - Plotting](@ref measurement_data_plotting)) for further information.
"""
function du_options(n::Integer; colors =colormap("RdBu", n)[end:-1:1], opacities = [1], args...)
	return (colors = colors, opacities = opacities, args...)
end

























# Plotting
####################################################################################################

export bin_analysis_plot, peak_analysis_plot, uncertainty_plot


"""
	bin_analysis_plot(results::Union{AdaptiveResult,Nothing},
		data = nothing,
		replicates = nothing; 
		keywords...
	)
Create and return basic plots `(dr_plot, density_plot)` for the [`DoseResponseResult`](@ref) and the K_τ-density `gird`. 

If `results` is an [`AdaptiveResult`](@ref), both the fitted gird is plotted into the density plot and the corresponding curve is plotted into the dose-response plot.

If `data` is a [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object the data points are plotted into the dose-response plot. Similarly, if `replicates` is an array of [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects, the data points are plotted as replicates in the dose-response plot.

**Keywords**

* `dr_plot = `[`dr_base_plot()`](@ref): The base plot onto which the `AdaptiveResult.result` and the [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects (`data` and `replicates`) are plotted. It can be any [`Plots.jl`](https://docs.juliaplots.org/stable/) plot (e.g. another dose-response plot). 
* `density_plot = `[`density_base_plot()`](@ref): The base plot onto which the `AdaptiveResult.grid` is plotted. It can be any [`Plots.jl`](https://docs.juliaplots.org/stable/) plot (e.g. another `K_τ` density plot).
* `fit_arguments = `[`fit_options()`](@ref): Keyword argument tuple for the `AdaptiveResult.result` data-series.
* `data_arguments = `[`data_options()`](@ref): Keyword argument tuple for the [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) data-series.
* `replicate_arguments = `[`replicate_options()`](@ref): Keyword argument tuple for the replicate data-series.
* `density_arguments = `[`density_options()`](@ref): Keyword argument tuple for the `AdaptiveResult.grid` data-series.
* `annotation_arguments = NamedTuple()`: Keyword argument tuple for bin-annotations in the density plot. See below for further information.

**Annotation bins**

Annotation bins allow to annotate the density plot with the number of epitopes (in units of the density values) in the respective bin. The following keywords can be used for the `annotation_arguments`:

* `annotation_bins = []`: The bins as grid domain ranges, e.g. `[[1e-10,1e-9], [1e-9,1e-8], [1e-5,1e-2]]`.
* `annotation_size = 10`: Size of the annotation font.
* `annotation_offset = 0.05`: Relative vertical offset from the top of the plot for the annotation labels. If a single number is provided, every other annotation is offsetted.  If a length-matched vector of numbers is provided, each annotation is individually offsetted accordingly.
* `annotation_color = :black`: Color of the annotations.
* `hover_points = false`: If true, adds scatter-points with tool-tips for the [Plotly/PlotlyJs backend](https://docs.juliaplots.org/latest/backends/).

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
Create and return plots to analyze the effect of peaks in the K_τ density on the corresponding dose-response curve.

Returns `(individual_dr_plot, cumulative_dr_plot, density_plot)`, where

* `density_plot` contains a plot of the K_τ-density with different colors for the different peaks (specified by the bins).
* `individual_dr_plot` contains the individual dose-response curves (color matched) that originate from the different peaks.
* `cumulative_dr_plot` contains the cumulative dose-response curves, i.e. dose-responses include the response increases caused by peaks with smaller K_τ. Again the curves are color matched with the peaks.

If `data` is a [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object the data points are plotted in the `cumulative_dr_plot`.

**Keywords**

* `individual_dr_plot = `[`dr_base_plot()`](@ref): The base plot onto which the individual dose-response curves are plotted. It can be any [`Plots.jl`](https://docs.juliaplots.org/stable/) plot (e.g. another dose-response plot). 
* `cumulative_dr_plot = `[`dr_base_plot()`](@ref): The base plot onto which the cumulative dose-response curves are plotted. It can be any [`Plots.jl`](https://docs.juliaplots.org/stable/) plot (e.g. another dose-response plot). 
* `density_plot = `[`density_base_plot()`](@ref): The base plot onto which the `K_τ`-peaks are plotted. It can be any [`Plots.jl`](https://docs.juliaplots.org/stable/) plot (e.g. another `K_τ`-density plot). 
* `fit_arguments = `[`fit_options()`](@ref): Keyword argument tuple for the `AdaptiveResult.result` data-series. Both `color` and `label` get overwritten by the `colors` keyword and the peak number (automatically determined).
* `density_arguments = `[`density_options()`](@ref): Keyword argument tuple for the `AdaptiveResult.grid` data-series. Both `color` and `label` get overwritten by the `colors` keyword and the peak number (automatically determined).
* `bins = `[`peak_detection`](@ref)`(results.grid, fill = false)[2]`: The bins that define the K_τ-peaks as grid domain ranges, e.g. `[[1e-10,1e-9], [1e-9,1e-8], [1e-5,1e-2]]`.
* `colors = collect(1:length(bins))`: The colors for the different bins.
* `join_bins = true`: If true, extends the bins if needed, s.t. there remains no gap between the bins.
"""
function peak_analysis_plot(results::AdaptiveResult,data::Union{FittingData, Nothing} = nothing;
	individual_dr_plot = dr_base_plot(),
	cumulative_dr_plot = dr_base_plot(),
	fit_arguments = fit_options(),
	density_arguments = density_options(),
	density_plot = density_base_plot(),
	bins = peak_detection(results.grid, fill = false)[2],
	colors = collect(1:length(bins)),
	join_bins = true
	)
	
	
	if !isnothing(data)
		cumulative_dr_plot = plot(cumulative_dr_plot)
		scatter!(data, color = :black, label = "data")
	end

	local_fit_keys = [k for k in keys(fit_arguments) if k !== :color && k !== :label]
	local_fit_arguments = NamedTuple(zip(local_fit_keys,[fit_arguments[k] for k in local_fit_keys]))
	local_density_keys = [k for k in keys(density_arguments) if k !== :color && k !== :label]
	local_density_arguments = NamedTuple(zip(local_density_keys,[density_arguments[k] for k in local_density_keys]))

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
		plot!(DensityPlot(single_peak_grid); color = colors[i],  label = "peak $i", local_density_arguments...)

		if length(results.optimizer) > length(results.grid)
			dose_response = DoseResponseResult(single_peak_grid,results.result.concentrations, offset = results.optimizer[end])
		else
			dose_response = DoseResponseResult(single_peak_grid,results.result.concentrations)
		end
		individual_dr_plot = plot(individual_dr_plot)
		plot!(dose_response; color = colors[i], label = "peak $i", local_fit_arguments...)

		cumulative_peak_grid = restrict_domain!(deepcopy(results.grid), upper = sorted_bins[i][2], weight_distribution = :log)
		if length(results.optimizer) > length(results.grid)
			dose_response = DoseResponseResult(cumulative_peak_grid,results.result.concentrations, offset = results.optimizer[end])
		else
			dose_response = DoseResponseResult(cumulative_peak_grid,results.result.concentrations)
		end
		cumulative_dr_plot = plot(cumulative_dr_plot)
		plot!(dose_response; color = colors[i], label = "peak $i", local_fit_arguments...)
	end
	return individual_dr_plot, cumulative_dr_plot, density_plot
end












"""
	uncertainty_plot(e_uncertainty::EpitopeUncertainty,
		d_uncertainty::DoseResponseUncertainty,
		grid::OneDimGrid; 
		keywords...
	)

Create and return uncertainty visualizations `(dr_uncertainty_plot, density_uncertainty_plot)`.

The estimated bounds of the [`DoseResponseUncertainty`](@ref) and [`EpitopeUncertainty`](@ref) are plotted as color-matched ribbons.

**Keywords**

* `dr_plot = `[`dr_base_plot()`](@ref): The base plot onto which the `AdaptiveResult.result` and uncertainty ribbons are plotted. It can be any [`Plots.jl`](https://docs.juliaplots.org/stable/) plot (e.g. another dose-response plot). 
* `density_plot = `[`density_base_plot()`](@ref): The base plot onto which the `AdaptiveResult.grid` and uncertainty ribbons are plotted. It can be any [`Plots.jl`](https://docs.juliaplots.org/stable/) plot (e.g. another `K_τ` density plot). 
* `eu_arguments = `[`eu_options`](@ref)`(length(e_uncertainty.levels))`: Keyword argument tuple for the [`EpitopeUncertainty`](@ref) data series. See [`EpitopeUncertainty` - plotting](@ref EpitopeUncertainty-plotting) for further information.
* `du_arguments = `[`du_options`](@ref)`(length(d_uncertainty.levels))`: Keyword argument tuple for the [`DoseResponseUncertainty`](@ref) data series. See [`DoseResponseUncertainty` - plotting](@ref DoseResponseUncertainty-plotting) for further information.

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