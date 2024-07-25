####################################################################################################
# Convenience methods to fit data.
####################################################################################################





# Internal functions
####################################################################################################

# Obtain the visual volumes in a logarithmic scale.
function log_volumes(centers,volumes)
	return log10.(centers .+ volumes/2) .- log10.(centers .- volumes/2)
end



# Determine grid range (for fitting) from the concentration range.
# Shifts (of the exponent) allow to increase/decrease the grid range.
# round_exponent determines if the exponents are rounded to integer values.
function grid_range(concentrations,lower_shift::Real,upper_shift::Real,round_exponent::Bool)

	lower, upper = extrema(concentrations)
	
	if iszero(lower)
		@warn("Automatic grid range estimation encountered concentration = 0. To avoid -Inf in the log-scale, eps() was used instead. Nevertheless, too large gird ranges should be avoided. It is recommended to provide an explicitly created grid.")
		lower = eps()
	end
	
	log_lower = log10(lower)
	log_upper = log10(upper)
	if round_exponent
		log_lower = floor(log_lower)
		log_upper = ceil(log_upper)
	end
	return 10.0^(log_lower+lower_shift), 10.0^(log_upper+upper_shift)
end











# Exported functions
####################################################################################################

export scaled_log_volume_prior, minimizer_generator, FittingCondition, fit_condition, fit_conditions



@doc raw"""
	scaled_log_volume_prior(scale::Real = 1)
Create a prior generator (see [`AdaptiveOptions`](@ref)) that generates the following prior:

```math
\text{prior}(\lambda) = - \frac{\text{scale}}{\text{length}(λ)^2} \cdot \left( \text{offset}^2 + \sum_{i=2}^{\text{length}(\lambda)} \left(\frac{\lambda_{i-1}}{\log_{10}(r_{i-1})-\log_{10}(l_{i-1})} - \frac{\lambda_{i}}{\log_{10}(r_{i})-\log_{10}(l_{i})}\right)^2 \right)
```

where ``[l_i,r_i]`` are the intervals corresponding to ``\lambda_i``.
"""
function scaled_log_volume_prior(scale::Real = 1)
	return function(centers, volumes, offset) 
		lV = log_volumes(centers,volumes)
		if isnothing(offset)
			return λ ->  -scale*(sum((λ[i]/lV[i]-λ[i+1]/lV[i+1])^2 for i in 1:length(λ)-1))/length(λ)^2
		else
			return λ -> -scale*(sum((λ[i]/lV[i]-λ[i+1]/lV[i+1])^2 for i in 1:length(λ)-2) + λ[end]^2)/length(λ)^2
		end
	end
end


"""
	minimizer_generator(optim_minimizer; 
		options = Optim.Options(g_tol = 1e-12, iterations = 2000), 
		gradient::Bool = false
	)
	
Create a minimization function `(f,∇f,initial_point) -> minimizing_point` as specified in [`adaptive_dose_response_fit`](@ref), using minimizers from [`Optim.jl`](https://julianlsolvers.github.io/Optim.jl/stable/) (e.g. `NelderMead()` or `LBFGS()`).

If `gradient = false`, the gradient function `∇f` is ignored, useful e.g. for [`adaptive_dose_response_fit`](@ref) where `∇f = nothing` is internally passed as argument in some cases.

**Examples**

```julia
minimizer_generator(NelderMead())
```
```julia
minimizer_generator(LBFGS(), 
	options = Optim.Options(g_tol = 1e-6, iterations = 400), 
	gradient = true
)
```
"""
function minimizer_generator(optim_minimizer; options = Optim.Options(g_tol = 1e-12, iterations = 2000), gradient::Bool = false)
	if !gradient
		return function(f,∇f,init)
			lower = zeros(length(init))
			upper = fill(Inf, length(init))
			return optimize(f,lower,upper, init, Fminbox(optim_minimizer),options).minimizer
		end
	else
		return function(f,∇f,init)
			lower = zeros(length(init))
			upper = fill(Inf, length(init))
			return optimize(f,∇f,lower,upper, init, Fminbox(optim_minimizer),options).minimizer
		end
	end
end




"""
	mutable struct FittingCondition

Data type storing the necessary information for the common workflow of [`adaptive_dose_response_fit`](@ref). Convenience constructors with recommended options are implemented. Note that all objects are "deepcopied" to avoid unwanted mutation.

**Fields**

* `data`: The [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object containing the dose-response data.
* `replicates`: Vector of [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects that constitute the replicates of measured dose-response curves.
* `grid`: Initial `K_τ` grid that is adaptively refined.
* `path`: Fitting results are saved to `path` if `path != ""`.
* `options_1`: [`AdaptiveOptions`](@ref) for the first run of [`adaptive_dose_response_fit`](@ref).
* `options_2`: [`AdaptiveOptions`](@ref) for a second run of [`adaptive_dose_response_fit`](@ref), e.g. to use a final, non-adaptive, gradient-based fit.
* `minimizer_1`: Minimization function for the first run of [`adaptive_dose_response_fit`](@ref).
* `minimizer_2`: Minimization function for the second run of [`adaptive_dose_response_fit`](@ref).
* `result_concentrations`: Concentrations to be used for the dose-response curve calculated from the fit result. This allows to obtain smooth result curves. If `result_concentrations = nothing`, the concentrations of `data` are used.

**Convenience constructors**

	FittingCondition(data::FittingData, replicates = nothing; keywords...)
Manual specification of the [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object and the optional replicates. Recommended fitting options are predefined and `path=""` is set to avoid accidental creation of files. 


	FittingCondition(concentrations::AbstractVector, response_replicates::AbstractVector...; 
		keywords...)

Construct a FittingCondition from a concentration vector and response vectors (variable number of arguments). This automatically creates the main [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object, and the vector of [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects for the replicates.

The main [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object uses the mean values of the responses together with the standard deviations as uncertainties. The uncertainty distributions are unnormalized logarithmic normal distributions:

	(y,m,Δy)-> -(y-m)^2/Δy^2

**Keywords**

The keywords correspond to the struct fields, except for the additional `scale` keyword. Setting `scale` overwrites the `objective` field to `:log_posterior` and the `prior_generator` field to [`scaled_log_volume_prior`](@ref) for both options: `options_1` and `options_2`.

The keyword defaults are:

* `grid = create_grid(LogRange(extrema(concentrations)...,3))`
* `path = ""`
* `options_1 = AdaptiveOptions(objective = :lsq, offset = eps(), iterations = 30)`
* `options_2 = AdaptiveOptions(objective = :lsq, offset = eps())`
* `minimizer_1 = minimizer_generator(NelderMead())`
* `minimizer_2 = minimizer_generator(LBFGS())`
* `result_concentrations = nothing`
* `scale = nothing`
"""
mutable struct FittingCondition
	data::FittingData
	replicates::Union{AbstractVector{FittingData},Nothing}
	grid::AdaptiveDensityApproximation.OneDimGrid
	path::AbstractString
	options_1::AdaptiveOptions
	options_2::AdaptiveOptions
	minimizer_1::Function
	minimizer_2::Union{Function,Nothing}
	result_concentrations

	function FittingCondition(data::FittingData,
		replicates = nothing;
		path="", 
		grid = create_grid(LogRange(grid_range(data.independent,-1,0,true)...,3)),
		options_1::AdaptiveOptions = AdaptiveOptions(objective = :lsq, offset = eps(), iterations = 30),
		options_2::AdaptiveOptions = AdaptiveOptions(objective = :lsq, offset = eps()),
		minimizer_1::Function = minimizer_generator(NelderMead()),
		minimizer_2 = minimizer_generator(LBFGS()) , 
		result_concentrations::Union{Nothing,AbstractArray{T,N}} = nothing,
		scale::Union{Nothing,Real}=nothing) where {N, T <: Real}

		# Default constructors should check for proper dose-response data. Improper dose-response data can be used by explicit mutation if need be.
		dose_response_check(data)
		if !isnothing(replicates)
			for replicate in replicates
				dose_response_check(replicate)
			end
		end
		
		if !isnothing(scale)
			options_1.prior_generator = scaled_log_volume_prior(scale)
			options_1.objective = :log_posterior
			options_2.prior_generator = scaled_log_volume_prior(scale)
			options_2.objective = :log_posterior
		end

		if !isnothing(result_concentrations)
			AntibodyMethodsDoseResponse.regular_positive_numbers_check(result_concentrations)
		end

		return new(	deepcopy(data),
					deepcopy(replicates),
					deepcopy(grid),
					deepcopy(path),
					deepcopy(options_1),
					deepcopy(options_2),
					deepcopy(minimizer_1),
					deepcopy(minimizer_2), 
					deepcopy(result_concentrations)
				)
	end
end






# Documentation in struct docstring.
function FittingCondition(concentrations::AbstractVector{T}, response_replicates::AbstractVector{S}...; args...) where {T <: Real, S <: Real}
	if length(response_replicates) == 1
		if length(concentrations) != length(response_replicates[1])
			throw(DimensionMismatch("Number of concentrations and number of responses do not match!"))
		end
		return FittingCondition(FittingData(concentrations,response_replicates[1], distributions = (y,m,Δy)-> -(y-m)^2/Δy^2); args...)
	end

	for (i,responses) in enumerate(response_replicates)
		if length(concentrations) != length(responses)
			throw(DimensionMismatch("Number of concentrations and number of responses for the $i-th replicate do not match!"))
		end
	end
	
	mean_responses = [mean([response_replicates[i][j] for i in eachindex(response_replicates)]) for j in eachindex(concentrations)]
	errors = [std([response_replicates[i][j] for i in eachindex(response_replicates)]) for j in eachindex(concentrations)]
	replicates = [FittingData(concentrations,response_replicates[i]) for i in eachindex(response_replicates)]

	# If all samples are the same, the std is 0, which causes problems during the fitting.
	for i in eachindex(errors)
		if iszero(errors[i])
			errors[i] = eps()
		end
	end

	return FittingCondition(FittingData(concentrations,mean_responses,errors, distributions = (y,m,Δy)-> -(y-m)^2/Δy^2), replicates; args...)
end







"""
	fit_condition(condition::FittingCondition)
Obtain the results for a `FittingCondition` object (i.e. fitting the data).

Returns the [`AdaptiveResult`](@ref) object and saves the data ([`FittingCondition`](@ref) and [`AdaptiveResult`](@ref) objects) to `FittingCondition.path` if not `FittingCondition.path = ""`.
"""
function fit_condition(condition::FittingCondition)
	fitting_grid = deepcopy(condition.grid)
	result = adaptive_dose_response_fit(condition.grid,condition.data, condition.minimizer_1,options = condition.options_1)
	fitting_time = result.time

	if !isnothing(condition.minimizer_2)
		try
			if isnothing(condition.options_1.offset)
				result = adaptive_dose_response_fit(result.grid,condition.data, condition.minimizer_2,options = condition.options_2)
			else
				# Overwrite offset in options_2 to use the estimated offset from the previous adaptive fit.
				condition.options_2.offset = result.optimizer[end]
				result = adaptive_dose_response_fit(result.grid,condition.data, condition.minimizer_2,options = condition.options_2)
			end
			fitting_time += result.time
			result.time = fitting_time


		catch
			@warn("Second optimization was skipped. $(condition.path)")

		end
	end

	# Re-construct dose-response curve with result_concentrations if they are not nothing.
	if !isnothing(condition.result_concentrations)
		if !isnothing(condition.options_1.offset)
			dr_result = DoseResponseResult(result.grid,condition.result_concentrations, model = condition.options_1.model, offset = result.optimizer[end])
		else
			dr_result = DoseResponseResult(result.grid,condition.result_concentrations, model = condition.options_1.model)
		end
		result.result = dr_result
	end


	if !isempty(condition.path)
		mkpath(condition.path)
		serialize(joinpath(condition.path,"data.jld"), condition.data)
		serialize(joinpath(condition.path,"result.jld"), result.result)
		serialize(joinpath(condition.path,"grid.jld"), result.grid)
		serialize(joinpath(condition.path,"optimizer.jld"), result.optimizer)
		serialize(joinpath(condition.path,"objective_value.jld"), result.objective_value)
		serialize(joinpath(condition.path,"time.jld"), result.time)
		if !isnothing(condition.replicates)
			serialize(joinpath(condition.path,"replicates.jld"), condition.replicates)
		end
	end

	return result
end


"""
	fit_conditions(conditions)
Multithreaded application of [`fit_condition`](@ref) to a collection of [`FittingCondition`](@ref) objects (`conditions`).
"""
function fit_conditions(conditions)
	fitting_progress = Progress(length(conditions))
	Threads.@threads for i in eachindex(conditions)
		fit_condition(conditions[i])
		next!(fitting_progress)
	end
	finish!(fitting_progress)
	return nothing
end
