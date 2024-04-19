using AntibodyMethodsDoseResponseConvenience
using Test, Statistics

cd(@__DIR__)

@testset "AntibodyMethodsDoseResponseConvenience.jl" begin

    @testset "Auxiliary functions" begin

        # log_volumes(centers, volumes) | should calculate the volume of the intervals with logarithmic boundaries [log(left), log(right)].
        log_volumes = AntibodyMethodsDoseResponseConvenience.log_volumes
        @test log_volumes([10,20,30],[2,4,6]) == [log10(11)-log10(9), log10(22)-log10(18), log10(33)-log10(27)]

        # grid_range(concentrations, lower_shift, upper_shift, exponent_rounding)
        # Should return the range of the concentrations, adjusted by the shifts and rounded exponents (if exponent_rounding = true).
        grid_range = AntibodyMethodsDoseResponseConvenience.grid_range
        @test grid_range([2e-10,2e-8,2e-6,1e-4,2e-2],-1,+1,true) == (1e-11,1e-0)
        @test prod(grid_range([2e-10,2e-8,2e-6,1e-4,2e-2],0,0,false) .≈ (2e-10,2e-2))
        # Zero concentration should raise a warning:
        @test_warn "Automatic" grid_range([0,2e-10,2e-8,2e-6,1e-4,2e-2], 0,0,false)
        @test grid_range([0,2e-10,2e-8,2e-6,1e-4,2e-2],0,0,false)[1] ≈ eps()




    
    end



    @testset "scaled_log_volume_prior" begin
        # scaled_log_volume_prior(scale) creates a prior_generator which itself returns a prior function.
        # prior should be -scale *(offset^2 +  ∑_i (λ[i] / log_vol[i] - λ[i+1]/log_vol[i+1])^2 ).
        log_volumes = AntibodyMethodsDoseResponseConvenience.log_volumes
        centers = [10,20,30]
        volumes = [5,6,8]
        lv = log_volumes(centers,volumes)

        # Without offset (3-element parameter).
        P(λ) = (-(λ[1]/lv[1]-λ[2]/lv[2])^2 - (λ[3]/lv[3]-λ[2]/lv[2])^2)/9
        # With offset (4-element parameter).
        Q(λ) = (-(λ[1]/lv[1]-λ[2]/lv[2])^2 - (λ[3]/lv[3]-λ[2]/lv[2])^2 - λ[4]^2)/16

        prior_generator = scaled_log_volume_prior(4)
        p = prior_generator(centers,volumes,nothing)
        q = prior_generator(centers,volumes,1) # offset != nothing

        # P and Q are unscaled -> multiply with scale = 4.
        @test p([1,2,3]) .≈ 4*P([1,2,3])
        @test q([1,2,3,4]) .≈ 4*Q([1,2,3,4])
    end


    @testset "minimizer_generator" begin
        # Test minimizer with simple function. Check if minimization is within the unit square of the theoretical minimum.

        # Gradient free minimizer.
        minimizer = minimizer_generator(NelderMead())
        @test prod(minimizer(x-> sum(x .^2), nothing,[2.0,2.0,2.0]) .< [1,1,1])

        # Gradient minimizer.
        minimizer = minimizer_generator(BFGS(), gradient = true)
        g! = function(grad,x) grad .= 2 .* x end
        @test prod(minimizer(x-> sum(x .^2), g! ,[2.0,2.0,2.0]) .< [1,1,1])

        # Test if gradient option picks the correct minimizer (Optim fails if nothing is passed as gradient).
        @test_throws MethodError minimizer(x-> sum(x .^2), nothing ,[2.0,2.0,2.0]) 
    end






    @testset "FittingCondition" begin
        
        # Check that dose_response_check is used.

        @test_throws DomainError FittingCondition(FittingData([1,2,3],[-1,1,2]))
        @test_throws DomainError FittingCondition([1,2,3],[-1,1,2],[3,2,4])

        concentrations = [2e-10,2e-8,2e-6,1e-4,2e-2]
        responses = [2e-10,2e-8,2e-6,1e-4,2e-2]
        data = FittingData(concentrations,responses)
        

        # Test different default constructors.
            # Test options by checking objective == :log_posterior (default would be :lsq).
            # Test minimizer with simple function (compare to minimizer_generator which was tested above).

        # Explicit constructor and immutability tests.
        m_data = FittingData([1,2,3],[1,2,3])
        m_path = "path"
        m_grid = create_grid([1,2,3])
        m_options = AdaptiveOptions()
        m_minimizer = x -> x^2
        m_result_concentrations = [1,1.5,2,2.5,3]
        # scale is not used as it alters the options (tested below) which is no issue if input data is properly copied (detached from original objects).

        condition = FittingCondition(m_data, path = m_path, grid = m_grid, options_1 = m_options, options_2 = m_options, minimizer_1  = m_minimizer, minimizer_2 = m_minimizer, result_concentrations = m_result_concentrations)
        
        # Mutate the original objects.
        m_data.independent = [4,5,6]
        m_path = "other path"
        refine!(m_grid)
        m_options.objective = :posterior
        m_minimizer = x -> -x
        push!(m_result_concentrations,5)

        # Testing object equality not possible, as recreating the original objects are new, different instances.
        @test condition.data.independent == [1,2,3]
        @test condition.path == "path"
        @test export_all(condition.grid) == export_all(create_grid([1,2,3]))
        @test condition.options_1.objective == :lsq
        @test condition.options_2.objective == :lsq
        @test condition.minimizer_1(3) == 9
        @test condition.minimizer_2(3) == 9
        @test condition.result_concentrations == [1,1.5,2,2.5,3]




        # Constructor using FittingData.
        condition = FittingCondition(data)
        @test condition.data.independent == data.independent
        @test condition.options_1.objective == :lsq
        @test condition.options_2.objective == :lsq
        @test isnothing(condition.replicates)
        @test length(condition.grid) == 2
        @test condition.minimizer_1(x-> sum(x),nothing,[1.0,2.0,3.0]) ≈ minimizer_generator(NelderMead())(x-> sum(x), nothing, [1.0,2.0,3.0])
        @test condition.minimizer_2(x-> sum(x),nothing,[1.0,2.0,3.0]) ≈ minimizer_generator(BFGS())(x-> sum(x), nothing, [1.0,2.0,3.0])
        @test isnothing(condition.result_concentrations)

        # Constructor using replicates | only single replicate.
        condition = FittingCondition(concentrations,responses)
        @test condition.data.independent == concentrations
        @test condition.options_1.objective == :lsq
        @test condition.options_2.objective == :lsq
        @test isnothing(condition.replicates)
        @test length(condition.grid) == 2
        @test condition.minimizer_1(x-> sum(x),nothing,[1.0,2.0,3.0]) ≈ minimizer_generator(NelderMead())(x-> sum(x), nothing, [1.0,2.0,3.0])
        @test condition.minimizer_2(x-> sum(x),nothing,[1.0,2.0,3.0]) ≈ minimizer_generator(BFGS())(x-> sum(x), nothing, [1.0,2.0,3.0])
        @test isnothing(condition.result_concentrations)

        # Constructor using replicates | multiple replicates.
        condition = FittingCondition(concentrations,responses,responses .+ 1)
        @test condition.data.independent == concentrations
        @test length(condition.replicates) == 2
        @test condition.options_1.objective == :lsq
        @test condition.options_2.objective == :lsq
        @test prod(condition.data.dependent .≈ responses .+ 0.5)
        @test prod(condition.data.errors .≈ std.([[r, r+1] for r in responses]))
        @test condition.minimizer_1(x-> sum(x),nothing,[1.0,2.0,3.0]) ≈ minimizer_generator(NelderMead())(x-> sum(x), nothing, [1.0,2.0,3.0])
        @test condition.minimizer_2(x-> sum(x),nothing,[1.0,2.0,3.0]) ≈ minimizer_generator(BFGS())(x-> sum(x), nothing, [1.0,2.0,3.0])
        @test isnothing(condition.result_concentrations)


        # Test effect of keywords
            # Use default AdaptiveOptions() -> objective now :lsq.
            # Use trivial optimizer x->x.
            # Use scale (compare to scaled_log_volume_prior).

        # Constructor using FittingData is already tested (mutability test above).

        # Constructor using replicates.
        condition = FittingCondition(concentrations,responses,options_1 = AdaptiveOptions(), options_2 = AdaptiveOptions(), minimizer_1 = x->x, minimizer_2 = x->x, scale = 3, result_concentrations = [1,2,3])
        @test condition.options_1.objective == :log_posterior
        @test condition.options_2.objective == :log_posterior
        @test condition.options_1.prior_generator([10,20,30],[1,2,3],nothing)([2,3,4]) == scaled_log_volume_prior(3)([10,20,30],[1,2,3],nothing)([2,3,4])
        @test condition.minimizer_1(10) == 10
        @test condition.minimizer_2(10) == 10
        @test condition.result_concentrations == [1,2,3]

        # Test mean and error calculation.

        # Single response.
        condition = FittingCondition(concentrations,responses)
        @test condition.data.dependent == responses
        @test condition.data.errors == ones(length(concentrations))

        # Well-behaved replicates.
        condition = FittingCondition(concentrations,responses, responses .+ 0.5)
        @test prod(condition.data.dependent .≈ responses .+ 0.25)
        @test prod(condition.data.errors .≈ std.([[r, r + 0.5] for r in responses]))

        # Tiny differences.
        condition = FittingCondition(concentrations,responses, responses .+ 1e-17)
        @test prod(condition.data.dependent .≈ responses .+ 1e-17/2)
        @test prod(condition.data.errors .≈ std.([[r, r + 1e-17] for r in responses]) )
        @test prod(condition.data.errors .!= zeros(length(concentrations)))

        # Identical replicates.
        condition = FittingCondition(concentrations,responses, responses)
        @test prod(condition.data.dependent .≈ responses)
        @test prod(condition.data.errors .≈ eps() .* ones(length(concentrations)) )
        @test prod(condition.data.errors .≈ eps() .* ones(length(concentrations)) )
        @test prod(condition.data.errors .!=  zeros(length(concentrations)) )
    end


    @testset "fit_condition" begin
        # Test return type (not result details), adaptive fitting by increased grid size, and saved files.

        density_grid = create_grid(LogRange(1e-8,1e-2,20))
        dr = DoseResponseResult(density_grid,LogRange(1e-8,1e-2,20))

        # Condition without path and without result_concentrations.
        condition = FittingCondition(dr.concentrations,dr.responses, options_1 = AdaptiveOptions(iterations = 2))

        result = fit_condition(condition)
        @test typeof(result.result) <: DoseResponseResult
        @test result.result.concentrations == dr.concentrations
        @test 2 < length(result.grid) < length(density_grid)
        # Without path, no data should be written.
        @test !("test" in readdir())

        # Condition without path but with result_concentrations.
        condition = FittingCondition(dr.concentrations,dr.responses, options_1 = AdaptiveOptions(iterations = 2), result_concentrations = LogRange(1e-10,1e-2,40))

        result = fit_condition(condition)
        @test typeof(result.result) <: DoseResponseResult
        @test result.result.concentrations == LogRange(1e-10,1e-2,40)
        @test 2 < length(result.grid) < length(density_grid)
        # Without path, no data should be written.
        @test !("test" in readdir())

        # Condition with path.
        condition = FittingCondition(dr.concentrations,dr.responses, options_1 = AdaptiveOptions(iterations = 2), path ="test")

        result = fit_condition(condition)
        @test typeof(result.result) <: DoseResponseResult
        @test 2 < length(result.grid) < length(density_grid)
        @test "test" in readdir()
        @test prod(file in readdir("test") for file in ["data.jld", "result.jld", "grid.jld", "optimizer.jld","objective_value.jld", "time.jld"])

        # Remove data to avoid state-dependent tests.
        rm("test", recursive = true)


        # Test that fitting multiple conditions does not throw an error.
        condition = FittingCondition(dr.concentrations,dr.responses, options_1 = AdaptiveOptions(iterations = 2))
        conditions = fill(condition,10)
        @test isnothing(fit_conditions(conditions))
        
    end




    @testset "load_results" begin
        density_grid = create_grid(LogRange(1e-8,1e-2,20))
        dr = DoseResponseResult(density_grid,LogRange(1e-8,1e-2,20))

        # Use fit_condition (tested above) to create correct files.
        condition = FittingCondition(dr.concentrations,dr.responses, options_1 = AdaptiveOptions(iterations = 2), path ="test")
        result = fit_condition(condition)
        
        loaded_result, data, replicates = load_results("test")

        @test loaded_result.result.concentrations == result.result.concentrations
        @test loaded_result.result.responses == result.result.responses
        @test loaded_result.time == result.time
        @test data.independent == dr.concentrations
        @test isnothing(replicates)
        
        rm("test", recursive = true)
    end


end
