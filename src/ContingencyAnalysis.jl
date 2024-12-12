module ContingencyAnalysis

using HypothesisTests, Statistics, Random

export monte_carlo_fisher,
    calculate_chi_statistic,
    simulate_contingency_table,
    analyse_matrix,
    analyse_contingency

"""
    monte_carlo_fisher(matrix::Matrix{T}, n_simulations::Int=10000, alternative::Symbol=:two) where {T<:Number}

Performs Monte Carlo simulation of Fisher's exact test with specified alternative hypothesis.

# Arguments
- `matrix`: Input contingency table
- `n_simulations`: Number of Monte Carlo simulations for Fisher's test
- `alternative`: Type of test (:two for two-sided, :less or :greater for one-sided)

# Returns
- Estimated p-value from the Monte Carlo simulation
"""
function monte_carlo_fisher(matrix::Matrix{T}, n_simulations::Int=10000,
    alternative::Symbol=:two) where {T<:Number}

    if !(alternative in [:two, :less, :greater])
        throw(ArgumentError("alternative must be :two, :less, or :greater"))
    end

    observed_stat = calculate_chi_statistic(matrix)

    if alternative == :two || size(matrix) != (2, 2)
        # Use chi-square statistic for two-sided tests or tables larger than 2x2
        n_greater = count(1:n_simulations) do _
            sim_matrix = simulate_contingency_table(sum(matrix, dims=2), sum(matrix, dims=1))
            sim_stat = calculate_chi_statistic(sim_matrix)
            sim_stat >= observed_stat
        end
    else
        # Use odds ratio for one-sided tests with 2x2 tables
        obs_or = (matrix[1, 1] * matrix[2, 2]) / (matrix[1, 2] * matrix[2, 1])
        n_greater = count(1:n_simulations) do _
            sim_matrix = simulate_contingency_table(sum(matrix, dims=2), sum(matrix, dims=1))
            sim_or = (sim_matrix[1, 1] * sim_matrix[2, 2]) / (sim_matrix[1, 2] * sim_matrix[2, 1])
            alternative == :greater ? sim_or >= obs_or : sim_or <= obs_or
        end
    end

    return n_greater / n_simulations
end

"""
    calculate_chi_statistic(matrix::Matrix{T}) where {T<:Number}

Calculates chi-square statistic for a contingency table.

# Returns
- Chi-square statistic value
"""
function calculate_chi_statistic(matrix::Matrix{T}) where {T<:Number}
    row_sums = sum(matrix, dims=2)
    col_sums = sum(matrix, dims=1)
    expected = (row_sums * col_sums) ./ sum(matrix)
    return sum((matrix .- expected) .^ 2 ./ expected)
end

"""
    simulate_contingency_table(row_sums::AbstractMatrix, col_sums::AbstractMatrix)

Simulates a contingency table preserving given marginal sums using probability-based allocation.

# Returns
- Simulated contingency table maintaining row and column sums
"""
function simulate_contingency_table(row_sums::AbstractMatrix, col_sums::AbstractMatrix)
    simulated = zeros(eltype(row_sums), size(row_sums, 1), size(col_sums, 2))
    remaining_rows = copy(vec(row_sums))
    remaining_cols = copy(vec(col_sums))

    while sum(remaining_rows) > 0
        for i in 1:size(simulated, 1), j in 1:size(simulated, 2)
            if remaining_rows[i] > 0 && remaining_cols[j] > 0
                prob = (remaining_rows[i] * remaining_cols[j]) /
                       (sum(remaining_rows) * sum(remaining_cols))
                value = rand() < prob ? 1 : 0
                if value == 1
                    simulated[i, j] += 1
                    remaining_rows[i] -= 1
                    remaining_cols[j] -= 1
                end
            end
        end
    end
    return simulated
end

"""
    analyse_matrix(matrix::Matrix{T}) where {T<:Number}

Analyses contingency table characteristics and removes zero marginals if present.

# Returns
- Named tuple containing cleaned matrix and its properties
"""
function analyse_matrix(matrix::Matrix{T}) where {T<:Number}
    non_zero_rows = vec(sum(matrix, dims=2) .> 0)
    non_zero_cols = vec(sum(matrix, dims=1) .> 0)

    cleaned_matrix = matrix[non_zero_rows, non_zero_cols]
    was_cleaned = !all(non_zero_rows) || !all(non_zero_cols)

    if was_cleaned
        println("\nRemoved $(count(.!non_zero_rows)) rows and $(count(.!non_zero_cols)) columns with zero marginals")
    end

    expected = (sum(cleaned_matrix, dims=2) * sum(cleaned_matrix, dims=1)) ./ sum(cleaned_matrix)
    prop_small = sum(expected .< 5) / prod(size(cleaned_matrix))

    return (cleaned_matrix=cleaned_matrix, prop_small=prop_small, total=sum(cleaned_matrix))
end

"""
    analyse_contingency(matrix::Matrix{T}, n_simulations::Int=10000; alternative::Symbol=:two) where {T<:Number}

Performs appropriate statistical test (Fisher's or Chi-square) based on matrix characteristics.

# Arguments
- `matrix`: Input contingency table
- `n_simulations`: Number of Monte Carlo simulations for Fisher's test
- `alternative`: Type of test (:two for two-sided, :less or :greater for one-sided)

# Returns
- Nothing (prints results to stdout)
"""
function analyse_contingency(matrix::Matrix{T}, n_simulations::Int=10000;
    alternative::Symbol=:two) where {T<:Number}
    if !(alternative in [:two, :less, :greater])
        throw(ArgumentError("alternative must be :two, :less, or :greater"))
    end

    characteristics = analyse_matrix(matrix)
    matrix = characteristics.cleaned_matrix

    # Add input validation
    if any(x -> x < 0, matrix)
        throw(ArgumentError("Contingency table cannot contain negative values"))
    end

    # Define threshold for exact vs Monte Carlo Fisher's test
    fisher_exact_threshold = 2000  # Maximum size for exact computation
    total_cells = prod(size(matrix))

    if alternative != :two && size(matrix) != (2, 2)
        @warn "One-sided tests are only meaningful for 2x2 tables. Using two-sided test instead."
        alternative = :two
    end

    if characteristics.total < 40 || (characteristics.prop_small > 0.2 && characteristics.total < 100)
        if total_cells <= fisher_exact_threshold && size(matrix) == (2, 2)
            test_result = perform_fishers_exact_2x2(matrix, alternative)
            if test_result !== nothing
                println("\nFisher's Exact Test ($(alternative)-sided): p-value = $(round(pvalue(test_result), digits=4))")
            else
                p_value = monte_carlo_fisher(matrix, n_simulations, alternative)
                println("\nFisher's Exact Test (Monte Carlo, $(alternative)-sided): p-value = $(round(p_value, digits=4))")
            end
        else
            p_value = monte_carlo_fisher(matrix, n_simulations, alternative)
            println("\nFisher's Exact Test (Monte Carlo, $(alternative)-sided): p-value = $(round(p_value, digits=4))")
        end
    else
        # Chi-square test is always two-sided
        test_result = ChisqTest(matrix)
        if characteristics.prop_small > 0.2
            if total_cells <= fisher_exact_threshold && size(matrix) == (2, 2)
                fisher_result = perform_fishers_exact_2x2(matrix, alternative)
                println("\nChi-square test: χ² = $(round(test_result.stat, digits=2)), df = $(test_result.df)")
                println("Asymptotic p-value (two-sided): $(round(pvalue(test_result), digits=4))")
                if fisher_result !== nothing
                    println("Fisher's Exact p-value ($(alternative)-sided): $(round(pvalue(fisher_result), digits=4))")
                else
                    mc_pvalue = monte_carlo_fisher(matrix, n_simulations, alternative)
                    println("Monte Carlo p-value ($(alternative)-sided): $(round(mc_pvalue, digits=4))")
                end
            else
                mc_pvalue = monte_carlo_fisher(matrix, n_simulations, alternative)
                println("\nChi-square test: χ² = $(round(test_result.stat, digits=2)), df = $(test_result.df)")
                println("Asymptotic p-value (two-sided): $(round(pvalue(test_result), digits=4))")
                println("Monte Carlo p-value ($(alternative)-sided): $(round(mc_pvalue, digits=4))")
            end
        else
            println("\nChi-square test results (two-sided):")
            println(test_result)
        end
    end
end

"""
    perform_fishers_exact_2x2(matrix::Matrix{T}, alternative::Symbol) where {T<:Number}

Helper function for 2x2 Fisher's exact test with specified alternative hypothesis.
"""
function perform_fishers_exact_2x2(matrix::Matrix{T}, alternative::Symbol) where {T<:Number}
    try
        return FisherExactTest(Int(matrix[1, 1]), Int(matrix[1, 2]),
            Int(matrix[2, 1]), Int(matrix[2, 2]),
            alternative=alternative)
    catch e
        @warn "Failed to perform exact Fisher's test, falling back to Monte Carlo simulation"
        return nothing
    end
end

end # module
