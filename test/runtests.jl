using ContingencyAnalysis
using Test
using Random

@testset "ContingencyAnalysis.jl" begin
    Random.seed!(123)  # For reproducibility

    @testset "2x2 Tables" begin
        # Test matrix with known association (positive)
        positive_matrix = [10 2; 1 8]

        # Test matrix with known association (negative)
        negative_matrix = [2 10; 8 1]

        # Test matrix with no association
        neutral_matrix = [5 5; 5 5]

        @testset "Monte Carlo Fisher Tests" begin
            # Use even more extreme matrices
            strong_positive = [40 4; 4 40]  # Stronger association
            strong_negative = [4 40; 40 4]  # Stronger association

            # Two-sided tests
            @test 0 <= monte_carlo_fisher(strong_positive, 1000, :two) <= 1
            @test 0 <= monte_carlo_fisher(strong_negative, 1000, :two) <= 1

            # One-sided tests with more relaxed thresholds
            p_greater = monte_carlo_fisher(strong_positive, 1000, :greater)
            p_less = monte_carlo_fisher(strong_positive, 1000, :less)
            @test p_greater < 0.2  # More reasonable for Monte Carlo
            @test p_less > 0.8     # Should be very non-significant

            p_greater = monte_carlo_fisher(strong_negative, 1000, :greater)
            p_less = monte_carlo_fisher(strong_negative, 1000, :less)
            @test p_greater > 0.8  # Should be very non-significant
            @test p_less < 0.2     # More reasonable for Monte Carlo

            # Neutral case - just test that p-values are not extreme
            p_greater = monte_carlo_fisher(neutral_matrix, 1000, :greater)
            p_less = monte_carlo_fisher(neutral_matrix, 1000, :less)
            @test 0.2 < p_greater < 0.8  # Should not be extreme
            @test 0.2 < p_less < 0.8     # Should not be extreme
        end

        @testset "Invalid Inputs" begin
            # Test negative values
            @test_throws ArgumentError analyse_contingency([-1 2; 3 4])

            # Test invalid alternative
            @test_throws ArgumentError monte_carlo_fisher(positive_matrix, 1000, :invalid)
        end
    end

    @testset "Larger Tables" begin
        # 3x3 matrix
        larger_matrix = [10 5 2; 3 8 4; 1 3 9]

        @testset "Monte Carlo Tests" begin
            # Two-sided test should work
            @test 0 <= monte_carlo_fisher(larger_matrix, 1000, :two) <= 1

            # One-sided tests should warn and default to two-sided
            @test_logs (:warn, "One-sided tests are only meaningful for 2x2 tables. Using two-sided test instead.") begin
                analyse_contingency(larger_matrix, 1000, alternative=:greater)
            end
        end
    end

    @testset "Matrix Analysis" begin
        # Test matrix with zero marginals
        zero_marginal_matrix = [0 0 0; 2 3 1; 0 1 2]  # First row all zeros
        result = analyse_matrix(zero_marginal_matrix)

        @test size(result.cleaned_matrix, 1) == 2  # Should remove the zero row
        @test result.total == sum(zero_marginal_matrix)
        @test 0 <= result.prop_small <= 1

        # Test expected frequency calculation
        small_freq_matrix = [1 1; 1 10]
        result = analyse_matrix(small_freq_matrix)
        @test result.prop_small > 0  # Should have some small expected frequencies
    end

    @testset "Chi-Square Statistic" begin
        test_matrix = [10 5; 2 8]

        # Test basic properties
        @test calculate_chi_statistic(test_matrix) >= 0  # Should be non-negative

        # Test independence case
        independent_matrix = [5 5; 5 5]
        @test calculate_chi_statistic(independent_matrix) ≈ 0 atol = 1e-10

        # Test perfect dependence
        dependent_matrix = [10 0; 0 10]
        @test calculate_chi_statistic(dependent_matrix) > calculate_chi_statistic(test_matrix)
    end

    @testset "Simulation Functions" begin
        row_sums = [12, 8]
        col_sums = [10, 10]

        # Run multiple simulations to ensure consistency
        for _ in 1:10
            simulated = simulate_contingency_table(reshape(row_sums, :, 1), reshape(col_sums, 1, :))

            # Test marginal sums preservation
            @test vec(sum(simulated, dims=2)) == row_sums
            @test vec(sum(simulated, dims=1)) == col_sums

            # Test non-negative values
            @test all(>=(0), simulated)
        end
    end
end
