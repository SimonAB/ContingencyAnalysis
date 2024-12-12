Pkg.develop(path="/Users/s_a_b/Documents/Work/Coding/ContingencyAnalysis")
using ContingencyAnalysis
using DataFrames

# Example 1: Tea Tasting Test
# A British woman claimed to be able to distinguish whether milk or tea was added first
# Testing null hypothesis of no association between true order and guess
tea_matrix = [3 1; 1 3]
println("\nAnalysing Tea Tasting Matrix:")
analyse_contingency(tea_matrix)

# Example 2: Criminal convictions of like-sex twins
convictions = [2 10; 15 3]
println("\nAnalysing Twin Convictions Matrix:")
analyse_contingency(convictions)

# Example 3: Job Satisfaction by Income
job_matrix = [
    1 2 1 0;
    3 3 6 1;
    10 10 14 9;
    6 7 12 11
]
println("\nAnalysing Job Satisfaction Matrix:")
analyse_contingency(job_matrix)

# Example 4: Mehta & Patel example
mp6_matrix = [
    1 2 2 1 1 0 1;
    2 0 0 2 3 0 0;
    0 1 1 1 2 7 3;
    1 1 2 0 0 0 1;
    0 1 1 1 1 0 0
]
println("\nAnalysing MP6 Matrix:")
analyse_contingency(mp6_matrix)
