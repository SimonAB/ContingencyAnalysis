using ContingencyAnalysis
using RDatasets, DataFrames

# Load Iris dataset
iris = dataset("datasets", "iris");

# Add a new column 'size' based on the condition
iris.size = ifelse.(iris.SepalLength .< 6, "small", "large");

# Create a contingency table
contingency_df = unstack(combine(groupby(iris, [:Species, :size]), nrow => :count),
    :Species, :size, :count, fill=0)

# Convert to matrix, dropping the Species column names
iris_matrix = Matrix{Int64}(contingency_df[:, 2:end])

# Analyse the contingency table
println("\nAnalysing Iris Contingency Table:")
analyse_contingency(iris_matrix)
