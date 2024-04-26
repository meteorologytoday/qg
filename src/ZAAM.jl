

module ZAAM

    using LinearAlgebra
    using SparseArrays

    include("constants.jl")
    

    include("Grid.jl")
    include("BasicMatrixOperators.jl")
    include("AdvancedMatrixOperators.jl")
    include("PhyParams.jl")

    include("psi_solver.jl")
    include("Env.jl")
    include("Core.jl")
    include("State.jl")
    include("Model.jl")
    

    include("dynamics.jl")
end
