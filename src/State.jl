mutable struct State

    B   :: AbstractArray{Float64, 1}
    Γ   :: AbstractArray{Float64, 1}
    Ψ   :: AbstractArray{Float64, 1}
    SST :: AbstractArray{Float64, 1}

    function State(
        ev :: Env,
    )
        Ny = ev.gd.Ny
 
        B   = zeros( Float64, Ny ) 
        Γ   = zeros( Float64, Ny) 
        Ψ   = zeros( Float64, Ny+1) 
        SST = zeros( Float64, Ny ) 

        return new(
            B,
            Γ,
            Ψ,
            SST,
        )    
    end
end

