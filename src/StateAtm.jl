mutable struct StateAtm

    h   :: AbstractArray{Float64, 3}
    u   :: AbstractArray{Float64, 3}
    v   :: AbstractArray{Float64, 3}
    q   :: AbstractArray{Float64, 3}
    zb  :: AbstractArray{Float64, 2}
    
    function StateAtm(
        ev :: Env,
    )
        Nx = ev.gd.Nx
        Ny = ev.gd.Ny
        Nz = ev.gd.Nz
         
        h = zeros( Float64, Nx, Ny,   Nz) 
        u = zeros( Float64, Nx, Ny,   Nz) 
        v = zeros( Float64, Nx, Ny+1, Nz) 
        q = zeros( Float64, Nx, Ny,   Nz) 
        zb = zeros( Float64, Nx, Ny,) 

        return new(
            h, u, v, q, zb,
        )    
    end
end

