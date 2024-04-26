struct Grid

    Ω   :: Float64
    R   :: Float64    # Radius of planet
    H   :: Float64    # Atmospheric

    ϕn :: Float64
    ϕs :: Float64

    Nx :: Int64
    Ny :: Int64
    Nz :: Int64
 
    mask :: AbstractArray{Int64, 3}

    Δx_T  :: AbstractArray{Float64, 3}
    Δy_T  :: AbstractArray{Float64, 3}
    Δz_T  :: AbstractArray{Float64, 3}

    Δy_U  :: AbstractArray{Float64, 3}
    Δz_U  :: AbstractArray{Float64, 3}
 
    Δx_V  :: AbstractArray{Float64, 3}
    Δy_V  :: AbstractArray{Float64, 3}
    Δz_V  :: AbstractArray{Float64, 3}
 
    Δx_W  :: AbstractArray{Float64, 3}
    Δy_W  :: AbstractArray{Float64, 3}
    Δz_W  :: AbstractArray{Float64, 3}

    Δx_VW :: AbstractArray{Float64, 3}
    Δy_VW :: AbstractArray{Float64, 3}
    Δz_VW :: AbstractArray{Float64, 3}
   
    ϕ_T   :: AbstractArray{Float64, 3}
    ϕ_U   :: AbstractArray{Float64, 3}
    ϕ_V   :: AbstractArray{Float64, 3}
    ϕ_W   :: AbstractArray{Float64, 3}
    ϕ_VW   :: AbstractArray{Float64, 3}

    z_T   :: AbstractArray{Float64, 3}
    z_U   :: AbstractArray{Float64, 3}
    z_V   :: AbstractArray{Float64, 3}
    z_W   :: AbstractArray{Float64, 3}
    z_VW   :: AbstractArray{Float64, 3}

    Δv_T  :: AbstractArray{Float64, 3}
    ∫Δv   :: Float64
    
    Δλ_W  :: AbstractArray{Float64, 3}

    function Grid(;
        Ny    :: Int64,
        Ω     :: Float64,
        ϕn    :: Float64, 
        ϕs    :: Float64,
        H     :: Float64,
        R     :: Float64,
    )

        Nx = 1
        Nz = 1

        mask = ones(Int64, Nx, Ny, Nz)
        
        _z_T, _z_W, _Δz_T, _Δz_W = genVerticalGrid(Nz=Nz, H=H)
        _ϕ_T, _ϕ_V, _Δϕ_T, _Δϕ_V = genHorizontalGrid(Nϕ=Ny, ϕs=ϕs, ϕn=ϕn)

        _z_U  = copy(_z_T)
        _Δz_U = copy(_Δz_T)
        _ϕ_U  = copy(_ϕ_T)
        _Δϕ_U = copy(_Δϕ_T)

        z_makeMesh = (a, nx, ny) -> repeat( reshape(a, 1, 1, :), outer=(nx, ny,  1) )
        y_makeMesh = (a, nx, nz) -> repeat( reshape(a, 1, :, 1), outer=(nx, 1,  nz) )
        x_makeMesh = (a, ny, nz) -> repeat( reshape(a, :, 1, 1), outer=(1,  ny, nz) )

        z_T   = z_makeMesh(_z_T,  Nx,   Ny)
        z_U   = z_makeMesh(_z_T,  Nx+1, Ny)
        z_V   = z_makeMesh(_z_T,  Nx,   Ny+1)
        z_W   = z_makeMesh(_z_W,  Nx,   Ny)
        z_VW  = z_makeMesh(_z_W,  Nx,   Ny+1)
    
        Δz_T  = z_makeMesh(_Δz_T, Nx,   Ny)
        Δz_U  = z_makeMesh(_Δz_U, Nx+1, Ny)
        Δz_V  = z_makeMesh(_Δz_T, Nx,   Ny+1)
        Δz_W  = z_makeMesh(_Δz_W, Nx,   Ny)
        Δz_VW = z_makeMesh(_Δz_W, Nx,   Ny+1)


        ϕ_T   = y_makeMesh(_ϕ_T,  Nx,   Nz  )
        ϕ_U   = y_makeMesh(_ϕ_U,  Nx+1, Nz  )
        ϕ_V   = y_makeMesh(_ϕ_V,  Nx,   Nz  )
        ϕ_W   = y_makeMesh(_ϕ_T,  Nx,   Nz+1)
        ϕ_VW  = y_makeMesh(_ϕ_V,  Nx,   Nz+1)

        Δϕ_T   = y_makeMesh(_Δϕ_T, Nx,   Nz  )
        Δϕ_U   = y_makeMesh(_Δϕ_T, Nx+1, Nz  )
        Δϕ_V   = y_makeMesh(_Δϕ_V, Nx,   Nz  )
        Δϕ_W   = y_makeMesh(_Δϕ_T, Nx,   Nz+1)
        Δϕ_VW  = y_makeMesh(_Δϕ_V, Nx,   Nz+1)

        Δλ = [2π, ]
        Δλ_T   = x_makeMesh(Δλ, Ny,   Nz)
        Δλ_V   = x_makeMesh(Δλ, Ny+1, Nz)
        Δλ_W   = x_makeMesh(Δλ, Ny,   Nz+1)
        Δλ_VW  = x_makeMesh(Δλ, Ny+1, Nz+1)

        # Uncomment the following line to remove curvature effect
        #cos(x) = 1.0

        # Calculating horizontal grid edges
        Δx_T = R * cos.(ϕ_T) .* Δλ_T;
        Δy_T = R * Δϕ_T;

        Δy_U = R * Δϕ_U;

        Δx_V = R * cos.(ϕ_V) .* Δλ_V;
        Δy_V = R * Δϕ_V;
 
        Δx_W = R * cos.(ϕ_W) .* Δλ_W;
        Δy_W = R * Δϕ_W;

        Δx_VW = R * cos.(ϕ_VW) .* Δλ_VW;
        Δy_VW = R * Δϕ_VW;


        Δv_T = Δx_T .* Δy_T .* Δz_T
        ∫Δv = sum(Δv_T)

        
        return new(
            
            Ω,
            R,
            H,

            ϕn,
            ϕs,

            Nx,
            Ny,
            Nz,
         
            mask,

            Δx_T,
            Δy_T,
            Δz_T,

            Δy_U,
            Δz_U,

            Δx_V,
            Δy_V,
            Δz_V,
 
            Δx_W,
            Δy_W,
            Δz_W,
        
            Δx_VW,
            Δy_VW,
            Δz_VW,

            ϕ_T,
            ϕ_U,
            ϕ_V,
            ϕ_W,
            ϕ_VW,

            z_T,
            z_U,
            z_V,
            z_W,
            z_VW,

            Δv_T,
            ∫Δv,

            Δλ_W,
        ) 
        
    end
end


function genHorizontalGrid(;
    Nϕ :: Int64,
    ϕs :: Float64,
    ϕn :: Float64,
)

    ϕ_V = collect(Float64, range(ϕs, ϕn, length=Nϕ+1))
    ϕ_T = (ϕ_V[1:end-1] + ϕ_V[2:end]) / 2.0

    Δϕ_T = similar(ϕ_T)
    Δϕ_V = similar(ϕ_V)

    δϕ = ϕ_V[2] - ϕ_V[1]

    Δϕ_T .= δϕ
    Δϕ_V .= δϕ
    
    return ϕ_T, ϕ_V, Δϕ_T, Δϕ_V

end


function genVerticalGrid(;
    Nz    :: Int64,
    H     :: Float64,
)

    η = collect(Float64, range(0, 1, length=Nz+1))
    z_W  = H * η

    Δz_W = similar(z_W)

    Δz_T = z_W[1:end-1] - z_W[2:end]
    Δz_W[2:end-1] = ( Δz_T[2:end] + Δz_T[1:end-1] ) / 2.0
    Δz_W[1] = Δz_W[2]
    Δz_W[end] = Δz_W[end-1]

    z_T = (z_W[1:end-1] + z_W[2:end]) / 2.0

    return z_T, z_W, Δz_T, Δz_W
end
