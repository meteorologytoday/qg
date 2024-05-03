struct Grid

    Nx :: Int64
    Ny :: Int64
    Nz :: Int64
    
    Lx   :: Float64    # atmospheric depth
    Ly   :: Float64    # atmospheric depth
    Lz   :: Float64    # atmospheric depth
 
    mask :: AbstractArray{Int64, 3}

    Δx_T  :: AbstractArray{Float64, 3}
    Δy_T  :: AbstractArray{Float64, 3}
    Δz_T  :: AbstractArray{Float64, 3}

    Δx_U  :: AbstractArray{Float64, 3}
    Δy_U  :: AbstractArray{Float64, 3}
    Δz_U  :: AbstractArray{Float64, 3}
 
    Δx_V  :: AbstractArray{Float64, 3}
    Δy_V  :: AbstractArray{Float64, 3}
    Δz_V  :: AbstractArray{Float64, 3}
 
    Δx_W  :: AbstractArray{Float64, 3}
    Δy_W  :: AbstractArray{Float64, 3}
    Δz_W  :: AbstractArray{Float64, 3}

    Δx_UV :: AbstractArray{Float64, 3}
    Δy_UV :: AbstractArray{Float64, 3}
    Δz_UV :: AbstractArray{Float64, 3}
 
    x_T   :: AbstractArray{Float64, 3}
    x_U   :: AbstractArray{Float64, 3}
    x_V   :: AbstractArray{Float64, 3}
    x_W   :: AbstractArray{Float64, 3}
    x_UV  :: AbstractArray{Float64, 3}

    y_T   :: AbstractArray{Float64, 3}
    y_U   :: AbstractArray{Float64, 3}
    y_V   :: AbstractArray{Float64, 3}
    y_W   :: AbstractArray{Float64, 3}
    y_UV  :: AbstractArray{Float64, 3}

    z_T   :: AbstractArray{Float64, 3}
    z_U   :: AbstractArray{Float64, 3}
    z_V   :: AbstractArray{Float64, 3}
    z_W   :: AbstractArray{Float64, 3}
    z_UV  :: AbstractArray{Float64, 3}

    Δv_T  :: AbstractArray{Float64, 3}
    ∫Δv   :: Float64
    
    function Grid(;
        Nx    :: Int64,
        Ny    :: Int64,
        Nz    :: Int64,
        Lx    :: Float64,
        Ly    :: Float64,
        Lz    :: Float64,
    )

        mask = ones(Int64, Nx, Ny, Nz)
        
        _x_T, _x_U, _Δx_T, _Δx_U = genGrid(Nϕ=Nx, ϕs=0, ϕn=Lx, periodic=true )
        _y_T, _y_V, _Δy_T, _Δy_V = genGrid(Nϕ=Ny, ϕs=0, ϕn=Ly, periodic=false)
        _z_T, _z_W, _Δz_T, _Δz_W = genVerticalGrid(Nz=Nz, H=Lz)

        _z_U  = copy(  _z_T )
        _Δz_U = copy( _Δz_T )
        
        _y_U  = copy(  _y_T )
        _Δy_U = copy( _Δy_T )


        z_makeMesh = (a, nx, ny) -> repeat( reshape(a, 1, 1, :), outer=(nx, ny,  1) )
        y_makeMesh = (a, nx, nz) -> repeat( reshape(a, 1, :, 1), outer=(nx, 1,  nz) )
        x_makeMesh = (a, ny, nz) -> repeat( reshape(a, :, 1, 1), outer=(1,  ny, nz) )

        z_T   = z_makeMesh(_z_T,  Nx, Ny)
        z_U   = z_makeMesh(_z_T,  Nx, Ny)
        z_V   = z_makeMesh(_z_T,  Nx, Ny+1)
        z_W   = z_makeMesh(_z_W,  Nx, Ny)
        z_UV  = z_makeMesh(_z_W,  Nx, Ny+1)
    
        Δz_T  = z_makeMesh(_Δz_T, Nx,   Ny)
        Δz_U  = z_makeMesh(_Δz_U, Nx,   Ny)
        Δz_V  = z_makeMesh(_Δz_T, Nx,   Ny+1)
        Δz_W  = z_makeMesh(_Δz_W, Nx,   Ny)
        Δz_UV = z_makeMesh(_Δz_W, Nx,   Ny+1)

        y_T   = y_makeMesh(_y_T,  Nx,   Nz  )
        y_U   = y_makeMesh(_y_T,  Nx,   Nz  )
        y_V   = y_makeMesh(_y_V,  Nx,   Nz  )
        y_W   = y_makeMesh(_y_T,  Nx,   Nz+1)
        y_UV  = y_makeMesh(_y_V,  Nx,   Nz+1)

        Δy_T   = y_makeMesh(_Δy_T, Nx,   Nz  )
        Δy_U   = y_makeMesh(_Δy_T, Nx,   Nz  )
        Δy_V   = y_makeMesh(_Δy_V, Nx,   Nz  )
        Δy_W   = y_makeMesh(_Δy_T, Nx,   Nz+1)
        Δy_UV  = y_makeMesh(_Δy_V, Nx,   Nz+1)

        x_T   = x_makeMesh(_x_T,  Ny,   Nz  )
        x_U   = x_makeMesh(_x_U,  Ny,   Nz  )
        x_V   = x_makeMesh(_x_T,  Ny+1, Nz  )
        x_W   = x_makeMesh(_x_T,  Ny,   Nz+1)
        x_UV  = x_makeMesh(_x_U,  Ny+1, Nz+1)

        Δx_T   = x_makeMesh(_Δx_T, Ny,  Nz  )
        Δx_U   = x_makeMesh(_Δx_U, Ny,  Nz  )
        Δx_V   = x_makeMesh(_Δx_T, Ny+1,Nz  )
        Δx_W   = x_makeMesh(_Δx_T, Ny,  Nz+1)
        Δx_UV  = x_makeMesh(_Δx_U, Ny+1,Nz+1)


        Δv_T = Δx_T .* Δy_T .* Δz_T
        ∫Δv = sum(Δv_T)

        
        return new(

            Nx,
            Ny,
            Nz,

            Lx,
            Ly,
            Lz,
         
            mask,

            Δx_T,
            Δy_T,
            Δz_T,

            Δx_U,
            Δy_U,
            Δz_U,

            Δx_V,
            Δy_V,
            Δz_V,
 
            Δx_W,
            Δy_W,
            Δz_W,
        
            Δx_UV,
            Δy_UV,
            Δz_UV,

            x_T,
            x_U,
            x_V,
            x_W,
            x_UV,

            y_T,
            y_U,
            y_V,
            y_W,
            y_UV,

            z_T,
            z_U,
            z_V,
            z_W,
            z_UV,

            Δv_T,
            ∫Δv,

        ) 
        
    end
end


function genGrid(;
    Nϕ :: Int64,
    ϕs :: Float64,
    ϕn :: Float64,
    periodic :: Bool = false,
)

    ϕ_V = collect(Float64, range(ϕs, ϕn, length=Nϕ+1))
    ϕ_T = (ϕ_V[1:end-1] + ϕ_V[2:end]) / 2.0

    Δϕ_T = similar(ϕ_T)
    Δϕ_V = similar(ϕ_V)

    δϕ = ϕ_V[2] - ϕ_V[1]

    Δϕ_T .= δϕ
    Δϕ_V .= δϕ
   
    if periodic
        ϕ_V = ϕ_V[1:end-1]
        Δϕ_V = Δϕ_V[1:end-1]
    end
 
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
