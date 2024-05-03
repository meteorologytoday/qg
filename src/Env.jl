mutable struct Env


    pp :: PhyParams
    gd :: Grid

    function Env(;
        Nz :: Int64,
        Nx :: Int64,
        Ny :: Int64,
        f0 :: Float64,
        β  :: AbstractArray{Float64, 1},
    )

        gd = Grid(
            Nx=Nx,
            Ny=Ny,
            Nz=Nz,
        )

        pp = PhyParams(f0, β)

        return new(
            pp, gd,
        )
        
    end
end
