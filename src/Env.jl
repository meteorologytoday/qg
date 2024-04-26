mutable struct Env


    pp :: PhyParams
    gd :: Grid

    function Env(;
        ϕs :: Float64,
        ϕn :: Float64,
        npts :: Integer,
        J :: Float64,
        L :: Union{Float64, AbstractArray{Float64, 1}, Function},
        K :: Float64,
        H :: Float64,
        N :: Float64,
        C :: Float64,
    )

        gd = Grid(
            Ny=npts,
            Ω=Ω,
            ϕn = ϕn,
            ϕs = ϕs,
            H = H,
            R = R,
        )

        if typeof(L) <: Float64

            _L = zeros(Float64, npts+1)
            _L .= L

            L = _L

        elseif typeof(L) <: Function

            L = L(gd)

            if length(L) != npts+1
                throw(ErrorException("Length of L is wrong."))
            end

        end

        

        pp = PhyParams(J, K, N, C, π/H, L)

        return new(
            pp, gd,
        )
        
    end
end
