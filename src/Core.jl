mutable struct Core

    amo :: AdvancedMatrixOperators
    ops :: Dict

    function Core(
        ev :: Env,
    )

        amo = AdvancedMatrixOperators(gd=ev.gd)

        ops = Dict(
        )

        return new(
            amo,
            ops,
        )    
    end
end
