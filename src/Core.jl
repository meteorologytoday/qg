mutable struct Core

    amo :: AdvancedMatrixOperators
    psi_solver :: StreamfunctionSolver
    ops :: Dict

    function Core(
        ev :: Env,
    )

        amo = AdvancedMatrixOperators(gd=ev.gd)
        psi_solver = StreamfunctionSolver(amo=amo, pp=ev.pp)

        ops = Dict(
            :ydiff    => amo.T_DIVy_V * spdiagm(0=>ev.pp.L) * amo.V_∂y_T,
            #:ydiff    => amo.T_DIVy_V * ev.pp.K * amo.V_∂y_T,
        )

        return new(
            amo,
            psi_solver,
            ops,
        )    
    end
end
