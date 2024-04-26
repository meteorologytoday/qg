function step!(
    m :: Model,
    dt :: Float64,
)

    st = m.st

    B = st.B

    diagΨ!(m)    

    A_DIVΨ, F_DIVΨ   = getHeatTransportMatrix(m)
    A_F_TOA, F_F_TOA = getF_TOAMatrix(m)
    A_F_AO,  F_F_AO  = getF_AOMatrix(m)
    A_YDiff, F_YDiff = getYDiffMatrix(m)


    A = A_F_TOA + A_F_AO + A_YDiff + A_DIVΨ
    F = F_F_TOA + F_F_AO + F_YDiff + F_DIVΨ

    #A = A_F_AO# + A_YDiff
    #F = F_F_AO# + F_YDiff



    B[:] = ( I - dt * A ) \ ( B + dt * F )
    
end

function diagΨ!(
    m :: Model,
)

    psi_solver = m.co.psi_solver
    m.st.Γ[:] = psi_solver.T_solveΓ_fromB_T * m.st.B
    m.st.Ψ[:] = psi_solver.V_solveΨ_fromB_T * m.st.B

end

function getHeatTransportMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)
    
    A = - m.ev.pp.N * m.co.amo.T_DIVy_V * m.co.psi_solver.V_solveΨ_fromB_T
    F = m.st.B * 0

    jacobian = A
    if return_jacobian

        return A, F, jacobian
    else

        return A, F

    end

end


function getYDiffMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)

    A = m.co.ops[:ydiff]
    F = m.st.B * 0

    jacobian = m.co.ops[:ydiff]


    if return_jacobian

        return A, F, jacobian
    else

        return A, F

    end

end

function computeTS(
    m :: Model,
    B :: AbstractArray,
)

    pp = m.ev.pp
    gd = m.ev.gd


    return θ0 * (1 .+ (B .- (pp.N * gd.H) / 2) / g0)

end



function computeT_TOA(
    m :: Model,
    B :: AbstractArray,
)

    pp = m.ev.pp
    gd = m.ev.gd


    return θ0 * (1 .+ (B .+ (pp.N * gd.H) / 2) / g0)

end


function getF_TOAMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)

    ev = m.ev
    pp = ev.pp
    gd = ev.gd

    Npts = gd.Ny

    factor = g0 / ( θ0 * c_p * ρ0 * gd.H ) * (gd.H * pp.k / 2) 

    T_TOA = computeT_TOA(m, m.st.B)
    tmp = σ_boltz * T_TOA.^3
    A   = spdiagm(0 => (- factor * tmp * θ0 / g0))
    F   = - factor * tmp * ( θ0 * ( 1 + pp.N * gd.H / (2 * g0)) )  

    jacobian_A = A
    jacobian_F = factor * spdiagm( 0 => - 3 * σ_boltz * T_TOA.^2 * ( θ0 * ( 1 + pp.N * gd.H / (2 * g0) ) ) * (θ0/g0) )


    if return_jacobian
        
        return A, F, jacobian_A + jacobian_F
        
    else
        return A, F
    end
end

function getF_AOMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)

    ev = m.ev
    pp = ev.pp
    gd = ev.gd
    st = m.st

    Npts = gd.Ny
    
    factor = g0 / ( θ0 * c_p * ρ0 * gd.H ) * (gd.H * pp.k / 2) 

    # F_ao
    A  = - factor * pp.C * θ0 / g0 * sparse(I, Npts, Npts)
    F = factor * ( st.SST .- θ0 * ( 1 - pp.N * gd.H / (2 * g0)) ) * pp.C 
    #println(Matrix(A))
    jacobian_A = A
    jacobian_F = 0 * A


    if return_jacobian
        
        return A, F, jacobian_A + jacobian_F
        
    else
        return A, F
    end
    
end
