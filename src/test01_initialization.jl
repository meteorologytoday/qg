include("ZAAM.jl")




m = ZAAM.Model(ZAAM.Env(
    ϕn =  deg2rad(90),
    ϕs =  deg2rad(-90),
    npts = 51,
    J = 1e3,
    L = 1e5 * 0,
    K = 1e5,
    H = 10e3, 
    N = 1e-6 * 0,
    C = ZAAM.ρ0 * ZAAM.c_p * 1e-3 * 10.0,
))


m.st.SST .= ZAAM.θ0 .+ exp.( - (m.ev.gd.ϕ_T[:] / deg2rad(10)).^2 / 2 ) * 30
#m.st.B[5] = ZAAM.g0

dt = 86400.0
steps = 300

save_B   = zeros(Float64, length(m.st.B),   steps)
save_SST = zeros(Float64, length(m.st.SST), steps)
save_TS  = zeros(Float64, length(m.st.B), steps)


for t=1:steps
   
    save_B[:, t]   = m.st.B
    save_SST[:, t] = m.st.SST
    save_TS[:, t]  = ZAAM.computeTS(m, m.st.B)
    
    println("Step Model: ", t) 

    ZAAM.step!(m, dt)


end


println("Try to output data...")

println("Loading libraries...")
using YAXArrays, NetCDF
using Dates
println("Done")

axlist = [
    RangeAxis("y", collect(1:length(m.ev.gd.ϕ_T))),
    RangeAxis("time", collect(1:steps)),
]

da_B   = YAXArray(axlist, save_B)
da_SST = YAXArray(axlist, save_SST)
da_TS  = YAXArray(axlist, save_TS)

ds = Dataset(B=da_B, SST=da_SST, TS=da_TS)


savedataset(ds, path="test.nc", driver=:netcdf, overwrite=true)
