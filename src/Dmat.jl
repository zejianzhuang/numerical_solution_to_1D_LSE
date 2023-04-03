using CSV
using DataFrames
using LinearAlgebra
using NLsolve


include("./functions.jl")
include("./xw.jl")



const N = 298
const V0 = -20                  # depth of well Unit: eV
const L = 4054.184              # Width of well. Unit: MeV^(-1) ≈ 8A
const up = L / 2 + 500          # Upper limit of integration. Unit: MeV^(-1)
const low = -up                 # Lower limit of integration. Unit: MeV^(-1)
node, weight = xw[N-5+1]

function variable()
    v = Dict(
    :node => node,
    :weight => weight,
    :low => low,
    :up => up,
    :N => N,
    :V0 => V0,
    :L => L)
   #if new_var != missing
    #    return merge(v, new_var)
   #else
        return v
   #end
end


function Dmat(E, var::Dict; return_det=true)
    #=
    The syntax for vars must be look like as follows:
        var = Dict(
            :node => Vector,
            :weight => Vector,
            :low => low,
            :up => up,
            :N => N,
            :V0 => V0,
            :L => L,
            :En => (E1, E2) # E1: real part of En, E2: imaginary part of En
            )
    =#
    eps = 1e-18
    x, wp = change_of_var(var[:node], var[:weight], var[:low], var[:up], var[:N])
    W = diagm(wp)
    V = square_poten_well(x, var[:N], var[:L], var[:V0])
    G = Green_func(ki(E), x, x, var[:N] )
    D = diagm(ones(N) ) - G * V * W

    if return_det == true
        deter = det(D)
        #detvec = [] # Used for storing the det. Split det into real part and imag part.
        if abs(real(deter) ) <= eps && abs(imag(deter) ) <= eps
            #append!(detvec, [0.0, 0.0])
            detvec = 0.0, 0.0
        elseif abs(imag(deter) ) <= eps
            #append!(detvec, [real(deter), 0.0])
            detvec = real(deter), 0.0
        else
            #append!(detvec, [real(deter), imag(deter) ] )
            detvec = real(deter), imag(deter)
        end
        return D, detvec
    else
        return D
    end
end

function pole_position(init_guess; file_path=nothing)
    tol = 1e-18
    res = []
    ims = []
    for init in init_guess
        r = nlsolve(E -> Dmat(E, variable())[2], init, ftol=tol)
        if converged(r) == true
            append!(res, r.zero[1])
            append!(ims, r.zero[2])
        end
    end
    df = DataFrame(repole=res, impole=ims)
    if file_path !== nothing
        CSV.write(file_path, df)
    else
        println(df)
    end
end;
# 🌠🌠🌠🌠 Investigate pole to bound state 🌠🌠🌠🌠
#init_guess = [[-3.54, 0.0], [-8.33, 0.0], [-12.52, 0.0], [-15.8, 0.0], [-18.14, 0.0], [-19.61, 0.0] ]

#pole_position(init_guess)
