



using CSV
using DataFrames
using LinearAlgebra
using NLsolve


include("./src/functions.jl")
include("./src/xw.jl")

const N = 298
const V0 = -20                  # depth of well Unit: eV
const L = 4054.184              # Width of well. Unit: MeV^(-1) ≈ 8A
const up = L / 2 + 500          # Upper limit of integration. Unit: MeV^(-1)
const low = -up                 # Lower limit of integration. Unit: MeV^(-1)



node, weight = xw[N-5+1]

x, wp = change_of_var(node, weight, low, up, N)
xp = x
W = diagm(wp)
V = square_poten_well(x, N, L, V0)

# 🌵🌵🌵🌵🌵🌵🌵🌵🌵 Scattering 🌵🌵🌵🌵🌵🌵🌵🌵🌵
function investigate_pole_position!(E; two_dims=false, dim=nothing, file_path::String="./output/invest_pole_position_to_gs.csv")
    eps = 1e-20
    df = DataFrame()
    if two_dims == false
        redet = []
        imdet = []
        for w in E
            G = Green_func(ki(w), x, xp, N)
            D = diagm(ones(N)) - G * V * W
            deter = det(D)
            if abs(real(deter)) <= eps && abs(imag(deter)) <= eps
                append!(redet, 0.0)
                append!(imdet, 0.0)

            elseif abs(imag(deter)) <= eps
                append!(redet, real(deter))
                append!(imdet, 0.0)
            else
                append!(redet, real(deter))
                append!(imdet, imag(deter))
            end
        end
        df[!, "E"] = E
        df[!, "redet"] = redet
        df[!, "imdet"] = imdet
        ###
    elseif two_dims == true && dim !== nothing
        reE = []
        imE = []
        inv_abs_det = []
        for rew in E[1]
            append!(reE, ones(dim).*rew)
            append!(imE, E[2] )
            append!(inv_abs_det, [1/(abs(det(diagm(ones(N)) - Green_func(ki(rew+imw*1.0im), x, xp, N) * V * W) ) + eps) for imw in E[2] ] )
            
        end
        
        df[!, "reE"] = reE
        df[!, "imE"] = imE
        df[!, "inv_abs_det"] = inv_abs_det
    end
    CSV.write(file_path, df)
end
# 🌠🌠🌠🌠 Investigate pole to bound state 🌠🌠🌠🌠
#E = -20.0:0.01:-3.0 # eV
#investigate_pole_position!(E)


# 🌠🌠🌠🌠 Investigate pole to resonance state 🌠🌠🌠🌠
#dim = 100
#E = [range(0.01, 70, dim), range(-15.0, 0.0, dim)]
#investigate_pole_position!(E, two_dims=true, dim=dim, file_path="output/invest_pole_position_to_resonance.csv")

# 🌠🌠🌠🌠 Extract the poles for bound states 🌠🌠🌠🌠
init_guess_gs = [[-3.5], [-8.5], [-12.0], [-15.5], [-18.0], [-19.8]] # Getting these from the figure

det_Tmat(E) = [real(det(diagm(ones(N)) - Green_func(ki(E[1]), x, xp, N) * V * W))]

function pole_position_gs!(init_gauss)

    s = []
    for guess in init_gauss
        r = nlsolve(det_Tmat, guess)
        if converged(r) == true
            append!(s, round(r.zero[1], digits=2) )
        end
    end
    df = DataFrame(gs=["\$$(s[i])\$" for i in 1:6])
    #println(df)
    CSV.write("./to_latex_tab/ground.csv", df)
end
#pole_position_gs!(init_guess_gs)

# 🌠🌠🌠🌠 Extract the poles for resonance 🌠🌠🌠🌠
init_guess_reson = [[64.0, -12.0], [51.0, -10.9], [38.0, -8.5], [25.1, -7.0], [16.0, -5.0], [9.0, -3.0], [1.01, -0.6] ] # Getting these from figure of 1 / |det(D)|

function det_Tmat1(E)
    w = E[1] + E[2]*1im
    deter = det(diagm(ones(N)) - Green_func(ki(w), x, xp, N) * V * W)

    return [real(deter), imag(deter) ]
end

function pole_position_resonance(init_guess)
    res = []
    ims = []
    for guess in init_guess
        r = nlsolve(det_Tmat1, guess)
        if converged(r) == true
            append!(res, round(r.zero[1], digits=2 ) )
            append!(ims, round(r.zero[2], digits=2 ) )
        end
    end
    df = DataFrame(Resonance=["\$$(res[i])$(ims[i])i\$" for i in 1:7] )
    #println(df)
    CSV.write("./to_latex_tab/pole_of_resonace.csv", df)
end
#pole_position_resonance(init_guess_reson)























println("💖💖💖💖💖💖💖💖💖 This programing has run successfully!!! 💖💖💖💖💖💖💖💖💖")







