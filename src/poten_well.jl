

using CSV
using DataFrames

const V0 = -20                  # depth of well Unit: eV
const L = 4054.184              # Width of well. Unit: MeV^(-1) ≈ 8A
const up = L/2 + 500                # Upper limit of integration. Unit: MeV^(-1)
const low = -up  

function poten_well!(x)
    V_vec = zeros(length(x) )
    for (i, xx) in enumerate(x)
        if abs(xx) <= L/2
            V_vec[i] = V0
        end
    end
    return V_vec
end

x = low:1:up

df = DataFrame()
df[!, :x] = x
df[!, :V] = poten_well!(x)
CSV.write("../output/poten_well.csv", df)