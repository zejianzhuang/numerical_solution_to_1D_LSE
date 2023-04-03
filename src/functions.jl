
const me = 0.511      # Electronic mass. Unit: MeV


function square_poten_well(x::Vector, N::Int, L, V0)
    mat_V = zeros(N, N)
    for (i, xx) in enumerate(x)
        if abs(xx) <= L/2
            mat_V[i, i] = V0 * 1e-6 # Uhit of V is eV. Here we convert it to MeV.
        end
    end
    return mat_V
end



function Green_func(k, x::Vector, xp::Vector, N::Int)
    G = ones(ComplexF64, N, N)
    for i in 1:N
        G[i, :] = [-1.0im*me / k * exp(1.0im*k*abs(x[i]-xp[j]) ) for j in 1:N]
    
    end
    return G
end

φ(k, x::Vector, N::Int) = [exp(1.0im*k*x[i]) for i in 1:N]

function change_of_var(node, weight, a, b, N)
    nop = [(b-a) * node[i] / 2.0 + (a+b) / 2.0 for i in 1:N]
    wp = [(b-a) / 2.0 * weight[i] for i in 1:N]
    return nop, wp
end

ki(E) = sqrt(2*me*E*1e-6 + 0.0im) # Unit of E: eV

