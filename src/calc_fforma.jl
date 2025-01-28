function calc_fforma(ξ, elem, deriv = true) # funções de forma gerais para m elementos
    ξs = elem.ξs
    tipo = length(ξs)
    N = ones(typeof(ξ), tipo) #vetor linha das funções de forma N_1, N_2, ..., N_m
    cte = zeros(typeof(ξ), tipo)
    for i = 1:tipo
        for j = 1:tipo  # m(numero de nós) = tipo
            if j != i
                # @infiltrate
                if (ξ - ξs[j]) == 0
                    ξ = ξ + 1e-9
                end
                N[i] = @views N[i] * (ξ - ξs[j]) / (ξs[i] - ξs[j]) #Eq. interpolador de lagrange
                cte[i] += @views 1 / (ξ - ξs[j])
            end
        end
    end
    if deriv
        N, cte .* N
    else
        N
    end
end

function calc_fforma_gen(ξ, ξs, deriv = true) # funções de forma gerais para m elementos
    tipo = length(ξs)
    N = ones(typeof(ξ), tipo) #vetor linha das funções de forma N_1, N_2, ..., N_m
    cte = zeros(typeof(ξ), tipo)
    for i = 1:tipo
        for j = 1:tipo  # m(numero de nós) = tipo
            if j != i
                # @infiltrate
                N[i] = N[i] * (ξ - ξs[j]) / (ξs[i] - ξs[j]) #Eq. interpolador de lagrange
                cte[i] += 1 / (ξ - ξs[j])
            end
        end
    end
    if deriv
        return N, cte .* N
    else
        return N
    end
end

function calc_dfforma_gen(qsi, qsis) # funções de forma gerais para m elementos
    ForwardDiff.derivative(x -> calc_fforma_gen(x, qsis), qsi)
end

function Nlinear(qsi)
    # Calcula as funções de forma lineares contínuas N1 e N2
    N1 = 1 / 2 * (1.0 .- qsi) # Função de forma N1 => linear contínua
    N2 = 1 / 2 * (1.0 .+ qsi) # Função de forma N2 => linear contínua
    return N1, N2
end

"pg ponto interpolado
x ponto interpolador"
function lagrange(pg, x, n)
    ni = length(pg)
    L = ones(ni, n)
    for j = 1:n
        for i = 1:n
            if (i != j)
                L[:, j] = L[:, j] .* (pg .- x[i]) / (x[j] - x[i])
            end
        end
    end
    return L
end

function lagrange(pg, x1, n1, x2, n2)
    l1 = lagrange(pg[:, 1], x1, n1)
    l2 = lagrange(pg[:, 2], x2, n2)


    ni = size(pg, 1)
    L = zeros(ni, n1 * n2)
    for i = 1:ni
        L[i, :] = (l1[i, :]*l2[i, :]')[:]
    end
    L
end
function lagrange(pg, x1, n1, x2, n2, x3, n3)
    l1 = lagrange(pg[:, 1:2], x1, n1, x2, n2)
    l2 = lagrange(pg[:, 3], x3, n3)


    ni = size(pg, 1)
    L = zeros(ni, n1 * n2 * n3)
    for i = 1:ni
        L[i, :] = (l1[i, :]*l2[i, :]')[:]
    end
    L
end
function criapontosinterp(n, t = 1)
    if t == 1
        x = cos.((2 * (1:n) .- 1) * pi / 2 / n)
    elseif t == 2
        x = cos.((1:n) * pi / (n + 1))
    elseif t == 3
        x = cos.(pi * (0:n-1) / (n - 1))
    end
    x
end


function recompressão(L, H)
    w1, r1 = qr(L)
    w1 = Array(w1)
    w2, r2 = qr(H')
    w2 = Array(w2)

    U, S, V = svd(r1 * r2')
    ind1 = S .> ϵ * S[1]#novo posto
    # ind2=S1.>ϵ*S1[1]#novo posto

    H1 = w1[:, ind1] * U[ind1, ind1] * Diagonal(S[ind1])

    H2 = V[ind1, ind1]' * w2[:, ind1]'
    H1 * H2

    # G1=w1[:,ind2]*U1[ind2,ind2]*diagm(S1[ind2])
    # G2=V1[ind2,ind2]'*w3[:,ind2]'
    # G1*G2
    H1, H2
end



function Bernsteins(p, t)
    h = t + 1.0
    u = 1.0 - t
    B = zeros(p + 1)
    if t < 0
        for k = 0:p
            B[k+1] = (u / 2)^p * binomial(p, k) * (h / u)^k
        end
    else
        for k = 0:p
            B[k+1] = (h / 2)^p * binomial(p, k) * (u / h)^(p - k)
        end
    end
    return B
end

function dBernsteins(p, t)
    b1 = Bernsteins(p - 1, t)
    dB = ([0; b1] - [b1; 0]) * p
end

function greville(knots, p, β = 0.0)

    n = size(knots, 1) - p - 1
    result = []
    for i = 1:n
        result = append!(result, sum(knots[i+1:i+p]) / p)
    end
    if β > 0.0
        aux = result[1] + β * (result[2] - result[1])
        result[end] = result[end] + β * (result[end-1] - result[end])
        result[1] = aux
    end
    result
end

function calc_fforma(t, elem_j::bezier, w)
    Be = Bernsteins(elem_j.p, t)
    dBedxi = dBernsteins(elem_j.p, t)
    Ce = elem_j.C
    Wb = elem_j.Wb # Bézier weights (nodal values)
    # @infiltrate
    wb = dot(Be, Wb)   #  Bézier weight function
    dwbdxi = dot(dBedxi, Wb)  # derivada de w com relaçao a xi
    # Shape function and derivatives
    R = w .* Ce * Be / wb # Equação (65) do artigo do Borden R tem dim [(p+1) x 1]
    dRdxi = w .* Ce * (dBedxi[:, 1] / wb - dwbdxi * Be / (wb * wb)) # Eq. (66)
    return R, dRdxi
end
