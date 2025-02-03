"""
    calc_fforma(ξ, elem, deriv = true)

Calcula as funções de forma polinomiais gerais.

# Parâmetros
- `ξ`: Ponto de avaliação.
- `elem`: Estrutura que contém os pontos nodais `ξs`.
- `deriv`: Booleano que indica se a derivada das funções de forma deve ser calculada. O padrão é `true`.

# Retorno
- Se `deriv` for `true`, retorna uma tupla `(N, dN)` onde `N` é o vetor das funções de forma e `dN` é o vetor das derivadas das funções de forma.
- Se `deriv` for `false`, retorna apenas `N`.

# Notas
- As funções de forma são calculadas usando o interpolador de Lagrange.
- Se `ξ` for igual a algum dos pontos nodais `ξs`, um pequeno valor é adicionado a `ξ` para evitar divisão por zero.

# Exemplo
"""
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
                    ξ = ξ + 1e-12
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
"""
    calc_fforma_gen(ξ, ξs, deriv = true)

Calcula as funções de forma gerais para elementos de um método numérico.

# Parâmetros
- `ξ`: Ponto de avaliação.
- `ξs`: Vetor contendo as coordenadas dos nós.
- `deriv`: Booleano opcional que indica se a derivada das funções de forma deve ser calculada. O padrão é `true`.

# Retorno
- Se `deriv` for `true`, retorna uma tupla `(N, dN)`, onde `N` é o vetor das funções de forma e `dN` é o vetor das derivadas das funções de forma.
- Se `deriv` for `false`, retorna apenas `N`.

"""
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

""" Calcula as funções de forma lineares contínuas N1 e N2
"""
function Nlinear(qsi)
    N1 = 1 / 2 * (1.0 .- qsi) # Função de forma N1 => linear contínua
    N2 = 1 / 2 * (1.0 .+ qsi) # Função de forma N2 => linear contínua
    return N1, N2
end

"""
    lagrange(pg, x, n)

Calcula o polinômio interpolador de Lagrange.

# Parâmetros
- `pg`: Vetor de pontos de interpolação.
- `x`: Ponto onde o polinômio será avaliado.
- `n`: Número de pontos de interpolação.

# Retorno
- Valor do polinômio interpolador de Lagrange avaliado em `x`.
"""
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
"""
    lagrange(pg, x1, n1, x2, n2)

Calcula a matriz de interpolação de Lagrange para os pontos fornecidos.

# Parâmetros
- `pg`: Matriz de pontos de grade, onde cada linha representa um ponto e cada coluna representa uma dimensão.
- `x1`: Vetor de coordenadas na primeira dimensão.
- `n1`: Número de pontos na primeira dimensão.
- `x2`: Vetor de coordenadas na segunda dimensão.
- `n2`: Número de pontos na segunda dimensão.

# Retorno
- `L`: Matriz de interpolação de Lagrange de tamanho `(ni, n1 * n2)`, onde `ni` é o número de pontos na grade `pg`.
"""
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
"""
    lagrange(pg, x1, n1, x2, n2, x3, n3)

Calcula a matriz de interpolação de Lagrange para os pontos fornecidos.

# Parâmetros
- `pg`: Matriz de pontos de grade, onde cada linha representa um ponto e cada coluna representa uma dimensão.
- `x1`: Vetor de coordenadas na primeira dimensão.
- `n1`: Número de pontos na primeira dimensão.
- `x2`: Vetor de coordenadas na segunda dimensão.
- `n2`: Número de pontos na segunda dimensão.
- `x3`: Vetor de coordenadas na terceira dimensão.
- `n3`: Número de pontos na terceira dimensão.

# Retorno
- `L`: Matriz de interpolação de Lagrange de tamanho `(ni, n1 * n2 * n3)`, onde `ni` é o número de pontos na grade `pg`.
"""
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
"""
    recompressão(L, H)

Função que realiza a recompressão de matrizes utilizando decomposição QR e SVD.

# Parâmetros
- `L`: Matriz de entrada.
- `H`: Matriz de entrada.

# Retorno
- `H1`: Primeira matriz recomprimida.
- `H2`: Segunda matriz recomprimida.

# Descrição
A função `recompressão` realiza os seguintes passos:
1. Decomposição QR das matrizes `L` e `H'`.
2. Conversão das matrizes resultantes `w1` e `w2` para arrays.
3. Decomposição SVD do produto das matrizes `r1` e `r2'`.
4. Seleção dos índices onde os valores singulares são maiores que um limiar `ϵ` multiplicado pelo maior valor singular.
5. Cálculo das matrizes recomprimidas `H1` e `H2` utilizando os índices selecionados.
6. Retorno das matrizes `H1` e `H2`.


"""

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
    # H1 * H2

    # G1=w1[:,ind2]*U1[ind2,ind2]*diagm(S1[ind2])
    # G2=V1[ind2,ind2]'*w3[:,ind2]'
    # G1*G2
    H1, H2
end


"""
    Bernsteins(p, t)

Calcula os polinômios de Bernstein para um dado grau `p` e parâmetro `t`.

# Parâmetros
- `p::Int`: Grau do polinômio de Bernstein.
- `t::Float64`: Parâmetro no qual o polinômio de Bernstein será avaliado.

# Retorno
- `h::Float64`: Resultado da avaliação do polinômio de Bernstein.

"""
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
"""
    dBernsteins(p, t)

Calcula a derivada das funções de Bernstein de grau `p` no ponto `t`.

# Parâmetros
- `p::Int`: O grau das funções de Bernstein.
- `t::Float64`: O ponto no qual a derivada será calculada.

# Retorno
- `dB::Vector{Float64}`: Um vetor contendo os valores das derivadas das funções de Bernstein de grau `p` no ponto `t`.

# Exemplo
"""
function dBernsteins(p, t)
    b1 = Bernsteins(p - 1, t)
    dB = ([0; b1] - [b1; 0]) * p
end
"""
    greville(knots, p, β = 0.0)

Calcula os nós de Greville para uma dada sequência de nós e grau da B-spline.

# Parâmetros
- `knots::Vector{Float64}`: Vetor contendo a sequência de nós.
- `p::Int`: Grau da B-spline.
- `β::Float64`: Parâmetro opcional para ajuste dos nós de Greville. O valor padrão é 0.0.

# Retorno
- `Vector{Float64}`: Vetor contendo os nós de Greville calculados.

# Exemplo
"""
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
"""
    calc_fforma(t, elem_j::bezier, w)

Calcula a forma de função para um elemento de Bézier.

# Parâmetros
- `t`: Ponto de avaliação.
- `elem_j::bezier`: Elemento de Bézier.
- `w`: Peso associado ao ponto de avaliação.

# Retorno
Retorna a forma de função calculada no ponto `t` para o elemento `elem_j` com o peso `w`.
"""
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
