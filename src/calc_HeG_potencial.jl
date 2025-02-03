"""
    calc_HeG(dad::potencial, npg=8; Pint=false)

Calcula as matrizes H e G para o problema de potencial usando o Método dos Elementos de Contorno (BEM).

# Argumentos
- `dad::potencial`: Os dados do problema de potencial, que incluem os elementos (ELEM) e nós (NOS).
- `npg::Int=8`: O número de pontos de quadratura de Gauss-Legendre a serem usados para a integração numérica.
- `Pint::Bool=false`: Uma flag para indicar se deve usar pontos internos (o padrão é falso).

# Retornos
- `H::Matrix{Float64}`: A matriz H.
- `G::Matrix{Float64}`: A matriz G.

# Descrição
Esta função calcula as matrizes H e G para um dado problema de potencial usando o Método dos Elementos de Contorno (BEM). 
Ela itera sobre os pontos fontes e elementos, realiza a integração numérica usando a quadratura de Gauss-Legendre, 
e preenche as matrizes H e G com os valores computados.

# Exemplo
```julia
H, G = calc_HeG(dad, 8)"""
function calc_HeG(dad::potencial, npg = 8; Pint = false)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    H = zeros(n, n)
    G = zeros(n, n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    contafonte = 1

    for elem_i in dad.ELEM  #Laço dos pontos fontes
        for ind_elem in elem_i.indices
            pf = dad.NOS[ind_elem, :]   # Coordenada (x,y)  dos pontos fonte
            for elem_j in dad.ELEM  #Laço dos elementos
                x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
                Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
                eet =
                    (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                    norm(Δelem)^2 + elem_j.ξs[1]
                N_geo, ~ = calc_fforma(eet, elem_j)
                ps = N_geo' * x
                b = norm(ps' - pf)#/norm(Δelem)
                # eta, Jt = telles(qsi, eet, b)
                eta, Jt = sinhtrans(qsi, eet, b)
                # @show eet,b
                # eta, wt = pontosintegra(dad.NOS, elem_j, ind_elem, qsi, w)
                # @show norm(eta-eta1)

                h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
                # h, g = integraelem(pf, x, qsi, w, elem_j, dad)
                # @infiltrate contafonte == 2
                H[contafonte, elem_j.indices] = h
                G[contafonte, elem_j.indices] = g

            end
            contafonte += 1
        end
    end
    for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
        # H[i, i] += -0.5
        # G[i, i] = 0
        H[i, i] = 0
        H[i, i] = -sum(H[i, :])
    end
    # somaH = H * (dad.NOS[:, 1] + dad.NOS[:, 2])
    # somaG = G * (normal_fonte[:, 1] + normal_fonte[:, 2]) * dad.k
    # for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    #   G[i, i] = (-somaH[i] - somaG[i]) / (normal_fonte[i, 1] + normal_fonte[i, 2])
    # end
    # corrigediag!(H, G, dad)
    if Pint
        ni = size(dad.pontos_internos, 1)
        Hi = zeros(ni, n)
        Gi = zeros(ni, n)
        for i = 1:ni
            pf = dad.pontos_internos[i, :]   # Coordenada (x,y)  dos pontos fonte
            for elem_j in dad.ELEM  #Laço dos elementos
                x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
                Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
                eet =
                    (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                    norm(Δelem)^2 + elem_j.ξs[1]
                N_geo = calc_fforma(eet, elem_j, false)
                ps = N_geo' * x
                b = norm(ps' - pf) / norm(Δelem)
                eta, Jt = sinhtrans(qsi, eet, b)
                # eta,Jt=telles(qsi,eet)
                h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
                cols = elem_j.indices
                # @infiltrate
                Hi[i, cols] = h
                Gi[i, cols] = g
            end
        end
        H = [H zeros(n, ni); Hi -I]
        G = [G; Gi]
    end
    H, G
end

"""
    calc_HeG_hiper(dad::potencial, npg=8)

Calcula as matrizes H e G hiper-singulares para o problema de potencial usando o Método dos Elementos de Contorno (BEM).

# Argumentos
- `dad::potencial`: Os dados do problema de potencial, que incluem os elementos (ELEM) e nós (NOS).
- `npg::Int=8`: O número de pontos de quadratura de Gauss-Legendre a serem usados para a integração numérica.

# Retornos
- `H::Matrix{Float64}`: A matriz H hiper-singular.
- `G::Matrix{Float64}`: A matriz G hiper-singular.

# Descrição
Esta função calcula as matrizes H e G hiper-singulares para um dado problema de potencial usando o Método dos Elementos de Contorno (BEM). 
Ela itera sobre os pontos fontes e elementos, realiza a integração numérica usando a quadratura de Gauss-Legendre, 
e preenche as matrizes H e G com os valores computados.


# Exemplo
```julia
H, G = calc_HeG(dad, 8)"""
function calc_HeG_hiper(dad::potencial, npg = 8)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    H = zeros(n, n)
    G = zeros(n, n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    contafonte = 1

    for elem_i in dad.ELEM  #Laço dos pontos fontes
        for ind_elem in elem_i.indices
            pf = dad.NOS[ind_elem, :]   # Coordenada (x,y)  dos pontos fonte
            for elem_j in dad.ELEM  #Laço dos elementos
                x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos

                nosing = elem_j.indices .== contafonte
                if sum(nosing) == 1
                    no_pf = findfirst(nosing)
                    xi0 = elem_j.ξs[no_pf]
                    # h, g = integraelemsing(pf, nf, x, qsi2, w2, elem_j, dad, xi0)
                    # hn, gn = integraelemsing_num(pf, nf, x, elem_j, dad, pre, xi0, 30)
                    h, g = integraelem_hiper_sing(
                        pf,
                        dad.normal[contafonte, :],
                        x,
                        xi0,
                        elem_j,
                        dad,
                    )
                else
                    Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
                    eet =
                        (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                        norm(Δelem)^2 + elem_j.ξs[1]
                    N_geo, ~ = calc_fforma(eet, elem_j)
                    ps = N_geo' * x
                    b = norm(ps' - pf)#/norm(Δelem)
                    # eta, Jt = telles(qsi, eet, b)
                    eta, Jt = sinhtrans(qsi, eet, b)
                    # eta, Jt2 = sinhtrans(eta, eet, b)
                    # @show eet,b
                    # eta, wt = pontosintegra(dad.NOS, elem_j, ind_elem, qsi, w)
                    # @show norm(eta-eta1)
                    # @infiltrate


                    h, g = integraelem_hiper(
                        pf,
                        dad.normal[contafonte, :],
                        x,
                        eta,
                        w .* Jt,
                        elem_j,
                        dad,
                    )
                end
                # h, g = integraelem_hiper(pf, dad.normal[contafonte, :], x, qsi, w, elem_j, dad)
                # @infiltrate contafonte == 2
                H[contafonte, elem_j.indices] = h
                G[contafonte, elem_j.indices] = g

            end
            contafonte += 1
        end
    end
    for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
        G[i, i] += -0.5

        # G[i, i] = 0
        # H[i, i] = 0
        # H[i, i] = -sum(H[i, :])
    end
    # somaH = H * (dad.NOS[:, 1] + dad.NOS[:, 2])
    # somaG = G * (dad.normal[:, 1] + dad.normal[:, 2]) * dad.k
    # for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    #   G[i, i] = (-somaH[i] - somaG[i]) / (dad.normal[i, 1] + dad.normal[i, 2])
    # end
    # corrigediag!(H, G, dad)

    H, G
end

"""
Função para integrar elementos.

# Parâmetros
- `pf`: Coordenadas do ponto de fonte.
- `x`: Coordenadas dos pontos pertencentes ao elemento.
- `eta`: Coordenadas eta dos pontos de integração.
- `w`: Pesos dos pontos de integração.
- `elem`: Elemento a ser integrado. 
- `dad:: estrutura: Estrutura contendo os dados do problema.

# Retorno
- `h`: Valor da integral de h no elemento.
- `g`: Valor da integral de g no elemento.
# Descrição
Esta função realiza a integração dos elementos fornecidos utilizando os pontos de integração e pesos especificados.
"""
function integraelem(pf, x, eta, w, elem, dad::Union{potencial,helmholtz})
    h = zeros(Float64, size(elem))
    g = zeros(Float64, size(elem))

    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x    # Ponto de gauss interpolador
        r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = dN' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        Qast, Tast = calsolfund(r, [sy, -sx], dad)
        # @infiltrate
        # @show pg,pf,r,[sy,-sx],Qast
        h += N * Qast * dgamadqsi * w[k]
        g += N * Tast * dgamadqsi * w[k]

    end
    h, g
end
"""
Função para integrar elementos considerando a solução fundamental hipersingular.

# Parâmetros
- `pf`: Coordenadas do ponto de fonte.
- `x`: Coordenadas dos pontos pertencentes ao elemento.
- `eta`: Coordenadas eta dos pontos de integração.
- `w`: Pesos dos pontos de integração.
- `elem`: Elemento a ser integrado. 
- `dad:: estrutura: Estrutura contendo os dados do problema.

# Retorno
- `h`: Valor da integral de h no elemento.
- `g`: Valor da integral de g no elemento.
# Descrição
Esta função realiza a integração dos elementos fornecidos utilizando os pontos de integração e pesos especificados.
"""
function integraelem_hiper(pf, nf, x, eta, w, elem, dad::Union{potencial,helmholtz})
    h = zeros(Float64, size(elem))
    g = zeros(Float64, size(elem))

    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x    # Ponto de gauss interpolador
        r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = dN' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        Qast, Tast = calsolfund_hiper(r, [sy, -sx], nf, dad)
        # @infiltrate
        # @show pg,pf,r,[sy,-sx],Qast
        h += N * Qast * dgamadqsi * w[k]
        g += N * Tast * dgamadqsi * w[k]

    end
    h, g
end
"""
Função para integrar elementos singulares considerando a solução fundamental hipersingular.

# Parâmetros
- `pf`: Coordenadas do ponto de fonte.
- `nf`: Normal do ponto de fonte.
- `x`: Coordenadas dos pontos pertencentes ao elemento.
- `xi0`: Coordenada eta dos ponto singular.
- `elem`: Elemento a ser integrado. 
- `dad:: estrutura: Estrutura contendo os dados do problema.
- `npg::Int=20`: Número de pontos de quadratura PTVSI a serem usados para a integração numérica.
# Retorno
- `h`: Valor da integral de h no elemento.
- `g`: Valor da integral de g no elemento.
# Descrição
Esta função realiza a integração dos elementos singulares.
"""
function integraelem_hiper_sing(
    pf,
    nf,
    x,
    xi0,
    elem,
    dad::Union{potencial,helmholtz},
    npg = 20,
)
    h = zeros(Float64, size(elem))
    g = zeros(Float64, size(elem))
    eta, w = novelquad(2, xi0, 20 * npg)

    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x    # Ponto de gauss interpolador
        r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = dN' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        Qast, Tast = calsolfund_hiper(r, [sy, -sx], nf, dad)
        # @infiltrate
        # @show pg,pf,r,[sy,-sx],Qast
        h += N * Qast * dgamadqsi * w[k]
        g += N * Tast * dgamadqsi * w[k]

    end
    h, g
end


"""
    calc_Ti(dad::Union{potencial, helmholtz}, T, q, npg=8)

Calcula a matriz `Ti` com as temperaturas nos pontos internos para um dado problema de potencial ou Helmholtz.

# Parâmetros
- `dad::Union{potencial, helmholtz}`: Estrutura contendo os dados do problema, incluindo elementos, nós e pontos internos.
- `T`: Vetor de temperaturas nos nós.
- `q`: Vetor de fluxos de calor nos nós.
- `npg`: Número de pontos de Gauss para a quadratura (padrão é 8).

# Retorno
- `Ti`: Vetor resultante após a integração dos elementos.

# Descrição
A função percorre todos os pontos internos e elementos do contorno, calculando a contribuição de cada elemento para o vetor `Ti` utilizando a quadratura de Gauss. Para cada ponto interno, a função:
1. Obtém a coordenada do ponto fonte.
2. Percorre todos os elementos do contorno.
3. Calcula as coordenadas dos nós geométricos do elemento.
4. Calcula a transformação de coordenadas e os fatores de forma.
5. Realiza a integração dos elementos utilizando a quadratura de Gauss.
6. Atualiza o vetor `Ti` com as contribuições de cada elemento.
"""
function calc_Ti(dad::Union{potencial,helmholtz}, T, q, npg = 8)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.pontos_internos, 1)
    Ti = zeros(n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss

    for i = 1:n
        pf = dad.pontos_internos[i, :]   # Coordenada (x,y) dos pontos fonte
        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
            eet =
                (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                norm(Δelem)^2 + elem_j.ξs[1]
            N_geo, ~ = calc_fforma(eet, elem_j)
            ps = N_geo' * x
            b = norm(ps' - pf) / norm(Δelem)
            eta, Jt = sinhtrans(qsi, eet, b)
            h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
            # h, g = integraelem(pf, x, qsi, w, elem_j, dad)
            Ti[i] += h' * T[elem_j.indices] - g' * q[elem_j.indices]
        end
    end
    Ti
end
"""
    calc_Ti(dad::potencial_iga, T, q, npg = 8)

Calcula a matriz `Ti` com as temperaturas nos pontos internos para um dado problema de potencial com elementos de contorno isogeométricos.

# Parâmetros
- `dad::Union{potencial, helmholtz}`: Estrutura contendo os dados do problema, incluindo elementos, nós e pontos internos.
- `T`: Vetor de temperaturas nos nós.
- `q`: Vetor de fluxos de calor nos nós.
- `npg`: Número de pontos de Gauss para a quadratura (padrão é 8).

# Retorno
- `Ti`: Vetor resultante após a integração dos elementos.

# Descrição
A função percorre todos os pontos internos e elementos do contorno, calculando a contribuição de cada elemento para o vetor `Ti` utilizando a quadratura de Gauss. Para cada ponto interno, a função:
1. Obtém a coordenada do ponto fonte.
2. Percorre todos os elementos do contorno.
3. Calcula as coordenadas dos nós geométricos do elemento.
4. Calcula a transformação de coordenadas e os fatores de forma.
5. Realiza a integração dos elementos utilizando a quadratura de Gauss.
6. Atualiza o vetor `Ti` com as contribuições de cada elemento.
"""
function calc_Ti(dad::potencial_iga, T, q, npg = 8)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.pontos_internos, 1)
    Ti = zeros(n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss

    for i = 1:n
        pf = dad.pontos_internos[i, :]   # Coordenada (x,y) dos pontos fonte
        for elem_j in dad.ELEM  #Laço dos elementos
            xf = elem_j.limites[:, 2]
            x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
            Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
            eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
            pc = dad.pontos_controle[:, elem_j.indices]
            # @infiltrate
            cf = pc[1:2, :] ./ pc[4, :]'
            N, dN = calc_fforma(eet, elem_j, pc[4, :])
            ps = cf * N
            b = norm(ps - pf) / norm(Δelem)

            eta, Jt = sinhtrans(qsi, eet, b)
            h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)

            Ti[i] += h' * T[elem_j.indices] - g' * q[elem_j.indices]
        end
    end
    Ti
end

"""
    calc_Aeb(dad::Union{potencial, helmholtz}, npg=8)

Calcula a matriz A e o velor b para os dados fornecidos.

# Parâmetros
- `dad::Union{potencial, helmholtz}`: Tipo de dado que pode ser `potencial` ou `helmholtz`.
- `npg`: Número de pontos de Gauss (opcional, padrão é 8).

# Retorna
- A matriz A e o vetor b calculados com base nos parâmetros fornecidos.

"""
function calc_Aeb(dad::Union{potencial,helmholtz}, npg = 8)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    ni = size(dad.pontos_internos, 1)
    A = zeros(n + ni, n)
    B = zeros(n + ni)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss

    for i = 1:n+ni  #Laço dos pontos fontes
        pf = [dad.NOS; dad.pontos_internos][i, :]   # Coordenada (x,y)  dos pontos fonte
        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
            eet =
                (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                norm(Δelem)^2 + elem_j.ξs[1]
            N_geo, ~ = calc_fforma(eet, elem_j)
            ps = N_geo' * x
            b = norm(ps' - pf) / norm(Δelem)
            eta, Jt = sinhtrans(qsi, eet, b)
            h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
            h[elem_j.indices.==i] = h[elem_j.indices.==i] .- 0.5
            if elem_j.tipoCDC == 1
                A[i, elem_j.indices] = h
                B[i] += dot(g, elem_j.valorCDC)
            else
                A[i, elem_j.indices] = -g
                B[i] += -dot(h, elem_j.valorCDC)
            end

        end
    end
    [A [zeros(n, ni); -diagm(ones(ni))]], B
end



"""
    calsolfund_hiper(r, n, nf, prob::Union{potencial,potencial_iga})

Calcula a solução fundamental hipersingular.

# Parâmetros
- `r`: Distância radial.
- `n`: normal do ponto de integração.
- `nf`: normal do ponto fonte.
- `prob`: Tipo de problema, pode ser `potencial` ou `potencial_iga`.

# Retorno
Retorna a solução fundamental hipersingular calculada com base nos parâmetros fornecidos.
"""
function calsolfund_hiper(r, n, nf, prob::Union{potencial,potencial_iga})
    R = norm(r)
    # @infiltrate
    Qast = -(dot(nf, n) - 2 * dot(r, n) / R * dot(r, nf) / R) / R^2 / (2 * π)
    Tast = dot(r, nf) / R^2 / (2 * π * prob.k)
    Qast, Tast
end
"""
    calsolfund(r, n prob::Union{potencial,potencial_iga})

Calcula a solução fundamental.

# Parâmetros
- `r`: Distância radial.
- `n`: normal do ponto de integração.
- `prob`: Tipo de problema, pode ser `potencial` ou `potencial_iga`.

# Retorno
Retorna a solução fundamental calculada com base nos parâmetros fornecidos.
"""
function calsolfund(r, n, prob::Union{potencial,potencial_iga})
    R = norm(r)
    # @infiltrate
    Qast = dot(r, n) / R^2 / (2 * π)       # Equação 4.36
    Tast = -log(R) / (2 * π * prob.k)
    Qast, Tast
end

"""
    Monta_M_RIMd(dad::potencial, npg)

Função que monta a matriz M utilizando o DIBEM (Direct interpolation method).

# Parâmetros
- `dad::potencial`: Estrutura de dados contendo as informações do problema potencial.
- `npg`: Número de pontos de Gauss para integração.

# Retorno
- Matriz `A` resultante da montagem utilizando o método RIMd.

# Descrição
A função realiza os seguintes passos:
1. Calcula o número de nós (`n_nos`), elementos (`nelem`) e nós internos (`n_noi`).
5. Calcula as matrizes de funções radiais `F`  e das soluções fundamentais `D` utilizando a função `FeD`.
6. Calcula as matrizes `M` e `M1` utilizando a função `calcMs`.
7. Monta a matriz `A` utilizando as matrizes `M`, `F` e `D`.
8. Ajusta os elementos da diagonal principal de `A`.
9. Retorna a matriz `A` somada com a matriz diagonal `M1`.

"""
function Monta_M_RIMd(dad::potencial, npg)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi, w = gausslegendre(npg)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(n_pontos)
    M1 = zeros(n_pontos)
    F, D = FeD(dad, nodes)
    M, M1 = calcMs(dad, npg)
    # @show size(M)
    # @show length(M)
    A = ones(length(M)) * M' / F .* D
    for i = 1:n_pontos #Laço dos pontos radiais
        A[i, i] = 0
        A[i, i] = -sum(A[i, :])
    end
    A + diagm(0 => M1)
    # M, M1, F, D
end
"""
A função `FeD` calcula duas matrizes, `F` e `D`, com base nas coordenadas dos nós fornecidos.

# Parâmetros
- `dad`: Estrutura de dados `k`.
- `nodes`: Matriz `n x 2` onde cada linha representa as coordenadas `(x, y)` de um nó.

# Retorno
- `F`: Matriz `n x n` onde cada elemento `F[i, j]` é o resultado da função `interpola` aplicada à distância entre os nós `i` e `j`.
- `D`: Matriz `n x n` onde cada elemento `D[i, j]` é o valor `-log(r) / (2 * π * dad.k)`, sendo `r` a distância entre os nós `i` e `j`.

# Notas
- A função ignora a diagonal principal das matrizes `F` e `D` (onde `i == j`).
- A função `interpola` deve ser definida em outro lugar no código.
"""
function FeD(dad, nodes)
    n = size(nodes, 1)
    F = zeros(n, n)
    D = zeros(n, n)

    for i = 1:n
        xi = nodes[i, 1]
        yi = nodes[i, 2]
        for j = 1:n
            if i == j
                continue
            end
            xj = nodes[j, 1]
            yj = nodes[j, 2]
            r = sqrt((xi - xj)^2 + (yi - yj)^2)
            F[i, j] = interpola(r)
            D[i, j] = -log(r) / (2 * π * dad.k)
        end
    end
    F, D
end
"""
    calcMs(dad::potencial, npg)

Calcula os valores das matrizes `M` e `M1` para um dado potencial `dad` utilizando a quadratura de Gauss-Legendre com `npg` pontos.

# Parâmetros
- `dad::potencial`: Estrutura contendo os dados do problema, incluindo nós (`NOS`), pontos internos (`pontos_internos`), elementos (`ELEM`) e constante `k`.
- `npg`: Número de pontos de Gauss-Legendre a serem utilizados na quadratura.

# Retorno
- `M`: Vetor com as integrais das funções radiais calculados para cada ponto fonte.
- `M1`: Vetor com as integrais das soluções fundamentais calculados para cada ponto fonte.

# Descrição
A função percorre todos os pontos radiais e elementos do problema, calculando os valores das matrizes `M` e `M1` através da função `calc_md`, que utiliza a quadratura de Gauss-Legendre para integração numérica. Os resultados são acumulados nos vetores `M` e `M1` e retornados ao final da execução.
"""
function calcMs(dad::potencial, npg)
    nodes = [dad.NOS; dad.pontos_internos]
    n_pontos = size(nodes, 1)
    M = zeros(n_pontos)
    M1 = zeros(n_pontos)
    qsi, w = gausslegendre(npg)

    for i = 1:n_pontos #Laço dos pontos radiais
        pf = nodes[i, :]
        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            m_el, m_el1 = calc_md(x, pf, dad.k, qsi, w, elem_j)
            M[i] = M[i] + m_el
            M1[i] = M1[i] + m_el1
        end
    end
    M, M1
end
"""
    calc_md(x, pf, k, qsi, w, elem)

Calcula o potencial e a sua derivada normal em um ponto fonte `pf` devido a um elemento `elem`.

# Parâmetros
- `x::Vector{Float64}`: Coordenadas dos nós do elemento.
- `pf::Vector{Float64}`: Coordenadas do ponto fonte.
- `k::Float64`: Coeficiente de condutividade térmica.
- `qsi::Vector{Float64}`: Pontos de Gauss para integração.
- `w::Vector{Float64}`: Pesos de Gauss para integração.
- `elem`: Dados do elemento 

# Retorna
- `m_el::Float64`: Integral no elemento da função de base radial.
- `m_el1::Float64`: Integral no elemento da soluçao fundamental.

# Descrição
A função `calc_md` realiza a integração numérica utilizando a técnica dos pontos de Gauss para calcular o potencial e sua derivada normal em um ponto fonte `pf` devido a um elemento `elem`. A função utiliza as funções de forma `N` e suas derivadas `dN_geo` para interpolar os pontos de Gauss e calcular as distâncias e vetores normais necessários para a integração. O resultado é o potencial `m_el` e sua derivada normal `m_el1` no ponto fonte `pf`.
"""
function calc_md(x, pf, k, qsi, w, elem)
    npg = length(w)
    m_el, m_el1 = 0, 0

    for i = 1:npg
        N, dN_geo = calc_fforma(qsi[i], elem)
        pg = N' * x    # Ponto de gauss interpolador
        r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = dN_geo' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

        nx = sy # Componente x do vetor normal unit�rio
        ny = -sx # Componente y do vetor normal unit�rio
        # @infiltrate
        R = norm(r)
        m = int_interpolaρdρ(R)
        m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)
        # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
        m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
    end
    return m_el, m_el1
end
