

function interpola(r)
    if r == 0
        return 0
    end
    r^2 * log(r)
    # r + 1

end
function int_interpolaρdρ(r)
    (4 * r^4 * log(r) - r^4) / 16
    # r^3/3
end


function Monta_M_RIMd(dad::potencial, npg)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi, w = gausslegendre(npg)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(n_pontos)
    M1 = zeros(n_pontos)
    F, D = FeD(dad)
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


function Finterp(dad, Tree1, Tree2, block; ninterp=3, compressão=true, ϵ=1e-3)
    # arg = [NOS1,NOS_GEO1,tipoCDC,valorCDC,normal,ELEM1,k]
    #         1      2        3      4        5     6   7
    n = size(block, 1)               # Quantidade de Submatrizes
    Faca = Array{Any}(undef, n, 2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    Daca = Array{Any}(undef, n)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    nodes = [dad.NOS; dad.pontos_internos]
    M, M1 = calcMs(dad, 8)
    for i = 1:n                       # Para cada Submatriz
        # @timeit to "Para cada Submatriz" begin
        b1 = Tree1[block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # Submatriz = Produto cartesiano I x J
        cols = vcat([dad.ELEM[ii].indices for ii in b2]...)

        if block[i, 3] == 0                # Se esses blocos não são admissiveis
            Faca[i, 1], Daca[i, 1] = FeD(dad, b1, cols)
            # Salva na linha i da 1º coluna da matriz Aaca a matriz H e salva G na matriz B
        else                              # Caso contrario (Se blocos são admissiveis)
            Faca[i, 1], Faca[i, 2], Daca[i] = FeDinterp(dad, b1, cols, 8, ninterp)
        end
        # end

    end
    return Faca, Daca
end

function FeD(dad, b1=0, b2=0)
    nodes = [dad.NOS; dad.pontos_internos]
    if b1 == 0
        n_pontos = size(nodes, 1)
        b1 = 1:n_pontos
        b2 = 1:n_pontos
    end
    n1, n2 = size(b1, 1), size(b2, 1)
    F = zeros(n1, n2)
    D = zeros(n1, n2)

    for i = 1:n1
        xi = nodes[b1[i], 1]
        yi = nodes[b1[i], 2]
        for j = 1:n2
            if b1[i] == b2[j]
                continue
            end
            xj = nodes[b2[j], 1]
            yj = nodes[b2[j], 2]
            r = sqrt((xi - xj)^2 + (yi - yj)^2)
            F[i, j] = interpola(r)
            D[i, j] = -log(r) / (2 * π * dad.k)
        end
    end
    F, D
end




function FeDinterp(dad::potencial, b1, b2, npg=8, ninterp=3)
    nodes = [dad.NOS; dad.pontos_internos][b1, :]

    nodes2 = [dad.NOS; dad.pontos_internos][b2, :]
    xmax = maximum(nodes, dims=1)
    xmin = minimum(nodes, dims=1)

    xs = criapontosinterp(ninterp)
    fontes, L, ninterp1, ninterp2 = gera_interpolação(ninterp, nodes, xmax, xmin, xs)

    F = zeros(ninterp1 * ninterp2, 0)
    D = zeros(ninterp1 * ninterp2, 0)
    n1, n2 = Nlinear(xs)
    xks = n1 * xmin + n2 * xmax
    nb2 = size(b2, 1)

    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    for j = 1:nb2

        f1 = zeros(ninterp1 * ninterp2)
        d1 = zeros(ninterp1 * ninterp2)
        ci = 0
        for i2 = 1:ninterp1
            for i1 = 1:ninterp2
                ci += 1
                pf = [xks[i1, 1], xks[i2, 2]]   # Coordenada (x,y)  dos pontos fonte
                r = norm(pf - nodes2[j, :])
                f = interpola(r)
                d = -log(r) / (2 * π * dad.k)

                f1[ci] += f
                d1[ci] += d
            end
        end
        F = [F f1]
        D = [D d1]
    end
    L, F, D
end




function Monta_M_RIM(dad::placa_fina, npg1=10, npg2=10)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    n_canto = size(dad.k.cantos, 1)

    n2 = 2n_nos + n_noi + n_canto
    qsi1, w1 = gausslegendre(npg1)
    qsi2, w2 = gausslegendre(npg2)
    n_pontos = n_nos + n_noi + n_canto
    nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    M = zeros(2n_nos + n_noi + n_canto, n_pontos)
    F = zeros(n_pontos, n_pontos)
    normal_fonte = calc_normais(dad)

    # Cálculo da matriz [F]
    @showprogress "Montando F" for i = 1:n_pontos
        xi = nodes[i, 1]
        yi = nodes[i, 2]
        for j = 1:n_pontos
            xj = nodes[j, 1]
            yj = nodes[j, 2]
            r = sqrt((xi - xj)^2 + (yi - yj)^2)
            F[i, j] = interpola(r)
        end
    end
    @showprogress "Montando M" for i = 1:n_pontos
        if i <= n_nos
            pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
            nf = normal_fonte[i, :]
            caso = "contorno"
        elseif i <= n_nos + n_noi
            pf = dad.pontos_internos[i-n_nos, :] # Coordenada (x,y)
            nf = zeros(2)
            caso = "interno"
        else
            pf = dad.k.cantos[i-n_nos-n_noi, 2:3] # Coordenada (x,y)
            nf = zeros(2)
            caso = "canto"
        end
        pr = nodes[i, :]
        for x = [:R0, :S0, :dRdx0, :dRdy0, :dSdx0, :dSdy0]
            @eval $x = zeros(2)
        end
        pre = [zeros(2) for idx in 1:6]
        for j = 1:n_pontos #Laço dos pontos fontes
            pr = nodes[j, :]
            for el = 1:nelem
                elem_j = dad.ELEM[el]#Laço dos elementos
                x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos

                m_el = calc_m(x, pf, nf, pr, qsi1, w1, elem_j, dad, pre)
                if caso == "contorno"
                    M[2i-1:2i, j] += m_el
                else
                    M[n_nos+i, j] += m_el[1]
                end
            end
        end
    end

    # @infiltrate
    M = M / F
    M2 = zeros(n2, n2)
    for i = 1:n_pontos
        if i <= n_nos
            M2[:, 2*i-1] = M[:, i]
        else
            M2[:, n_nos+i] = M[:, i]
        end
    end
    M2
end

function calc_m(x, pf, nf, pr, qsi1, w1, elem, dad, pre)

    npg = length(w1)
    m_el = zeros(2)

    for i = 1:npg
        N, dN = calc_fforma(qsi1[i], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        ny = -dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        nx = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        r = pg' - pf
        m = calcula_F(pr, pf, nf, pg, [nx, ny], qsi1, w1, dad, pre)
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w1[i]
    end
    return m_el
end

function calcula_F(pr, pf, nf, pg, n, qsi, w, dad, pre) #???????????????/
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(2) # Inicializa a integral de F_area
    theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
        rline = norm(xc - pr)
        ro = r / 2 * (qsi[i] + 1)
        term = compute_domainterms(ro, theta, nf, dad, pre)
        f = interpola(rline)
        # @infiltrate
        F = F + term * f * ro * drodqsi * w[i]# Integral de F_area
    end
    return F
end

function calcula_F1(pr, pf, nf, pg, n, qsi, w, dad, pre) #???????????????/
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(2) # Inicializa a integral de F_area
    theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
        rline = norm(xc - pr)
        ro = r / 2 * (qsi[i] + 1)
        term = compute_domainterms(ro, theta, nf, dad, pre)
        f = interpola(rline)
        # @infiltrate
        F = F + term * ro * drodqsi * w[i]# Integral de F_area
    end
    return F
end
function compute_domainterms(r, theta, nf, dad::Union{placa_fina}, pre)


    m1 = nf[1]
    m2 = nf[2]
    # Distance of source and field points
    r1 = r .* cos(theta)
    r2 = r .* sin(theta)

    R, S, dRdx, dRdy, dSdx, dSdy = pre

    G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2

    H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2

    C1 = ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[1])

    C2 = ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[2])

    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)

    a = 1


    for i = 1:2
        R[i] = r^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) * (log(r^2 / a^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)) - 3) - 4 * r^2 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        S[i] = r^2 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) * (log(r^2 / a^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)) - 3) + r^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dRdx[i] = 2 * r * (cos(theta) + dad.k.d[i] * sin(theta)) * (log(r^2 / a^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)) - 2) - 4 * r * dad.k.e[i] * sin(theta) * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dRdy[i] = 2 * r * (dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) - dad.k.e[i]^2 * sin(theta)) * (log(r^2 / a^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)) - 2) - 4 * r * dad.k.e[i] * (cos(theta) + 2 * dad.k.d[i] * sin(theta)) * atan((dad.k.e[i] * sin(theta)), (cos(theta) + dad.k.d[i] * sin(theta)))

        dSdx[i] = r * dad.k.e[i] * sin(theta) * (log(r^2 / a^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)) - 2) + 2 * r * (cos(theta) + dad.k.d[i] * sin(theta)) * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dSdy[i] = r * dad.k.e[i] * (cos(theta) + 2 * dad.k.d[i] * sin(theta)) * (log(r^2 / a^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)) - 2) + 2 * r * (dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) - dad.k.e[i]^2 * sin(theta)) * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
    end


    w = 1 / (8 * pi * dad.k.D22) * (C1 * R[1] + C2 * R[2] + C3 * (S[1] - S[2]))

    dwdx = 1 / (8 * pi) * (C1 * dRdx[1] + C2 * dRdx[2] + C3 * (dSdx[1] - dSdx[2]))
    dwdy = 1 / (8 * pi) * (C1 * dRdy[1] + C2 * dRdy[2] + C3 * (dSdy[1] - dSdy[2]))

    dwdm = -(dwdx * m1 + dwdy * m2) / dad.k.D22

    [w, dwdm]
end
