function calsolfund(pg, pf, n, nf, Casca::casca, pre)
    w, V = calsolfund(pg, pf, n, nf, Casca.dadplaca, pre)
    u, t = calsolfund(pg, pf, n, Casca.dadpe)
    k11 = 0
    k22 = 0
    k12 = 0
    R12 = 0
    A = Casca.dadpe.k.A3

    k11 = 0
    k22 = 0
    if (Casca.R11 != 0)
        k11 = 1 / Casca.R11
        k22 = nu / Casca.R11
    end
    if (Casca.R22 != 0)
        k11 = k11 + nu / Casca.R22
        k22 = k22 + 1 / Casca.R22
    end

    kn = dot([k11, k12], n)

    uc = [
        u zeros(2, 2)
        zeros(2, 2) w
    ]

    pc = [
        t zeros(2, 2)
        w.*kn V
    ]
    uc, pc
end


function calsolfund(pg, pf, n, nf, Casca::casca_aniso, pre)
    w, V = calsolfund(pg, pf, n, nf, Casca.dadplaca, pre)
    u, t = calsolfund(pg, pf, n, Casca.dadpe)
    k11 = 0
    k22 = 0
    k12 = 0
    R12 = 0
    A = dadpe.k.A3
    k11 = 0
    k22 = 0

    if (Casca.R11 != 0)
        k11 = A[1, 1] / Casca.R11
        k22 = A[1, 2] / Casca.R11
        k12 = A[1, 3] / Casca.R11
    end
    if (Casca.R22 != 0)
        k11 = k11 + A[1, 2] / Casca.R22
        k22 = k22 + A[2, 2] / Casca.R22
        k12 = k12 + A[2, 3] / Casca.R22
    end
    if (Casca.R12 != 0)
        k11 = k11 + (2 * A[1, 3]) / Casca.R12
        k22 = k22 + (2 * A[2, 3]) / Casca.R12
        k12 = k12 + (2 * A[3, 3]) / Casca.R12
    end

    kn = dot([k11, k12], n)

    uc = [
        u zeros(2, 2)
        zeros(2, 2) w
    ]

    pc = [
        t zeros(2, 2)
        w.*kn V
    ]
    uc, pc
end

function compute_domain_terms(r, theta, nf, n, Casca::casca, pre)
    m1 = nf[1]
    m2 = nf[2]
    D = Casca.dadplaca.k.D

    E, v = Casca.dadpe.k.E, Casca.dadpe.k.nu
    r = pg - pf      # Distancia entre ponto de gauss e ponto fonte

    GE = E / (2 * (1 + v))
    C = E * h / (1 - nu^2)

    rx = cos(theta)
    ry = sin(theta)

    k11 = 0
    k22 = 0
    if (Casca.R11 != 0)
        k11 = 1 / Casca.R11
        k22 = nu / Casca.R11
    end
    if (Casca.R22 != 0)
        k11 = k11 + nu / Casca.R22
        k22 = k22 + 1 / Casca.R22
    end
    prod1 = 4 * pi * (1 - v)
    prod2 = (3 - 4 * v) * log(1 / R)

    u11 = (prod2 + r1[1]^2) / (2 * prod1 * GE)
    u22 = (prod2 + r1[2]^2) / (2 * prod1 * GE)
    u12 = (r1[1] * r1[2]) / (2 * prod1 * GE)
    u21 = u12

    u11x =
        (rx .* ((3 - 4 * nu) * rx .^ 2 + (1 - 4 * nu) * ry .^ 2)) ./
        (8 * GE * (-1 + nu) * pi * r)
    u22y =
        (ry .* ((1 - 4 * nu) * rx .^ 2 + (3 - 4 * nu) * ry .^ 2)) ./
        (8 * GE * (-1 + nu) * pi * r)
    u12x = ((rx - ry) .* ry .* (rx + ry)) ./ (8 * GE * (-1 + nu) * pi * r)
    u12y = (-rx .^ 3 + rx .* ry .^ 2) ./ (8 * GE * (-1 + nu) * pi * r)

    term1 = C * (k11 * u11x + k22 * u12y)#P1
    term2 = C * (k11 * u12x + k22 * u22y)#P1

    # Distance of source and field points
    w = r^2 / (8 * pi * D) * (log(r) - 1 / 2)
    dwdx = -(r * cos(theta) * (-log(r^2))) / (8 * pi * D)
    dwdy = -(r * (-log(r^2)) * sin(theta)) / (8 * pi * D)

    dwdm = -(dwdx * m1 + dwdy * m2)

    d2wdx2 = (1 + cos(2 * theta) + log(r .^ 2)) / (8 * pi * D)
    d2wdxdy = sin(2 * theta) / (8 * pi * D)
    d2wdy2 = -(-1 + cos(2 * theta) - log(r .^ 2)) / (8 * pi * D)


    d2wdmdx = -(d2wdx2 * m1 + d2wdxdy * m2)
    d2wdmdy = -(d2wdxdy * m1 + d2wdy2 * m2)

    term3 = -C * k11 * dwdx#P3
    term4 = -C * k22 * dwdy#P3
    term5 = -C * k11 * d2wdmdx#P5
    term6 = -C * k22 * d2wdmdy#P5
    term7 = zeros(size(w))#P2
    term8 = zeros(size(dwdm))#P4
    if (Casca.R11 != 0)
        term7 = C * (k11 / Casca.R11) * w
        term8 = C * (k11 / Casca.R11) * dwdm
    end
    if (Casca.R22 != 0)
        term7 = term7 + C * (k22 / Casca.R22) * w
        term8 = term8 + C * (k22 / Casca.R22) * dwdm
    end


    termcasca = [
        0 0 term1 0
        0 0 term2 0
        term3 term4 term7 0
        term5 term6 term8 0
    ]
    termdinamico = [
        u11 u12 0 0
        u21 u22 0 0
        0 0 w 0
        0 0 dwdm 0
    ]
    termdinamico, termcasca
end



function Monta_M_RIMd(Casca::casca, npg)
    dad = Casca.dadplaca

    n_nos = size(dad.NOS, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    n_canto = size(dad.k.cantos, 1)

    pre = [zeros(2) for idx = 1:20]


    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    F = zeros(n_pontos, n_pontos)
    DP = zeros(4n_pontos, 4n_pontos)
    DM = zeros(4n_pontos, 4n_pontos)
    # Dx = zeros(2n_pontos, 2n_pontos)
    # Dy = zeros(2n_pontos, 2n_pontos)
    # dA = zeros(2n_pontos, 4n_pontos)
    M1 = zeros(n_pontos)
    M2 = zeros(2n_pontos, 2)
    # M2x = zeros(2n_pontos, 2)
    # M2y = zeros(2n_pontos, 2)
    normal_fonte = calc_normais(dad)
    for x in [:R0, :S0, :dRdx0, :dRdy0, :dSdx0, :dSdy0]
        @eval $x = zeros(2)
    end
    pre = [zeros(2) for idx = 1:6]
    # Cálculo da matriz [F]
    @showprogress "Montando F e D" for i = 1:n_pontos
        if i <= n_nos
            pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
            nf = normal_fonte[i, :]
            caso = "contorno"
        elseif i <= n_nos + n_noi
            pf = dad.pontos_internos[i-n_nos, :] # Coordenada (x,y)
            nf = zeros(2)
            caso = "interno"
        end
        for j = 1:n_pontos
            # @show i, j
            pr = nodes[j, :]
            r = pr - pf
            R = norm(r)
            # @infiltrate
            F[i, j] = interpola(R)
            if i == j
                continue
                if j > n_nos
                    nc = zeros(2)
                else
                    nc = normal_fonte[j, :]
                end
                dm, dp = compute_domain_terms(R, atan(r[2], r[1]), nf, nc, dad, pre)

                if caso == "contorno"
                    DM[4i-3:4i, 4j-3:4j] += dm
                    DP[4i-3:4i, 4j-3:4j] += dp
                elseif caso == "interno"
                    DM[4n_nos+3(i-n_nos)-2:n_nos+3(i-n_nos), 4j-3:4j] += dm[1:3, :]
                    DP[4n_nos+3(i-n_nos)-2:n_nos+3(i-n_nos), 4j-3:4j] += dp[1:3, :]
                elseif caso == "canto"
                    DM[3n_nos+2n_noi+i, 4j-3:4j] += dm[3, :]
                    DP[3n_nos+2n_noi+i, 4j-3:4j] += dp[3, :]
                end
            end
            qsi, w = gausslegendre(npg)

            for elem_j in dad.ELEM  #Laço dos elementos
                x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
                # @infiltrate
                m_el, m_elm, m_elp = calc_md(x, pf, nf, qsi, w, elem_j, dad, pre)
                # m_el, m_el1, m_el1x, m_el1y = calc_md(x, pf, nf, qsi, w, elem_j, dad, pre)
                if caso == "contorno"
                    M1[i] = M1[i] + m_el
                    Mm[4i-3:4i] += m_elm
                    Mp[4i-3:4i] += m_elp
                elseif caso == "interno"
                    M1[i] = M1[i] + m_el
                    Mm[4n_nos+3(i-n_nos)-2:n_nos+3(i-n_nos), :] += m_elm[1:3, :]
                    Mp[4n_nos+3(i-n_nos)-2:n_nos+3(i-n_nos), :] += m_elp[1:3, :]
                elseif caso == "canto"
                    M1[i] = M1[i] + m_el
                    Mm[3n_nos+2n_noi+i, :] += m_elm[1:3, :]
                    Mp[3n_nos+2n_noi+i, :] += m_elp[1:3, :]
                end

            end
        end
        # @show size(M)
        # @show length(M)
        aux = M1' / F
        aux = [aux; aux; aux; aux][:]'
        # aux = [aux aux]'
        # @infiltrate
        AM = aux .* DM
        AP = aux .* DP
        # Ax = aux .* Dx
        # Ay = aux .* Dy
        for i = 1:n_nos #Laço dos pontos radiais
            AM[4i-3:4i, 4i-3:4i] .= 0
            AM[4i-3:4i, 4i-3:4i] =
                -[sum(AM[4i-3:4i, 1:4:end], dims = 2) sum(AM[4i-3:4i, 2:4:end], dims = 2) sum(
                    AM[4i-3:4i, 3:4:end],
                    dims = 2,
                ) sum(AM[4i-3:4i, 4:4:end], dims = 2)] + Mm[4i-3:4i, :]
            AP[4i-3:4i, 4i-3:4i] .= 0
            AP[4i-3:4i, 4i-3:4i] =
                -[sum(AP[4i-3:4i, 1:4:end], dims = 2) sum(AP[4i-3:4i, 2:4:end], dims = 2) sum(
                    AP[4i-3:4i, 3:4:end],
                    dims = 2,
                ) sum(AP[4i-3:4i, 4:4:end], dims = 2)] + Mp[4i-3:4i, :]
        end
        AM, AP
    end

    function calc_md(x, pf, nf, qsi, w, elem, dad::casca, pre)
        npg = length(w)
        m_el, m_elm, m_elp = 0, zeros(4), zeros(4)
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
            r = pg' - pf
            R = norm(r)
            m = int_interpolaρdρ(R)
            mm, mp = calcula_F1(pf, nf, pg, [nx, ny], qsi, w, dad, pre)
            # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
            # @infiltrate
            m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
            m_elm += dot([nx, ny], r) / norm(r)^2 * mm * dgamadqsi * w[i]
            m_elp += dot([nx, ny], r) / norm(r)^2 * mp * dgamadqsi * w[i]
        end
        return m_el, m_el1
    end

    function calcula_F1(pf, nf, pg, n, qsi, w, dad::casca, pre) #
        npg = length(w)
        R = (pg' - pf)
        r = norm(R)
        drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
        #    para qsi (dr/dqsi)
        Fm, Fp = zeros(4), zeros(4) # Inicializa a integral de F_area
        theta = atan(R[2], R[1])
        for i = 1:npg # Percorre os pontos de integra��o
            ro = r / 2 * (qsi[i] + 1)
            termm, termp = compute_domainterms(ro, theta, nf, dad, pre)
            Fm = Fm + termm * ro * drodqsi * w[i]# Integral de F_area
            Fp = Fp + termp * ro * drodqsi * w[i]# Integral de F_area
        end
        return Fm, Fp
    end

    function calc_HeG(Casca::casca, npg = 8)
        dad = Casca.dadplaca
        nelem = size(dad.ELEM, 1)# Quantidade de elementos discretizados no contorno
        n_fis = size(dad.NOS, 1)
        n_internos = size(dad.pontos_internos, 1)
        n_cantos = size(dad.k.cantos, 1)
        H = zeros(4 * n_fis + n_cantos + 3n_internos, 4 * n_fis + n_cantos + 3n_internos)
        G = zeros(4 * n_fis + n_cantos + 3n_internos, 4 * n_fis + n_cantos + 3n_internos)
        q = zeros(4 * n_fis + n_cantos + 3n_internos)
        qsi, w = gausslegendre(npg)# Quadratura de gauss
        # qsi2, w2 = gausslegendre(2npg)# Quadratura de gauss
        normal_fonte = calc_normais(dad)
        pre = [zeros(2) for idx = 1:29]


        @showprogress "Montando H e G" for i = 1:n_fis+n_cantos+n_internos
            if i <= n_fis
                pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
                nf = normal_fonte[i, :]
                caso = "contorno"
            elseif i <= n_fis + n_internos
                pf = dad.pontos_internos[i-n_fis, :] # Coordenada (x,y)
                nf = zeros(2)
                caso = "interno"
            else
                pf = dad.k.cantos[i-n_fis-n_internos, 2:3] # Coordenada (x,y)
                nf = zeros(2)
                caso = "canto"
            end
            for j = 1:length(dad.ELEM)#Laço dos elementos
                elem_j = dad.ELEM[j]#Laço dos elementos
                x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
                Δelem = (x[end, :] - x[1, :]) # Δx e Δy entre o primeiro e ultimo nó geometrico
                eet =
                    (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                    norm(Δelem)^2 + elem_j.ξs[1]
                N_geo, dN = calc_fforma(eet, elem_j)
                ps = N_geo' * x
                b = norm(ps' - pf) / norm(Δelem)
                eta, Jt = sinhtrans(qsi, eet, b)
                # h, g = integraelem(pf, nf, x, qsi, w, elem_j, dad, pre)
                # @infiltrate
                h, g = integraelem(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                # @infiltrate
                # @show h
                # @infiltrate
                if caso == "contorno"
                    nosing = elem_j.indices .== i
                    if sum(nosing) == 1
                        no_pf = findfirst(nosing)
                        xi0 = elem_j.ξs[no_pf]
                        # h, g = integraelemsing(pf, nf, x, qsi2, w2, elem_j, dad, xi0)
                        h, g = integraelemsing(x, dad, xi0)
                    end
                    cols =
                        [4elem_j.indices .- 3 4elem_j.indices .- 2 4elem_j.indices .- 1 4elem_j.indices]'[:]
                    # @infiltrate
                    H[4*i-3:4*i, cols] = h
                    G[4*i-3:4*i, cols] = g
                    q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                    q[4*i-1:4*i] += q_el
                elseif caso == "canto"
                    if j == dad.k.cantos[i-n_fis-n_internos, 8]
                        xi0 = 1.0
                        h, g = integraelemsing(x, dad, xi0)

                    elseif j == dad.k.cantos[i-n_fis-n_internos, 9]
                        xi0 = -1.0
                        h, g = integraelemsing(x, dad, xi0)

                    end
                    cols =
                        [4elem_j.indices .- 3 4elem_j.indices .- 2 4elem_j.indices .- 1 4elem_j.indices]'[:]
                    ind = i + 3n_fis + 2n_noi
                    H[ind, cols] = h[1, :]
                    G[ind, cols] = g[1, :]
                    q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                    q[ind] += q_el[1]
                else
                    cols =
                        [4elem_j.indices .- 3 4elem_j.indices .- 2 4elem_j.indices .- 1 4elem_j.indices]'[:]
                    ind = (i - n_fis) * 3 + 4n_fis - 2
                    H[ind:ind+2, cols] = h[1, :]
                    G[ind:ind+2, cols] = g[1, :]
                    q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                    q[ind+2] += q_el[1]
                end
            end
            if caso == "contorno"
                R, W = compute_Rw(pf, nf, dad, pre)
                cols = (n_fis * 4 + 3n_internos) .+ (1:n_cantos)
                ind = i + 3n_fis + 2n_noi
                H[4*i-1:4*i, cols] = R
                G[4*i-1:4*i, cols] = W
            elseif caso == "canto"
                R, W = compute_Rw(pf, nf, dad, pre)
                ind = i + 3n_fis + 2n_noi
                cols = (n_fis * 4 + 3n_internos) .+ (1:n_cantos)
                H[ind, cols] = R[1, :]
                G[ind, cols] = W[1, :]
            else
                R, W = compute_Rw(pf, nf, dad, pre)
                ind = (i - n_fis) * 3 + 4n_fis - 2
                cols = (n_fis * 4 + 3n_internos) .+ (1:n_cantos)
                H[ind+2, cols] = R[1, :]
                G[ind+2, cols] = W[1, :]

            end
        end
        # for i = 1:n#i=1:size(dad.NOS,1) #Laço dos pontos fontes
        # H[2i-1:2i,2i-1:2i].=0
        # H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
        # end
        # @infiltrate
        # H[1:2n_fis, 1:2n_fis] += I / 2
        H[(4n_fis).+(1:3n_internos), (4n_fis).+(1:3n_internos)] += I
        # H[2n_fis+n_internos+1:end, 2n_fis+n_internos+1:end] += I / 4
        H, G, q
    end
end
