

function interpola(r; tipo = "tps")
    if tipo == "tps"
        if r == 0
            return 0
        end
        return r^2 * log(r)
    elseif tipo == "r"
        return r
    elseif tipo == "d4"
        C = 1
        return (2C - r) / (r + C)^4
    end
end
function int_interpolaρdρ(r; tipo = "tps")
    if tipo == "tps"
        if r == 0
            return 0
        end
        return (4 * r^4 * log(r) - r^4) / 16 # int r^3 * log(r)
    elseif tipo == "r"
        return r^3 / 3
    elseif tipo == "d4"
        C = 1
        return r^2 / (r + C)^3
    end
end
function int_interpola(r; tipo = "tps")
    if tipo == "tps"
        if r == 0
            return 0
        end
        return (1 / 3) * (r^3) * log(r) - (1 / 9) * (r^3)#  integrate(r^2 * log(r))
    elseif tipo == "r"
        return r^2 / 2
    elseif tipo == "d4"
        C = 1
        return (-C + r) / (2 * (C + r)^3)
    end
end

function Monta_M_RIMd(dad::Union{placa_fina,placa_fina_isotropica}, npg)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    n_canto = size(dad.k.cantos, 1)


    n_pontos = n_nos + n_noi + n_canto
    nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    F = zeros(n_pontos, n_pontos)
    D = zeros(n_pontos + n_nos, n_pontos)
    M1 = zeros(n_pontos)
    M2 = zeros(n_pontos + n_nos)
    normal_fonte = calc_normais(dad)
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
        else
            pf = dad.k.cantos[i-n_nos-n_noi, 2:3] # Coordenada (x,y)
            nf = zeros(2)
            caso = "canto"
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
            end
            if caso == "contorno"
                D[2i-1:2i, j] = compute_domainterms(R, atan(r[2], r[1]), nf, dad, pre)
            else
                D[i+n_nos, j] = compute_domainterms(R, atan(r[2], r[1]), nf, dad, pre)[1]
            end
        end
        qsi, w = gausslegendre(npg)

        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            m_el, m_el1 = calc_md(x, pf, nf, qsi, w, elem_j, dad, pre)
            # @infiltrate
            if caso == "contorno"
                # M1[2i-1] = M1[2i-1] + m_el
                # M1[2i] = M1[2i] + m_el
                M1[i] = M1[i] + m_el
                M2[2i-1:2i] = M2[2i-1:2i] + m_el1
            else
                M1[i] = M1[i] + m_el
                M2[i+n_nos] = M2[i+n_nos] + m_el1[1]
            end
        end
    end
    # @show size(M)
    # @show length(M)
    # @infiltrate
    A = M1' / F .* D
    for i = 1:n_nos #Laço dos pontos radiais
        A[2i-1:2i, i] .= 0
        A[2i-1:2i, i] = -sum(A[2i-1:2i, :], dims = 2) + M2[2i-1:2i]
    end
    for i = n_nos+1:n_pontos#Laço dos pontos radiais
        A[i+n_nos, i] = 0
        A[i+n_nos, i] = -sum(A[i+n_nos, :]) + M2[i+n_nos]
    end
    M = zeros(n_pontos + n_nos, n_pontos + n_nos)
    for i = 1:n_nos #Laço dos pontos radiais
        M[:, 2*i-1] = A[:, i]
    end
    for i = n_nos+1:n_pontos#Laço dos pontos radiais
        M[:, n_nos+i] = A[:, i]
    end
    M
end

function calc_md(x, pf, nf, qsi, w, elem, dad::Union{placa_fina,placa_fina_isotropica}, pre)
    npg = length(w)
    m_el, m_el1 = 0, zeros(2)
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
        m1 = calcula_F1(pf, nf, pg, [nx, ny], qsi, w, dad, pre)
        # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
        m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
    end
    return m_el, m_el1
end

function Monta_M_RIM(dad::Union{placa_fina,placa_fina_isotropica}, npg1 = 10, npg2 = 10)
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
        for x in [:R0, :S0, :dRdx0, :dRdy0, :dSdx0, :dSdy0]
            @eval $x = zeros(2)
        end
        pre = [zeros(2) for idx = 1:6]
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

function calcula_F1(pf, nf, pg, n, qsi, w, dad, pre) #
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(2) # Inicializa a integral de F_area
    theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        ro = r / 2 * (qsi[i] + 1)
        term = compute_domainterms(ro, theta, nf, dad, pre)
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

    C1 =
        ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[1])

    C2 =
        ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[2])

    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)

    a = 1


    for i = 1:2
        R[i] =
            r^2 *
            ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 3
            ) -
            4 *
            r^2 *
            dad.k.e[i] *
            sin(theta) *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        S[i] =
            r^2 *
            dad.k.e[i] *
            sin(theta) *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 3
            ) +
            r^2 *
            ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dRdx[i] =
            2 *
            r *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) -
            4 *
            r *
            dad.k.e[i] *
            sin(theta) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dRdy[i] =
            2 *
            r *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                dad.k.e[i]^2 * sin(theta)
            ) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) -
            4 *
            r *
            dad.k.e[i] *
            (cos(theta) + 2 * dad.k.d[i] * sin(theta)) *
            atan((dad.k.e[i] * sin(theta)), (cos(theta) + dad.k.d[i] * sin(theta)))

        dSdx[i] =
            r *
            dad.k.e[i] *
            sin(theta) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) +
            2 *
            r *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dSdy[i] =
            r *
            dad.k.e[i] *
            (cos(theta) + 2 * dad.k.d[i] * sin(theta)) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) +
            2 *
            r *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                dad.k.e[i]^2 * sin(theta)
            ) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
    end


    w = 1 / (8 * pi * dad.k.D22) * (C1 * R[1] + C2 * R[2] + C3 * (S[1] - S[2]))

    dwdx = 1 / (8 * pi) * (C1 * dRdx[1] + C2 * dRdx[2] + C3 * (dSdx[1] - dSdx[2]))
    dwdy = 1 / (8 * pi) * (C1 * dRdy[1] + C2 * dRdy[2] + C3 * (dSdy[1] - dSdy[2]))

    dwdm = -(dwdx * m1 + dwdy * m2) / dad.k.D22

    [w, dwdm]
end

function compute_domainterms(r, theta, nf, dad::Union{placa_fina_isotropica}, pre)
    m1 = nf[1]
    m2 = nf[2]
    D = dad.k.D

    # Distance of source and field points
    w = r^2 / (8 * pi * D) * (log(r) - 1 / 2)
    dwdx = -(r * cos(theta) * (-log(r^2))) / (8 * pi * D)
    dwdy = -(r * (-log(r^2)) * sin(theta)) / (8 * pi * D)

    dwdm = -(dwdx * m1 + dwdy * m2)

    [w, dwdm]
end
function Monta_M_RIMd_alt(dad, npg)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    n_canto = size(dad.k.cantos, 1)


    n_pontos = n_nos + n_noi + n_canto
    nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    F = zeros(n_pontos, n_pontos)
    D = zeros(n_pontos + n_nos, n_pontos)
    M1 = zeros(n_pontos)
    M2 = zeros(n_pontos + n_nos)
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
        else
            pf = dad.k.cantos[i-n_nos-n_noi, 2:3] # Coordenada (x,y)
            nf = zeros(2)
            caso = "canto"
        end
        for j = 1:n_pontos
            # @show i, j
            pr = nodes[j, :]
            r = pr - pf
            R = norm(r)
            if R == 0
                R = 1e-10
            end
            # @infiltrate
            F[i, j] = interpola(R)

            if caso == "contorno"
                D[2i-1:2i, j] = compute_domainterms(R, atan(r[2], r[1]), nf, dad, pre) .* R
            else
                D[i+n_nos, j] =
                    compute_domainterms(R, atan(r[2], r[1]), nf, dad, pre)[1] .* R
            end
        end
        qsi, w = gausslegendre(npg)

        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            m_el = calc_md_alt(x, pf, nf, qsi, w, elem_j, dad, pre)
            # @infiltrate
            if caso == "contorno"
                # M1[2i-1] = M1[2i-1] + m_el
                # M1[2i] = M1[2i] + m_el
                M1[i] = M1[i] + m_el
                # M2[2i-1:2i] = M2[2i-1:2i] + m_el1
            else
                M1[i] = M1[i] + m_el
                # M2[i+n_nos] = M2[i+n_nos] + m_el1[1]
            end
        end
    end
    # @show size(M)
    # @show length(M)
    @infiltrate
    # A = M1' / F .* D
    A = M1' .* (D / F)

    M = zeros(n_pontos + n_nos, n_pontos + n_nos)
    for i = 1:n_nos #Laço dos pontos radiais
        M[:, 2*i-1] = A[:, i]
    end
    for i = n_nos+1:n_pontos#Laço dos pontos radiais
        M[:, n_nos+i] = A[:, i]
    end
    M
end
function calc_md_alt(x, pf, nf, qsi, w, elem, dad, pre)
    npg = length(w)
    m_el, m_el1 = 0, zeros(2)

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
        m = int_interpola(R)
        # m1 = calcula_F1(pf, nf, pg, [nx, ny], qsi, w, dad, pre)
        # @infiltrate
        m_el += dot([nx, ny], r) / R^2 * m * dgamadqsi * w[i]
        # m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
    end
    return m_el
end


function DRM(dad, Hg, Gg)
    D = dad.k.D11
    n_nos = size(dad.NOS, 1)
    n_el = nelem = size(dad.ELEM, 1)# N. total de elementos
    n_cantos = size(dad.k.cantos, 1) # N. total de cantos
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    c1 = 1
    c2 = -1  #-1/(225*D)
    c3 = 0
    c4 = 0


    # Carregamento constante no domínio
    carga_dom = 1

    n_pontos = n_nos + n_noi + n_cantos
    nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]

    F = zeros(n_pontos, n_pontos)


    # Cálculo da matriz [F]
    for i = 1:n_pontos
        pf = nodes[i, :]
        for j = 1:n_pontos
            pr = nodes[j, :]
            R = pr - pf
            r = norm(R)
            if i == j
                rx = 0
                ry = 0

            else
                rx = R[1] / r
                ry = R[2] / r
            end
            drdx = rx
            drdy = ry

            # /@infiltrate
            d4wdx4 =
                3 * (
                    8 * c1 - 5 * c2 * (-3 - 6 * drdx^2 + drdx^4) * r +
                    24 * (c3 + 4 * c3 * drdx^2) * r^2 +
                    35 * c4 * (1 + 6 * drdx^2 + drdx^4) * r^3
                )

            d4wdx3dy =
                3 *
                drdx *
                drdy *
                r *
                (-5 * c2 * (-3 + drdx^2) + 48 * c3 * r + 35 * c4 * (3 + drdx^2) * r^2)

            d4wdx2dy2 =
                8 * c1 - 15 * c2 * (-1 - drdx^2 + (-1 + drdx^2) * drdy^2) * r +
                24 * c3 * (1 + 2 * drdx^2 + 2 * drdy^2) * r^2 +
                35 * c4 * (1 + 3 * drdx^2 + 3 * (1 + drdx^2) * drdy^2) * r^3

            d4wdxdy3 =
                3 *
                drdx *
                drdy *
                r *
                (-5 * c2 * (-3 + drdy^2) + 48 * c3 * r + 35 * c4 * (3 + drdy^2) * r^2)

            d4wdy4 =
                3 * (
                    8 * c1 - 5 * c2 * (-3 - 6 * drdy^2 + drdy^4) * r +
                    24 * (c3 + 4 * c3 * drdy^2) * r^2 +
                    35 * c4 * (1 + 6 * drdy^2 + drdy^4) * r^3
                )
            #A força F será representada pela equação diferencial da placa.
            F[i, j] =
                dad.k.D11 * d4wdx4 +
                4 * dad.k.D16 * d4wdx3dy +
                2 * (dad.k.D12 + 2 * dad.k.D66) * d4wdx2dy2 +
                4 * dad.k.D26 * d4wdxdy3 +
                dad.k.D22 * d4wdy4

        end

    end
    # @infiltrate

    # Calcula wc e Vnc
    # ================

    # Inicialização das matrizes H e G
    WC = zeros(2 * n_nos + n_cantos + n_noi, n_pontos)
    VnC = zeros(2 * n_nos + n_cantos + n_noi, n_pontos)
    normal_fonte = calc_normais(dad)

    # loop sobre os pontos de reciprocidade dual
    # formam as colunas da matriz WC
    # imaginar como se fossem pontos fontes

    @showprogress "Montando DR" for j = 1:n_pontos
        pf = nodes[j, :] # Coordenada (x,y)


        # loop sobre os elementos
        # formam as linhas da matriz WC
        # imaginar como se fossem os pontos campos()
        for i = 1:n_pontos
            if i <= n_nos
                pr = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
                nf = normal_fonte[i, :]
                caso = "contorno"
            elseif i <= n_nos + n_noi
                pr = dad.pontos_internos[i-n_nos, :] # Coordenada (x,y)
                nf = zeros(2)
                caso = "interno"
            else
                pr = dad.k.cantos[i-n_nos-n_noi, 2:3] # Coordenada (x,y)
                nf = zeros(2)
                caso = "canto"
            end
            R = pr - pf
            r = norm(R)
            if r == 0
                rx = 0
                ry = 0

            else
                rx = R[1] / r
                ry = R[2] / r
            end

            drdx = rx
            drdy = ry
            if caso == "contorno"

                nx = nf[1]
                ny = nf[2]

                ##############
                wc = c1 * r^4 + c2 * r^5 + c3 * r^6 + c4 * r^7
                ##############

                dwdx = drdx * r^3 * (4 * c1 + 5 * c2 * r + 6 * c3 * r^2 + 7 * c4 * r^3)
                dwdy = drdy * r^3 * (4 * c1 + 5 * c2 * r + 6 * c3 * r^2 + 7 * c4 * r^3)


                d2wdx2 =
                    r^2 * (
                        4 * (c1 + 2 * c1 * drdx^2) +
                        5 * (c2 + 3 * c2 * drdx^2) * r +
                        6 * (c3 + 4 * c3 * drdx^2) * r^2 +
                        7 * (c4 + 5 * c4 * drdx^2) * r^3
                    )
                d2wdy2 =
                    r^2 * (
                        4 * (c1 + 2 * c1 * drdy^2) +
                        5 * (c2 + 3 * c2 * drdy^2) * r +
                        6 * (c3 + 4 * c3 * drdy^2) * r^2 +
                        7 * (c4 + 5 * c4 * drdy^2) * r^3
                    )


                d2wdxdy =
                    drdx *
                    drdy *
                    r^2 *
                    (8 * c1 + 15 * c2 * r + 24 * c3 * r^2 + 35 * c4 * r^3)


                d3wdx3 =
                    3 *
                    drdx *
                    r *
                    (
                        8 * c1 +
                        5 * c2 * (3 + drdx^2) * r +
                        8 * c3 * (3 + 2 * drdx^2) * r^2 +
                        35 * c4 * (1 + drdx^2) * r^3
                    )
                d3wdx2dy =
                    drdy *
                    r *
                    (
                        8 * c1 +
                        15 * c2 * (1 + drdx^2) * r +
                        24 * (c3 + 2 * c3 * drdx^2) * r^2 +
                        35 * (c4 + 3 * c4 * drdx^2) * r^3
                    )
                d3wdxdy2 =
                    drdx *
                    r *
                    (
                        8 * c1 +
                        15 * c2 * (1 + drdy^2) * r +
                        24 * (c3 + 2 * c3 * drdy^2) * r^2 +
                        35 * (c4 + 3 * c4 * drdy^2) * r^3
                    )
                d3wdy3 =
                    3 *
                    drdy *
                    r *
                    (
                        8 * c1 +
                        5 * c2 * (3 + drdy^2) * r +
                        8 * c3 * (3 + 2 * drdy^2) * r^2 +
                        35 * c4 * (1 + drdy^2) * r^3
                    )

                f1 = dad.k.D11 * nx^2 + 2 * dad.k.D16 * nx * ny + dad.k.D12 * ny^2
                f2 = 2 * (dad.k.D16 * nx^2 + 2 * dad.k.D66 * nx * ny + dad.k.D26 * ny^2)
                f3 = dad.k.D12 * nx^2 + 2 * dad.k.D26 * nx * ny + dad.k.D22 * ny^2

                h1 =
                    dad.k.D11 * nx * (1 + ny^2) + 2 * dad.k.D16 * ny^3 -
                    dad.k.D12 * nx * ny^2
                h2 =
                    4 * dad.k.D16 * nx +
                    dad.k.D12 * ny * (1 + nx^2) +
                    4 * dad.k.D66 * ny^3 - dad.k.D11 * nx^2 * ny - 2 * dad.k.D26 * nx * ny^2
                h3 =
                    4 * dad.k.D26 * ny +
                    dad.k.D12 * nx * (1 + ny^2) +
                    4 * dad.k.D66 * nx^3 - dad.k.D22 * nx * ny^2 - 2 * dad.k.D16 * nx^2 * ny
                h4 =
                    dad.k.D22 * ny * (1 + nx^2) + 2 * dad.k.D26 * nx^3 -
                    dad.k.D12 * nx^2 * ny

                dwdn = (dwdx * nx + dwdy * ny) #rotação
                vn = -(h1 * d3wdx3 + h2 * d3wdx2dy + h3 * d3wdxdy2 + h4 * d3wdy3)
                mn = -(f1 * d2wdx2 + f2 * d2wdxdy + f3 * d2wdy2)

                #dwcdn
                dwcdn = dwdn

                #Vnc
                Vnc = vn

                #Mnc
                Mnc = mn
                WC[2i-1:2i, j] = [wc dwcdn]'
                VnC[2i-1:2i, j] = [Vnc Mnc]'
            elseif caso == "interno"
                wci = c1 * r^4 + c2 * r^5 + c3 * r^6 + c4 * r^7

                WC[n_nos+i, j] = wci
                VnC[n_nos+i, j] = 0
            else
                na1 = dad.k.cantos[i-n_noi-n_nos, 4]
                na2 = dad.k.cantos[i-n_noi-n_nos, 5]
                nd1 = dad.k.cantos[i-n_noi-n_nos, 6]
                nd2 = dad.k.cantos[i-n_noi-n_nos, 7]
                wci = c1 * r^4 + c2 * r^5 + c3 * r^6 + c4 * r^7



                d2wdx2 =
                    r^2 * (
                        4 * (c1 + 2 * c1 * drdx^2) +
                        5 * (c2 + 3 * c2 * drdx^2) * r +
                        6 * (c3 + 4 * c3 * drdx^2) * r^2 +
                        7 * (c4 + 5 * c4 * drdx^2) * r^3
                    )
                d2wdy2 =
                    r^2 * (
                        4 * (c1 + 2 * c1 * drdy^2) +
                        5 * (c2 + 3 * c2 * drdy^2) * r +
                        6 * (c3 + 4 * c3 * drdy^2) * r^2 +
                        7 * (c4 + 5 * c4 * drdy^2) * r^3
                    )


                d2wdxdy =
                    drdx *
                    drdy *
                    r^2 *
                    (8 * c1 + 15 * c2 * r + 24 * c3 * r^2 + 35 * c4 * r^3)


                g1a = (dad.k.D12 - dad.k.D11) * na1 * na2 + dad.k.D16 * (na1^2 - na2^2)
                g2a =
                    2 * (dad.k.D26 - dad.k.D16) * na1 * na2 +
                    2 * dad.k.D66 * (na1^2 - na2^2)
                g3a = (dad.k.D22 - dad.k.D12) * na1 * na2 + dad.k.D26 * (na1^2 - na2^2)


                g1d = (dad.k.D12 - dad.k.D11) * nd1 * nd2 + dad.k.D16 * (nd1^2 - nd2^2)
                g2d =
                    2 * (dad.k.D26 - dad.k.D16) * nd1 * nd2 +
                    2 * dad.k.D66 * (nd1^2 - nd2^2)
                g3d = (dad.k.D22 - dad.k.D12) * nd1 * nd2 + dad.k.D26 * (nd1^2 - nd2^2)


                tna = -(g1a * d2wdx2 + g2a * d2wdxdy + g3a * d2wdy2)
                tnd = -(g1d * d2wdx2 + g2d * d2wdxdy + g3d * d2wdy2)

                #   RciC
                RciC = (tnd - tna)
                # @infiltrate
                WC[n_nos+i, j] = wci
                VnC[n_nos+i, j] = RciC
            end
        end
    end
    # @infiltrate
    M = (Hg * WC - Gg * VnC) * inv(F)

    M2 = zeros(n_pontos + n_nos, n_pontos + n_nos)
    for i = 1:n_pontos
        if i <= n_nos
            M2[:, 2*i-1] = M[:, i]
        else
            M2[:, n_nos+i] = M[:, i]
        end
    end
    M2
end


function Monta_dM_RIM(dad::Union{placa_fina,placa_fina_isotropica}, npg1 = 10, npg2 = 10)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    n_canto = size(dad.k.cantos, 1)

    qsi1, w1 = gausslegendre(npg1)
    n_pontos = n_nos + n_noi + n_canto
    nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    Mx = zeros(2n_nos + n_noi + n_canto, n_pontos)
    My = zeros(2n_nos + n_noi + n_canto, n_pontos)
    Mxx = zeros(2n_nos + n_noi + n_canto, n_pontos)
    Myy = zeros(2n_nos + n_noi + n_canto, n_pontos)
    Mxy = zeros(2n_nos + n_noi + n_canto, n_pontos)
    dMx = zeros(2n_nos + 2n_noi, n_pontos)
    dMy = zeros(2n_nos + 2n_noi, n_pontos)
    # dMdtdm = zeros(2n_nos, n_pontos)
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
    @showprogress "Montando dMs" for i = 1:n_pontos
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

        pre = [zeros(2) for idx = 1:18]
        for j = 1:n_pontos #Laço dos pontos fontes
            pr = nodes[j, :]
            for el = 1:nelem
                elem_j = dad.ELEM[el]#Laço dos elementos
                x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos

                m_el, m_el3 = calc_dm(x, pf, nf, pr, qsi1, w1, elem_j, dad, pre)
                # @infiltrate
                if caso == "contorno"
                    Mx[2i-1:2i, j] += m_el[:, 1]
                    My[2i-1:2i, j] += m_el[:, 2]

                    Mxx[2i-1:2i, j] += m_el3[:, 1]
                    Myy[2i-1:2i, j] += m_el3[:, 2]
                    Mxy[2i-1:2i, j] += m_el3[:, 3]

                    dMx[2i-1:2i, j] += m_el[:, 3]
                    dMy[2i-1:2i, j] += m_el[:, 4]
                    # dMdtdm[2i-1:2i, j] += m_el[:, 5]
                elseif caso == "interno"
                    Mx[n_nos+i, j] += m_el[1, 1]
                    My[n_nos+i, j] += m_el[1, 2]

                    Mxx[n_nos+i, j] += m_el3[1, 1]
                    Myy[n_nos+i, j] += m_el3[1, 2]
                    Mxy[n_nos+i, j] += m_el3[1, 3]

                    dMx[2i-1:2i, j] += m_el[:, 3]
                    dMy[2i-1:2i, j] += m_el[:, 4]
                elseif caso == "canto"
                    Mx[n_nos+i, j] += m_el[1, 1]
                    My[n_nos+i, j] += m_el[1, 2]


                    Mxx[n_nos+i, j] += m_el3[1, 1]
                    Myy[n_nos+i, j] += m_el3[1, 2]
                    Mxy[n_nos+i, j] += m_el3[1, 3]
                end
            end
        end
    end

    # @infiltrate
    Mx = Mx / F
    My = My / F

    Mxx = Mxx / F
    Myy = Myy / F
    Mxy = Mxy / F

    dMx = dMx / F
    dMy = dMy / F
    # dMdtdm = dMdtdm / F

    lineM2 = 2 * n_nos + n_noi + n_canto
    columnM2 = 2 * n_nos + 2 * n_noi + 2 * n_canto
    linesdM2dtdm = 2 * n_nos
    M2 = zeros(lineM2, columnM2)
    M3 = zeros(2 * n_nos + n_noi + n_canto, 3 * n_nos + 3 * n_noi + 3 * n_canto)
    dM2dtdm = zeros(linesdM2dtdm, lineM2)
    dM2dx2 = zeros(n_noi, lineM2)
    dM2dxdy = zeros(n_noi, lineM2)
    dM2dy2 = zeros(n_noi, lineM2)
    for i = 1:n_pontos
        M2[:, 2*i-1] = Mx[:, i]
        M2[:, 2*i] = My[:, i]

        M3[:, 3*i-2] = Mxx[:, i]
        M3[:, 3*i-1] = Myy[:, i]
        M3[:, 3*i] = 2Mxy[:, i]
        if i <= n_nos
            # dM2dtdm[:, 2*i-1] = dMdtdm[1:linesdM2dtdm, i]
            dM2dtdm[1:2:end, 2*i-1] = -(
                Mx[1:2:2n_nos, i] .* -normal_fonte[:, 2] +
                My[1:2:2n_nos, i] .* normal_fonte[:, 1]
            )
            dM2dtdm[2:2:end, 2*i-1] = -(
                Mx[1:2:2n_nos, i] .* normal_fonte[:, 1] +
                My[1:2:2n_nos, i] .* normal_fonte[:, 2]
            )
            dM2dx2[1:n_noi, 2*i-1] = dMx[2*n_nos+1:2:2*(n_nos+n_noi)-1, i]
            dM2dxdy[1:n_noi, 2*i-1] = dMx[2*n_nos+2:2:2*(n_nos+n_noi), i]
            dM2dy2[1:n_noi, 2*i-1] = dMy[2*n_nos+2:2:2*(n_nos+n_noi), i]
            # @infiltrate
        else
            # dM2dtdm[:, n_nos+i] = dMdtdm[1:linesdM2dtdm, i]
            dM2dtdm[1:2:end, n_nos+i] = -(
                Mx[1:2:2n_nos, i] .* -normal_fonte[:, 2] +
                My[1:2:2n_nos, i] .* normal_fonte[:, 1]
            )
            dM2dtdm[2:2:end, n_nos+i] = -(
                Mx[1:2:2n_nos, i] .* normal_fonte[:, 1] +
                My[1:2:2n_nos, i] .* normal_fonte[:, 2]
            )
            dM2dx2[1:n_noi, n_nos+i] = dMx[2*n_nos+1:2:2*(n_nos+n_noi)-1, i]
            dM2dxdy[1:n_noi, n_nos+i] = dMx[2*n_nos+2:2:2*(n_nos+n_noi), i]
            dM2dy2[1:n_noi, n_nos+i] = dMy[2*n_nos+2:2:2*(n_nos+n_noi), i]
        end
    end

    M2, dM2dx2, dM2dy2, dM2dxdy, dM2dtdm, M3
end


function calc_dm(x, pf, nf, pr, qsi1, w1, elem, dad, pre)

    npg = length(w1)
    m_el = zeros(2, 4)
    m_el3 = zeros(2, 3)


    for i = 1:npg
        N, dN = calc_fforma(qsi1[i], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        ny = -dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        nx = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        r = pg' - pf
        m, m3 = calcula_dF(pr, pf, nf, pg, [nx, ny], qsi1, w1, dad, pre)
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w1[i]
        m_el3 += dot([nx, ny], r) / norm(r)^2 * m3 * dgamadqsi * w1[i]
    end
    return m_el, m_el3
end

function calcula_dF(pr, pf, nf, pg, n, qsi, w, dad, pre) #???????????????/
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(2, 4) # Inicializa a integral de F_area
    F3 = zeros(2, 3) # Inicializa a integral de F_area
    theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
        rline = norm(xc - pr)
        ro = r / 2 * (qsi[i] + 1)
        term, term3 = compute_domain_d_terms(ro, theta, nf, n, dad, pre)
        f = interpola(rline)
        F = F + term * f * ro * drodqsi * w[i]# Integral de F_area
        F3 = F3 + term3 * f * ro * drodqsi * w[i]# Integral de F_area
    end
    # @infiltrate
    return F, F3
end

function compute_domain_d_terms(r, theta, nf, n, dad::Union{placa_fina}, pre)


    m1 = nf[1]
    m2 = nf[2]
    t1 = -nf[2]
    t2 = nf[1]   # Distance of source and field points
    r1 = r .* cos(theta)
    r2 = r .* sin(theta)
    # @infiltrate
    dRdx,
    dRdy,
    dSdx,
    dSdy,
    d2Rdx2,
    d2Rdxdy,
    d2Rdy2,
    d2Sdx2,
    d2Sdxdy,
    d2Sdy2,
    d3Rdx2dy,
    d3Rdx3,
    d3Rdxdy2,
    d3Rdy3,
    d3Sdx2dy,
    d3Sdx3,
    d3Sdxdy2,
    d3Sdy3 = pre

    G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2

    H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2

    C1 =
        ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[1])

    C2 =
        ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[2])

    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)

    a = 1


    for i = 1:2
        # R[i] = r^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) * (log(r^2 / a^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)) - 3) - 4 * r^2 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        # S[i] = r^2 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) * (log(r^2 / a^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)) - 3) + r^2 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dRdx[i] =
            2 *
            r *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) -
            4 *
            r *
            dad.k.e[i] *
            sin(theta) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dRdy[i] =
            2 *
            r *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                dad.k.e[i]^2 * sin(theta)
            ) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) -
            4 *
            r *
            dad.k.e[i] *
            (cos(theta) + 2 * dad.k.d[i] * sin(theta)) *
            atan((dad.k.e[i] * sin(theta)), (cos(theta) + dad.k.d[i] * sin(theta)))

        dSdx[i] =
            r *
            dad.k.e[i] *
            sin(theta) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) +
            2 *
            r *
            (cos(theta) + dad.k.d[i] * sin(theta)) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        dSdy[i] =
            r *
            dad.k.e[i] *
            (cos(theta) + 2 * dad.k.d[i] * sin(theta)) *
            (
                log(
                    r^2 / a^2 * (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ),
                ) - 2
            ) +
            2 *
            r *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                dad.k.e[i]^2 * sin(theta)
            ) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))


        d2Rdx2[i] =
            2 * log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            )

        d2Rdxdy[i] =
            2 *
            dad.k.d[i] *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) -
            4 *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Rdy2[i] =
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) -
            8 *
            dad.k.d[i] *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdx2[i] =
            2 * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdxdy[i] =
            dad.k.e[i] * (log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            )) +
            2 *
            dad.k.d[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdy2[i] =
            2 *
            dad.k.d[i] *
            dad.k.e[i] *
            log(
                r .^ 2 / a^2 .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                ),
            ) +
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d3Rdx3[i] =
            4 * (cos(theta) + dad.k.d[i] * sin(theta)) ./ (
                r .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                )
            )

        d3Rdx2dy[i] =
            4 * (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) ./ (
                r .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                )
            )

        d3Rdxdy2[i] =
            4 * (
                (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
            ) ./ (
                r .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                )
            )

        d3Rdy3[i] =
            4 * (
                dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
            ) ./ (
                r .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                )
            )

        d3Sdx3[i] =
            -2 * dad.k.e[i] * sin(theta) ./ (
                r .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                )
            )


        d3Sdx2dy[i] =
            2 * dad.k.e[i] * cos(theta) ./ (
                r .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                )
            )

        d3Sdxdy2[i] =
            2 *
            dad.k.e[i] *
            (
                2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
            ) ./ (
                r .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                )
            )

        d3Sdy3[i] =
            2 *
            dad.k.e[i] *
            (
                (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
            ) ./ (
                r .* (
                    (cos(theta) + dad.k.d[i] * sin(theta)) .^ 2 +
                    dad.k.e[i]^2 * sin(theta) .^ 2
                )
            )
    end



    dwdx =
        1 / (8 * pi * dad.k.D22) * (C1 * dRdx[1] + C2 * dRdx[2] + C3 * (dSdx[1] - dSdx[2]))
    dwdy =
        1 / (8 * pi * dad.k.D22) * (C1 * dRdy[1] + C2 * dRdy[2] + C3 * (dSdy[1] - dSdy[2]))

    d2wdx2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
    d2wdxdy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
    d2wdy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))

    d2wdmdx = -(d2wdx2 * m1 + d2wdxdy * m2)
    d2wdmdy = -(d2wdxdy * m1 + d2wdy2 * m2)
    d3wdx3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdx3[1] + C2 * d3Rdx3[2] + C3 * (d3Sdx3[1] - d3Sdx3[2]))
    d3wdx2dy =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdx2dy[1] + C2 * d3Rdx2dy[2] + C3 * (d3Sdx2dy[1] - d3Sdx2dy[2]))
    d3wdxdy2 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdxdy2[1] + C2 * d3Rdxdy2[2] + C3 * (d3Sdxdy2[1] - d3Sdxdy2[2]))
    d3wdy3 =
        1 / (8 * pi * dad.k.D22) *
        (C1 * d3Rdy3[1] + C2 * d3Rdy3[2] + C3 * (d3Sdy3[1] - d3Sdy3[2]))


    d3wdx2dm = -(d3wdx3 * m1 + d3wdx2dy * m2)
    d3wdy2dm = -(d3wdxdy2 * m1 + d3wdy3 * m2)
    d3wdxdydm = -(d3wdx2dy * m1 + d3wdxdy2 * m2)
    # dwdt = -(dwdx * t1 + dwdy * t2)
    # dwdm = -(dwdx * m1 + dwdy * m2)
    # @infiltrate
    [dwdx dwdy d2wdx2 d2wdxdy; d2wdmdx d2wdmdy d2wdxdy d2wdy2],
    [d2wdx2 d2wdy2 d2wdxdy; d3wdx2dm d3wdy2dm d3wdxdydm]
end


function Monta_dM_RIMd(dad::Union{placa_fina,placa_fina_isotropica}, npg)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    n_canto = size(dad.k.cantos, 1)


    n_pontos = n_nos + n_noi + n_canto
    nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    F = zeros(n_pontos, n_pontos)
    # D = zeros(n_pontos + n_nos, n_pontos)

    Dx = zeros(2n_nos + n_noi + n_canto, n_pontos)
    Dy = zeros(2n_nos + n_noi + n_canto, n_pontos)
    Dxx = zeros(2n_nos + n_noi + n_canto, n_pontos)
    Dyy = zeros(2n_nos + n_noi + n_canto, n_pontos)
    Dxy = zeros(2n_nos + n_noi + n_canto, n_pontos)

    dDx = zeros(2n_nos + 2n_noi, n_pontos)
    dDy = zeros(2n_nos + 2n_noi, n_pontos)
    # dDdtdm = zeros(2n_nos, n_pontos)

    M1 = zeros(n_pontos)
    Mx = zeros(n_pontos + n_nos)
    My = zeros(n_pontos + n_nos)
    Mxx = zeros(n_pontos + n_nos)
    Myy = zeros(n_pontos + n_nos)
    Mxy = zeros(n_pontos + n_nos)
    dMx = zeros(2n_nos + 2n_noi)
    dMy = zeros(2n_nos + 2n_noi)
    # dMdtdm = zeros(n_pontos + n_nos)

    normal_fonte = calc_normais(dad)

    pre = [zeros(2) for idx = 1:20]

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
        else
            pf = dad.k.cantos[i-n_nos-n_noi, 2:3] # Coordenada (x,y)
            nf = zeros(2)
            caso = "canto"
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
            end
            if j > n_nos
                nc = zeros(2)
            else
                nc = normal_fonte[j, :]
            end
            d, d3 = compute_domain_d_terms(R, atan(r[2], r[1]), nf, nc, dad, pre)
            if caso == "contorno"
                Dx[2i-1:2i, j] += d[:, 1]
                Dy[2i-1:2i, j] += d[:, 2]

                Dxx[2i-1:2i, j] += d3[:, 1]
                Dyy[2i-1:2i, j] += d3[:, 2]
                Dxy[2i-1:2i, j] += d3[:, 3]

                dDx[2i-1:2i, j] += d[:, 3]
                dDy[2i-1:2i, j] += d[:, 4]
                # dDdtdm[2i-1:2i, j] += d[:, 5]
                # dDdtdm[2i-1, j] += -(-d[1, 1] * nf[2] + d[1, 2] * nf[1])
                # dDdtdm[2i, j] += -(d[1, 1] * nf[1] + d[1, 2] * nf[2])

            elseif caso == "interno"
                Dx[n_nos+i, j] += d[1, 1]
                Dy[n_nos+i, j] += d[1, 2]
                Dxx[n_nos+i, j] += d3[1, 1]
                Dyy[n_nos+i, j] += d3[1, 2]
                Dxy[n_nos+i, j] += d3[1, 3]

                dDx[2i-1:2i, j] += d[:, 3]
                dDy[2i-1:2i, j] += d[:, 4]
            elseif caso == "canto"
                Dx[n_nos+i, j] += d[1, 1]
                Dy[n_nos+i, j] += d[1, 2]
                Dxx[n_nos+i, j] += d3[1, 1]
                Dyy[n_nos+i, j] += d3[1, 2]
                Dxy[n_nos+i, j] += d3[1, 3]
            end
        end
        qsi, w = gausslegendre(npg)

        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            m_el, m_el1, m_el3 = calc_dmd(x, pf, nf, qsi, w, elem_j, dad, pre)
            # @infiltrate
            if caso == "contorno"
                M1[i] = M1[i] + m_el
                Mx[2i-1:2i] += m_el1[:, 1]
                My[2i-1:2i] += m_el1[:, 2]
                dMx[2i-1:2i] += m_el1[:, 3]
                dMy[2i-1:2i] += m_el1[:, 4]

                Mxx[2i-1:2i] += m_el3[:, 1]
                Myy[2i-1:2i] += m_el3[:, 2]
                Mxy[2i-1:2i] += m_el3[:, 3]
                # dMdtdm[2i-1:2i] += m_el1[:, 5]
            elseif caso == "interno"
                M1[i] = M1[i] + m_el
                Mx[n_nos+i] += m_el1[1, 1]
                My[n_nos+i] += m_el1[1, 2]
                # dMx[n_nos+i] += m_el1[:, 3]
                # dMy[n_nos+i] += m_el1[:, 4]
                dMx[2i-1:2i] += m_el1[:, 3]
                dMy[2i-1:2i] += m_el1[:, 4]

                Mxx[n_nos+i] += m_el3[1, 1]
                Myy[n_nos+i] += m_el3[1, 2]
                Mxy[n_nos+i] += m_el3[1, 3]
                # dMdtdm[n_nos+i] += m_el1[1, 5]
            elseif caso == "canto"
                M1[i] = M1[i] + m_el
                Mx[n_nos+i] += m_el1[1, 1]
                My[n_nos+i] += m_el1[1, 2]

                Mxx[n_nos+i] += m_el3[1, 1]
                Myy[n_nos+i] += m_el3[1, 2]
                Mxy[n_nos+i] += m_el3[1, 3]
                # dMdtdm[n_nos+i] += m_el1[1, 5]
            end
            #     M1[i] = M1[i] + m_el
            #     M2[2i-1:2i] = M2[2i-1:2i] + m_el1
            # else
            #     M1[i] = M1[i] + m_el
            #     M2[i+n_nos] = M2[i+n_nos] + m_el1[1]
            # end
        end
    end
    # @show size(M)
    # @show length(M)
    # @infiltrate
    # A = M1' / F .* D
    aux = M1' / F

    Ax = aux .* Dx
    Ay = aux .* Dy
    Axx = aux .* Dxx
    Ayy = aux .* Dyy
    Axy = aux .* Dxy
    dAx = aux .* dDx
    dAy = aux .* dDy
    # dAdtdm = M1' / F .* dDdtdm
    for i = 1:n_nos #Laço dos pontos radiais
        Ax[2i-1:2i, i] .= 0
        Ax[2i-1:2i, i] = -sum(Ax[2i-1:2i, :], dims = 2) + Mx[2i-1:2i]

        Ay[2i-1:2i, i] .= 0
        Ay[2i-1:2i, i] = -sum(Ay[2i-1:2i, :], dims = 2) + My[2i-1:2i]

        Axx[2i-1:2i, i] .= 0
        Axx[2i-1:2i, i] = -sum(Axx[2i-1:2i, :], dims = 2) + Mxx[2i-1:2i]

        Ayy[2i-1:2i, i] .= 0
        Ayy[2i-1:2i, i] = -sum(Ayy[2i-1:2i, :], dims = 2) + Myy[2i-1:2i]

        Axy[2i-1:2i, i] .= 0
        Axy[2i-1:2i, i] = -sum(Axy[2i-1:2i, :], dims = 2) + Mxy[2i-1:2i]

        dAx[2i-1:2i, i] .= 0
        dAx[2i-1:2i, i] = -sum(dAx[2i-1:2i, :], dims = 2) + dMx[2i-1:2i]

        dAy[2i-1:2i, i] .= 0
        dAy[2i-1:2i, i] = -sum(dAy[2i-1:2i, :], dims = 2) + dMy[2i-1:2i]

        # dAdtdm[2i-1:2i, i] .= 0
        # dAdtdm[2i-1:2i, i] = -sum(dAdtdm[2i-1:2i, :], dims=2) + dMdtdm[2i-1:2i]
    end
    for i = n_nos+1:n_nos+n_noi#Laço dos pontos radiais
        Ax[i+n_nos, i] = 0
        # @infiltrate
        Ax[i+n_nos, i] = -sum(Ax[i+n_nos, :]) + Mx[i+n_nos]

        Ay[i+n_nos, i] = 0
        Ay[i+n_nos, i] = -sum(Ay[i+n_nos, :]) + My[i+n_nos]

        Axx[i+n_nos, i] = 0
        # @infiltrate
        Axx[i+n_nos, i] = -sum(Axx[i+n_nos, :]) + Mxx[i+n_nos]

        Ayy[i+n_nos, i] = 0
        Ayy[i+n_nos, i] = -sum(Ayy[i+n_nos, :]) + Myy[i+n_nos]

        Axy[i+n_nos, i] = 0
        Axy[i+n_nos, i] = -sum(Axy[i+n_nos, :]) + Mxy[i+n_nos]


        dAx[2i-1:2i, i] .= 0
        dAx[2i-1:2i, i] = -sum(dAx[2i-1:2i, :], dims = 2) + dMx[2i-1:2i]

        dAy[2i-1:2i, i] .= 0
        dAy[2i-1:2i, i] = -sum(dAy[2i-1:2i, :], dims = 2) + dMy[2i-1:2i]


    end

    for i = n_nos+n_noi+1:n_pontos#Laço dos pontos radiais
        Ax[i+n_nos, i] = 0
        # @infiltrate
        Ax[i+n_nos, i] = -sum(Ax[i+n_nos, :]) + Mx[i+n_nos]

        Ay[i+n_nos, i] = 0
        Ay[i+n_nos, i] = -sum(Ay[i+n_nos, :]) + My[i+n_nos]

        Axx[i+n_nos, i] = 0
        # @infiltrate
        Axx[i+n_nos, i] = -sum(Axx[i+n_nos, :]) + Mxx[i+n_nos]

        Ayy[i+n_nos, i] = 0
        Ayy[i+n_nos, i] = -sum(Ayy[i+n_nos, :]) + Myy[i+n_nos]

        Axy[i+n_nos, i] = 0
        Axy[i+n_nos, i] = -sum(Axy[i+n_nos, :]) + Mxy[i+n_nos]
    end

    M2 = zeros(2 * n_nos + n_noi + n_canto, 2 * n_nos + 2 * n_noi + 2 * n_canto)
    dM2dtdm = zeros(2n_nos, 2 * n_nos + n_noi + n_canto)
    dM2dx2 = zeros(n_noi, 2 * n_nos + n_noi + n_canto)
    dM2dxdy = zeros(n_noi, 2 * n_nos + n_noi + n_canto)
    dM2dy2 = zeros(n_noi, 2 * n_nos + n_noi + n_canto)
    M3 = zeros(2 * n_nos + n_noi + n_canto, 3 * (n_nos + n_noi + n_canto))

    for i = 1:n_pontos
        M2[:, 2*i-1] = Ax[:, i]
        M2[:, 2*i] = Ay[:, i]
        M3[:, 3*i-2] = Axx[:, i]
        M3[:, 3*i-1] = Ayy[:, i]
        M3[:, 3*i] = 2Axy[:, i]
        if i <= n_nos
            dM2dtdm[1:2:end, 2*i-1] = -(
                Ax[1:2:2n_nos, i] .* -normal_fonte[:, 2] +
                Ay[1:2:2n_nos, i] .* normal_fonte[:, 1]
            )
            dM2dtdm[2:2:end, 2*i-1] = -(
                Ax[1:2:2n_nos, i] .* normal_fonte[:, 1] +
                Ay[1:2:2n_nos, i] .* normal_fonte[:, 2]
            )
            dM2dx2[1:n_noi, 2*i-1] = dAx[2*n_nos+1:2:2*(n_nos+n_noi)-1, i]
            dM2dxdy[1:n_noi, 2*i-1] = dAx[2*n_nos+2:2:2*(n_nos+n_noi), i]
            dM2dy2[1:n_noi, 2*i-1] = dAy[2*n_nos+2:2:2*(n_nos+n_noi), i]
        else
            dM2dtdm[1:2:end, n_nos+i] = -(
                Ax[1:2:2n_nos, i] .* -normal_fonte[:, 2] +
                Ay[1:2:2n_nos, i] .* normal_fonte[:, 1]
            )
            dM2dtdm[2:2:end, n_nos+i] = -(
                Ax[1:2:2n_nos, i] .* normal_fonte[:, 1] +
                Ay[1:2:2n_nos, i] .* normal_fonte[:, 2]
            )
            dM2dx2[1:n_noi, n_nos+i] = dAx[2*n_nos+1:2:2*(n_nos+n_noi)-1, i]
            dM2dxdy[1:n_noi, n_nos+i] = dAx[2*n_nos+2:2:2*(n_nos+n_noi), i]
            dM2dy2[1:n_noi, n_nos+i] = dAy[2*n_nos+2:2:2*(n_nos+n_noi), i]
        end
    end
    M2, dM2dx2, dM2dy2, dM2dxdy, dM2dtdm, M3
end


function calc_dmd(x, pf, nf, qsi, w, elem, dad, pre)
    npg = length(w)
    m_el, m_el1 = 0, zeros(2, 4)
    m_el3 = zeros(2, 3)

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
        m1, m3 = calcula_dF1(pf, nf, pg, [nx, ny], qsi, w, dad, pre)
        # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
        m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
        m_el3 += dot([nx, ny], r) / norm(r)^2 * m3 * dgamadqsi * w[i]
    end
    return m_el, m_el1, m_el3
end

function calcula_dF1(pf, nf, pg, n, qsi, w, dad, pre) #
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(2, 4) # Inicializa a integral de F_area
    F3 = zeros(2, 3) # Inicializa a integral de F_area
    theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        ro = r / 2 * (qsi[i] + 1)
        term, term3 = compute_domain_d_terms(ro, theta, nf, n, dad, pre)
        # @infiltrate
        F = F + term * ro * drodqsi * w[i]# Integral de F_area
        F3 = F3 + term3 * ro * drodqsi * w[i]# Integral de F_area
    end
    return F, F3
end




function compute_domain_d_terms(r, theta, nf, n, dad::placa_fina_isotropica, pre)

    D = dad.k.D
    ni = dad.k.nu

    m1 = nf[1]
    m2 = nf[2]
    t1 = -nf[2]
    t2 = nf[1]   # Distance of source and field points
    r1 = r .* cos(theta)
    r2 = r .* sin(theta)


    dwdx = -(r * cos(theta) * (-log(r^2))) / (8 * pi * D)
    dwdy = -(r * (-log(r^2)) * sin(theta)) / (8 * pi * D)

    d2wdx2 = (1 + cos(2 * theta) + log(r .^ 2)) / (8 * pi * D)
    d2wdxdy = sin(2 * theta) / (8 * pi * D)
    d2wdy2 = -(-1 + cos(2 * theta) - log(r .^ 2)) / (8 * pi * D)

    d3wdx3 = 1 / (4 * pi * r * D) .* cos(theta) .* (3 - 2 * cos(theta) .^ 2)
    d3wdx2dy = 1 / (4 * pi * r * D) .* sin(theta) .* (1 - 2 * cos(theta) .^ 2)
    d3wdxdy2 = 1 / (4 * pi * r * D) .* cos(theta) .* (1 - 2 * sin(theta) .^ 2)
    d3wdy3 = 1 / (4 * pi * r * D) .* sin(theta) .* (3 - 2 * sin(theta) .^ 2)

    d2wdmdx = -(d2wdx2 * m1 + d2wdxdy * m2)
    d2wdmdy = -(d2wdxdy * m1 + d2wdy2 * m2)

    d3wdx2dm = -(d3wdx3 * m1 + d3wdx2dy * m2)
    d3wdy2dm = -(d3wdxdy2 * m1 + d3wdy3 * m2)
    d3wdxdydm = -(d3wdx2dy * m1 + d3wdxdy2 * m2)
    # dwdt = -(dwdx * t1 + dwdy * t2)
    # dwdm = -(dwdx * m1 + dwdy * m2)
    # @infiltrate
    [dwdx dwdy d2wdx2 d2wdxdy; d2wdmdx d2wdmdy d2wdxdy d2wdy2],
    [d2wdx2 d2wdy2 d2wdxdy; d3wdx2dm d3wdy2dm d3wdxdydm]
end
