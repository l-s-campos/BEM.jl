function calc_HeG(dad::Union{placa_espessa,placa_espessa_isotropica}, npg = 8)
    nelem = size(dad.ELEM, 1)# Quantidade de elementos discretizados no contorno
    n_fis = size(dad.NOS, 1)
    n_internos = size(dad.pontos_internos, 1)
    H = zeros(3 * n_fis + 3 * n_internos, 3 * n_fis + 3 * n_internos)
    G = zeros(3 * n_fis + 3 * n_internos, 3 * n_fis + 3 * n_internos)
    q = zeros(3 * n_fis + 3 * n_internos)
    qsi, w = gausslegendre(npg)# Quadratura de gauss
    # qsi2, w2 = gausslegendre(2npg)# Quadratura de gauss
    normal_fonte = calc_normais(dad)

    @showprogress "Montando H e G" for i = 1:n_fis+n_internos
        if i <= n_fis
            pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
            nf = normal_fonte[i, :]
            caso = "contorno"
        elseif i <= n_fis + n_internos
            pf = dad.pontos_internos[i-n_fis, :] # Coordenada (x,y)
            nf = zeros(2)
            caso = "interno"
        end
        Cij = zeros(3, 3)
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

            # h, g, Cij = integraelem(pf, x, qsi, w, elem_j, dad)
            # @infiltrate
            h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
            # @infiltrate
            # @show h
            # @infiltrate
            if caso == "contorno"
                nosing = elem_j.indices .== i
                if sum(nosing) == 1
                    no_pf = findfirst(nosing)
                    xi0 = elem_j.ξs[no_pf]
                    h, g = integraelemsing_num(pf, x, elem_j, dad, xi0, 100)
                    # g[2, 2:2:end] = g2
                    # @show [hn h], [gm g]
                end
            end
            cols = [3elem_j.indices .- 2 3elem_j.indices .- 1 3elem_j.indices]'[:]
            # @infiltrate
            H[3i-2:3i, cols] = h
            G[3i-2:3i, cols] = g
            q_el = compute_q(pf, x, qsi, w, elem_j, dad)
            q[3i-2:3i] += q_el
            # H[3i-2:3i, 3i-2:3i] -= Cij
        end

    end
    H[1:3n_fis, 1:3n_fis] += I / 2
    H[(3n_fis).+(1:3n_internos), (3n_fis).+(1:3n_internos)] += I
    H, G, q
    # calsolfund([1,1],[0,0],[0,1],dad)
end
function calsolfund(pg, pf, n, dad::placa_espessa_isotropica)

    rs = pg - pf
    # Distance of source and field points
    RX = rs[1]
    RY = rs[2]
    R = norm(rs)
    v = dad.k.nu
    hpl = dad.k.h
    k1 = 128 * π * dad.k.D
    k2 = 8 * π * dad.k.D
    Lamd = sqrt(10) / hpl  # hpl = espessura da placa
    # Lamd = dad.k.λ  # hpl = espessura da placa
    z = Lamd * R

    # K0 = BEM.besselk(0, z)
    # K1 = BEM.besselk(1, z)
    # Az = K0 + 2 / z * (K1 - 1 / z)
    # Bz = K0 + 1 / z * (K1 - 1 / z)
    Az, Bz, K1z, K0 = Bess(z)
    # ----- CALCULA DERIVADAS
    DRx = RX / R
    DRy = RY / R
    DRn = DRx * n[1] + DRy * n[2]

    # ----- SOLUCION FUNDAMENTAL U*
    U11 =
        1 / (k2 * (1 - v)) *
        (8 * Bz - (1 - v) * (2 * log(z) - 1) - (8 * Az + 2 * (1 - v)) * DRx^2)
    U12 = -1 / (k2 * (1 - v)) * (8 * Az + 2 * (1 - v)) * DRx * DRy
    U13 = 1 / k2 * (2 * log(z) - 1) * R * DRx
    U21 = U12
    U22 =
        1 / (k2 * (1 - v)) *
        (8 * Bz - (1 - v) * (2 * log(z) - 1) - (8 * Az + 2 * (1 - v)) * DRy^2)
    U23 = 1 / k2 * (2 * log(z) - 1) * R * DRy
    U31 = -1 / k2 * (2 * log(z) - 1) * R * DRx
    U32 = -1 / k2 * (2 * log(z) - 1) * R * DRy
    U33 = 1 / (k2 * (1 - v) * Lamd^2) * ((1 - v) * z^2 * (log(z) - 1) - 8 * log(z))

    # ----- SOLUCION FUNDAMENTAL P*
    AA = 4 * Az + 2 * z * K1z + 1 - v
    BB = 4 * Az + 1 + v
    CC = 2 * (8 * Az + 2 * z * K1z + 1 - v)
    DD = Lamd^2 / (2 * π)
    GG = -(1 - v) / (8 * π)
    FF = 2 * (1 + v) / (1 - v)
    k3 = -1 / (4 * π * R)
    P11 = k3 * (AA * (DRn + DRx * n[1]) + BB * DRx * n[1] - CC * DRx^2 * DRn)
    P12 = k3 * (AA * DRy * n[1] + BB * DRx * n[2] - CC * DRy * DRx * DRn)
    P21 = k3 * (AA * DRx * n[2] + BB * DRy * n[1] - CC * DRx * DRy * DRn)
    P22 = k3 * (AA * (DRn + DRy * n[2]) + BB * DRy * n[2] - CC * DRy^2 * DRn)
    P13 = DD * (Bz * n[1] - Az * DRx * DRn)
    P23 = DD * (Bz * n[2] - Az * DRy * DRn)
    P31 = GG * ((FF * log(z) - 1) * n[1] + 2 * DRx * DRn)
    P32 = GG * ((FF * log(z) - 1) * n[2] + 2 * DRy * DRn)
    P33 = -1 / (2 * π * R) * DRn

    # ----- ARGUMENTO DE LAS INTEGRALES He y Ge
    Ge = [U11 U12 U13; U21 U22 U23; U31 U32 U33]
    He = [P11 P12 P13; P21 P22 P23; P31 P32 P33]

    # # ----- ARGUMENTO DE LA INTEGRAL Cij
    # HS11 = P11 - RX * P13
    # HS12 = P12 - RY * P13
    # HS21 = P21 - RX * P23
    # HS22 = P22 - RY * P23
    # HS31 = P31 - RX * P33
    # HS32 = P32 - RY * P33
    # HS13 = P13
    # HS23 = P23
    # HS33 = P33
    # FCij = [HS11 HS12 HS13; HS21 HS22 HS23; HS31 HS32 HS33]

    # @infiltrate
    return He, Ge#, FCij  # @infiltrate
end
function integraelem(
    pf,
    x,
    eta,
    w,
    elem,
    dad::Union{placa_espessa,placa_espessa_isotropica},
)
    h = zeros(Float64, 3, 3 * size(elem))
    g = zeros(Float64, 3, 3 * size(elem))
    # Cij = zeros(Float64, 3, 3)
    Nm = zeros(Float64, 3, 3 * size(elem))

    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        Past, Wast = calsolfund(pg', pf, [sy, -sx], dad)
        Nm[1, 1:3:end] = N
        Nm[2, 2:3:end] = N
        Nm[3, 3:3:end] = N
        h += Past * Nm * dgamadqsi * w[k]
        g += Wast * Nm * dgamadqsi * w[k]
        # Cij += FCij * dgamadqsi * w[k]
        # @infiltrate
    end
    h, g#, Cij
end
function integraelemsing_num(
    pf,
    x,
    elem,
    dad::Union{placa_espessa_isotropica},
    xi0,
    npg = 30,
)
    h = zeros(Float64, 3, 3 * size(elem))
    g = zeros(Float64, 3, 3 * size(elem))
    # Cij = zeros(Float64, 3, 3)

    Nm = zeros(Float64, 3, 3 * size(elem))
    eta, w = novelquad(2, xi0, npg)
    # @infiltrate

    for k = 1:size(w, 1)

        # @infiltrate
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        Past, Wast = calsolfund(pg', pf, [sy, -sx], dad)
        Nm[1, 1:3:end] = N
        Nm[2, 2:3:end] = N
        Nm[3, 3:3:end] = N
        h += Past * Nm * dgamadqsi * w[k]
        g += Wast * Nm * dgamadqsi * w[k]
        # Cij += FCij * dgamadqsi * w[k]

        # @infiltrate
    end
    h, g#, Cij
end

function compute_q(
    pf,
    x,
    eta,
    Gauss_w,
    elem,
    dad::Union{placa_espessa,placa_espessa_isotropica},
)

    # Inicialização da variável q_el
    q_el = zeros(3)
    # Início da integração numérica sobre o elemento
    for k = 1:size(Gauss_w, 1)# Integração da carga distribuída sobre o domínio
        N, dN = calc_fforma(eta[k], elem)
        pc = N' * x # Ponto de gauss interpolador

        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

        n = [sy, -sx]
        rs = pc' - pf
        r = norm(rs)
        theta = atan(rs[2], rs[1])
        int = zeros(3)
        jacob = r / 2

        nr = (n[1] * rs[1] + n[2] * rs[2]) / r

        for kr = 1:size(Gauss_w, 1)# Integração da carga distribuída sobre o domínio
            rho = (eta[kr] / 2 + 1 / 2) * r
            xc = pf[1] + rho * cos(theta)
            yc = pf[2] + rho * sin(theta)

            ~, Wast = calsolfund([xc, yc], pf, n, dad)
            C = dad.k.carga[3]
            int = int + Wast * [0, 0, C] * rho * jacob * Gauss_w[kr]
        end
        # Cálculo dos componentes de q_el [2 X 2]
        q_el += int * nr / r * dgamadqsi * Gauss_w[k]
        # @infiltrate
    end
    #--------------------------------------------------------------------------------
    return q_el
end

function separa(dad::Union{placa_espessa,placa_espessa_isotropica}, x)
    n_canto = size(dad.k.cantos, 1) # Number of corners
    if typeof(dad) == placa_espessa_isotropica
        scale = dad.k.D
    else
        scale = dad.k.D22
    end
    n_ipoints = length(dad.pontos_internos[:, 1]) # Number of internal nodes
    tra = zeros()
    n = size(dad.NOS, 1)
    T = zeros(3n)
    q = zeros(3n)
    for elem_i in dad.ELEM, i = 1:3   # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[i] == 0
            T[3ind_elem.+(i-3)] = elem_i.valorCDC[i, :] # A temperatura é a condição de contorno
            q[3ind_elem.+(i-3)] = x[3ind_elem.+(i-3)] * scale# O fluxo é o valor calculado
        elseif elem_i.tipoCDC[i] == 1
            T[3ind_elem.+(i-3)] = x[3ind_elem.+(i-3)] #
            q[3ind_elem.+(i-3)] = elem_i.valorCDC[i, :]
        end
    end
    Tint = x[3n.+(1:3n_ipoints)]
    [T; Tint], q
end
function Bess(z)
    K0 = BEM.besselk(0, z)
    K1 = BEM.besselk(1, z)
    Az = K0 + 2 / z * (K1 - 1 / z)
    Bz = K0 + 1 / z * (K1 - 1 / z)
    Az, Bz, K1, K0
end
function Monta_M_RIM(dad::Union{placa_espessa,placa_espessa_isotropica}, npg1 = 10)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi1, w1 = gausslegendre(npg1)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(3n_pontos, 3n_pontos)
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
        end
        pr = nodes[i, :]
        for j = 1:n_pontos #Laço dos pontos fontes
            pr = nodes[j, :]
            for el = 1:nelem
                elem_j = dad.ELEM[el]#Laço dos elementos
                x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
                m_el = calc_m(x, pf, pr, qsi1, w1, elem_j, dad)
                M[3i-2:3i, 3j-2:3j] += m_el
            end
        end
    end
    F3 = zeros(3n_pontos, 3n_pontos)
    for i = 1:n_pontos
        for j = 1:n_pontos
            F3[3i-2, 3j-2] = F[i, j]
            F3[3i-1, 3j-1] = F[i, j]
            F3[3i, 3j] = F[i, j]
        end
    end
    # @infiltrate
    M = M / F3
    M
end
function calc_m(
    x,
    pf,
    pr,
    qsi1,
    w1,
    elem,
    dad::Union{placa_espessa,placa_espessa_isotropica},
)

    npg = length(w1)
    m_el = zeros(3, 3)

    for i = 1:npg
        N, dN = calc_fforma(qsi1[i], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        ny = -dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        nx = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        r = pg' - pf
        m = calcula_F(pr, pf, pg, [nx, ny], qsi1, w1, dad)
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w1[i]
    end
    return m_el
end
function calcula_F(
    pr,
    pf,
    pg,
    n,
    qsi,
    w,
    dad::Union{placa_espessa,placa_espessa_isotropica},
) #
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(3, 3) # Inicializa a integral de F_area
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
        rline = norm(xc - pr)
        ro = r / 2 * (qsi[i] + 1)
        ~, term = calsolfund(xc, pf, n, dad)
        f = interpola(rline)
        # @infiltrate
        F = F + term * f * ro * drodqsi * w[i]# Integral de F_area
    end
    return F
end
function autovalor(H, G, M, dad::Union{placa_espessa,placa_espessa_isotropica}) #

    ind_deslconhecido = zeros(Int, 0)
    nc = size(dad.NOS, 1)
    ncanto = size(dad.k.bc_canto, 1)
    nh = size(H, 2)
    for elem_i in dad.ELEM, i = 1:3  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        # @show elem_i.tipoCDC[i]
        if elem_i.tipoCDC[i] == 0
            ind_deslconhecido = [ind_deslconhecido; 3ind_elem .+ (i - 3)]
        end
    end

    sort!(ind_deslconhecido)
    ind_forcaconhecida = setdiff(1:nh, ind_deslconhecido)
    H22 = H[ind_forcaconhecida, ind_forcaconhecida]
    H12 = H[ind_deslconhecido, ind_forcaconhecida]
    M22 = M[ind_forcaconhecida, ind_forcaconhecida]
    M12 = M[ind_deslconhecido, ind_forcaconhecida]
    G11 = G[ind_deslconhecido, ind_deslconhecido]
    G21 = G[ind_forcaconhecida, ind_deslconhecido]

    Hb = H22 - G21 * (G11 \ H12)
    Mb = M22 - G21 * (G11 \ M12)
    u = zeros(nh)

    a, v = eigen(Hb, Mb)
    ~, ind = findmin(abs, a)
    u[ind_forcaconhecida] = real.(v[:, ind])
    # @infiltrate

    real.(sort(abs.(a))[1:3]), [u[1:2:2nc]; u[2nc+1:end-ncanto]]
end

function aplicaT(dad::Union{placa_espessa,placa_espessa_isotropica}, M, N)
    n = size(N, 1)
    NOS = [dad.NOS; dad.pontos_internos]
    K = zeros(2 * size(NOS, 1), size(NOS, 1))
    F, Fx, Fy = BEM.montaFs(NOS)
    # @infiltrate
    Fyx = Fy * Fx
    Fxy = Fx * Fy
    Fyy = Fy * Fy
    Fxx = Fx * Fx
    Nx = Fx * N
    Ny = Fy * N

    K =
        (Nx[:, 1] + Ny[:, 3])' .* Fx +
        (Nx[:, 3] + Ny[:, 2])' .* Fy +
        (N[:, 1]' .* Fxx + N[:, 3]' .* Fxy + N[:, 2]' .* Fyy + N[:, 3]' .* Fyx)
    k3 = zeros(3 * n, 3 * n)
    # @infiltrate
    for i = 1:n
        for j = 1:n
            k3[3*i, 3*j] = K[i, j]
        end
    end
    M = M * k3
end

function Monta_M_RIMd(dad::Union{placa_espessa,placa_espessa_isotropica}, npg1 = 10)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi1, w1 = gausslegendre(npg1)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(3n_pontos, 3)
    M1 = zeros(n_pontos)
    D = zeros(3n_pontos, 3n_pontos)

    F = zeros(n_pontos, n_pontos)
    # Cálculo da matriz [F]
    @showprogress "Montando F" for i = 1:n_pontos
        p1 = nodes[i, :]
        for j = 1:n_pontos
            p2 = nodes[j, :]
            r = norm(p1 - p2)
            F[i, j] = interpola(r)
            if i == j
                continue
            end
            ~, D[3i-2:3i, 3j-2:3j] = calsolfund(nodes[i, :], nodes[j, :], [0, 0], dad)
        end
    end
    @showprogress "Montando Md" for i = 1:n_pontos
        pf = nodes[i, :] # Coordenada (x,y)dos pontos fonte
        for el = 1:nelem
            elem_j = dad.ELEM[el]#Laço dos elementos
            x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
            m_el, m_el1 = calc_md(x, pf, qsi1, w1, elem_j, dad)
            M1[i] += m_el
            M[3i-2:3i, :] += m_el1
        end
    end
    aux = M1' / F
    aux = [aux; aux; aux][:]'
    A = aux .* D
    # @infiltrate
    for i = 1:n_pontos #Laço dos pontos radiais
        A[3i-2:3i, 3i-2:3i] .= 0
        A[3i-2:3i, 3i-2:3i] =
            -[sum(A[3i-2:3i, 1:3:end], dims = 2) sum(A[3i-2:3i, 2:3:end], dims = 2) sum(
                A[3i-2:3i, 3:3:end],
                dims = 2,
            )] + M[3i-2:3i, :]
    end
    A
end
function calc_md(x, pf, qsi, w, elem, dad::Union{placa_espessa,placa_espessa_isotropica})
    npg = length(w)
    m_el, m_el1 = 0, zeros(3, 3)
    for i = 1:npg
        N, dN = calc_fforma(qsi[i], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        ny = -dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        nx = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        r = pg' - pf
        R = norm(r)
        m = int_interpolaρdρ(R)
        m1 = calcula_F1(pf, pg, [nx, ny], qsi, w, dad)
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
        m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
    end
    return m_el, m_el1
end

function calcula_F1(pf, pg, n, qsi, w, dad::Union{placa_espessa,placa_espessa_isotropica}) #
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(3, 3) # Inicializa a integral de F_area
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
        ro = r / 2 * (qsi[i] + 1)
        ~, term = calsolfund(xc, pf, n, dad)
        # f = interpola(rline)
        # @infiltrate
        F = F + term * ro * drodqsi * w[i]# Integral de F_area
    end
    return F
end
