function Compute_Material_Placa(Material)
    n_laminae = size(Material, 1)
    D11 = 0
    D22 = 0
    D66 = 0
    D12 = 0
    D16 = 0
    D26 = 0
    hant = 0
    for i = 1:n_laminae
        # i = 1
        theta = Material[i, 7]
        theta = theta * pi / 180
        E1 = Material[i, 2]
        E2 = Material[i, 3]
        G12 = Material[i, 4]
        ni12 = Material[i, 5]
        h = Material[i, 6]
        ni21 = ni12 * E2 / E1

        Q11 = E1 / (1 - ni12 * ni21)
        Q22 = E2 / (1 - ni12 * ni21)
        Q66 = G12
        Q16 = 0
        Q26 = 0
        Q12 = ni21 * E1 / (1 - ni12 * ni21)
        Q = [
            Q11 Q12 Q16
            Q12 Q22 Q26
            Q16 Q26 Q66
        ]
        m = cos(theta)
        n = sin(theta)
        T = [
            m^2 n^2 2*m*n
            n^2 m^2 -2*m*n
            -m*n m*n m^2-n^2
        ]
        Q = inv(T) * Q * inv(T')
        compliance = inv(Q)
        a11 = compliance[1, 1]
        a12 = compliance[1, 2]
        a16 = compliance[1, 3]
        a22 = compliance[2, 2]
        a26 = compliance[2, 3]
        a66 = compliance[3, 3]
        A = [
            a11 a12 a16
            a12 a22 a26
            a16 a26 a66
        ]
        del = 1 / det(A)
        B11 = del * (a22 * a66 - a26^2)
        B22 = del * (a11 * a66 - a16^2)
        B66 = del * (a11 * a22 - a12^2)
        B12 = del * (a16 * a26 - a12 * a66)
        B16 = del * (a12 * a26 - a22 * a16)
        B26 = del * (a12 * a16 - a11 * a26)
        D11 = D11 + 2 * ((h + hant)^3 - hant^3) / 3 * B11
        D22 = D22 + 2 * ((h + hant)^3 - hant^3) / 3 * B22
        D66 = D66 + 2 * ((h + hant)^3 - hant^3) / 3 * B66
        D12 = D12 + 2 * ((h + hant)^3 - hant^3) / 3 * B12
        D16 = D16 + 2 * ((h + hant)^3 - hant^3) / 3 * B16
        D26 = D26 + 2 * ((h + hant)^3 - hant^3) / 3 * B26
        hant = hant + h
    end
    # using PolynomialRoots
    charac_poly = [(D11), (4 * D16), (2 * D12 + 4 * D66), (4 * D26), (D22)] # Characteristic polynomial
    # charac_poly=[(D22),(4 * D26),(2 * D12 + 4 * D66),(4 * D16),(dad.k.D11)]; # Characteristic polynomial
    roots_poly = roots(charac_poly)
    # @show roots_poly
    # Conforme Lekhnitskii, p�gina 28,
    # parauma placa isotr�pica:mi1 = mi2 =iemi1c = mi2c = -i
    # Pega somente as ra�zes que tem a parte imagin�ria positiva.
    # Exclui o conjugado complexo.
    aux = zeros(Complex, 0)
    for i = 1:4
        if (imag(roots_poly[i]) > 0)
            push!(aux, roots_poly[i])
        end
    end
    mi = zeros(Complex, 2)
    # Coloca os n�meros em ordem crescente de acordo com a parte real.
    if (real(aux[1]) > real(aux[2]))
        mi[1] = aux[2]
        mi[2] = aux[1]
    else
        mi[1] = aux[1]
        mi[2] = aux[2]
    end
    d = zeros(2)
    e = zeros(2)
    d[1] = real(mi[1])
    e[1] = imag(mi[1])
    d[2] = real(mi[2])
    e[2] = imag(mi[2])
    # @infiltrate
    (
        d = d,
        e = e,
        D11 = D11,
        D22 = D22,
        D66 = D66,
        D12 = D12,
        D16 = D16,
        D26 = D26,
        Material = Material,
    )
end

function calc_HeG(dad::Union{placa_fina,placa_fina_isotropica}, npg = 8)
    nelem = size(dad.ELEM, 1)# Quantidade de elementos discretizados no contorno
    n_fis = size(dad.NOS, 1)
    n_internos = size(dad.pontos_internos, 1)
    n_cantos = size(dad.k.cantos, 1)
    H = zeros(2 * n_fis + n_cantos + n_internos, 2 * n_fis + n_cantos + n_internos)
    G = zeros(2 * n_fis + n_cantos + n_internos, 2 * n_fis + n_cantos + n_internos)
    q = zeros(2 * n_fis + n_cantos + n_internos)
    qsi, w = gausslegendre(npg)# Quadratura de gauss
    # qsi2, w2 = gausslegendre(2npg)# Quadratura de gauss
    normal_fonte = dad.normal
    pre = [zeros(2) for idx = 1:30]


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
                    # hn, gn = integraelemsing_num(pf, nf, x, elem_j, dad, pre, xi0, 30)
                    h, g2 = integraelemsing(x, dad, xi0)
                    g[2, 2:2:end] = g2
                    # @show [hn h], [gm g]
                end
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                # @infiltrate
                H[2i-1:2i, cols] = h
                G[2i-1:2i, cols] = g
                q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                q[2i-1:2i] += q_el
            elseif caso == "canto"
                if j == dad.k.cantos[i-n_fis-n_internos, 8]
                    xi0 = 1.0
                    h, g2 = integraelemsing(x, dad, xi0)
                    g[2, 2:2:end] = g2
                elseif j == dad.k.cantos[i-n_fis-n_internos, 9]
                    xi0 = -1.0
                    h, g2 = integraelemsing(x, dad, xi0)
                    g[2, 2:2:end] = g2
                end
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                ind = i - n_fis + 2n_fis
                H[ind, cols] = h[1, :]
                G[ind, cols] = g[1, :]
                q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                q[ind] += q_el[1]
            else
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                ind = i - n_fis + 2n_fis
                H[ind, cols] = h[1, :]
                G[ind, cols] = g[1, :]
                q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                q[ind] += q_el[1]
            end
        end
        if caso == "contorno"
            R, W = compute_Rw(pf, nf, dad, pre)
            cols = (n_fis * 2 + n_internos) .+ (1:n_cantos)
            H[2i-1:2i, cols] = R
            G[2i-1:2i, cols] = W
        else
            R, W = compute_Rw(pf, nf, dad, pre)
            ind = i + n_fis
            cols = (n_fis * 2 + n_internos) .+ (1:n_cantos)
            H[ind, cols] = R[1, :]
            G[ind, cols] = W[1, :]
        end
    end
    # for i = 1:n#i=1:size(dad.NOS,1) #Laço dos pontos fontes
    # H[2i-1:2i,2i-1:2i].=0
    # H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
    # end
    # @infiltrate
    H[1:2n_fis, 1:2n_fis] += I / 2
    H[(2n_fis).+(1:n_internos), (2n_fis).+(1:n_internos)] += I
    # H[2n_fis+n_internos+1:end, 2n_fis+n_internos+1:end] += I / 4
    H, G, q
end
function calsolfund(pg, pf, n, nf, dad::Union{placa_fina}, pre)

    R,
    S,
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
    d3Sdy3,
    d4Rdx2dy2,
    d4Rdx3dy,
    d4Rdx4,
    d4Rdxdy3,
    d4Rdy4,
    d4Sdx2dy2,
    d4Sdx3dy,
    d4Sdx4,
    d4Sdxdy3,
    d4Sdy4,
    dRdx,
    dRdy,
    dSdx,
    dSdy = pre
    rs = pg - pf
    nx = n[1]
    ny = n[2]
    m1 = nf[1]
    m2 = nf[2]
    # Distance of source and field points
    r1 = rs[1]
    r2 = rs[2]
    r = norm(rs)
    # Thin plate fundamental solutions
    theta = atan(r2, r1)

    G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
    H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
    C1 =
        ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[1])
    C2 =
        ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[2])
    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
    a = 1
    # a=constFS;

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
        d2Rdx2[i] =
            2 * log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            )
        d2Rdxdy[i] =
            2 *
            dad.k.d[i] *
            log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            ) -
            4 *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
        d2Rdy2[i] =
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            ) -
            8 *
            dad.k.d[i] *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
        d3Rdx3[i] =
            4 * (cos(theta) + dad.k.d[i] * sin(theta)) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Rdx2dy[i] =
            4 * (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Rdxdy2[i] =
            4 * (
                (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Rdy3[i] =
            4 * (
                dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d4Rdx4[i] =
            -4 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) /
            (
                r^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Rdx3dy[i] =
            -4 / r^2 * (
                dad.k.d[i] /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) +
                2 * dad.k.e[i]^2 * sin(theta) * cos(theta) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Rdx2dy2[i] =
            -4 / r^2 * (
                (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 * dad.k.e[i]^2 * cos(theta)^2 /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Rdxdy3[i] =
            -4 / r^2 * (
                dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 *
                dad.k.e[i]^2 *
                cos(theta) *
                (2 * dad.k.d[i] * cos(theta) + (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Rdy4[i] =
            -4 / r^2 * (
                (dad.k.d[i]^4 - dad.k.e[i]^4) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 *
                dad.k.e[i]^2 *
                cos(theta) *
                (
                    (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                ) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )

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
        d2Sdx2[i] =
            2 * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
        d2Sdxdy[i] =
            dad.k.e[i] * (log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            )) +
            2 *
            dad.k.d[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdy2[i] =
            2 *
            dad.k.d[i] *
            dad.k.e[i] *
            log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            ) +
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
        d3Sdx3[i] =
            -2 * dad.k.e[i] * sin(theta) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))

        d3Sdx2dy[i] =
            2 * dad.k.e[i] * cos(theta) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Sdxdy2[i] =
            2 *
            dad.k.e[i] *
            (
                2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Sdy3[i] =
            2 *
            dad.k.e[i] *
            (
                (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d4Sdx4[i] =
            4 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) / (
                r^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Sdx3dy[i] =
            2 * dad.k.e[i] / r^2 * (
                1 /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 * cos(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Sdx2dy2[i] =
            -4 *
            dad.k.e[i] *
            cos(theta) *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) / (
                r^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Sdxdy3[i] =
            -2 * dad.k.e[i] / r^2 * (
                (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) +
                (
                    2 *
                    (dad.k.d[i]^2 + dad.k.e[i]^2) *
                    cos(theta) *
                    (cos(theta) + dad.k.d[i] * sin(theta)) -
                    4 * dad.k.e[i]^2 * cos(theta)^2
                ) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )

        d4Sdy4[i] =
            -4 * dad.k.e[i] / r^2 * (
                dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) +
                cos(theta) * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                ) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
    end

    w = 1 / (8 * pi * dad.k.D22) * (C1 * R[1] + C2 * R[2] + C3 * (S[1] - S[2]))
    dwdx = 1 / (8 * pi) * (C1 * dRdx[1] + C2 * dRdx[2] + C3 * (dSdx[1] - dSdx[2]))
    dwdy = 1 / (8 * pi) * (C1 * dRdy[1] + C2 * dRdy[2] + C3 * (dSdy[1] - dSdy[2]))
    d2wdx2 = 1 / (8 * pi) * (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
    d2wdxdy =
        1 / (8 * pi) * (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
    d2wdy2 = 1 / (8 * pi) * (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))
    d3wdx3 = 1 / (8 * pi) * (C1 * d3Rdx3[1] + C2 * d3Rdx3[2] + C3 * (d3Sdx3[1] - d3Sdx3[2]))
    d3wdx2dy =
        1 / (8 * pi) *
        (C1 * d3Rdx2dy[1] + C2 * d3Rdx2dy[2] + C3 * (d3Sdx2dy[1] - d3Sdx2dy[2]))
    d3wdxdy2 =
        1 / (8 * pi) *
        (C1 * d3Rdxdy2[1] + C2 * d3Rdxdy2[2] + C3 * (d3Sdxdy2[1] - d3Sdxdy2[2]))
    d3wdy3 = 1 / (8 * pi) * (C1 * d3Rdy3[1] + C2 * d3Rdy3[2] + C3 * (d3Sdy3[1] - d3Sdy3[2]))
    d4wdx4 = 1 / (8 * pi) * (C1 * d4Rdx4[1] + C2 * d4Rdx4[2] + C3 * (d4Sdx4[1] - d4Sdx4[2]))
    d4wdx3dy =
        1 / (8 * pi) *
        (C1 * d4Rdx3dy[1] + C2 * d4Rdx3dy[2] + C3 * (d4Sdx3dy[1] - d4Sdx3dy[2]))
    d4wdx2dy2 =
        1 / (8 * pi) *
        (C1 * d4Rdx2dy2[1] + C2 * d4Rdx2dy2[2] + C3 * (d4Sdx2dy2[1] - d4Sdx2dy2[2]))
    d4wdxdy3 =
        1 / (8 * pi) *
        (C1 * d4Rdxdy3[1] + C2 * d4Rdxdy3[2] + C3 * (d4Sdxdy3[1] - d4Sdxdy3[2]))
    d4wdy4 = 1 / (8 * pi) * (C1 * d4Rdy4[1] + C2 * d4Rdy4[2] + C3 * (d4Sdy4[1] - d4Sdy4[2]))
    f1 = dad.k.D11 * nx^2 + 2 * dad.k.D16 * nx * ny + dad.k.D12 * ny^2
    f2 = 2 * (dad.k.D16 * nx^2 + 2 * dad.k.D66 * nx * ny + dad.k.D26 * ny^2)
    f3 = dad.k.D12 * nx^2 + 2 * dad.k.D26 * nx * ny + dad.k.D22 * ny^2

    h1 = dad.k.D11 * nx * (1 + ny^2) + 2 * dad.k.D16 * ny^3 - dad.k.D12 * nx * ny^2
    h2 =
        4 * dad.k.D16 * nx + dad.k.D12 * ny * (1 + nx^2) + 4 * dad.k.D66 * ny^3 -
        dad.k.D11 * nx^2 * ny - 2 * dad.k.D26 * nx * ny^2
    h3 =
        4 * dad.k.D26 * ny + dad.k.D12 * nx * (1 + ny^2) + 4 * dad.k.D66 * nx^3 -
        dad.k.D22 * nx * ny^2 - 2 * dad.k.D16 * nx^2 * ny
    h4 = dad.k.D22 * ny * (1 + nx^2) + 2 * dad.k.D26 * nx^3 - dad.k.D12 * nx^2 * ny

    dwdn = (dwdx * nx + dwdy * ny) / dad.k.D22
    mn = -(f1 * d2wdx2 + f2 * d2wdxdy + f3 * d2wdy2) / dad.k.D22
    vn = -(h1 * d3wdx3 + h2 * d3wdx2dy + h3 * d3wdxdy2 + h4 * d3wdy3) / dad.k.D22

    dmndx = -(f1 * d3wdx3 + f2 * d3wdx2dy + f3 * d3wdxdy2)
    dmndy = -(f1 * d3wdx2dy + f2 * d3wdxdy2 + f3 * d3wdy3)
    dvndx = -(h1 * d4wdx4 + h2 * d4wdx3dy + h3 * d4wdx2dy2 + h4 * d4wdxdy3)
    dvndy = -(h1 * d4wdx3dy + h2 * d4wdx2dy2 + h3 * d4wdxdy3 + h4 * d4wdy4)
    dwdm = -(dwdx * m1 + dwdy * m2) / dad.k.D22
    d2wdndm =
        -(d2wdx2 * nx * m1 + d2wdxdy * (nx * m2 + ny * m1) + d2wdy2 * ny * m2) / dad.k.D22
    dmndm = -(dmndx * m1 + dmndy * m2) / dad.k.D22
    dvndm = -(dvndx * m1 + dvndy * m2) / dad.k.D22
    # Assembly of matrices that contain fundamental solutions.
    u = [
        w -dwdn
        dwdm -d2wdndm
    ]
    p = [
        vn -mn
        dvndm -dmndm
    ]
    # @infiltrate
    u, p
    # @infiltrate
end
function calsolfund2(pg, pf, n, nf, dad::Union{placa_fina}, pre)

    R,
    S,
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
    d3Sdy3,
    d4Rdx2dy2,
    d4Rdx3dy,
    d4Rdx4,
    d4Rdxdy3,
    d4Rdy4,
    d4Sdx2dy2,
    d4Sdx3dy,
    d4Sdx4,
    d4Sdxdy3,
    d4Sdy4,
    dRdx,
    dRdy,
    dSdx,
    dSdy = pre
    rs = pg - pf
    nx = n[1]
    ny = n[2]
    m1 = nf[1]
    m2 = nf[2]
    s1 = -ny
    s2 = nx
    # Distance of source and field points
    r1 = rs[1]
    r2 = rs[2]
    r = norm(rs)
    # Thin plate fundamental solutions
    theta = atan(r2, r1)

    G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
    H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
    C1 =
        ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[1])
    C2 =
        ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[2])
    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
    a = 1
    # a=constFS;

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
        d2Rdx2[i] =
            2 * log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            )
        d2Rdxdy[i] =
            2 *
            dad.k.d[i] *
            log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            ) -
            4 *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
        d2Rdy2[i] =
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            ) -
            8 *
            dad.k.d[i] *
            dad.k.e[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
        d3Rdx3[i] =
            4 * (cos(theta) + dad.k.d[i] * sin(theta)) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Rdx2dy[i] =
            4 * (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Rdxdy2[i] =
            4 * (
                (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Rdy3[i] =
            4 * (
                dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d4Rdx4[i] =
            -4 * ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) /
            (
                r^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Rdx3dy[i] =
            -4 / r^2 * (
                dad.k.d[i] /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) +
                2 * dad.k.e[i]^2 * sin(theta) * cos(theta) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Rdx2dy2[i] =
            -4 / r^2 * (
                (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 * dad.k.e[i]^2 * cos(theta)^2 /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Rdxdy3[i] =
            -4 / r^2 * (
                dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 *
                dad.k.e[i]^2 *
                cos(theta) *
                (2 * dad.k.d[i] * cos(theta) + (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Rdy4[i] =
            -4 / r^2 * (
                (dad.k.d[i]^4 - dad.k.e[i]^4) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 *
                dad.k.e[i]^2 *
                cos(theta) *
                (
                    (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                ) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )

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
        d2Sdx2[i] =
            2 * atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
        d2Sdxdy[i] =
            dad.k.e[i] * (log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            )) +
            2 *
            dad.k.d[i] *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

        d2Sdy2[i] =
            2 *
            dad.k.d[i] *
            dad.k.e[i] *
            log(
                r^2 / a^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2),
            ) +
            2 *
            (dad.k.d[i]^2 - dad.k.e[i]^2) *
            atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
        d3Sdx3[i] =
            -2 * dad.k.e[i] * sin(theta) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))

        d3Sdx2dy[i] =
            2 * dad.k.e[i] * cos(theta) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Sdxdy2[i] =
            2 *
            dad.k.e[i] *
            (
                2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d3Sdy3[i] =
            2 *
            dad.k.e[i] *
            (
                (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
            ) /
            (r * ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2))
        d4Sdx4[i] =
            4 * dad.k.e[i] * sin(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) / (
                r^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Sdx3dy[i] =
            2 * dad.k.e[i] / r^2 * (
                1 /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) -
                2 * cos(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Sdx2dy2[i] =
            -4 *
            dad.k.e[i] *
            cos(theta) *
            (
                dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                dad.k.e[i]^2 * sin(theta)
            ) / (
                r^2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
        d4Sdxdy3[i] =
            -2 * dad.k.e[i] / r^2 * (
                (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) +
                (
                    2 *
                    (dad.k.d[i]^2 + dad.k.e[i]^2) *
                    cos(theta) *
                    (cos(theta) + dad.k.d[i] * sin(theta)) -
                    4 * dad.k.e[i]^2 * cos(theta)^2
                ) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )

        d4Sdy4[i] =
            -4 * dad.k.e[i] / r^2 * (
                dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2) +
                cos(theta) * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                ) /
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2)^2
            )
    end

    w = 1 / (8 * pi * dad.k.D22) * (C1 * R[1] + C2 * R[2] + C3 * (S[1] - S[2]))
    dwdx = 1 / (8 * pi) * (C1 * dRdx[1] + C2 * dRdx[2] + C3 * (dSdx[1] - dSdx[2]))
    dwdy = 1 / (8 * pi) * (C1 * dRdy[1] + C2 * dRdy[2] + C3 * (dSdy[1] - dSdy[2]))
    d2wdx2 = 1 / (8 * pi) * (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
    d2wdxdy =
        1 / (8 * pi) * (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
    d2wdy2 = 1 / (8 * pi) * (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))
    d3wdx3 = 1 / (8 * pi) * (C1 * d3Rdx3[1] + C2 * d3Rdx3[2] + C3 * (d3Sdx3[1] - d3Sdx3[2]))
    d3wdx2dy =
        1 / (8 * pi) *
        (C1 * d3Rdx2dy[1] + C2 * d3Rdx2dy[2] + C3 * (d3Sdx2dy[1] - d3Sdx2dy[2]))
    d3wdxdy2 =
        1 / (8 * pi) *
        (C1 * d3Rdxdy2[1] + C2 * d3Rdxdy2[2] + C3 * (d3Sdxdy2[1] - d3Sdxdy2[2]))
    d3wdy3 = 1 / (8 * pi) * (C1 * d3Rdy3[1] + C2 * d3Rdy3[2] + C3 * (d3Sdy3[1] - d3Sdy3[2]))
    d4wdx4 = 1 / (8 * pi) * (C1 * d4Rdx4[1] + C2 * d4Rdx4[2] + C3 * (d4Sdx4[1] - d4Sdx4[2]))
    d4wdx3dy =
        1 / (8 * pi) *
        (C1 * d4Rdx3dy[1] + C2 * d4Rdx3dy[2] + C3 * (d4Sdx3dy[1] - d4Sdx3dy[2]))
    d4wdx2dy2 =
        1 / (8 * pi) *
        (C1 * d4Rdx2dy2[1] + C2 * d4Rdx2dy2[2] + C3 * (d4Sdx2dy2[1] - d4Sdx2dy2[2]))
    d4wdxdy3 =
        1 / (8 * pi) *
        (C1 * d4Rdxdy3[1] + C2 * d4Rdxdy3[2] + C3 * (d4Sdxdy3[1] - d4Sdxdy3[2]))
    d4wdy4 = 1 / (8 * pi) * (C1 * d4Rdy4[1] + C2 * d4Rdy4[2] + C3 * (d4Sdy4[1] - d4Sdy4[2]))
    f1 = dad.k.D11 * nx^2 + 2 * dad.k.D16 * nx * ny + dad.k.D12 * ny^2
    f2 = 2 * (dad.k.D16 * nx^2 + 2 * dad.k.D66 * nx * ny + dad.k.D26 * ny^2)
    f3 = dad.k.D12 * nx^2 + 2 * dad.k.D26 * nx * ny + dad.k.D22 * ny^2

    h1 = dad.k.D11 * nx * (1 + ny^2) + 2 * dad.k.D16 * ny^3 - dad.k.D12 * nx * ny^2
    h2 =
        4 * dad.k.D16 * nx + dad.k.D12 * ny * (1 + nx^2) + 4 * dad.k.D66 * ny^3 -
        dad.k.D11 * nx^2 * ny - 2 * dad.k.D26 * nx * ny^2
    h3 =
        4 * dad.k.D26 * ny + dad.k.D12 * nx * (1 + ny^2) + 4 * dad.k.D66 * nx^3 -
        dad.k.D22 * nx * ny^2 - 2 * dad.k.D16 * nx^2 * ny
    h4 = dad.k.D22 * ny * (1 + nx^2) + 2 * dad.k.D26 * nx^3 - dad.k.D12 * nx^2 * ny

    dwdn = (dwdx * nx + dwdy * ny) / dad.k.D22
    mn = -(f1 * d2wdx2 + f2 * d2wdxdy + f3 * d2wdy2) / dad.k.D22
    vn = -(h1 * d3wdx3 + h2 * d3wdx2dy + h3 * d3wdxdy2 + h4 * d3wdy3) / dad.k.D22

    dmndx = -(f1 * d3wdx3 + f2 * d3wdx2dy + f3 * d3wdxdy2)
    dmndy = -(f1 * d3wdx2dy + f2 * d3wdxdy2 + f3 * d3wdy3)
    dvndx = -(h1 * d4wdx4 + h2 * d4wdx3dy + h3 * d4wdx2dy2 + h4 * d4wdxdy3)
    dvndy = -(h1 * d4wdx3dy + h2 * d4wdx2dy2 + h3 * d4wdxdy3 + h4 * d4wdy4)
    dwdm = -(dwdx * m1 + dwdy * m2) / dad.k.D22
    d2wdndm =
        -(d2wdx2 * nx * m1 + d2wdxdy * (nx * m2 + ny * m1) + d2wdy2 * ny * m2) / dad.k.D22
    dmndm = -(dmndx * m1 + dmndy * m2) / dad.k.D22
    dvndm = -(dvndx * m1 + dvndy * m2) / dad.k.D22
    # Assembly of matrices that contain fundamental solutions.

    dwds = (dwdx * s1 + dwdy * s2) / dad.k.D22
    d2wdmds =
        (d2wdx2 * m1 * s1 + d2wdxdy * (m1 * s2 + m2 * s1) + d2wdy2 * m2 * s2) / dad.k.D22
    #  @infiltrate
    u = [
        w -dwdn
        dwdm -d2wdndm
    ]
    p = [
        vn -mn
        dvndm -dmndm
    ]
    # @infiltrate
    u2 = [
        -dwds
        -d2wdmds
    ]
    u, p, u2
    # @infiltrate
end
function integraelem(
    pf,
    nf,
    x,
    eta,
    w,
    elem,
    dad::Union{placa_fina,placa_fina_isotropica},
    pre,
)
    h = zeros(Float64, 2, 2 * size(elem))
    g = zeros(Float64, 2, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        wast, vast = calsolfund(pg', pf, [sy, -sx], nf, dad, pre)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        h += wast * Nm * dgamadqsi * w[k]
        g += vast * Nm * dgamadqsi * w[k]

        # @infiltrate
    end
    h, g
end
function integraelem2(
    pf,
    nf,
    x,
    eta,
    w,
    elem,
    dad::Union{placa_fina,placa_fina_isotropica},
    pre,
)
    h = zeros(Float64, 2, 2 * size(elem))
    g = zeros(Float64, 2, 2 * size(elem))
    g2 = zeros(Float64, 1, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        wast, vast, wast2 = calsolfund2(pg', pf, [sy, -sx], nf, dad, pre)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        h += vast * Nm * dgamadqsi * w[k]
        g += wast * Nm * dgamadqsi * w[k]
        # @infiltrate
        g2 += wast2' * Nm * dgamadqsi * w[k]

    end
    h, g, g2
end

function geraNiv(xi)
    n = size(xi, 1)
    N = zeros(n, n)
    for i = 1:n
        N[i, 1:n], ~ = calc_fforma_gen(xi[i], xi)
    end
    inv(N)
end
function compute_Rw(pf, nf, dad::placa_fina, pre)
    m1 = nf[1]
    m2 = nf[2]
    R,
    S,
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
    d3Sdy3,
    dRdx,
    dRdy,
    dSdx,
    dSdy,
    extra = pre
    RS = zeros(2, size(dad.k.cantos, 1))
    WS = zeros(2, size(dad.k.cantos, 1))
    for i = 1:size(dad.k.cantos, 1)
        pc = dad.k.cantos[i, 2:3]
        na1 = dad.k.cantos[i, 4]
        na2 = dad.k.cantos[i, 5]
        nd1 = dad.k.cantos[i, 6]
        nd2 = dad.k.cantos[i, 7]
        rs = pc - pf
        # Distance of source and field points
        r1 = rs[1]
        r2 = rs[2]
        r = norm(rs)
        if r != 0

            G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
            H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
            C1 =
                ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) /
                (G * H * dad.k.e[1])
            C2 =
                ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) /
                (G * H * dad.k.e[2])
            C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
            # a = constFS
            a = 1
            theta = atan(r2, r1)

            for i = 1:2
                R[i] =
                    r^2 *
                    (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 -
                        dad.k.e[i]^2 * sin(theta)^2
                    ) *
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
                    (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 -
                        dad.k.e[i]^2 * sin(theta)^2
                    ) *
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
                d2Rdx2[i] =
                    2 * log(
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        ),
                    )
                ####
                d2Rdxdy[i] =
                    2 *
                    dad.k.d[i] *
                    log(
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        ),
                    ) -
                    4 *
                    dad.k.e[i] *
                    atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
                ###
                d2Rdy2[i] =
                    2 *
                    (dad.k.d[i]^2 - dad.k.e[i]^2) *
                    log(
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        ),
                    ) -
                    8 *
                    dad.k.d[i] *
                    dad.k.e[i] *
                    atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
                d3Rdx3[i] =
                    4 * (cos(theta) + dad.k.d[i] * sin(theta)) / (
                        r * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        )
                    )
                d3Rdx2dy[i] =
                    4 * (
                        dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                        dad.k.e[i]^2 * sin(theta)
                    ) / (
                        r * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        )
                    )
                d3Rdxdy2[i] =
                    4 * (
                        (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                        (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
                    ) / (
                        r * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        )
                    )
                d3Rdy3[i] =
                    4 * (
                        dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                        (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                    ) / (
                        r * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        )
                    )


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

                d2Sdx2[i] =
                    2 *
                    atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
                d2Sdxdy[i] =
                    dad.k.e[i] * (log(
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
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
                        r^2 / a^2 * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        ),
                    ) +
                    2 *
                    (dad.k.d[i]^2 - dad.k.e[i]^2) *
                    atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))
                d3Sdx3[i] =
                    -2 * dad.k.e[i] * sin(theta) / (
                        r * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        )
                    )

                d3Sdx2dy[i] =
                    2 * dad.k.e[i] * cos(theta) / (
                        r * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        )
                    )
                d3Sdxdy2[i] =
                    2 *
                    dad.k.e[i] *
                    (
                        2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                        (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                    ) / (
                        r * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        )
                    )
                d3Sdy3[i] =
                    2 *
                    dad.k.e[i] *
                    (
                        (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                        2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                    ) / (
                        r * (
                            (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                            dad.k.e[i]^2 * sin(theta)^2
                        )
                    )
            end
            we = 1 / (8 * pi * dad.k.D22) * (C1 * R[1] + C2 * R[2] + C3 * (S[1] - S[2]))

            dwdx = 1 / (8 * pi) * (C1 * dRdx[1] + C2 * dRdx[2] + C3 * (dSdx[1] - dSdx[2]))
            dwdy = 1 / (8 * pi) * (C1 * dRdy[1] + C2 * dRdy[2] + C3 * (dSdy[1] - dSdy[2]))
            d2wdx2 =
                1 / (8 * pi) *
                (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
            d2wdxdy =
                1 / (8 * pi) *
                (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
            d2wdy2 =
                1 / (8 * pi) *
                (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))
            d3wdx3 =
                1 / (8 * pi) *
                (C1 * d3Rdx3[1] + C2 * d3Rdx3[2] + C3 * (d3Sdx3[1] - d3Sdx3[2]))
            d3wdx2dy =
                1 / (8 * pi) *
                (C1 * d3Rdx2dy[1] + C2 * d3Rdx2dy[2] + C3 * (d3Sdx2dy[1] - d3Sdx2dy[2]))
            d3wdxdy2 =
                1 / (8 * pi) *
                (C1 * d3Rdxdy2[1] + C2 * d3Rdxdy2[2] + C3 * (d3Sdxdy2[1] - d3Sdxdy2[2]))
            d3wdy3 =
                1 / (8 * pi) *
                (C1 * d3Rdy3[1] + C2 * d3Rdy3[2] + C3 * (d3Sdy3[1] - d3Sdy3[2]))
            g1a = (dad.k.D12 - dad.k.D11) * na1 * na2 + dad.k.D16 * (na1^2 - na2^2)
            g2a = 2 * (dad.k.D26 - dad.k.D16) * na1 * na2 + 2 * dad.k.D66 * (na1^2 - na2^2)
            g3a = (dad.k.D22 - dad.k.D12) * na1 * na2 + dad.k.D26 * (na1^2 - na2^2)

            g1d = (dad.k.D12 - dad.k.D11) * nd1 * nd2 + dad.k.D16 * (nd1^2 - nd2^2)
            g2d = 2 * (dad.k.D26 - dad.k.D16) * nd1 * nd2 + 2 * dad.k.D66 * (nd1^2 - nd2^2)
            g3d = (dad.k.D22 - dad.k.D12) * nd1 * nd2 + dad.k.D26 * (nd1^2 - nd2^2)

            tna = -(g1a * d2wdx2 + g2a * d2wdxdy + g3a * d2wdy2)
            tnd = -(g1d * d2wdx2 + g2d * d2wdxdy + g3d * d2wdy2)

            dtndxa = -(g1a * d3wdx3 + g2a * d3wdx2dy + g3a * d3wdxdy2)
            dtndya = -(g1a * d3wdx2dy + g2a * d3wdxdy2 + g3a * d3wdy3)
            dtndxd = -(g1d * d3wdx3 + g2d * d3wdx2dy + g3d * d3wdxdy2)
            dtndyd = -(g1d * d3wdx2dy + g2d * d3wdxdy2 + g3d * d3wdy3)

            dtndma = dtndxa * m1 + dtndya * m2
            dtndmd = dtndxd * m1 + dtndyd * m2

            Rci = (tnd - tna) / dad.k.D22
            dRci = -(dtndmd - dtndma) / dad.k.D22

            dwdm = -(dwdx * m1 + dwdy * m2) / dad.k.D22
            R = [Rci dRci]'
            w = [we dwdm]'
        else
            sa1 = -na2
            sa2 = na1
            sd1 = -nd2
            sd2 = nd1
            cbetac = -(sa1 * sd1 + sa2 * sd2)
            sbetac = -(na1 * sd1 + na2 * sd2)
            betac = atan(sbetac, cbetac)
            R = [betac / (2 * pi) 0]'
            w = [0 0]'
        end
        RS[:, i] = R
        WS[:, i] = w
    end
    return RS, WS
end
function compute_q(pf, nf, x, eta, Gauss_w, elem, dad::placa_fina, pre)

    # Inicialização da variável q_el
    q_el = zeros(2, 1)
    # Início da integração numérica sobre o elemento
    for k = 1:size(Gauss_w, 1)# Integração da carga distribuída sobre o domínio
        N, dN = calc_fforma(eta[k], elem)
        pc = N' * x # Ponto de gauss interpolador
        rs = pc' - pf
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        n1 = sy
        n2 = -sx
        xf = pf[1]
        yf = pf[2]
        # Distance of source and field points
        r1 = rs[1]
        r2 = rs[2]
        r = norm(rs)
        rd1 = r1 / r
        rd2 = r2 / r
        m1 = nf[1]
        m2 = nf[2]
        theta = atan(r2, r1)
        nr = n1 * rd1 + n2 * rd2

        G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
        H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
        C1 =
            ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[1])
        C2 =
            ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) /
            (G * H * dad.k.e[2])
        C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
        a = 1
        A = dad.k.carga[1]
        B = dad.k.carga[2]
        C = dad.k.carga[3]
        # Clocal = 1
        Clocal = A * xf + B * yf + C
        intRkrdr,
        intSkrdr,
        intRABrdr,
        intSABrdr,
        intdRdxABrdr,
        intdRdyABrdr,
        intdSdxABrdr,
        intdSdyABrdr,
        intdRdxrdr,
        intdRdyrdr,
        intdSdxrdr,
        intdSdyrdr,
        extra = pre

        for i = 1:2
            intRkrdr[i] =
                Clocal * (
                    r^3 * (
                        -16 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        sin(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) -
                        (
                            -7 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (
                            -1 - dad.k.d[i]^2 +
                            dad.k.e[i]^2 +
                            (-1 + dad.k.d[i]^2 - dad.k.e[i]^2) * cos(2 * theta) -
                            2 * dad.k.d[i] * sin(2 * theta)
                        )
                    )
                ) / 16
            intSkrdr[i] =
                Clocal * (
                    r^3 * (
                        2 *
                        dad.k.e[i] *
                        (
                            -7 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        sin(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) +
                        2 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (
                            1 + dad.k.d[i]^2 - dad.k.e[i]^2 +
                            (1 - dad.k.d[i]^2 + dad.k.e[i]^2) * cos(2 * theta) +
                            2 * dad.k.d[i] * sin(2 * theta)
                        )
                    )
                ) / 16
            intRABrdr[i] =
                (
                    r^4 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        -40 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        sin(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) -
                        (
                            -17 +
                            5 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (
                            -1 - dad.k.d[i]^2 +
                            dad.k.e[i]^2 +
                            (-1 + dad.k.d[i]^2 - dad.k.e[i]^2) * cos(2 * theta) -
                            2 * dad.k.d[i] * sin(2 * theta)
                        )
                    )
                ) / 50
            intSABrdr[i] =
                (
                    r^4 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        2 *
                        dad.k.e[i] *
                        (
                            -17 +
                            5 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        sin(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) +
                        5 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (
                            1 + dad.k.d[i]^2 - dad.k.e[i]^2 +
                            (1 - dad.k.d[i]^2 + dad.k.e[i]^2) * cos(2 * theta) +
                            2 * dad.k.d[i] * sin(2 * theta)
                        )
                    )
                ) / 50
            intdRdxABrdr[i] =
                (
                    r^3 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        -4 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        sin(theta) +
                        (
                            -5 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (cos(theta) + dad.k.d[i] * sin(theta))
                    )
                ) / 4
            intdRdyABrdr[i] =
                (
                    r^3 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        -4 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (cos(theta) + 2 * dad.k.d[i] * sin(theta)) +
                        (
                            -5 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (
                            dad.k.d[i] * cos(theta) +
                            (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                        )
                    )
                ) / 4
            intdSdxABrdr[i] =
                (
                    r^3 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        dad.k.e[i] *
                        (
                            -5 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        sin(theta) +
                        4 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (cos(theta) + dad.k.d[i] * sin(theta))
                    )
                ) / 8
            intdSdyABrdr[i] =
                (
                    r^3 *
                    (A * cos(theta) + B * sin(theta)) *
                    (
                        dad.k.e[i] *
                        (
                            -5 +
                            2 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        (cos(theta) + 2 * dad.k.d[i] * sin(theta)) +
                        4 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (
                            dad.k.d[i] * cos(theta) +
                            (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                        )
                    )
                ) / 8
            intdRdxrdr[i] =
                Clocal * (
                    2 *
                    r^2 *
                    (
                        -6 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        sin(theta) +
                        (
                            -8 +
                            3 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (cos(theta) + dad.k.d[i] * sin(theta))
                    )
                ) / 9
            intdRdyrdr[i] =
                Clocal * (
                    2 *
                    r^2 *
                    (
                        -6 *
                        dad.k.e[i] *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (cos(theta) + 2 * dad.k.d[i] * sin(theta)) +
                        (
                            -8 +
                            3 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) * (
                            dad.k.d[i] * cos(theta) +
                            (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                        )
                    )
                ) / 9
            intdSdxrdr[i] =
                Clocal * (
                    r^2 * (
                        dad.k.e[i] *
                        (
                            -8 +
                            3 * log(
                                (
                                    r^2 * (
                                        dad.k.e[i]^2 * sin(theta)^2 +
                                        (cos(theta) + dad.k.d[i] * sin(theta))^2
                                    )
                                ) / a^2,
                            )
                        ) *
                        sin(theta) +
                        6 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (cos(theta) + dad.k.d[i] * sin(theta))
                    )
                ) / 9
            intdSdyrdr[i] =
                -Clocal * (
                    r^2 * (
                        -(
                            dad.k.e[i] *
                            (
                                -8 +
                                3 * log(
                                    (
                                        r^2 * (
                                            dad.k.e[i]^2 * sin(theta)^2 +
                                            (cos(theta) + dad.k.d[i] * sin(theta))^2
                                        )
                                    ) / a^2,
                                )
                            ) *
                            (cos(theta) + 2 * dad.k.d[i] * sin(theta))
                        ) -
                        6 *
                        atan(
                            (dad.k.e[i] * sin(theta)),
                            (cos(theta) + dad.k.d[i] * sin(theta)),
                        ) *
                        (
                            dad.k.d[i] * cos(theta) +
                            (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                        )
                    )
                ) / 9
        end

        int1 =
            1 / (8 * pi * dad.k.D22) *
            (
                C1 * (intRkrdr[1] + intRABrdr[1]) +
                C2 * (intRkrdr[2] + intRABrdr[2]) +
                C3 * ((intSkrdr[1] + intSABrdr[1]) - (intSkrdr[2] + intSABrdr[2]))
            ) *
            nr
        intdwdxrdr =
            1 / (8 * pi) * (
                C1 * (intdRdxrdr[1] + intdRdxABrdr[1]) +
                C2 * (intdRdxrdr[2] + intdRdxABrdr[2]) +
                C3 * ((intdSdxrdr[1] + intdSdxABrdr[1]) - (intdSdxrdr[2] + intdSdxABrdr[2]))
            )
        intdwdyrdr =
            1 / (8 * pi) * (
                C1 * (intdRdyrdr[1] + intdRdyABrdr[1]) +
                C2 * (intdRdyrdr[2] + intdRdyABrdr[2]) +
                C3 * ((intdSdyrdr[1] + intdSdyABrdr[1]) - (intdSdyrdr[2] + intdSdyABrdr[2]))
            )
        int2 = -1 * (intdwdxrdr * m1 + intdwdyrdr * m2) * nr / dad.k.D22
        int_gw = [
            int1
            int2
        ]
        # Cálculo dos componentes de q_el [2 X 2]
        q_el = q_el + int_gw * dgamadqsi * Gauss_w[k]
        # @infiltrate
    end
    #--------------------------------------------------------------------------------
    return q_el
end

function separa(x, dad::Union{placa_fina,placa_fina_isotropica})
    n_canto = size(dad.k.cantos, 1) # Number of corners
    if typeof(dad) == placa_fina_isotropica
        scale = dad.k.D
    else
        scale = dad.k.D22
    end
    n_ipoints = length(dad.pontos_internos[:, 1]) # Number of internal nodes
    tra = zeros()
    n = size(dad.NOS, 1)
    T = zeros(2n)
    q = zeros(2n)
    for elem_i in dad.ELEM, i = 1:2   # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[i] == 0
            T[2ind_elem.+(i-2)] = elem_i.valorCDC[i, :] # A temperatura é a condição de contorno
            q[2ind_elem.+(i-2)] = x[2ind_elem.+(i-2)] * scale# O fluxo é o valor calculado
        elseif elem_i.tipoCDC[i] == 1
            T[2ind_elem.+(i-2)] = x[2ind_elem.+(i-2)] #
            q[2ind_elem.+(i-2)] = elem_i.valorCDC[i, :]
        end
    end
    Tint = x[2n.+(1:n_ipoints)]
    qint = zeros(n_ipoints)
    Tcanto = zeros(n_canto)
    qcanto = zeros(n_canto)
    for c = 1:n_canto
        if dad.k.bc_canto[c, 1] == 1
            Tcanto[c] = x[end-n_canto+c]
        else
            qcanto[c] = x[end-n_canto+c] * scale
        end
    end
    # @infiltrate
    [T; Tint; Tcanto], [q; qint; qcanto]
end
function passonotempo(dad, A, B, q, bc_val, M, dM, npg)
    M2, dM2dx2, dM2dy2, dM2dxdy, Mdtdm = dM
    nt = 2 * size(dad.NOS, 1) + size(dad.k.cantos, 1) + size(dad.pontos_internos, 1)
    dis = zeros(dad.k.n_dt, nt)
    tens = zeros(dad.k.n_dt, size(dad.NOS, 1), 3)
    tens_i = zeros(dad.k.n_dt, size(dad.pontos_internos, 1), 3)
    U_ant = zeros(nt, 3)
    prob = LinearProblem(A, q)
    linsolve = init(prob)
    sol = solve(linsolve)
    @showprogress "Resolvendo no tempo" for idt = 1:dad.k.n_dt
        b =
            B * bc_val +
            M * (1 / dad.k.dt^2 * (5 * U_ant[:, 3] - 4 * U_ant[:, 2] + U_ant[:, 1])) +
            q
        linsolve.b = b
        sol = solve!(linsolve)
        x = sol.u
        # x = A \ b
        listdisp, listtra = separa(x, dad)

        d2w_sta = calc_dw_int(dad, listtra, listdisp, npg)
        # @show d2w_sta
        dwdt_sta = calc_dwdt_cont(dad, listtra, listdisp, npg)
        d2wdt2 =
            1 / dad.k.dt^2 *
            (2 * listdisp - 5 * U_ant[:, 3] + 4 * U_ant[:, 2] - U_ant[:, 1])

        der2M =
            [dM2dx2 * d2wdt2 dM2dy2 * d2wdt2 dM2dxdy * d2wdt2] * dad.k.rho * dad.k.thickness

        # @infiltrate
        dMdt = Mdtdm[1:2:end, :] * d2wdt2 * dad.k.rho * dad.k.thickness
        d2w = d2w_sta - der2M

        dwdt_cont = dwdt_sta - 2 * dMdt
        d2wnt = calc_d2w_cont(dad, listtra, listdisp, dwdt_cont)


        tensaoc, tensaoi = compute_tens(d2w, d2wnt, dad)
        # @infiltrate

        dis[idt, :] = listdisp
        tens[idt, :, :] = tensaoc
        tens_i[idt, :, :] = tensaoi[:, :, 1]

        U_ant[:, 1] = U_ant[:, 2]
        U_ant[:, 2] = U_ant[:, 3]
        U_ant[:, 3] = listdisp
        # @infiltrate
    end
    dis, tens, tens_i
end
function integraelemsing(x, dad::placa_fina, xi0)
    # Cálculo da distância do ponto fonte [xf,yf] ao ponto campo
    dx = (x[3, 1] - x[1, 1]) / 2 * 3
    dy = (x[3, 2] - x[1, 2]) / 2 * 3
    L = sqrt(dx^2 + dy^2)
    sx = dx / L
    sy = dy / L
    nx = sy
    ny = -sx
    # Ângulo entre r2 e r1
    theta = atan(dy, dx)
    G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2
    H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2
    C1 =
        ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[1])
    C2 =
        ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[2]^2)) / (G * H * dad.k.e[2])
    C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H)
    a = 1
    h_el = zeros(2, 6)
    g_el = zeros(1, 3)
    # Integrate[N1d,(qsi,-1,1)]
    for x in [:NlogLqqo, :intNsr2, :intNsr, :intN]
        @eval $x = zeros(3)
    end
    intN[1] = 3 / 4
    intN[2] = 1 / 2
    intN[3] = intN[1]

    if xi0 ≈ 0
        #   2*Integrate[N1d[qsi] Log[L [Abs[qsi - qsio]]/2],(qsi,-1,1)]*Jac
        #  Jac=L/2
        NlogLqqo[1] = -(L * (1 + log(8) - 3 * log(L))) / 8
        NlogLqqo[2] = (L * (-3 + log(L / 2))) / 4
        NlogLqqo[3] = -(L * (1 + log(8) - 3 * log(L))) / 8
        intNsr[1] = -3 / 2
        intNsr[2] = 0
        intNsr[3] = 3 / 2
        intNsr2[1] = 9 / 4
        intNsr2[2] = -13 / 2
        intNsr2[3] = 9 / 4
    elseif abs(xi0) ≈ 2 / 3
        NlogLqqo[1] = (L * (-39 + 10 * log(5) - 27 * log(6) + 27 * log(L))) / 72
        NlogLqqo[2] = (L * (-5 + (25 * log(5)) / 6 - 3 * log(6) + 3 * log(L))) / 12
        NlogLqqo[3] = -(L * (3 + 25 * log(6 / 5) + log(36) - 27 * log(L))) / 72

        intNsr[1] = (3 * (-4 + (4 * log(5)) / 3)) / 4
        intNsr[2] = 3
        intNsr[3] = 0
        intNsr2[1] = (3 * (-9 / 5 - 3 * log(5))) / 4
        intNsr2[2] = (-9 + 6 * log(5)) / 2
        intNsr2[3] = (3 * (3 - log(5))) / 4
    else
        NlogLqqo[1] = (L * (-7 + 3 * log(L))) / 8
        NlogLqqo[2] = (L * log(L)) / 4
        NlogLqqo[3] = (L * (-1 + 3 * log(L))) / 8
        intNsr[1] = (3 * (-10 + 5 * log(2))) / 8
        intNsr[2] = 9 / 2 - (5 * log(8)) / 4
        intNsr[3] = (3 * (-2 + log(2))) / 8
        intNsr2[1] = 0
        intNsr2[2] = 0
        intNsr2[3] = 0
    end
    for x in [
        :d2Rdx2,
        :d2Rdx2ns,
        :d2Rdxdy,
        :d2Rdxdyns,
        :d2Rdxdys,
        :d2Rdy2,
        :d2Rdy2ns,
        :d2Rdy2s,
        :d2Sdx2,
        :d2Sdx2ns,
        :d2Sdx2s,
        :d2Sdxdy,
        :d2Sdxdyns,
        :d2Sdxdys,
        :d2Sdy2,
        :d2Sdy2ns,
        :d2Sdy2s,
        :d3Rdx2dy,
        :d3Rdx3,
        :d3Rdxdy2,
        :d3Rdy3,
        :d3Sdx2dy,
        :d3Sdx3,
        :d3Sdxdy2,
        :d3Sdy3,
        :d4Rdx2dy2,
        :d4Rdx3dy,
        :d4Rdx4,
        :d4Rdxdy3,
        :d4Rdy4,
        :d4Sdx2dy2,
        :d4Sdx3dy,
        :d4Sdx4,
        :d4Sdxdy3,
        :d4Sdy4,
        :d2Rdx2s,
    ]
        @eval $x = zeros(2)
    end
    for j = 1:3
        for i = 1:2
            constA = (cos(theta) + dad.k.d[i] * sin(theta))^2 + dad.k.e[i]^2 * sin(theta)^2
            constB = atan(dad.k.e[i] * sin(theta), (cos(theta) + dad.k.d[i] * sin(theta)))

            #----------------------------------------------------------------------
            # Cálculo das integrais das derivadas segundas de R
            #----------------------------------------------------------------------
            # Cálculo das integrais da derivada segunda de R1 e R2 em relação a x conforme Shi e
            # Bezine
            d2Rdx2s[i] = 4 * NlogLqqo[j]
            d2Rdx2ns[i] = 2 * log(constA / a^2) * intN[j] * L / 2
            d2Rdx2[i] = d2Rdx2s[i] + d2Rdx2ns[i]
            # Cálculo das integrais da derivada segunda de R1 e R2 em relação a x e y conforme Shi e Bezine
            d2Rdxdys[i] = 4 * dad.k.d[i] * NlogLqqo[j]
            d2Rdxdyns[i] =
                (2 * dad.k.d[i] * log(constA / a^2) - 4 * dad.k.e[i] * constB) *
                intN[j] *
                L / 2
            d2Rdxdy[i] = d2Rdxdys[i] + d2Rdxdyns[i]
            # Cálculo das integrais da derivada segunda de R1 e R2 em relação a y conforme Shi e Bezine
            d2Rdy2s[i] = 4 * (dad.k.d[i]^2 - dad.k.e[i]^2) * NlogLqqo[j]
            d2Rdy2ns[i] =
                (
                    2 * (dad.k.d[i]^2 - dad.k.e[i]^2) * log(a^2 * constA) -
                    8 * dad.k.d[i] * dad.k.e[i] * constB
                ) *
                intN[j] *
                L / 2
            d2Rdy2[i] = d2Rdy2s[i] + d2Rdy2ns[i]

            #--------------------------------------------------------------------------------
            # Cálculo das integrais das derivadas segundas de S
            #--------------------------------------------------------------------------------
            # Cálculo das integrais das derivadas segunda de S1 e S2 em relação a x conforme Shi e Bezine
            d2Sdx2s[i] = 0
            d2Sdx2ns[i] = 2 * constB * intN[j] * L / 2
            d2Sdx2[i] = d2Sdx2s[i] + d2Sdx2ns[i]
            # Cálculo das integrais da derivada segunda de S1 e S2 em relação a x e y conforme Shi e Bezine
            d2Sdxdys[i] = 2 * dad.k.e[i] * NlogLqqo[j]
            d2Sdxdyns[i] =
                (dad.k.e[i] * log(constA / a^2) + 2 * dad.k.d[i] * constB) * intN[j] * L / 2
            d2Sdxdy[i] = d2Sdxdys[i] + d2Sdxdyns[i]
            # Cálculo das integrais da derivada segunda de S1 e S2 em relação a y conforme Shi e Bezine
            d2Sdy2s[i] = 4 * dad.k.d[i] * dad.k.e[i] * NlogLqqo[j]
            d2Sdy2ns[i] =
                (
                    2 * dad.k.d[i] * dad.k.e[i] * log(constA / a^2) +
                    2 * (dad.k.d[i]^2 - dad.k.e[i]^2) * constB
                ) *
                intN[j] *
                L / 2
            d2Sdy2[i] = d2Sdy2s[i] + d2Sdy2ns[i]
            #--------------------------------------------------------------------------------
            # Cálculo das integrais das derivadas terceiras de R
            #--------------------------------------------------------------------------------
            # Cálculo das integrais da derivada terceira de R1 e R2 em relação a x conforme Shi e Bezine
            d3Rdx3[i] =
                4 * (2 * intNsr[j] / L) * L / 2 * (cos(theta) + dad.k.d[i] * sin(theta)) /
                (constA)
            # Cálculo das integrais da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
            d3Rdx2dy[i] =
                4 * (2 * intNsr[j] / L) * L / 2 * (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                    dad.k.e[i]^2 * sin(theta)
                ) / (constA)
            # Cálculo das integrais da derivada terceira de R1 e R2 em relação a x e y conforme Shi e Bezine
            d3Rdxdy2[i] =
                4 * (2 * intNsr[j] / L) * L / 2 * (
                    (dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^2 + dad.k.e[i]^2) * dad.k.d[i] * sin(theta)
                ) / (constA)
            # Cálculo das integrais da derivada terceira de R1 e R2 em relação a y conforme Shi e Bezine
            d3Rdy3[i] =
                4 * (2 * intNsr[j] / L) * L / 2 * (
                    dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                    (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                ) / (constA)
            #--------------------------------------------------------------------------------
            # Cálculo das derivadas quartas de R
            #--------------------------------------------------------------------------------
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a x conforme Shi e Bezine
            d4Rdx4[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 *
                ((cos(theta) + dad.k.d[i] * sin(theta))^2 - dad.k.e[i]^2 * sin(theta)^2) /
                (constA^2)
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdx3dy[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 * (
                    dad.k.d[i] / constA +
                    2 * dad.k.e[i]^2 * sin(theta) * cos(theta) /
                    (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    )^2
                )
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdx2dy2[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 * (
                    (dad.k.d[i]^2 + dad.k.e[i]^2) / (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ) - 2 * dad.k.e[i]^2 * cos(theta)^2 / constA^2
                )
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a x e y conforme Shi e Bezine
            d4Rdxdy3[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 * (
                    dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) / (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    ) -
                    2 *
                    dad.k.e[i]^2 *
                    cos(theta) *
                    (
                        2 * dad.k.d[i] * cos(theta) +
                        (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                    ) /
                    (
                        (cos(theta) + dad.k.d[i] * sin(theta))^2 +
                        dad.k.e[i]^2 * sin(theta)^2
                    )^2
                )
            # Cálculo das integrais da derivada quarta de R1 e R2 em relação a y conforme Shi e Bezine
            d4Rdy4[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 * (
                    (dad.k.d[i]^4 - dad.k.e[i]^4) / constA -
                    2 *
                    dad.k.e[i]^2 *
                    cos(theta) *
                    (
                        (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                        2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                    ) / constA^2
                )

            #--------------------------------------------------------------------------------
            # Cálculo das das integrais derivadas terceiras de S
            #--------------------------------------------------------------------------------
            # Cálculo das integrais da derivada terceira de S1 e S2 em relação a x conforme Shi e Bezine
            d3Sdx3[i] =
                -2 * (2 * intNsr[j] / L) * L / 2 * dad.k.e[i] * sin(theta) / (constA)
            # Cálculo das integrais da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
            d3Sdx2dy[i] =
                2 * (2 * intNsr[j] / L) * L / 2 * dad.k.e[i] * cos(theta) / (constA)
            # Cálculo das integrais da derivada terceira de S1 e S2 em relação a x e y conforme Shi e Bezine
            d3Sdxdy2[i] =
                2 * (2 * intNsr[j] / L) * L / 2 *
                dad.k.e[i] *
                (
                    2 * dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) -
                    (dad.k.d[i]^2 - dad.k.e[i]^2) * sin(theta)
                ) / (constA)
            # Cálculo das integrais da derivada terceira de S1 e S2 em relação a y conforme Shi e Bezine
            d3Sdy3[i] =
                2 * (2 * intNsr[j] / L) * L / 2 *
                dad.k.e[i] *
                (
                    (3 * dad.k.d[i]^2 - dad.k.e[i]^2) * cos(theta) +
                    2 * dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) * sin(theta)
                ) / (constA)
            #--------------------------------------------------------------------------------
            # Cálculo das derivadas quartas de S
            #--------------------------------------------------------------------------------
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a x conforme Shi e Bezine
            d4Sdx4[i] =
                4 * (4 * intNsr2[j] / L^2) * L / 2 *
                dad.k.e[i] *
                sin(theta) *
                (cos(theta) + dad.k.d[i] * sin(theta)) / (constA^2)
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdx3dy[i] =
                2 * (4 * intNsr2[j] / L^2) * L / 2 *
                dad.k.e[i] *
                (
                    1 / constA -
                    2 * cos(theta) * (cos(theta) + dad.k.d[i] * sin(theta)) / constA^2
                )
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdx2dy2[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 *
                dad.k.e[i] *
                cos(theta) *
                (
                    dad.k.d[i] * (cos(theta) + dad.k.d[i] * sin(theta)) +
                    dad.k.e[i]^2 * sin(theta)
                ) / (constA^2)
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a x e y conforme Shi e Bezine
            d4Sdxdy3[i] =
                -2 * dad.k.e[i] * (4 * intNsr2[j] / L^2) * L / 2 * (
                    (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                    (
                        2 *
                        (dad.k.d[i]^2 + dad.k.e[i]^2) *
                        cos(theta) *
                        (cos(theta) + dad.k.d[i] * sin(theta)) -
                        4 * dad.k.e[i]^2 * cos(theta)^2
                    ) / constA^2
                )
            # Cálculo das integrais da derivada quarta de S1 e S2 em relação a y conforme Shi e Bezine
            d4Sdy4[i] =
                -4 * (4 * intNsr2[j] / L^2) * L / 2 *
                dad.k.e[i] *
                (
                    dad.k.d[i] * (dad.k.d[i]^2 + dad.k.e[i]^2) / constA +
                    cos(theta) * (
                        dad.k.d[i] * (dad.k.d[i]^2 - 3 * dad.k.e[i]^2) * cos(theta) +
                        (dad.k.d[i]^4 - dad.k.e[i]^4) * sin(theta)
                    ) / constA^2
                )
        end
        # Integrais das derivadas segundas da solução fundamental
        d2wdx2 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d2Rdx2[1] + C2 * d2Rdx2[2] + C3 * (d2Sdx2[1] - d2Sdx2[2]))
        d2wdxdy =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d2Rdxdy[1] + C2 * d2Rdxdy[2] + C3 * (d2Sdxdy[1] - d2Sdxdy[2]))
        d2wdy2 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d2Rdy2[1] + C2 * d2Rdy2[2] + C3 * (d2Sdy2[1] - d2Sdy2[2]))
        # Integrais das derivadas terceira da solução fundamental
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
        # Integrais das derivadas quarta da solução fundamental
        d4wdx4 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdx4[1] + C2 * d4Rdx4[2] + C3 * (d4Sdx4[1] - d4Sdx4[2]))
        d4wdx3dy =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdx3dy[1] + C2 * d4Rdx3dy[2] + C3 * (d4Sdx3dy[1] - d4Sdx3dy[2]))
        d4wdx2dy2 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdx2dy2[1] + C2 * d4Rdx2dy2[2] + C3 * (d4Sdx2dy2[1] - d4Sdx2dy2[2]))
        d4wdxdy3 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdxdy3[1] + C2 * d4Rdxdy3[2] + C3 * (d4Sdxdy3[1] - d4Sdxdy3[2]))
        d4wdy4 =
            1 / (8 * pi * dad.k.D22) *
            (C1 * d4Rdy4[1] + C2 * d4Rdy4[2] + C3 * (d4Sdy4[1] - d4Sdy4[2]))
        # Cálculo de f1; f2 e f3 conforme Shi e Bezine
        f1 = dad.k.D11 * nx^2 + 2 * dad.k.D16 * nx * ny + dad.k.D12 * ny^2
        f2 = 2 * (dad.k.D16 * nx^2 + 2 * dad.k.D66 * nx * ny + dad.k.D26 * ny^2)
        f3 = dad.k.D12 * nx^2 + 2 * dad.k.D26 * nx * ny + dad.k.D22 * ny^2
        # Cálculo de h1; h2; h3 e h4 conforme Shi e Bezine
        h1 = dad.k.D11 * nx * (1 + ny^2) + 2 * dad.k.D16 * ny^3 - dad.k.D12 * nx * ny^2
        h2 =
            4 * dad.k.D16 * nx + dad.k.D12 * ny * (1 + nx^2) + 4 * dad.k.D66 * ny^3 -
            dad.k.D11 * nx^2 * ny - 2 * dad.k.D26 * nx * ny^2
        h3 =
            4 * dad.k.D26 * ny + dad.k.D12 * nx * (1 + ny^2) + 4 * dad.k.D66 * nx^3 -
            dad.k.D22 * nx * ny^2 - 2 * dad.k.D16 * nx^2 * ny
        h4 = dad.k.D22 * ny * (1 + nx^2) + 2 * dad.k.D26 * nx^3 - dad.k.D12 * nx^2 * ny
        #--------------------------------------------------------------------------------
        # Cálculo das demais variáveis da solução fundamental
        #--------------------------------------------------------------------------------
        mn = -(f1 * d2wdx2 + f2 * d2wdxdy + f3 * d2wdy2)
        vn = -(h1 * d3wdx3 + h2 * d3wdx2dy + h3 * d3wdxdy2 + h4 * d3wdy3)
        # Cálculo das derivadas primeiras de Mn
        dmndx = -(f1 * d3wdx3 + f2 * d3wdx2dy + f3 * d3wdxdy2)
        dmndy = -(f1 * d3wdx2dy + f2 * d3wdxdy2 + f3 * d3wdy3)
        # Cálculo das derivadas primeiras de Vn
        dvndx = -(h1 * d4wdx4 + h2 * d4wdx3dy + h3 * d4wdx2dy2 + h4 * d4wdxdy3)
        dvndy = -(h1 * d4wdx3dy + h2 * d4wdx2dy2 + h3 * d4wdxdy3 + h4 * d4wdy4)
        #--------------------------------------------------------------------------------
        # Cálculo da derivada da solução fundamental
        #--------------------------------------------------------------------------------
        #--------------------------------------------------------------------------------
        # Cálculo das demais variáveis da derivada da solução fundamental
        #--------------------------------------------------------------------------------
        m1 = nx
        m2 = ny
        d2wdndm = -(d2wdx2 * nx * m1 + d2wdxdy * (nx * m2 + ny * m1) + d2wdy2 * ny * m2)
        dmndm = -(dmndx * m1 + dmndy * m2)
        dvndm = -(dvndx * m1 + dvndy * m2)
        # Atribuição das variáveis da solução fundamental - atribuição de w_ij à variável
        # de saída w_est
        g22 = -d2wdndm
        h11 = vn
        h22 = -dmndm
        h12 = -mn
        h21 = dvndm
        # @infiltrate
        if xi0 ≈ -2 / 3
            h_el[:, 2*j-1:2*j] = [
                h11 h12
                h21 h22
            ]
            g_el[j] = g22
        elseif xi0 ≈ 0
            h_el[:, 2*j-1:2*j] = [
                h11 h12
                h21 h22
            ]
            g_el[j] = g22
        elseif xi0 ≈ -1
            h_el[1, 2*j-1:2*j] = [h11 h12]
        elseif xi0 ≈ 1
            h_el[1, 7-2*j:8-2*j] = [h11 h12]
        else
            h_el[:, 7-2*j:8-2*j] = [
                h11 h12
                h21 h22
            ]
            g_el[4-j] = g22
        end
    end
    return h_el, g_el
end

function simula_placa_tempo(d::Dict)
    reset_timer!()
    @unpack nelem, NPX, npg, problema, metodo, nt = d
    NPY = NPX
    tdad = @timed dad = format_dad(
        getfield(Main, Symbol(problema))(nelem, 3, nt),
        NPX,
        NPY,
        canto = true,
    ) # dados
    # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
    # println("2. Montando a matriz A e o vetor b")
    tHeG = @timed H, G, q = calc_HeG(dad, npg)  #importante
    # println("3. Montando a matriz M")
    # M = BEM.Monta_M_RIM(dad, 10, 20)
    if metodo == "DRM"
        tM = @timed M = getfield(Main, Symbol(metodo))(dad, H, G)
    else
        tM = @timed M = getfield(Main, Symbol(metodo))(dad, npg)
    end
    dM = BEM.Monta_dM_RIMd(dad, 10)
    # mats = matread("data\\exp_raw\\matrizes.mat")
    # M = mats["M2"]
    # H1 = mats["Hg"]
    # G1 = mats["Gg"]
    # @infiltrate
    M = M * dad.k.rho * dad.k.thickness
    H_td = 2 / dad.k.dt^2 * M + H

    A, B, bc_val = aplicaCDC(H_td, G, dad) # Calcula a matriz A e o vetor b
    tpasso = @timed dis, tens, tens_i = passonotempo(dad, A, B, q, bc_val, M, dM, npg)
    # dis / dad.k.w_st, range(0, step=dad.k.dt, length=dad.k.n_dt) / dad.k.to
    centro = 2 * size(dad.NOS, 1) + ceil(Int64, size(dad.pontos_internos, 1) / 2)
    centroi = ceil(Int64, size(dad.pontos_internos, 1) / 2)
    # ws = autovalor(H, G, M, dad)

    fulld = copy(d)
    fulld["dis"] = dis[:, centro] / dad.k.w_st
    fulld["t"] = range(0, step = dad.k.dt, length = dad.k.n_dt) / dad.k.to
    fulld["timerM"] = tM.time
    fulld["timerHeG"] = tHeG.time
    fulld["timerpasso"] = tpasso.time
    fulld["timerdad"] = tdad.time
    fulld["tens"] = tens[:, centroi, 1]
    # fulld["timerdad"] = t
    # fulld["frequencias"] = ws
    fulld
end

function autovalor(H, G, M, dad::Union{placa_fina,placa_fina_isotropica})

    ind_deslconhecido = zeros(Int, 0)
    nc = size(dad.NOS, 1)
    ncanto = size(dad.k.bc_canto, 1)
    nh = size(H, 2)
    for elem_i in dad.ELEM, i = 1:2  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        # @show elem_i.tipoCDC[i]
        if elem_i.tipoCDC[i] == 0
            ind_deslconhecido = [ind_deslconhecido; 2ind_elem .+ (i - 2)]
        end
    end
    for c = 1:ncanto
        if dad.k.bc_canto[c, 1] == 0
            ind_deslconhecido = [ind_deslconhecido; nh - ncanto + c]
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

    # if num

    #   r = lobpcg(Hb, Mb, false, 1)
    #   a = r.λ
    #   v = r.X
    #   u[ind_forcaconhecida] = real.(v)
    #   return a, u
    # end

    a, v = eigen(Hb, Mb)
    ~, ind = findmin(abs, a)
    u[ind_forcaconhecida] = real.(v[:, ind])
    # @infiltrate

    real.(sort(abs.(a))[1:3]), [u[1:2:2nc]; u[2nc+1:end-ncanto]]
end


function calc_HeGeIt(dad::Union{placa_fina,placa_fina_isotropica}, npg = 8)
    n_fis = size(dad.NOS, 1)
    n_internos = size(dad.pontos_internos, 1)
    n_cantos = size(dad.k.cantos, 1)
    H = zeros(2 * n_fis + n_cantos + n_internos, 2 * n_fis + n_cantos + n_internos)
    G = zeros(2 * n_fis + n_cantos + n_internos, 2 * n_fis + n_cantos + n_internos)
    It = zeros(2 * n_fis + n_cantos + n_internos, n_fis + n_cantos + n_internos)
    dNs = zeros(n_fis, n_fis)
    q = zeros(2 * n_fis + n_cantos + n_internos)
    qsi, w = gausslegendre(npg)# Quadratura de gauss
    # qsi2, w2 = gausslegendre(2npg)# Quadratura de gauss
    normal_fonte = dad.normal
    pre = [zeros(2) for idx = 1:30]

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
            Δelem = (x[end, :] - x[1, :])  # Δx e Δy entre o primeiro e ultimo nó geometrico
            eet =
                (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                norm(Δelem)^2 + elem_j.ξs[1]
            N_geo, dN = calc_fforma(eet, elem_j)
            ps = N_geo' * x
            b = norm(ps' - pf) / norm(Δelem)
            # eta, Jt = qsi, 1
            eta, Jt = sinhtrans(qsi, eet, b)
            # h, g = integraelem(pf, nf, x, qsi, w, elem_j, dad, pre)
            # @infiltrate
            h, g, it = integraelem2(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
            # @infiltrate
            # @show h
            # @infiltrate
            if caso == "contorno"
                nosing = elem_j.indices .== i
                if sum(nosing) == 1
                    # @infiltrate
                    dNs[i, elem_j.indices] = dN
                    no_pf = findfirst(nosing)
                    xi0 = elem_j.ξs[no_pf]
                    # h, g = integraelemsing(pf, nf, x, qsi2, w2, elem_j, dad, xi0)

                    # h, g2 = integraelemsing(x, dad, xi0)
                    # g[2, 2:2:end] = g2
                    h, g = integraelemsing_num(pf, nf, x, elem_j, dad, pre, xi0, 30)
                    # @infiltrate
                    # @show h, hn
                end
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                H[2i-1:2i, cols] = h
                G[2i-1:2i, cols] = g
                It[2i-1:2i, elem_j.indices] = it

                q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                # @infiltrate
                # q_el = compute_q(pf, nf, x, qsi, w, elem_j, dad, pre)
                q[2i-1:2i] += q_el
            elseif caso == "canto"
                if j == dad.k.cantos[i-n_fis-n_internos, 8]
                    # dNs[i, elem_j.indices] = dN / 2
                    xi0 = 1.0
                    h, g2 = integraelemsing(x, dad, xi0)
                    g[2, 2:2:end] = g2
                elseif j == dad.k.cantos[i-n_fis-n_internos, 9]
                    # dNs[i, elem_j.indices] = dN / 2
                    xi0 = -1.0
                    h, g2 = integraelemsing(x, dad, xi0)
                    g[2, 2:2:end] = g2
                end
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                ind = i - n_fis + 2n_fis
                H[ind, cols] = h[1, :]
                G[ind, cols] = g[1, :]
                q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                q[ind] += q_el[1]
            else
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                ind = i - n_fis + 2n_fis
                H[ind, cols] = h[1, :]
                G[ind, cols] = g[1, :]
                q_el = compute_q(pf, nf, x, eta, w .* Jt, elem_j, dad, pre)
                q[ind] += q_el[1]
            end
        end
        if caso == "contorno"
            R, W = compute_Rw(pf, nf, dad, pre)
            cols = (n_fis * 2 + n_internos) .+ (1:n_cantos)
            H[2i-1:2i, cols] = R
            G[2i-1:2i, cols] = W
        else
            R, W = compute_Rw(pf, nf, dad, pre)
            ind = i + n_fis
            cols = (n_fis * 2 + n_internos) .+ (1:n_cantos)
            H[ind, cols] = R[1, :]
            G[ind, cols] = W[1, :]
        end
    end
    # for i = 1:n#i=1:size(dad.NOS,1) #Laço dos pontos fontes
    # H[2i-1:2i,2i-1:2i].=0
    # H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
    # end
    # @infiltrate
    H[1:2n_fis, 1:2n_fis] += I / 2
    H[(2n_fis).+(1:n_internos), (2n_fis).+(1:n_internos)] += I
    # H[2n_fis+n_internos+1:end, 2n_fis+n_internos+1:end] += I / 4
    H, G, q, It, dNs
end

function aplicaT(G, iT, T, dNs, M3, N)
    n = size(T, 1)
    n2 = size(N, 1)
    Nit = 0 * G
    # @infiltrate
    Nit[:, 1:2:2n] =
        G[:, 1:2:2n] .* T[:, 2]' * dNs - G[:, 2:2:2n] .* T[:, 1]' - iT[:, 1:n] .* T[:, 2]' +
        M3[:, 1:3:3n] .* N[1:n, 1]' +
        M3[:, 2:3:3n] .* N[1:n, 2]' +
        M3[:, 3:3:3n] .* N[1:n, 3]'
    Nit[:, 2n+1:n+n2] =
        M3[:, 3n+1:3:3n2] .* N[n+1:n2, 1]' +
        M3[:, 3n+2:3:3n2] .* N[n+1:n2, 2]' +
        M3[:, 3n+3:3:3n2] .* N[n+1:n2, 3]'
    Nit[:, 2:2:2n] = G[:, 1:2:2n] .* T[:, 1]'
    Nit
end
function aplicaT(dad::Union{placa_fina,placa_fina_isotropica}, G, T, dNs, M2, N)
    n = size(T, 1)
    NOS = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    Nit = 0 * G
    K = zeros(2 * size(NOS, 1), size(NOS, 1))
    F, Fx1, Fy1 = BEM.montaFs(
        [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]],
        [dad.NOS; dad.pontos_internos],
    )
    F1, Fx, Fy = BEM.montaFs(
        [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]],
        [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]],
    )
    # @infiltrate
    Ni = F * N

    K[1:2:end, :] = Ni[:, 1] .* Fx + Ni[:, 3] .* Fy
    K[2:2:end, :] = Ni[:, 2] .* Fy + Ni[:, 3] .* Fx

    k2 = M2 * K

    # @infiltrate
    Nit[:, 1:2:2n] = G[:, 1:2:2n] .* T[:, 2]' * dNs - k2[:, 1:n]
    Nit[:, 2n+1:end] = -k2[:, n+1:end]
    Nit[:, 2:2:2n] = G[:, 1:2:2n] .* T[:, 1]'
    Nit
end
function aplicaT(dad::Union{placa_fina,placa_fina_isotropica}, G, T, dNs, M2, Ni, Fx, Fy)
    n = size(T, 1)
    NOS = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    Nit = 0 * G
    K = zeros(2 * size(NOS, 1), size(NOS, 1))
    # @infiltrate
    # Ni = F * N

    K[1:2:end, :] = Ni[:, 1] .* Fx + Ni[:, 3] .* Fy
    K[2:2:end, :] = Ni[:, 2] .* Fy + Ni[:, 3] .* Fx

    k2 = M2 * K

    # @infiltrate
    Nit[:, 1:2:2n] = G[:, 1:2:2n] .* T[:, 2]' * dNs - k2[:, 1:n]
    Nit[:, 2n+1:end] = -k2[:, n+1:end]
    Nit[:, 2:2:2n] = G[:, 1:2:2n] .* T[:, 1]'
    Nit
end
function aplicaT(dad::Union{placa_fina,placa_fina_isotropica}, M, N)
    # @infiltrate
    n2 = size(N, 1)
    F1, Fx1, Fy1 = BEM.montaFs(
        [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]],
        [dad.NOS; dad.pontos_internos],
    )
    F, Fx, Fy = BEM.montaFs(
        [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]],
        [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]],
    )
    Fyx = Fy * Fx
    Fxy = Fx * Fy
    Fyy = Fy * Fy
    Fxx = Fx * Fx
    Nx = Fx1 * N
    Ny = Fy1 * N
    # K = (Nx[1] + Ny[3]) .* Fx + (Nx[3] + Ny[2]) .* Fy + (N[1] .* Fxx + N[3] .* Fxy + N[2] .* Fyy + N[3] .* Fyx)
    K =
        (Nx[:, 1] + Ny[:, 3])' .* Fx +
        (Nx[:, 3] + Ny[:, 2])' .* Fy +
        (
            (F1 * N[:, 1])' .* Fxx +
            (F1 * N[:, 3])' .* Fxy +
            (F1 * N[:, 2])' .* Fyy +
            (F1 * N[:, 3])' .* Fyx
        )

    k2 = zeros(
        2nc(dad) + ni(dad) + size(dad.k.cantos, 1),
        2nc(dad) + ni(dad) + size(dad.k.cantos, 1),
    )
    k2[1:2:2nc(dad), 1:2:2nc(dad)] = K[1:nc(dad), 1:nc(dad)]
    k2[2:2:2nc(dad), 1:2:2nc(dad)] = K[1:nc(dad), 1:nc(dad)]
    k2[1:2:2nc(dad), 2nc(dad)+1:end] = K[1:nc(dad), nc(dad)+1:end]
    k2[2:2:2nc(dad), 2nc(dad)+1:end] = K[1:nc(dad), nc(dad)+1:end]
    k2[2nc(dad)+1:end, 1:2:2nc(dad)] = K[nc(dad)+1:end, 1:nc(dad)]
    k2[2nc(dad)+1:end, 2nc(dad)+1:end] = K[nc(dad)+1:end, nc(dad)+1:end]
    # @infiltrate
    M * k2
end
function simula_placa_flambagem(d::Dict, ana = [], Nexato = [])
    reset_timer!()
    @unpack nelem, NPX, npg, problema, bc = d
    NPY = NPX
    entrada = getfield(Main, Symbol(problema))(nelem, 3, bc)

    dad = format_dad(entrada[1], NPX, NPY, canto = true) # dados
    dadpe = format_dad(entrada[2], NPX, NPY) # dados
    # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
    # println("2. Montando a matriz A e o vetor b")
    H, G, q, It, dNs = calc_HeGeIt(dad, npg)  #importante
    tdMd = @timed dMd = BEM.Monta_dM_RIMd(dad, npg)
    tMd = @timed Md = BEM.Monta_M_RIMd(dad, npg)
    tdM = @timed dM = BEM.Monta_dM_RIM(dad, npg)
    tM = @timed M = BEM.Monta_M_RIM(dad, npg)
    fulld = copy(d)


    # @infiltrate
    # @show dad.pontos_internos
    if isempty(Nexato)
        tHeGpl = @timed Hpe, Gpe = calc_HeG(dadpe, npg)  #importante
        fulld["timerHeGpl"] = tHeGpl.time
        # nosrestritos = [floor(Int, nelem / 2)+2 1
        # floor(Int, nelem / 2)+2+nelem*3 2
        # floor(Int, nelem / 2)+2+nelem*6 1
        # floor(Int, nelem / 2)+2+nelem*9 2]
        nosrestritos = [1 1]
        A, b = aplicaCDC(Hpe, Gpe, dadpe, nosrestritos) # Calcula a matriz A e o vetor b
        # @infiltrate
        x = A \ b
        u, t = separa(dadpe, x, nosrestritos) #importante
        tens_cont, tens_nt = calc_tens_cont(dadpe, u, t)
        tens_int = calc_tens_int(dadpe, u, t)
        tens = [tens_cont; tens_int]

    else
        tens_nt, tens = criatensoes(dad, Nexato)
    end
    Mit2 = BEM.aplicaT(G, It, tens_nt, dNs, dM[6], tens)
    Mit1 = BEM.aplicaT(dad, G, tens_nt, dNs, dMd[1], tens)
    Mit = BEM.aplicaT(dad, M, tens)

    Mit2d = BEM.aplicaT(G, It, tens_nt, dNs, dMd[6], tens)
    Mit1d = BEM.aplicaT(dad, G, tens_nt, dNs, dMd[1], tens)
    Mitd = BEM.aplicaT(dad, Md, tens)
    # @infiltrate
    fulld = copy(d)

    entrada = getfield(Main, Symbol(problema))(nelem, 3, bc)

    dad = format_dad(entrada[1], NPX, NPY, canto = true)

    a, v = BEM.autovalor(H, G, Mit, dad)
    a1, v = BEM.autovalor(H, G, Mit1, dad)
    a2, v = BEM.autovalor(H, G, Mit2, dad)

    ad, v = BEM.autovalor(H, G, Mitd, dad)
    a1d, v = BEM.autovalor(H, G, Mit1d, dad)
    a2d, v = BEM.autovalor(H, G, Mit2d, dad)
    # @infiltrate


    # D3 = dad.k.D12 + 2 * dad.k.D66
    # k = lambda / D3 / pi^2
    k = a / dad.k.D / pi^2
    k1 = a1 / dad.k.D / pi^2
    k2 = a2 / dad.k.D / pi^2
    kd = ad / dad.k.D / pi^2
    k1d = a1d / dad.k.D / pi^2
    k2d = a2d / dad.k.D / pi^2


    @show k, k1, k2, kd, k1d, k2d
    # fulld["bc"] = bc

    # @infiltrate
    fulld["k_Purbolaksono"] = k
    fulld["k_Wrobel"] = k1
    fulld["k_doval"] = k2
    fulld["kd_Purbolaksono"] = kd
    fulld["kd_Wrobel"] = k1d
    fulld["kd_doval"] = k2d

    fulld["e_k_Purbolaksono"] = k - ana[bc]
    fulld["e_k_Wrobel"] = k1 - ana[bc]
    fulld["e_k_doval"] = k2 - ana[bc]
    fulld["e_kd_Purbolaksono"] = kd - ana[bc]
    fulld["e_kd_Wrobel"] = k1d - ana[bc]
    fulld["e_kd_doval"] = k2d - ana[bc]

    fulld["timerMd"] = tMd.time
    fulld["timerdMd"] = tdMd.time
    fulld["timerM"] = tM.time
    fulld["timerdM"] = tdM.time
    # fulld["timerHeG"] = tHeG.time
    # fulld["timerdad"] = tdad.time + tdadpe.time
    fulld
end

function flambagem(fulld, d, H, G, It, dNs, dMd, Md, dM, M, Nexato = [], ana = [])
    reset_timer!()
    @unpack nelem, NPX, npg, problema, bc = d
    NPY = NPX
    entrada = getfield(Main, Symbol(problema))(nelem, 3, bc)
    dadpe = format_dad(entrada[2], NPX, NPY) # dados

    dad = format_dad(entrada[1], NPX, NPY, canto = true) # dados
    # @show dad.pontos_internos
    if isempty(Nexato)
        tHeGpl = @timed Hpe, Gpe = calc_HeG(dadpe, npg)  #importante
        fulld["timerHeGpl"] = tHeGpl.time
        # nosrestritos = [floor(Int, nelem / 2)+2 1
        # floor(Int, nelem / 2)+2+nelem*3 2
        # floor(Int, nelem / 2)+2+nelem*6 1
        # floor(Int, nelem / 2)+2+nelem*9 2]
        nosrestritos = [1 1]
        A, b = aplicaCDC(Hpe, Gpe, dadpe, nosrestritos) # Calcula a matriz A e o vetor b
        # @infiltrate
        x = A \ b
        u, t = separa(dadpe, x, nosrestritos) #importante
        tens_cont, tens_nt = calc_tens_cont(dadpe, u, t)
        tens_int = calc_tens_int(dadpe, u, t)
        tens = [tens_cont; tens_int]

    else
        tens_nt, tens = criatensoes(dad, Nexato)
    end
    Mit2 = BEM.aplicaT(G, It, tens_nt, dNs, dM[6], tens)
    Mit1 = BEM.aplicaT(dad, G, tens_nt, dNs, dM[1], tens)
    Mit = BEM.aplicaT(dad, M, tens)

    Mit2d = BEM.aplicaT(G, It, tens_nt, dNs, dMd[6], tens)
    Mit1d = BEM.aplicaT(dad, G, tens_nt, dNs, dMd[1], tens)
    Mitd = BEM.aplicaT(dad, Md, tens)
    # @infiltrate

    entrada = getfield(Main, Symbol(problema))(nelem, 3, bc)

    dad = format_dad(entrada[1], NPX, NPY, canto = true)

    a, v = BEM.autovalor(H, G, Mit, dad)
    a1, v = BEM.autovalor(H, G, Mit1, dad)
    a2, v = BEM.autovalor(H, G, Mit2, dad)

    ad, v = BEM.autovalor(H, G, Mitd, dad)
    a1d, v = BEM.autovalor(H, G, Mit1d, dad)
    a2d, v = BEM.autovalor(H, G, Mit2d, dad)
    # @infiltrate


    # D3 = dad.k.D12 + 2 * dad.k.D66
    # k = lambda / D3 / pi^2
    k = a / dad.k.D / pi^2
    k1 = a1 / dad.k.D / pi^2
    k2 = a2 / dad.k.D / pi^2
    kd = ad / dad.k.D / pi^2
    k1d = a1d / dad.k.D / pi^2
    k2d = a2d / dad.k.D / pi^2


    # @show k, k1, k2, kd, k1d, k2d
    # fulld["bc"] = bc

    # @infiltrate
    fulld["k_Purbolaksono"] = k
    fulld["bc"] = bc
    fulld["k_Wrobel"] = k1
    fulld["k_doval"] = k2
    fulld["kd_Purbolaksono"] = kd
    fulld["kd_Wrobel"] = k1d
    fulld["kd_doval"] = k2d

    # fulld["e_k_Purbolaksono"] = k - ana[bc]
    # fulld["e_k_Wrobel"] = k1 - ana[bc]
    # fulld["e_k_doval"] = k2 - ana[bc]
    # fulld["e_kd_Purbolaksono"] = kd - ana[bc]
    # fulld["e_kd_Wrobel"] = k1d - ana[bc]
    # fulld["e_kd_doval"] = k2d - ana[bc]

    # fulld["timerHeG"] = tHeG.time
    # fulld["timerdad"] = tdad.time + tdadpe.time
    fulld
end

function Matrizes_placa_flambagem(d::Dict, ana = [], Nexato = [])
    reset_timer!()
    @unpack nelem, NPX, npg, problema, bc = d
    NPY = NPX
    entrada = getfield(Main, Symbol(problema))(nelem, 3, bc)

    dad = format_dad(entrada[1], NPX, NPY, canto = true) # dados
    dadpe = format_dad(entrada[2], NPX, NPY) # dados
    # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
    # println("2. Montando a matriz A e o vetor b")
    H, G, q, It, dNs = calc_HeGeIt(dad, npg)  #importante
    tdMd = @timed dMd = BEM.Monta_dM_RIMd(dad, npg)
    tMd = @timed Md = BEM.Monta_M_RIMd(dad, npg)
    # tdM = @timed dM = BEM.Monta_dM_RIM(dad, npg)
    # tM = @timed M = BEM.Monta_M_RIM(dad, npg)
    tdM = tdMd
    tM = tMd
    dM = dMd
    M = Md

    fulld = copy(d)
    fulld["timerMd"] = tMd.time
    fulld["timerdMd"] = tdMd.time
    fulld["timerM"] = tM.time
    fulld["timerdM"] = tdM.time

    fulld, H, G, q, It, dNs, dMd, Md, dM, M
end
"
tens_nt, tens= criatensoes(dad, [-1,0,0])
"
function criatensoes(dad, Nexato)
    n_fis = size(dad.NOS, 1)
    n_internos = size(dad.pontos_internos, 1)
    n_cantos = size(dad.k.cantos, 1)
    tens_nt = zeros(n_fis, 3)
    tens = zeros(n_fis + n_internos, 3)
    tens[:, 1] .= Nexato[1]
    tens[:, 2] .= Nexato[2]
    tens[:, 3] .= Nexato[3]
    for elem_j in dad.ELEM  #Laço dos elementos

        x = dad.NOS[elem_j.indices, :]

        for i = 1:length(elem_j.indices)
            N, dN = calc_fforma(elem_j.ξs[i], elem_j)
            dxdqsi = dN' * x
            dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
            sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
            sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

            normal = [sy, -sx]

            # lij = [normal'
            #   -normal[2] normal[1]]

            # s = [Nexato[1] Nexato[3]; Nexato[3] Nexato[2]]
            # sigma = lij * s * lij'
            # sigma = lij'* s * lij
            # tx=valuex*nx;
            #      ty=valuey*ny;
            #      txy=valuex*ny;
            # @infiltrate
            tens_nt[elem_j.indices[i], :] = [
                Nexato[1] * normal[1] + Nexato[3] * normal[2],
                Nexato[2] * normal[2] + Nexato[3] * normal[1],
                Nexato[3] * normal[2],
            ]  #normal,  tangente
        end
    end
    tens_nt, tens
end
function criatensoesnt(dad, tens)
    n_fis = size(dad.NOS, 1)
    tens_nt = zeros(n_fis, 2)

    for elem_j in dad.ELEM  #Laço dos elementos

        x = dad.NOS[elem_j.indices, :]

        for i = 1:length(elem_j.indices)
            N, dN = calc_fforma(elem_j.ξs[i], elem_j)
            dxdqsi = dN' * x
            dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
            sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
            sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

            normal = [sy, -sx]

            # @infiltrate
            tens_nt[elem_j.indices[i], :] = [
                tens[i, 1] * normal[1] + tens[i, 3] * normal[2],
                tens[i, 2] * normal[2] + tens[i, 3] * normal[1],
            ]  #normal,  tangente
        end
    end
    tens_nt
end
function calsolfund2(pg, pf, n, nf, dad::Union{placa_fina_isotropica}, pre)
    D = dad.k.D
    ni = dad.k.nu
    rs = pg - pf
    n1 = n[1]
    n2 = n[2]
    m1 = nf[1]
    m2 = nf[2]
    s1 = -n2
    s2 = n1
    # Distance of source and field points
    r1 = rs[1]
    r2 = rs[2]
    r = norm(rs)
    # Thin plate fundamental solutions
    theta = atan(r2, r1)
    # Tangential vector
    rd1 = r1 / r
    rd2 = r2 / r

    # Inner products between vectors
    mr = m1 * rd1 + m2 * rd2
    nr = n1 * rd1 + n2 * rd2
    mn = n1 * m1 + n2 * m2
    ms = m1 * s1 + m2 * s2
    sr = s1 * rd1 + s2 * rd2

    # Fundamental solution of the first integral equation
    w = r^2 / (8 * pi * D) * (log(r) - 1 / 2)
    dwdn = r / (4 * pi * D) * log(r) * nr
    mnn = -1 / (4 * pi) * ((1 + ni) * log(r) + (1 - ni) * nr^2 + ni)
    vn = nr / (4 * pi * r) * (2 * (1 - ni) * sr^2 - 3 + ni)

    # Fundamental solution of the second integral equation
    dwdm = -r / (4 * pi * D) * log(r) * mr
    d2wdndm = -1 / (4 * pi * D) * (mr * nr + mn * log(r))
    dvndm =
        -1 / (4 * pi * r^2) * (
            2 * (1 - ni) * sr * (sr * mn + 2 * nr * ms - 4 * sr * mr * nr) +
            (3 - ni) * (-mn + 2 * mr * nr)
        )
    dmndm = 1 / (4 * pi * r) * ((1 + ni) * mr - 2 * (1 - ni) * nr * (mr * nr - mn))

    theta = atan(r2, r1)

    dwdx = -(r * cos(theta) * (-log(r^2))) / (8 * pi * D)
    dwdy = -(r * (-log(r^2)) * sin(theta)) / (8 * pi * D)

    d2wdx2 = (1 + cos(2 * theta) + log(r^2)) / (8 * pi * D)
    d2wdxdy = sin(2 * theta) / (8 * pi * D)
    d2wdy2 = -(-1 + cos(2 * theta) - log(r^2)) / (8 * pi * D)

    dwds = [dwdx dwdy] * [s1; s2]
    # d2wdmds = [d2wdx2 d2wdxdy; d2wdxdy d2wdy2] * [n1; n2] * [s1 s2]
    d2wdmds = (d2wdx2 * m1 * s1 + d2wdxdy * (m1 * s2 + m2 * s1) + d2wdy2 * m2 * s2)
    # @infiltrate
    u = [
        w -dwdn
        dwdm -d2wdndm
    ]
    p = [
        vn -mnn
        dvndm -dmndm
    ]
    # @infiltrate
    u2 = [
        -dwds
        -d2wdmds
    ]
    u, p, u2
    # @infiltrate
end

function integraelemsing(x, dad::placa_fina_isotropica, xi0)
    # Cálculo da distância do ponto fonte [xf,yf] ao ponto campo
    dx = (x[3, 1] - x[1, 1]) / 2 * 3
    dy = (x[3, 2] - x[1, 2]) / 2 * 3
    L = sqrt(dx^2 + dy^2)
    sx = dx / L
    sy = dy / L
    nx = sy
    ny = -sx
    # Ângulo entre r2 e r1
    theta = atan(dy, dx)
    D = dad.k.D
    nu = dad.k.nu
    D11 = D
    D22 = D
    D12 = nu * D
    D16 = 0
    D26 = 0
    D66 = (1 - nu^2) * D / (2 * (1 + nu))

    h_el = zeros(2, 6)
    g_el = zeros(1, 3)
    # Integrate[N1d,(qsi,-1,1)]
    for x in [:NlogLqqo, :intNsr2, :intNsr, :intN]
        @eval $x = zeros(3)
    end
    intN[1] = 3 / 4
    intN[2] = 1 / 2
    intN[3] = intN[1]

    if xi0 ≈ 0
        #   2*Integrate[N1d[qsi] Log[L [Abs[qsi - qsio]]/2],(qsi,-1,1)]*Jac
        #  Jac=L/2
        NlogLqqo[1] = -(L * (1 + log(8) - 3 * log(L))) / 8
        NlogLqqo[2] = (L * (-3 + log(L / 2))) / 4
        NlogLqqo[3] = -(L * (1 + log(8) - 3 * log(L))) / 8
        intNsr[1] = -3 / 2
        intNsr[2] = 0
        intNsr[3] = 3 / 2
        intNsr2[1] = 9 / 4
        intNsr2[2] = -13 / 2
        intNsr2[3] = 9 / 4
    elseif abs(xi0) ≈ 2 / 3
        NlogLqqo[1] = (L * (-39 + 10 * log(5) - 27 * log(6) + 27 * log(L))) / 72
        NlogLqqo[2] = (L * (-5 + (25 * log(5)) / 6 - 3 * log(6) + 3 * log(L))) / 12
        NlogLqqo[3] = -(L * (3 + 25 * log(6 / 5) + log(36) - 27 * log(L))) / 72

        intNsr[1] = (3 * (-4 + (4 * log(5)) / 3)) / 4
        intNsr[2] = 3
        intNsr[3] = 0
        intNsr2[1] = (3 * (-9 / 5 - 3 * log(5))) / 4
        intNsr2[2] = (-9 + 6 * log(5)) / 2
        intNsr2[3] = (3 * (3 - log(5))) / 4
    else
        NlogLqqo[1] = (L * (-7 + 3 * log(L))) / 8
        NlogLqqo[2] = (L * log(L)) / 4
        NlogLqqo[3] = (L * (-1 + 3 * log(L))) / 8
        intNsr[1] = (3 * (-10 + 5 * log(2))) / 8
        intNsr[2] = 9 / 2 - (5 * log(8)) / 4
        intNsr[3] = (3 * (-2 + log(2))) / 8
        intNsr2[1] = 0
        intNsr2[2] = 0
        intNsr2[3] = 0
    end

    for j = 1:3
        # Integrais das derivadas segundas da solução fundamental
        d2wdx2 =
            1 / (8 * pi * D) * ((1 + cos(2 * theta)) * intN[j] * L / 2 + 2 * NlogLqqo[j])
        d2wdxdy = 1 / (8 * pi * D) * sin(2 * theta) * intN[j] * L / 2
        d2wdy2 = -1 / (8 * pi * D) * ((-1 + cos(2 * theta)) * intN[j] * L - 2 * NlogLqqo[j])

        # Integrais das derivadas terceira da solução fundamental
        d3wdx3 =
            1 / (4 * pi * D) * cos(theta) * (3 - 2 * cos(theta)^2) * 2 * intNsr[j] / L * L /
            2
        d3wdx2dy =
            1 / (4 * pi * D) * sin(theta) * (1 - 2 * cos(theta)^2) * 2 * intNsr[j] / L * L /
            2
        d3wdxdy2 =
            1 / (4 * pi * D) * cos(theta) * (1 - 2 * sin(theta)^2) * 2 * intNsr[j] / L * L /
            2
        d3wdy3 =
            1 / (4 * pi * D) * sin(theta) * (3 - 2 * sin(theta)^2) * 2 * intNsr[j] / L * L /
            2

        # Integrais das derivadas quarta da solução fundamental
        d4wdx4 =
            1 / (4 * pi * D) *
            (8 * cos(theta)^4 - 12 * cos(theta)^2 + 3) *
            (4 * intNsr2[j] / L^2) *
            L / 2
        d4wdx3dy =
            1 / (4 * pi * D) *
            (
                -6 * cos(theta) * sin(theta) +
                4 * sin(theta) * cos(theta)^3 +
                4 * cos(theta)^3 * sin(theta)
            ) *
            (4 * intNsr2[j] / L^2) *
            L / 2
        d4wdx2dy2 =
            1 / (4 * pi * D) *
            (8 * (cos(theta) * sin(theta))^2 - 1) *
            (4 * intNsr2[j] / L^2) *
            L / 2
        d4wdxdy3 =
            1 / (4 * pi * D) *
            (
                -6 * cos(theta) * sin(theta) +
                4 * cos(theta) * sin(theta)^3 +
                4 * cos(theta) * sin(theta)^3
            ) *
            (4 * intNsr2[j] / L^2) *
            L / 2
        d4wdy4 =
            1 / (4 * pi * D) *
            (8 * sin(theta)^4 - 12 * sin(theta)^2 + 3) *
            (4 * intNsr2[j] / L^2) *
            L / 2

        # Cálculo de f1, f2 e f3 conforme Shi e Bezine
        f1 = D11 * nx^2 + 2 * D16 * nx * ny + D12 * ny^2
        f2 = 2 * (D16 * nx^2 + 2 * D66 * nx * ny + D26 * ny^2)
        f3 = D12 * nx^2 + 2 * D26 * nx * ny + D22 * ny^2

        # Cálculo de h1, h2, h3 e h4 conforme Shi e Bezine
        h1 = D11 * nx * (1 + ny^2) + 2 * D16 * ny^3 - D12 * nx * ny^2
        h2 =
            4 * D16 * nx + D12 * ny * (1 + nx^2) + 4 * D66 * ny^3 - D11 * nx^2 * ny -
            2 * D26 * nx * ny^2
        h3 =
            4 * D26 * ny + D12 * nx * (1 + ny^2) + 4 * D66 * nx^3 - D22 * nx * ny^2 -
            2 * D16 * nx^2 * ny
        h4 = D22 * ny * (1 + nx^2) + 2 * D26 * nx^3 - D12 * nx^2 * ny

        #--------------------------------------------------------------------------------
        # Cálculo das demais variáveis da solução fundamental
        #--------------------------------------------------------------------------------
        mn = -(f1 * d2wdx2 + f2 * d2wdxdy + f3 * d2wdy2)
        vn = -(h1 * d3wdx3 + h2 * d3wdx2dy + h3 * d3wdxdy2 + h4 * d3wdy3)

        # Cálculo das derivadas primeiras de Mn
        dmndx = -(f1 * d3wdx3 + f2 * d3wdx2dy + f3 * d3wdxdy2)
        dmndy = -(f1 * d3wdx2dy + f2 * d3wdxdy2 + f3 * d3wdy3)

        # Cálculo das derivadas primeiras de Vn
        dvndx = -(h1 * d4wdx4 + h2 * d4wdx3dy + h3 * d4wdx2dy2 + h4 * d4wdxdy3)
        dvndy = -(h1 * d4wdx3dy + h2 * d4wdx2dy2 + h3 * d4wdxdy3 + h4 * d4wdy4)
        m1 = nx
        m2 = ny
        d2wdndm = -(d2wdx2 * nx * m1 + d2wdxdy * (nx * m2 + ny * m1) + d2wdy2 * ny * m2)
        dmndm = -(dmndx * m1 + dmndy * m2)
        dvndm = -(dvndx * m1 + dvndy * m2)
        # Atribuição das variáveis da solução fundamental - atribuição de w_ij à variável
        # de saída w_est
        g22 = -d2wdndm
        h11 = vn
        h22 = -dmndm
        h12 = -mn
        h21 = dvndm

        # @infiltrate
        if xi0 ≈ -2 / 3
            h_el[:, 2*j-1:2*j] = [
                h11 h12
                h21 h22
            ]
            g_el[j] = g22
        elseif xi0 ≈ -1
            h_el[1, 2*j-1:2*j] = [h11 h12]
        elseif xi0 ≈ 1
            h_el[1, 7-2*j:8-2*j] = [h11 h12]
        else
            h_el[:, 7-2*j:8-2*j] = [
                h11 h12
                h21 h22
            ]
            g_el[4-j] = g22
        end
    end
    return h_el, g_el
end

function compute_q(pf, nf, x, eta, Gauss_w, elem, dad::placa_fina_isotropica, pre)

    # Inicialização da variável q_el
    q_el = zeros(2, 1)
    D = dad.k.D
    # Início da integração numérica sobre o elemento
    for k = 1:size(Gauss_w, 1)# Integração da carga distribuída sobre o domínio
        N, dN = calc_fforma(eta[k], elem)
        pc = N' * x # Ponto de gauss interpolador
        rs = pc' - pf
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        n1 = sy
        n2 = -sx
        xf = pf[1]
        yf = pf[2]
        # Distance of source and field points
        r1 = rs[1]
        r2 = rs[2]
        r = norm(rs)
        rd1 = r1 / r
        rd2 = r2 / r
        m1 = nf[1]
        m2 = nf[2]
        theta = atan(r2, r1)
        nr = n1 * rd1 + n2 * rd2

        a = 1
        A = dad.k.carga[1]
        B = dad.k.carga[2]
        C = dad.k.carga[3]
        Clocal = A * xf + B * yf + C
        # mr = m1 * rd1 + m2 * rd2
        nr = n1 * rd1 + n2 * rd2
        # mn = n1 * m1 + n2 * m2

        # (integral (C x fundamental solution x nr))/r
        # int1 = Clocal / (128 * pi * D) * r^3 * (4 * log(r) - 3) * nr
        # int2 = -Clocal / (128 * pi * D) * r^2 * ((4 * log(r) - 3) * (2 * mr * nr + mn) + 4 * mr * nr)



        intrdwdx = -(r^2 * cos(theta) * (2 - 6 * log(r))) / (72 * D * pi)

        intrdwdy = -(r^2 * (2 - 6 * log(r)) * sin(theta)) / (72 * D * pi)

        # int1 = -Clocal * nr * (r^3 * (1 * cos(theta)^2 - 4 * log(r))) / (128 * D * pi)
        # int2 = -Clocal * nr * (intrdwdx * m1 + intrdwdy * m2)

        int1 = -Clocal * nr * (r^3 * (3 - 4 * log(r))) / (128 * D * pi)
        int2 = -Clocal * nr * (intrdwdx * m1 + intrdwdy * m2)



        int_gw = [
            int1
            int2
        ]
        # Cálculo dos componentes de q_el [2 X 2]
        q_el = q_el + int_gw * dgamadqsi * Gauss_w[k]
        # @infiltrate
    end
    #--------------------------------------------------------------------------------
    return q_el
end

function compute_Rw(pf, nf, dad::placa_fina_isotropica, pre)
    m1 = nf[1]
    m2 = nf[2]
    D = dad.k.D
    ni = dad.k.nu
    RS = zeros(2, size(dad.k.cantos, 1))
    WS = zeros(2, size(dad.k.cantos, 1))
    for i = 1:size(dad.k.cantos, 1)
        pc = dad.k.cantos[i, 2:3]
        na1 = dad.k.cantos[i, 4]
        na2 = dad.k.cantos[i, 5]
        nd1 = dad.k.cantos[i, 6]
        nd2 = dad.k.cantos[i, 7]

        sa1 = -na2
        sa2 = na1
        sd1 = -nd2
        sd2 = nd1

        rs = pc - pf
        # Distance of source and field points
        r1 = rs[1]
        r2 = rs[2]
        r = norm(rs)
        if r != 0
            rd1 = r1 / r
            rd2 = r2 / r

            rsa = sa1 * rd1 + sa2 * rd2
            rsd = sd1 * rd1 + sd2 * rd2
            rna = na1 * rd1 + na2 * rd2
            rnd = nd1 * rd1 + nd2 * rd2
            mna = m1 * na1 + m2 * na2
            mnd = m1 * nd1 + m2 * nd2
            mr = m1 * rd1 + m2 * rd2
            msa = m1 * sa1 + m2 * sa2
            msd = m1 * sd1 + m2 * sd2

            sen2ba = -2 * rna * rsa
            sen2bd = -2 * rnd * rsd

            mnsa = (1 - ni) / (8 * pi) * sen2ba
            mnsd = (1 - ni) / (8 * pi) * sen2bd

            Rci = mnsd - mnsa

            dmnsa = -(1 - ni) / (4 * pi * r) * (2 * mr * rna * rsa - mna * rsa - msa * rna)
            dmnsd = -(1 - ni) / (4 * pi * r) * (2 * mr * rnd * rsd - mnd * rsd - msd * rnd)

            dRci = dmnsd - dmnsa

            we = 1 / (8 * pi * D) * r^2 * (log(r) - 0.5)
            dwdm = -r / (4 * pi * D) * mr * log(r)
            R = [Rci dRci]'
            w = [we dwdm]'
        else
            sa1 = -na2
            sa2 = na1
            sd1 = -nd2
            sd2 = nd1
            cbetac = -(sa1 * sd1 + sa2 * sd2)
            sbetac = -(na1 * sd1 + na2 * sd2)
            betac = atan(sbetac, cbetac) / (2 * pi)
            if betac < 0
                betac += 1
            end
            R = [betac 0]'
            w = [0 0]'
        end
        RS[:, i] = R
        WS[:, i] = w
    end
    return RS, WS
end


function autovalor_num(H, G, M, dad)

    ind_deslconhecido = zeros(Int, 0)
    nc = size(dad.NOS, 1)
    ncanto = size(dad.k.bc_canto, 1)
    nh = size(H, 2)
    for elem_i in dad.ELEM, i = 1:2  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        # @show elem_i.tipoCDC[i]
        if elem_i.tipoCDC[i] == 0
            ind_deslconhecido = [ind_deslconhecido; 2ind_elem .+ (i - 2)]
        end
    end
    for c = 1:ncanto
        if dad.k.bc_canto[c, 1] == 0
            ind_deslconhecido = [ind_deslconhecido; nh - ncanto + c]
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
    v = ones(size(Hb, 1))
    lold = 0
    l = 0
    for i = 1:100
        vnew = Hb \ (Mb * v)
        l = dot(v, vnew) / dot(vnew, vnew)
        if abs((l - lold) / l) < 1e-6
            continue
        end
        lold = l
        v = vnew / norm(vnew)
        # @show i
    end
    u[ind_forcaconhecida] = real.(v)
    # @infiltrate

    l, [u[1:2:2nc]; u[2nc+1:end-ncanto]]
end

function placa_grande1(dad::placa_fina_isotropica, dadpe, ntime, npg, e = 0.5)
    tHeG = @timed H, G, q, It, dNs = calc_HeGeIt(dad, npg)  #importante
    tHeGpe = @timed Hpe, Gpe = calc_HeG(dadpe, npg, interno = true)  #importante

    # Md = BEM.Monta_M_RIMd(dad, npg)
    dMd = BEM.Monta_dM_RIMd(dad, npg)
    Mpe = BEM.Monta_M_RIMd(dadpe, npg)
    A, B, bc_val = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b


    Ape, bpe = aplicaCDC(Hpe, Gpe, dadpe) # Calcula a matriz A e o vetor b

    CB = dadpe.k.E * dad.k.h / (1 - dadpe.k.nu^2)
    ~, Fx, Fy =
        BEM.montaFs([dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]], smooth = 1e-8)
    # F, ~, ~ = BEM.montaFs([dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]], [dad.NOS; dad.pontos_internos])

    tens = zeros(size(dad.NOS, 1) + size(dad.pontos_internos, 1) + size(dad.k.cantos, 1), 3)
    tens_nt = zeros(size(dad.NOS, 1), 3)
    ws = zeros(
        2 * size(dad.NOS, 1) + size(dad.pontos_internos, 1) + size(dad.k.cantos, 1),
        ntime,
    )
    wind = [
        1:2:2*nc(dad)
        2*nc(dad)+1:2*size(dad.NOS, 1)+size(dad.pontos_internos, 1)+size(dad.k.cantos, 1)
    ]
    @showprogress "Passo não-linear" for i = 1:ntime
        Mit = BEM.aplicaT(dad, G, tens_nt, dNs, dMd[1], tens, Fx, Fy)
        # Mit = BEM.aplicaT(G, It, tens_nt, dNs, dMd[6], tens)
        # Mit = BEM.aplicaT(dad, Md, tens[1:end-4, :])

        # Mit=(Mit1+Mit2)/2
        # @infiltrate

        # b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
        if i == 1
            b = B * bc_val + i * q / ntime
        elseif i == 2
            b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
        elseif i == 3
            b = B * bc_val + i * q / ntime + Mit * (e * ws[:, i-1] + (1 - e) * ws[:, i-2])
        else
            # b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
            b = B * bc_val + i * q / ntime + Mit * (e * ws[:, i-1] + (1 - e) * ws[:, i-2])
            # b = B * bc_val + i * q / ntime + Mit * (e * ws[:, i-1] + (1 - e) * (ws[:, i-1] + ws[:, i-2]) / 2)
        end
        x = A \ b
        listdisp, listtra = BEM.separa(x, dad)

        ws[:, i] = listdisp

        if i == 1
            w = ws[wind, i]
        elseif i == 2
            w = (e * ws[wind, i] + (1 - e) * ws[wind, i-1])
        else
            # w = (e * ws[wind, i] + (1 - e) * ws[wind, i-1])
            w = (e * ws[wind, i] + (1 - e) * (ws[wind, i-1] + ws[wind, i-2]) / 2)
        end
        # w = ws[wind, i]
        dwdx = Fx * w
        dwdy = Fy * w

        somadw = (dadpe.k.nu) / (1 - dadpe.k.nu) * (dwdx .^ 2 + dwdy .^ 2)
        N_nonlin =
            (1 - dadpe.k.nu) / 2 * CB * [dwdx .^ 2 + somadw dwdy .^ 2 + somadw dwdx .* dwdy]
        dNndx = Fx * N_nonlin
        dNndy = Fy * N_nonlin
        DN = [dNndx[:, 1] + dNndy[:, 3] dNndx[:, 3] + dNndy[:, 2]]'[1:2*(size(
            dad.NOS,
            1,
        )+size(dad.pontos_internos, 1))]
        qnl = Mpe * DN

        x = Ape \ (bpe + qnl)
        u, t, uint = separa(dadpe, x) #importante

        ut = [u; uint; zeros(4, 2)]
        dudx = Fx * ut
        dudy = Fy * ut
        somadu = 2(dadpe.k.nu) / (1 - dadpe.k.nu) * (dudx[:, 1] + dudy[:, 2])
        N_lin =
            (1 - dadpe.k.nu) / 2 *
            CB *
            [2dudx[:, 1] + somadu 2dudy[:, 2] + somadu dudy[:, 1] + dudx[:, 2]]
        # N_lin = dad.k.A3 * [duxdx duydy (duxdy + duydx)]
        # tens = F * (N_lin + N_nonlin)
        tens = (N_lin + N_nonlin)
        tens_nt = criatensoesnt(dad, tens)
        # @infiltrate


    end
    ws
end
function matrizes_grande(dad::Union{placa_fina,placa_fina_isotropica}, dadpe, npg)
    H, G, q, It, dNs = calc_HeGeIt(dad, npg)  #importante
    Hpe, Gpe = calc_HeG(dadpe, npg, interno = true)  #importante

    # Md = BEM.Monta_M_RIMd(dad, npg)
    dMd = BEM.Monta_dM_RIMd(dad, npg)
    Mpe = BEM.Monta_M_RIMd(dadpe, npg)
    A, B, bc_val = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
    H, G, q, It, dNs, Hpe, Gpe, dMd, Mpe, A, B, bc_val
end

function placa_grande(mats, ntime, dad::placa_fina, dadpe, e = 0.5)
    H, G, q, It, dNs, Hpe, Gpe, dMd, Mpe, A, B, bc_val = mats

    Ape, bpe = aplicaCDC(Hpe, Gpe, dadpe) # Calcula a matriz A e o vetor b

    ~, Fx, Fy =
        BEM.montaFs([dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]], smooth = 1e-8)
    # F, ~, ~ = BEM.montaFs([dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]], [dad.NOS; dad.pontos_internos])

    tens = zeros(size(dad.NOS, 1) + size(dad.pontos_internos, 1) + size(dad.k.cantos, 1), 3)
    tens_nt = zeros(size(dad.NOS, 1), 3)
    ws = zeros(
        2 * size(dad.NOS, 1) + size(dad.pontos_internos, 1) + size(dad.k.cantos, 1),
        ntime,
    )
    wind = [
        1:2:2*nc(dad)
        2*nc(dad)+1:2*size(dad.NOS, 1)+size(dad.pontos_internos, 1)+size(dad.k.cantos, 1)
    ]
    no_meio = ceil(Int, size(dad.NOS, 1) + size(dad.pontos_internos, 1) / 2)
    # tenso = zeros(ntime, 3)
    @showprogress "Passo não-linear" for i = 1:ntime
        Mit = BEM.aplicaT(dad, G, tens_nt, dNs, dMd[1], tens, Fx, Fy) * 1
        # Mit = BEM.aplicaT(G, It, tens_nt, dNs, dMd[6], tens)
        # Mit = BEM.aplicaT(dad, Md, tens[1:end-4, :])

        # Mit=(Mit1+Mit2)/2
        # @infiltrate
        # b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
        if i == 1
            b = B * bc_val + i * q / ntime
        elseif i == 2
            b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
        elseif i == 3
            b = B * bc_val + i * q / ntime + Mit * (e * ws[:, i-1] + (1 - e) * ws[:, i-2])
        else
            # b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
            b = B * bc_val + i * q / ntime + Mit * (e * ws[:, i-1] + (1 - e) * ws[:, i-2])

        end
        x = A \ b
        listdisp, listtra = BEM.separa(x, dad)
        # if norm(ws[no_meio, i] - listdisp[no_meio]) < 1e-2
        #   continue
        # end
        ws[:, i] = listdisp


        if i == 1
            w = ws[wind, i]
        elseif i == 2
            w = (e * ws[wind, i] + (1 - e) * ws[wind, i-1])
        else
            # w = (e * ws[wind, i] + (1 - e) * ws[wind, i-1])
            w = (e * ws[wind, i] + (1 - e) * (ws[wind, i-1] + ws[wind, i-2]) / 2)
        end
        # w = ws[wind, i]
        dwdx = Fx * w
        dwdy = Fy * w
        [dwdx dwdy]
        # @infiltrate
        # corrige_derivada(listdisp, dad, dwdx, dwdy)

        N_nonlin = (dadpe.k.A3 * [dwdx .^ 2 dwdy .^ 2 2dwdx .* dwdy]')' / 2
        dNndx = Fx * N_nonlin
        dNndy = Fy * N_nonlin
        DN = [dNndx[:, 1] + dNndy[:, 3] dNndx[:, 3] + dNndy[:, 2]]'[1:2*(size(
            dad.NOS,
            1,
        )+size(dad.pontos_internos, 1))]
        qnl = Mpe * DN

        x = Ape \ (bpe + qnl)
        u, t, uint = separa(dadpe, x) #importante
        ut = [u; uint; zeros(4, 2)]
        dudx = Fx * ut
        dudy = Fy * ut
        # somadu = 2(dadpe.k.nu) / (1 - dadpe.k.nu) * (dudx[:, 1] + dudy[:, 2])
        # N_lin = (1 - dadpe.k.nu) / 2 * CB * [2dudx[:, 1] + somadu 2dudy[:, 2] + somadu dudy[:, 1] + dudx[:, 2]]
        N_lin = (dadpe.k.A3 * [dudx[:, 1] dudy[:, 2] (dudy[:, 1] + dudx[:, 2])]')'
        # tens = F * (N_lin + N_nonlin)


        tens = (N_lin + N_nonlin)
        tens_nt = criatensoesnt(dad, tens)
        # tenso[i, :] = tens[no_meio, :]
        # @infiltrate
        # tens[no_meio, :]


    end
    ws
end
function corrige_derivada(w, dad, dwdx, dwdy)
    for elem_j in dad.ELEM  #Laço dos elementos

        x = dad.NOS[elem_j.indices, :]

        for i = 1:size(elem_j.indices, 1)
            N, dN = calc_fforma(elem_j.ξs[i], elem_j)
            dxdqsi = dN' * x
            dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
            sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
            sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

            normal = [sy -sx]

            dw = w[2elem_j.indices.-1]' * dN / dgamadqsi
            # @infiltrate
            dudnt = [-normal[2] normal[1]] * dw + normal * w[2elem_j.indices[i]]
            dwdx[elem_j.indices[i]] = dudnt[1]
            dwdy[elem_j.indices[i]] = dudnt[2]
        end
    end
end
function placa_grande(mats, ntime, dad::placa_fina_isotropica, dadpe, e = 0.5)
    #https://docs.sciml.ai/NonlinearSolve/stable/tutorials/nonlinear/#Using-Jacobian-Free-Newton-Krylov-(JNFK)-Methods
    H, G, q, It, dNs, Hpe, Gpe, dMd, Mpe, A, B, bc_val = mats

    Ape, bpe = aplicaCDC(Hpe, Gpe, dadpe) # Calcula a matriz A e o vetor b

    CB = dadpe.k.E * dad.k.h / (1 - dadpe.k.nu^2)

    ~, Fx, Fy =
        BEM.montaFs([dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]], smooth = 1e-8)
    # F, ~, ~ = BEM.montaFs([dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]], [dad.NOS; dad.pontos_internos])

    tens = zeros(size(dad.NOS, 1) + size(dad.pontos_internos, 1) + size(dad.k.cantos, 1), 3)
    tens_nt = zeros(size(dad.NOS, 1), 3)
    ws = zeros(
        2 * size(dad.NOS, 1) + size(dad.pontos_internos, 1) + size(dad.k.cantos, 1),
        ntime,
    )
    wind = [
        1:2:2*nc(dad)
        2*nc(dad)+1:2*size(dad.NOS, 1)+size(dad.pontos_internos, 1)+size(dad.k.cantos, 1)
    ]
    no_meio = ceil(Int, size(dad.NOS, 1) + size(dad.pontos_internos, 1) / 2)
    # tenso = zeros(ntime, 3)
    @showprogress "Passo não-linear" for i = 1:ntime
        Mit = BEM.aplicaT(dad, G, tens_nt, dNs, dMd[1], tens, Fx, Fy)
        # Mit = BEM.aplicaT(G, It, tens_nt, dNs, dMd[6], tens)
        # Mit = BEM.aplicaT(dad, Md, tens[1:end-4, :])

        # Mit=(Mit1+Mit2)/2
        # @infiltrate

        # b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
        if i == 1
            b = B * bc_val + i * q / ntime
        elseif i == 2
            b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
        elseif i == 3
            b = B * bc_val + i * q / ntime + Mit * (e * ws[:, i-1] + (1 - e) * ws[:, i-2])
        else
            # b = B * bc_val + i * q / ntime + Mit * ws[:, i-1]
            b = B * bc_val + i * q / ntime + Mit * (e * ws[:, i-1] + (1 - e) * ws[:, i-2])
            # b = B * bc_val + i * q / ntime + Mit * (e * ws[:, i-1] + (1 - e) * (ws[:, i-1] + ws[:, i-2]) / 2)
        end
        x = A \ b
        listdisp, listtra = BEM.separa(x, dad)

        ws[:, i] = listdisp

        if i == 1
            w = ws[wind, i]
        elseif i == 2
            w = (e * ws[wind, i] + (1 - e) * ws[wind, i-1])
        else
            # w = (e * ws[wind, i] + (1 - e) * ws[wind, i-1])
            w = (e * ws[wind, i] + (1 - e) * (ws[wind, i-1] + ws[wind, i-2]) / 2)
        end
        # w = ws[wind, i]
        dwdx = Fx * w
        dwdy = Fy * w

        somadw = (dadpe.k.nu) / (1 - dadpe.k.nu) * (dwdx .^ 2 + dwdy .^ 2)
        N_nonlin =
            (1 - dadpe.k.nu) / 2 * CB * [dwdx .^ 2 + somadw dwdy .^ 2 + somadw dwdx .* dwdy]
        dNndx = Fx * N_nonlin
        dNndy = Fy * N_nonlin
        DN = [dNndx[:, 1] + dNndy[:, 3] dNndx[:, 3] + dNndy[:, 2]]'[1:2*(size(
            dad.NOS,
            1,
        )+size(dad.pontos_internos, 1))]
        qnl = Mpe * DN

        x = Ape \ (bpe + qnl)
        u, t, uint = separa(dadpe, x) #importante

        ut = [u; uint; zeros(4, 2)]
        dudx = Fx * ut
        dudy = Fy * ut
        somadu = 2(dadpe.k.nu) / (1 - dadpe.k.nu) * (dudx[:, 1] + dudy[:, 2])
        N_lin =
            (1 - dadpe.k.nu) / 2 *
            CB *
            [2dudx[:, 1] + somadu 2dudy[:, 2] + somadu dudy[:, 1] + dudx[:, 2]]
        # N_lin = dad.k.A3 * [duxdx duydy (duxdy + duydx)]
        # tens = F * (N_lin + N_nonlin)
        tens = (N_lin + N_nonlin)
        tens_nt = criatensoesnt(dad, tens)
        # @infiltrate i == 300


    end
    ws
end
