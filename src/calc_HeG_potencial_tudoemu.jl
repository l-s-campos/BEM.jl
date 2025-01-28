function calc_HeGt(dad::potencial, npg = 8)

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
                eta, Jt = sinhtrans(qsi, eet, b)
                # @show eet,b
                # eta, wt = pontosintegra(dad.NOS, elem_j, ind_elem, qsi, w)
                # @show norm(eta-eta1)

                # h, g = integraelem(pf, x, eta, wt, elem_j, dad)
                h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
                # @infiltrate contafonte == 2
                H[contafonte, elem_j.indices] = h
                G[contafonte, elem_j.indices] = g

            end
            contafonte += 1
        end
    end
    for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
        H[i, i] = -0.5
    end

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

    Ht = [H zeros(n, ni); Hi -I]
    Gt = [G; Gi]
    # @infiltrate
    Ht, Gt
end
function aplicaCDC_u(Ht, Gt, dad)
    n = size(dad.NOS, 1)
    tipoCDC = zeros(Bool, n)
    valoresconhecidos = zeros(n)
    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[1] == 1
            tipoCDC[ind_elem] .= 1
            valoresconhecidos[ind_elem] = elem_i.valorCDC
        end
    end


    normal_fonte = calc_normais(dad)
    troca = tipoCDC[:] .== 0
    F, Fx, Fy = BEM.montaFs(dad.NOS[troca, :], [dad.NOS; dad.pontos_internos])

    b = zeros(n)
    Gu =
        -dad.k *
        Gt[:, troca] *
        (Fx .* normal_fonte[troca, 1] + Fy .* normal_fonte[troca, 2])
    # @infiltrate
    A = Ht - Gu
    b = Gt[:, (1:n)[tipoCDC]] * valoresconhecidos[tipoCDC]
    # b = (-H[:, troca] + Gu[:, troca]) * valoresconhecidos[troca] + Gt[:, tipoCDC] * valoresconhecidos[tipoCDC]

    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        # @infiltrate
        if elem_i.tipoCDC == 0
            for i in ind_elem
                A[i, :] .= 0
                A[i, i] = 1
                b[i] = valoresconhecidos[i]
            end
        end
    end
    A, b
end

function aplicaCDC_alpha(Ht, Gt, dad)
    normal_fonte = calc_normais(dad)
    F, dFdx, dFdy, dFdxx, dFdyy, dFdxy =
        BEM.montaF([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])
    dFdn = (
        dFdx[1:nc(dad), :] .* normal_fonte[:, 1] + dFdy[1:nc(dad), :] .* normal_fonte[:, 2]
    )
    Ga = -dad.k * Gt * dFdn
    Ha = Ht * F

    n = size(dad.NOS, 1)
    tipoCDC = zeros(Bool, n)
    valoresconhecidos = zeros(n)
    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[1] == 1
            tipoCDC[ind_elem] .= 1
            valoresconhecidos[ind_elem] = elem_i.valorCDC
        end
    end
    tipof = [tipoCDC; zeros(Bool, ni(dad))]
    A = [Ga + Ha; [F[tipof, :]; -dad.k * dFdn[tipoCDC.==0, :]]]
    b = zeros(2n + ni(dad))
    b[ni(dad)+n+1:end, :] = [valoresconhecidos[tipoCDC]; valoresconhecidos[tipoCDC.==0]]
    # @infiltrate

    a = A \ b
    T = F * a
    q = dFdn * a
    T, q
end


function corrige_autovalor(G, H, M, dad)
    n = nc(dad) + ni(dad)
    ind_ud = trues(n)
    ind_qd = falses(nc(dad))
    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[1] == 0
            ind_qd[ind_elem] .= 1
            ind_ud[ind_elem] .= 0
        end
    end
    guu = G[(1:nc(dad))[ind_qd], ind_qd]
    gqu = G[ind_ud, ind_qd]
    huq = H[(1:nc(dad))[ind_qd], ind_ud]
    hqq = H[ind_ud, ind_ud]
    muq = M[(1:nc(dad))[ind_qd], ind_ud]
    mqq = M[ind_ud, ind_ud]
    Hb = hqq - gqu * (guu \ huq)
    Mb = mqq - gqu * (guu \ muq)
    a, v = eigen(Hb, -Mb)
    v_wn = sort(sqrt.(abs.(real(a))))
end


function corrige_autovalor_u(Gt, Ht, M, dad)
    n = size(dad.NOS, 1)

    ud = trues(n + ni(dad))
    qd = falses(nc(dad))

    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[1] == 0
            ud[ind_elem] .= 0
            qd[ind_elem] .= 1
        end
    end

    normal_fonte = calc_normais(dad)
    F, Fx, Fy = BEM.montaFs(dad.NOS[qd, :], [dad.NOS; dad.pontos_internos])

    Gu = -dad.k * Gt[:, qd] * (Fx .* normal_fonte[qd, 1] + Fy .* normal_fonte[qd, 2])
    # @infiltrate
    Hu = Ht - Gu
    a, v = eigen(Hu[ud, ud], -M[ud, ud])
    v_wn = sort(sqrt.(abs.(real(a))))
end
