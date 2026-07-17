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


    normal_fonte = dad.normal
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
    normal_fonte = dad.normal
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
    A = [-Ga + Ha; [F[tipof, :]; -dad.k * dFdn[tipoCDC .== 0, :]]]
    b = zeros(2n + ni(dad))
    b[(ni(dad)+n+1):end, :] = [valoresconhecidos[tipoCDC]; valoresconhecidos[tipoCDC .== 0]]
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

    normal_fonte = dad.normal
    F, Fx, Fy = BEM.montaFs(dad.NOS[qd, :], [dad.NOS; dad.pontos_internos])

    Gu = -dad.k * Gt[:, qd] * (Fx .* normal_fonte[qd, 1] + Fy .* normal_fonte[qd, 2])
    # @infiltrate
    Hu = Ht - Gu
    a, v = eigen(Hu[ud, ud], -M[ud, ud])
    v_wn = sort(sqrt.(abs.(real(a))))
end
