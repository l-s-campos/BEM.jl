function termoelasticidade(dad, npg; θ=1, carga=0, metodo="dibem")
    H, G = calc_HeG(dad, npg, interno=true)  #importante

    A, b = BEM.aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b

    EE = dad.k.E
    v = dad.k.nu
    k = dad.k.k

    kchap = EE * k / (1 - 2 * v)

    np = nc(dad) + ni(dad)
    if θ isa Number
        theta = fill(θ, np)
    elseif θ isa Vector
        theta = θ
    elseif θ isa Function
        theta = [θ.(dad.NOS[:, 1], dad.NOS[:, 2]); θ.(dad.pontos_internos[:, 1], dad.pontos_internos[:, 2])]
    end

    F, Fx, Fy = BEM.montaFs([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])
    dF = [zeros(size(Fx)); zeros(size(Fy))]
    for i = 1:np
        dF[2i-1, :] = Fx[i, :]
        dF[2i, :] = Fy[i, :]
    end
    if θ isa Number && carga == 0
        q2 = 0
        dq2 = 0
        qc, dqc = 0, 0

    else
        if metodo == "DIBEM"
            Mpe = BEM.Monta_M_RIMd(dad, npg)
            dMpe = BEM.Monta_dM_RIMd(dad, npg)
            # Mpe = zeros(2np, 2np)
            # dMpe = zeros(3np, 2np)
        else
            Mpe = BEM.Monta_M_RIM(dad, npg)
            dMpe = BEM.Monta_dM_RIM(dad, npg)
        end
        q2 = Mpe * kchap * dF * theta
        dq2 = dMpe * kchap * dF * theta
        # @infiltrate
        # q2 = Mpe * kchap * dtemp'[:]
        # dq2 = dMpe * kchap * dtemp'[:]
    end
    if carga == 0
        qc = 0
        dqc = 0
    else
        carga_nodal = zeros(2np)
        for i = 1:nc(dad)
            carga_nodal[2i-1:2i] = carga(dad.NOS[i, 1], dad.NOS[i, 2])
        end
        for i = nc(dad)+1:np
            carga_nodal[2i-1:2i] = carga(dad.pontos_internos[i-nc(dad), 1], dad.pontos_internos[i-nc(dad), 2])
        end
        qc = Mpe * carga_nodal
        dqc = dMpe * carga_nodal
    end
    # qc1, dqc1 = calc_forçacorpo(dad, dta, 10)
    # qc, dqc = calc_forçacorpo(dad, carga, 10)

    # @infiltrate
    normal_fonte = BEM.calc_normais(dad)

    q1 = G * (normal_fonte .* theta[1:nc(dad)])'[:] * kchap

    # x = A \ (b + q1)
    x = A \ (b + q1 .- q2 .+ qc)

    u, t, uint = separa(dad, x) #importante


    S, D = calc_SeD(dad)
    dq1 = D * (normal_fonte .* theta[1:nc(dad)])'[:] * kchap
    dq3 = [kchap * theta kchap * theta 0 * theta]
    # σi = reshape(D * t'[:] - S * u'[:] + dq1  - dq3, 3, :)'
    fatorcontorno = [fill(2, nc(dad)); ones(ni(dad))]
    sigma = reshape(D * t'[:] - S * u'[:] + dq1 .- dq2 .+ dqc, 3, :)' .* fatorcontorno - dq3
    #  rbf
    # du = dF * [u; uint]
    # sigma1 = zeros(np, 3)
    # for i = 1:np
    #     # @infiltrate
    #     sigma1[i, 1] = 2GG * v / (1 - 2v) * (du[2i-1, 1] + du[2i, 2]) + 2GG * (du[2i-1, 1]) - kchap * theta[i]
    #     sigma1[i, 2] = 2GG * v / (1 - 2v) * (du[2i-1, 1] + du[2i, 2]) + 2GG * (du[2i, 2]) - kchap * theta[i]
    #     sigma1[i, 3] = GG * (du[2i-1, 2] + du[2i, 1])
    # end
    # @infiltrate
    #
    [u; uint], sigma
end

