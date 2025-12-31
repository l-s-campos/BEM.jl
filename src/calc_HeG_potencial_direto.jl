function calc_HeGd(dad::potencial, npg = 3, npg2 = 10; interno = false)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    ni = size(dad.pontos_internos, 1)
    intelems = zeros(n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    qsi2, w2 = gausslegendre(npg2)    # Quadratura de gauss
    normal_fonte = dad.normal
    k = dad.k

    ptos_contorno = dad.NOS
    ptos_interno = dad.pontos_internos
    ptos = [ptos_contorno; ptos_interno]

    for elem_j in dad.ELEM  #Laço dos elementos
        x = ptos_contorno[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelem = integraelemd(x, qsi, w, elem_j)
        # @infiltrate
        intelems[elem_j.indices] = intelem
    end

    if interno
        nt = n + ni
        H = zeros(nt, n)
        G = zeros(nt, n)
    else
        nt = n
        H = zeros(n, n)
        G = zeros(n, n)
    end

    for i = 1:nt
        pf = ptos[i, :]
        for elem_j in dad.ELEM
            xj = ptos_contorno[elem_j.indices[1], 1]
            yj = ptos_contorno[elem_j.indices[1], 2]

            r = sqrt((pf[1] - xj)^2 + (pf[2] - yj)^2)
            # @infiltrate i == 380
            if r <= 2 * elem_j.comprimento
                x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
                Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
                eet =
                    (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                    norm(Δelem)^2 + elem_j.ξs[1]
                N_geo, ~ = calc_fforma(eet, elem_j)
                ps = N_geo' * x
                # b = 1e-6
                b = norm(ps' - pf) / elem_j.comprimento

                eta, Jt = sinhtrans(qsi2, eet, b)

                h, g = integraelem(pf, x, eta, w2 .* Jt, elem_j, dad, k)

                H[i, elem_j.indices] = h
                G[i, elem_j.indices] = g
            else
                for j in elem_j.indices
                    if i == j
                        continue
                    end
                    xj = ptos_contorno[j, 1]
                    yj = ptos_contorno[j, 2]
                    # r = sqrt((xi - xj)^2 + (yi - yj)^2)

                    Qast, Tast =
                        calsolfund([xj - pf[1], yj - pf[2]], normal_fonte[j, :], dad, k)

                    H[i, j] = Qast * intelems[j]
                    G[i, j] = Tast * intelems[j]
                end
            end
        end
    end

    if interno
        H = [H [zeros(n, ni); -I]]
    end
    for i = 1:n                              #i=1:size(ptos_contorno,1) #Laço dos pontos fontes
        H[i, i] = 0.0
        # H[i, i] -= 0.5
        H[i, i] = -sum(H[i, :])
    end

    H, G
end

function calc_HeGd_simetria(dad::potencial, npg = 3; interno = false)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    intelems = zeros(n)
    H = zeros(n, n)
    G = zeros(n, n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    normal_fonte = dad.normal
    k = dad.k

    ptos_contorno = dad.NOS
    ptos_interno = dad.pontos_internos

    for elem_j in dad.ELEM  #Laço dos elementos
        x = ptos_contorno[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelem = integraelemd(x, qsi, w, elem_j)
        # @infiltrate
        intelems[elem_j.indices] = intelem
    end

    for i = 1:n
        xi = ptos_contorno[i, 1]
        yi = ptos_contorno[i, 2]
        for j = i:n
            if i == j
                continue
            end
            xj = ptos_contorno[j, 1]
            yj = ptos_contorno[j, 2]
            # r = sqrt((xi - xj)^2 + (yi - yj)^2)

            Qast, Tast = calsolfund([xj - xi, yj - yi], normal_fonte[j, :], dad, k)

            H[i, j] = Qast * intelems[j]
            G[i, j] = Tast * intelems[j]

            H[j, i] = Qast * intelems[i]
            G[j, i] = Tast * intelems[i]
        end
    end

    if interno
        ni = size(ptos_interno, 1)
        Hi = zeros(ni, n)
        Gi = zeros(ni, n)
        for i = 1:ni
            pf = ptos_interno[i, :]   # Coordenada (x,y)  dos pontos fonte
            for j = 1:n
                pc = ptos_contorno[j, :]

                Qast, Tast = calsolfund(pc - pf, normal_fonte[j, :], dad, k)
                # @infiltrate
                Hi[i, j] = Qast * intelems[j]
                Gi[i, j] = Tast * intelems[j]
            end
        end
        H = [H zeros(n, ni); Hi -I]
        G = [G; Gi]
    end

    for i = 1:n                              #i=1:size(ptos_contorno,1) #Laço dos pontos fontes
        # H[i, i] = -0.5
        H[i, i] = -sum(H[i, :])
    end

    if interno
        somaH =
            H * [
                ptos_contorno[:, 1] + ptos_contorno[:, 2]
                ptos_interno[:, 1] + ptos_interno[:, 2]
            ]
        somaG = G * (normal_fonte[:, 1] + normal_fonte[:, 2]) * k
    else
        somaH = H * (ptos_contorno[:, 1] + ptos_contorno[:, 2])
        somaG = G * (normal_fonte[:, 1] + normal_fonte[:, 2]) * k
    end

    for i = 1:n                              #i=1:size(ptos_contorno,1) #Laço dos pontos fontes
        G[i, i] = (-somaH[i] - somaG[i]) / (normal_fonte[i, 1] + normal_fonte[i, 2])
    end
    H, G
end

function calc_Hd(dad::potencial, npg = 3; interno = false)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    intelems = zeros(n)
    H = zeros(n, n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    normal_fonte = dad.normal
    k = dad.k

    ptos_contorno = dad.NOS
    ptos_interno = dad.pontos_internos

    for elem_j in dad.ELEM  #Laço dos elementos
        x = ptos_contorno[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelem = integraelemd(x, qsi, w, elem_j)
        # @infiltrate
        intelems[elem_j.indices] = intelem
    end

    for i = 1:n
        xi = ptos_contorno[i, 1]
        yi = ptos_contorno[i, 2]
        for j = 1:n
            if i == j
                continue
            end
            xj = ptos_contorno[j, 1]
            yj = ptos_contorno[j, 2]
            # r = sqrt((xi - xj)^2 + (yi - yj)^2)

            Qast, Tast = calsolfund([xj - xi, yj - yi], normal_fonte[j, :], dad, k)

            H[i, j] = Qast * intelems[j]
        end
    end

    if interno
        ni = size(ptos_interno, 1)
        Hi = zeros(ni, n)
        for i = 1:ni
            pf = ptos_interno[i, :]   # Coordenada (x,y)  dos pontos fonte
            for j = 1:n
                pc = ptos_contorno[j, :]

                Qast, Tast = calsolfund(pc - pf, normal_fonte[j, :], dad, k)
                # @infiltrate
                Hi[i, j] = Qast * intelems[j]
            end
        end
        H = [H zeros(n, ni); Hi -I]
    end


    for i = 1:n                              #i=1:size(ptos_contorno,1) #Laço dos pontos fontes
        # H[i, i] = -0.5
        H[i, i] = -sum(H[i, :])
    end

    if interno
        somaH =
            H * [
                ptos_contorno[:, 1] + ptos_contorno[:, 2]
                ptos_interno[:, 1] + ptos_interno[:, 2]
            ]
    else
        somaH = H * (ptos_contorno[:, 1] + ptos_contorno[:, 2])
    end

    H
end


function Monta_M_RIMdd(dad::potencial, npg; tiporadial = "tps", aug = false)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi, w = gausslegendre(npg)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(n_pontos)
    M1 = zeros(n_pontos)
    F, D = FeD(dad, nodes, tiporadial)

    nodes = [dad.NOS; dad.pontos_internos]
    normal = dad.normal
    k = dad.k

    intelems = zeros(n_nos)
    for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelem = integraelemd(x, qsi, w, elem_j)
        # @infiltrate
        intelems[elem_j.indices] = intelem
    end

    for i = 1:n_pontos #varia no ponto campo
        xi = nodes[i, 1]
        yi = nodes[i, 2]
        for j = 1:n_nos #varia nos pontos da integração de contorno

            xj = nodes[j, 1]
            yj = nodes[j, 2]

            r = [xj - xi, yj - yi]
            R = norm(r)

            if R == 0
                continue
            else
                m = int_interpolaρdρ(R, tipo = tiporadial)
                m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)

                M[i] = M[i] + dot(normal[j, :], r) / R^2 * m * intelems[j]
                M1[i] = M1[i] + dot(normal[j, :], r) / R^2 * m1 * intelems[j]
            end
        end
    end

    if aug
        P, Pint = calcPs(dad, npg)
        W = F \ P / (P' / F * P)
        Mf = M' / F
        A = (Mf - Mf * P * W' + Pint' * W') .* D
        for i = 1:n_pontos #Laço dos pontos radiais
            A[i, i] = 0
            A[i, i] = -sum(A[i, :])
        end
        return A + diagm(0 => M1)

    else
        #S = M' / F
        A = M' / F .* D
        for i = 1:n_pontos #Laço dos pontos radiais
            A[i, i] = 0
            A[i, i] = -sum(A[i, :])
        end
        return A + diagm(0 => M1)
    end

end



# function calc_HeG_Hd(dad::potencial; npg = 3, atol = 1e-4)
#     n = size(dad.NOS, 1)
#     intelems = zeros(n)
#     qsi, w = BEM.gausslegendre(npg)    # Quadratura de gauss
#     normal_fonte = BEM.calc_normais(dad)

#     for elem_j in dad.ELEM  #Laço dos elementos
#         x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
#         intelems[elem_j.indices] =  integraelemd(x, qsi, w, elem_j)
#         # @infiltrate
#         # intelem
#     end

#     # X = [[Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i in 1:nc(dad)]; [Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for i in 1:ni(dad)]]
#     # # X = [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i in 1:nc(dad)]

#     splitter = BEM.PrincipalComponentSplitter(; nmax = 10)
#     # Xclt = Yclt = ClusterTree(X, splitter)
#     Xclt, Yclt = ClusterTree(dad, splitter)

#     # @infiltrate
#     adm = StrongAdmissibilityStd(eta = 3)
#     comp = PartialACA(; atol)

#     KG = BEM.kernelG(dad, intelems)
#     KH = BEM.kernelH(dad, intelems)
#     HG = assemble_hmatrix(KG, Xclt, Yclt; adm, comp)
#     HH = assemble_hmatrix(KH, Xclt, Yclt; adm, comp)
#     # @infiltrate
#     corrige_diagonais!(dad, HH, HG, normal_fonte)
#     HH, HG
# end


function calc_HeG_Hd(dad::potencial; npg = 3, atol = 1e-6, nmax = 10, eta = 3)
    n = size(dad.NOS, 1)
    intelems = zeros(n)
    qsi, w = BEM.gausslegendre(npg)    # Quadratura de gauss
    normal_fonte = dad.normal

    for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelems[elem_j.indices] = integraelemd(x, qsi, w, elem_j)
    end

    splitter = BEM.PrincipalComponentSplitter(; nmax = nmax)
    Xclt, Yclt = ClusterTree(dad, splitter)

    V = [SVector{2,Float64}(row) for row in eachrow(dad.NOS)]
    Yclt_cont = ClusterTree(V, splitter)

    # @infiltrate
    adm = StrongAdmissibilityStd(eta = eta)
    comp = PartialACA(; atol)

    pontos = [SVector{2,Float64}(row) for row in eachrow([dad.NOS; dad.pontos_internos])]

    normal = [SVector{2,Float64}(row) for row in eachrow(normal_fonte)]

    k = dad.k

    KG = BEM.kernelG(pontos, n, k, intelems)
    KH = BEM.kernelH(pontos, n, normal, intelems)

    HG = assemble_hmatrix(KG, Xclt, Yclt_cont; adm, comp)
    HH = assemble_hmatrix(KH, Xclt, Yclt_cont; adm, comp)
    # @infiltrate
    corrige_diagonais!(dad, HH, HG, normal_fonte)
    HH, HG
end

#--------
function preencheH(
    MH::HMatrix{ClusterTree{2,Float64},Float64},
    dad::potencial;
    npg = 3,
    atol = 1e-6,
    rtol = 1e-5,
)
    normal_fonte = [Point2D(dad.normal[i, 1], dad.normal[i, 2]) for i = 1:nc(dad)]

    intelems = zeros(nc(dad))
    qsi, w = BEM.gausslegendre(npg)    # Quadratura de gauss

    for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelems[elem_j.indices] = BEM.integraelemd(x, qsi, w, elem_j)
    end

    X = [
        [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
        [Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for i = 1:ni(dad)]
    ]

    KG = BEM.kernelG(X, nc(dad), dad.k, intelems)
    KH = BEM.kernelH(X, nc(dad), normal_fonte, intelems)

    KG = BEM.PermutedMatrix(KG, BEM.loc2glob(MH.rowtree), BEM.loc2glob(MH.coltree))
    KH = BEM.PermutedMatrix(KH, BEM.loc2glob(MH.rowtree), BEM.loc2glob(MH.coltree))


    comp = PartialACA(; atol = atol, rtol = rtol)

    HH = deepcopy(MH)
    rmul!(HH, 0)
    GH = deepcopy(HH)

    BEM.assemble_contorno!(HH, KH, comp, BEM.ACABuffer(Float64), nc(dad))
    BEM.assemble_contorno!(GH, KG, comp, BEM.ACABuffer(Float64), nc(dad))

    corrige_diagonais!(dad, HH, GH, dad.normal)

    # compress!(HH, TSVD(atol = 1e-6, rtol = 1e-6))
    # compress!(GH, TSVD(atol = 1e-6, rtol = 1e-6))

    HH, GH
end

function assemble_contorno!(hmat, K, comp, bufs, nc)
    if isleaf(hmat) # base case
        if colperm(hmat)[colrange(hmat)][1] <= nc # se for contorno
            _process_leaf!(hmat, K, comp, bufs)
        end
    else
        # recurse on children
        for child in children(hmat)
            assemble_contorno!(child, K, comp, bufs, nc)
        end
    end
    #return compress!(hmat, TSVD(atol = 1e-6, rtol = 1e-6))
    hmat
end





function calc_H_Hd(dad::potencial; npg = 3, atol = 1e-6, nmax = 10, eta = 3)
    n = size(dad.NOS, 1)
    intelems = zeros(n)
    qsi, w = BEM.gausslegendre(npg)    # Quadratura de gauss
    normal_fonte = dad.normal


    for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelems[elem_j.indices] = integraelemd(x, qsi, w, elem_j)
    end


    splitter = BEM.PrincipalComponentSplitter(; nmax = nmax)
    Xclt, Yclt = ClusterTree(dad, splitter)

    # @infiltrate
    adm = StrongAdmissibilityStd(eta = eta)
    comp = PartialACA(; atol)

    pontos = [SVector{2,Float64}(row) for row in eachrow([dad.NOS; dad.pontos_internos])]

    normal = [SVector{2,Float64}(row) for row in eachrow(normal_fonte)]

    k = dad.k

    #KG = BEM.kernelG(pontos, n, k, intelems)
    KH = BEM.kernelH(pontos, n, normal, intelems)

    #HG = assemble_hmatrix(KG, Xclt, Yclt; adm, comp)
    HH = assemble_hmatrix(KH, Xclt, Yclt; adm, comp)
    # @infiltrate
    corrige_diagonaisH!(HH)

    #HH,HG
    HH
end


function corrige_diagonaisH!(Hmat::HMatrix)
    piv = pivot(Hmat)
    diag = Hmat * ones(size(Hmat, 2))

    for block in nodes(Hmat)
        hasdata(block) || continue
        !block.admissible || continue
        irange = rowrange(block) .- piv[1] .+ 1
        jrange = colrange(block) .- piv[2] .+ 1

        irangeg = rowperm(Hmat)[irange]
        jrangeg = colperm(Hmat)[jrange]
        inds = indexin(irangeg, jrangeg)
        for i = 1:size(irange, 1)
            inds[i] !== nothing || continue

            block.data[i, inds[i]] = -diag[irangeg[i]]
        end
    end
end

function calc_HeG_Hhelm(dad::potencial; npg = 3, atol = 1e-4)
    n = size(dad.NOS, 1)
    intelems = zeros(n)
    qsi, w = BEM.gausslegendre(npg)    # Quadratura de gauss
    normal_fonte = dad.normal

    for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelems[elem_j.indices] = integraelemd(x, qsi, w, elem_j)
        # @infiltrate
    end


    splitter = BEM.PrincipalComponentSplitter(; nmax = 10)
    Xclt, Yclt = ClusterTree(dad, splitter)

    # @infiltrate
    adm = StrongAdmissibilityStd(eta = 3)
    comp = PartialACA(; atol)

    KG = BEM.kernelG(dad, intelems)
    KH = BEM.kernelH(dad, intelems)
    HG = assemble_hmatrix(KG, Xclt, Yclt; adm, comp)
    HH = assemble_hmatrix(KH, Xclt, Yclt; adm, comp)
    # @infiltrate
    corrige_diagonais!(dad, HH, HG, normal_fonte)
    HH, HG
end

function corrige_diagonais!(
    dad::DadosBEM,
    Hmat::HMatrix,
    Gmat::HMatrix,
    normal_fonte::Array{Float64,2},
    #normal_fonte::Vector{SVector{2,Float64}},
)
    piv = pivot(Hmat)
    diag = Hmat * ones(size(Hmat, 2))

    for block in nodes(Hmat)
        hasdata(block) || continue
        !block.admissible || continue
        #block.admissible && continue
        irange = rowrange(block) .- piv[1] .+ 1
        jrange = colrange(block) .- piv[2] .+ 1

        irangeg = rowperm(Hmat)[irange]
        jrangeg = colperm(Hmat)[jrange]
        inds = indexin(irangeg, jrangeg)
        # @infiltrate
        # @show irangeg
        # @show jrangeg
        # @show inds
        for i = 1:size(irange, 1)
            inds[i] !== nothing || continue
            # @infiltrate block.admissible
            block.data[i, inds[i]] = -diag[irangeg[i]]
        end
    end

    diagG = (
        Hmat * [(dad.NOS[:, 1] + dad.NOS[:, 2]); zeros(ni(dad))] +
        Gmat * [(normal_fonte[:, 1] + normal_fonte[:, 2]); zeros(ni(dad))] * dad.k
    )
    diagG[end-ni(dad)+1:end] .= 0
    piv = pivot(Gmat)
    # @infiltrate
    for block in nodes(Gmat)
        hasdata(block) || continue
        !block.admissible || continue
        #block.admissible && continue
        irange = rowrange(block) .- piv[1] .+ 1
        jrange = colrange(block) .- piv[2] .+ 1
        irangeg = rowperm(Hmat)[irange]
        jrangeg = colperm(Hmat)[jrange]
        inds = indexin(irangeg, jrangeg)
        for i = 1:size(irange, 1)
            inds[i] !== nothing || continue
            irangeg[i] <= nc(dad) || continue
            block.data[i, inds[i]] =
                -diagG[irangeg[i]] /
                (normal_fonte[irangeg[i], 1] + normal_fonte[irangeg[i], 2])
        end
    end

end


function aplicaCDC(HH::HMatrix, HG::HMatrix, dad::Union{helmholtz,potencial})
    tipoCDC = BEM.tipoCDC(dad)
    valorCDC = [BEM.valorCDC(dad); zeros(ni(dad))]

    HH_aux = deepcopy(HH) #matrizes aux sao alteradas pela função
    HG_aux = deepcopy(HG)

    aux = Array
    piv = pivot(HH)
    nodesH = nodes(HH_aux)
    nodesG = nodes(HG_aux)
    # @infiltrate
    for block = 1:length(nodesH)
        # @show block, nodesH[block].admissible
        hasdata(nodesH[block]) || continue
        # @show block
        irange = rowrange(nodesH[block]) .- piv[1] .+ 1
        jrange = colrange(nodesH[block]) .- piv[2] .+ 1
        jrangeg = colperm(HH)[jrange]
        # @show tipoCDC[jrangeg]
        tipoCDC[jrangeg[1]] == 0 || continue

        aux = deepcopy(nodesH[block].data)

        nodesH[block].data = -1 * nodesG[block].data
        nodesG[block].data = -1 * (aux)
    end
    HG * valorCDC, HH_aux, HG_aux
    # @infiltrate
end

# function converte_hmat_tens(H::HMatrix)

#     for block in nodes(H)
#         # @show block, nodesH[block].admissible
#         hasdata(block) || continue
#         if !block.admissible
#              block.data=reduce(vcat, [reduce(hcat, row) for row in eachrow(block.data)])
#         else
#              block.data.A=reduce(vcat, [reduce(hcat, row) for row in eachrow(block.data.A)])
#              block.data.B=reduce(vcat, [reduce(hcat, row) for row in eachrow(block.data.B)])

#             end
#     end
#     # @infiltrate
# end
# function reduz_matdemat(A)
#     reduce(vcat, [reduce(hcat, row) for row in eachrow(A)])
# end

function calc_Tid(dad::potencial, T, q, npg = 3)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    nb = nc(dad)
    npi = ni(dad)
    intelems = zeros(nb)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    normal_fonte = dad.normal
    Ti = zeros(npi)
    k = dad.k
    for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelem = integraelemd(x, qsi, w, elem_j)
        # @infiltrate
        intelems[elem_j.indices] = intelem

    end
    for i = 1:npi
        xi = dad.pontos_internos[i, 1]
        yi = dad.pontos_internos[i, 2]
        for j = 1:nb
            xj = dad.NOS[j, 1]
            yj = dad.NOS[j, 2]
            # r = sqrt((xi - xj)^2 + (yi - yj)^2)
            Qast, Tast = calsolfund([xj - xi, yj - yi], normal_fonte[j, :], dad, k)
            Ti[i] += Qast * intelems[j] * T[j] - Tast * intelems[j] * q[j]
        end
    end
    Ti
end




#__________________________________________________________________________________________________________
"Funcao para calcular a das funções de forma"
function integraelemd(x, eta, w, elem)
    h = zeros(Float64, size(elem))

    for k in eachindex(w)
        N, dN = calc_fforma(eta[k], elem)
        # pg = N' * x    # Ponto de gauss interpolador
        dxdqsi = dN' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        # sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        # sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        # @infiltrate
        # @show pg,pf,r,[sy,-sx],Qast
        h += N * dgamadqsi * w[k]
    end
    h
end

function Monta_M_Hd(
    dad::potencial,
    npg;
    tiporadial = "tps",
    atol = 1e-6,
    nmax = 10,
    eta = 3,
    rtol = 1e-7,
)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi, w = gausslegendre(npg)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(n_pontos)
    M1 = zeros(n_pontos)

    splitter = BEM.PrincipalComponentSplitter(; nmax = nmax)
    Xclt, Yclt = ClusterTree(dad, splitter)

    adm = StrongAdmissibilityStd(eta = eta)
    comp = PartialACA(; atol)

    #p_internos = [SVector{2,Float64}(row) for row in eachrow(dad.pontos_internos)]
    #p_contorno = [SVector{2,Float64}(row) for row in eachrow(dad.NOS)]

    pontos = [SVector{2,Float64}(row) for row in eachrow([dad.NOS; dad.pontos_internos])]
    nc_dad = n_nos
    k = dad.k

    func_radial = radial_func(tiporadial)

    KF = BEM.kernelF(pontos, nc_dad, func_radial)

    HF = assemble_hmat(KF, Xclt, Yclt; adm, comp)
    M, M1 = calcMs(dad, npg, tiporadial)

    # @show size(M)
    # @show length(M)
    #@infiltrate
    # S = HF \ M
    S, _ = gmres(HF, M, rtol = rtol, itmax = 100) #ficou melhor com 10^-10
    KD = BEM.kernelD(pontos, nc_dad, k, S)
    HD = assemble_hmat(KD, Xclt, Yclt; adm, comp)
    corrige_diagonal_M!(HD, M1)

    # A = ones(length(M)) * M' / F .* D
    # for i = 1:n_pontos #Laço dos pontos radiais
    #   A[i, i] = 0
    #   A[i, i] = -sum(A[i, :])
    # end
    # A + diagm(0 => M1)
    HD
end


function corrige_diagonal_M!(Hmat::HMatrix, M1::Vector{Float64})
    piv = pivot(Hmat)
    diag = Hmat * ones(size(Hmat, 2)) - M1

    for block in nodes(Hmat)
        hasdata(block) || continue
        #!block.admissible || continue
        block.admissible && continue
        irange = rowrange(block) .- piv[1] .+ 1
        jrange = colrange(block) .- piv[2] .+ 1

        irangeg = rowperm(Hmat)[irange]
        jrangeg = colperm(Hmat)[jrange]

        inds = indexin(irangeg, jrangeg)
        for i = 1:size(irange, 1)
            inds[i] !== nothing || continue
            block.data[i, inds[i]] = -diag[irangeg[i]]
        end
    end
end



# function corrige_diagonal_M!(Hmat::HMatrix, M1::Vector{Float64})
#     piv = pivot(Hmat)
#     ncol = size(Hmat, 2)
#     diag = Hmat * ones(ncol) .- M1

#     for block in nodes(Hmat)
#         hasdata(block) || continue #pula se o bloco não tem dados

#         block.admissible && continue #pula se o bloco for admissível (pois estes não estarão nas diagonais)

#         # obter ranges locais — collect garante que dá um vetor indexável
#         ir_local = collect(rowrange(block)) .- piv[1] .+ 1
#         jr_local = collect(colrange(block)) .- piv[2] .+ 1

#         # converte indices locais para globais
#         ir_g = rowperm(Hmat)[ir_local]
#         jr_g = colperm(Hmat)[jr_local]

#         # colpos faz: colpos[indice global] = indice local. Ou seja, passando por jr_g, coloco correspondencia entre as colunas globais e locais
#         colpos = Dict{Int,Int}()
#         for (local_j, gidx) in enumerate(jr_g)
#             colpos[gidx] = local_j
#         end

#         # com o get, vejo se existe (!= 0) algum par de indices que representa i e j globais, que é = diagonal
#         for (local_i, gidx_i) in enumerate(ir_g)
#             local_j = get(colpos, gidx_i, 0)   # 0 significa "não existe" (evita nothing)
#             if local_j != 0
#                 # ajusta o elemento (linha local i, coluna local j) para -diag[índice_global]
#                 block.data[local_i, local_j] = -diag[gidx_i]
#             end
#         end
#     end

#     return nothing
# end


# function Monta_M_Hdd(dad::potencial, npg; tiporadial = "tps", atol = 1e-6, nmax = 10, eta = 3)
#     n_nos = size(dad.NOS, 1)
#     nelem = size(dad.ELEM, 1)
#     n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

#     qsi, w = gausslegendre(npg)
#     n_pontos = n_nos + n_noi
#     nodes = [dad.NOS; dad.pontos_internos]
#     M = zeros(n_pontos)
#     M1 = zeros(n_pontos)

#     normal = dad.normal
#     k = dad.k
#     intelems = zeros(n_nos)

#     for elem_j in dad.ELEM  #Laço dos elementos
#         x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
#         intelem = integraelemd(x, qsi, w, elem_j)
#         # @infiltrate
#         intelems[elem_j.indices] = intelem
#     end

#     for i = 1:n_pontos #varia no ponto campo
#       xi = nodes[i, 1]
#       yi = nodes[i, 2]

#       #for j =1:n_nos #varia nos pontos da integração de contorno
#       for elem_j in dad.ELEM #varia nos pontos da integração de contorno
#         for j in elem_j.indices
#           xj = nodes[j, 1]
#           yj = nodes[j, 2]

#           r = [xj-xi,yj-yi]
#           R = norm(r)

#           nosing = elem_j.indices .== i

#           if sum(nosing)==1

#             m_el, m_el1 = calc_md(nodes[elem_j.indices, :], [xi,yi], k, qsi, w, elem_j, tiporadial, npg)

#             # if m_el > 10^(-10) || m_el1 > 10^(-10)
#             #   @show m_el,m_el1
#             # end

#             M[i] = M[i] + m_el
#             M1[i] = M1[i] + m_el1
#           else
#             m = int_interpolaρdρ(R, tipo = tiporadial)
#             m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)

#             M[i] = M[i] + dot(normal[j,:], r) / R^2 * m * intelems[j]
#             M1[i] = M1[i] + dot(normal[j,:], r) / R^2 * m1 * intelems[j]
#           end
#         end
#       end
#     end


#     splitter = BEM.PrincipalComponentSplitter(; nmax = nmax)
#     Xclt, Yclt = ClusterTree(dad, splitter)

#     adm = StrongAdmissibilityStd(eta = eta)
#     comp = PartialACA(; atol)

#     pontos = [SVector{2,Float64}(row) for row in eachrow([dad.NOS;dad.pontos_internos])]

#     func_radial = radial_func(tiporadial)

#     KF = BEM.kernelF(pontos, n_nos, func_radial)

#     HF = assemble_hmat(KF, Xclt, Yclt; adm, comp)

#     # @show size(M)
#     # @show length(M)
#     #@infiltrate
#     # S = HF \ M
#     S,aux = gmres(HF, M, rtol = 1e-5, itmax = 100)
#     KD = BEM.kernelD(pontos,n_nos,k,S)
#     HD = assemble_hmat(KD, Xclt, Yclt; adm, comp)
#     corrige_diagonal_M!(HD, M1)

#     HD
# end


# function Monta_M_Hdd2(dad::potencial, npg; tiporadial = "tps", atol = 1e-6, nmax = 10, eta = 3)
#     n_nos = size(dad.NOS, 1)
#     nelem = size(dad.ELEM, 1)
#     n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

#     qsi, w = gausslegendre(npg)
#     n_pontos = n_nos + n_noi
#     nodes = [dad.NOS; dad.pontos_internos]
#     M = zeros(n_pontos)
#     M1 = zeros(n_pontos)

#     normal = dad.normal
#     k = dad.k
#     intelems = zeros(n_nos)
#     for elem_j in dad.ELEM  #Laço dos elementos
#         x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
#         intelem = integraelemd(x, qsi, w, elem_j)
#         # @infiltrate
#         intelems[elem_j.indices] = intelem
#     end

#     for j = 1:n_nos #varia nos pontos da integração de contorno
#       xj = nodes[j, 1]
#       yj = nodes[j, 2]

#       for i = 1:n_pontos #varia no ponto campo
#         if i<j
#           continue
#         else
#           if i<=n_nos
#             xi = nodes[i, 1]
#             yi = nodes[i, 2]

#             r = [xj-xi,yj-yi]
#             R = norm(r)

#             if R==0
#               continue
#             else
#               m = int_interpolaρdρ(R, tipo = tiporadial)
#               m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)

#               M[i] = M[i] + dot(normal[j,:], r) / R^2 * m * (intelems[j]+intelems[i])
#               M1[i] = M1[i] + dot(normal[j,:], r) / R^2 * m1 * (intelems[j]+intelems[i])
#             end
#           else
#             xi = nodes[i, 1]
#             yi = nodes[i, 2]

#             r = [xj-xi,yj-yi]
#             R = norm(r)

#             if R==0
#               continue
#             else
#               m = int_interpolaρdρ(R, tipo = tiporadial)
#               m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)

#               M[i] = M[i] + dot(normal[j,:], r) / R^2 * m * intelems[j]
#               M1[i] = M1[i] + dot(normal[j,:], r) / R^2 * m1 * intelems[j]
#             end
#           end
#         end
#       end
#     end

#     splitter = BEM.PrincipalComponentSplitter(; nmax = nmax)
#     Xclt, Yclt = ClusterTree(dad, splitter)

#     adm = StrongAdmissibilityStd(eta = eta)
#     comp = PartialACA(; atol)

#     pontos = [SVector{2,Float64}(row) for row in eachrow([dad.NOS;dad.pontos_internos])]

#     func_radial = radial_func(tiporadial)

#     KF = BEM.kernelF(pontos, n_nos, func_radial)

#     HF = assemble_hmat(KF, Xclt, Yclt; adm, comp)

#     # @show size(M)
#     # @show length(M)
#     #@infiltrate
#     # S = HF \ M
#     S,aux = gmres(HF, M, rtol = 1e-5, itmax = 100)
#     KD = BEM.kernelD(pontos,n_nos,k,S)
#     HD = assemble_hmat(KD, Xclt, Yclt; adm, comp)
#     corrige_diagonal_M!(HD, M1)

#     HD
# end


function Monta_M_HddH(
    dad::potencial,
    npg;
    tiporadial = "tps",
    atol = 1e-6,
    nmax = 10,
    eta = 3,
    rtol = 1e-5,
)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi, w = gausslegendre(npg)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    k = dad.k
    #vetor M e M1
    M = zeros(n_pontos)
    M1 = zeros(n_pontos)

    #intelems
    intelems = zeros(n_nos)
    for elem_j in dad.ELEM
        x = dad.NOS[elem_j.indices, :]
        intelem = integraelemd(x, qsi, w, elem_j)
        intelems[elem_j.indices] = intelem
    end

    #define variaiveis das hierarquicas
    splitter = BEM.PrincipalComponentSplitter(; nmax = nmax)
    Xclt, Yclt = ClusterTree(dad, splitter)

    V = [SVector{2,Float64}(row) for row in eachrow(dad.NOS)]
    Yclt_cont = ClusterTree(V, splitter)

    # Xrand = [rand(1000) rand(1000)]
    # Xrand = [SVector{2,Float64}(row) for row in eachrow(Xrand)]
    # Yrand = Xrand

    adm_aux = StrongAdmissibilityStd(eta = 10)
    adm = StrongAdmissibilityStd(eta = eta)

    comp = PartialACA(; atol = atol)

    pontos = [SVector{2,Float64}(row) for row in eachrow(nodes)]
    normal = [SVector{2,Float64}(row) for row in eachrow(dad.normal)]

    func_radial = radial_func(tiporadial)
    func_radial_M = int_interpolaρdρ_func(tiporadial)
    func_radial_M1 = R -> -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)

    #kernels e matrizes

    KM = BEM.kernelM(pontos, n_nos, func_radial_M, normal)
    KM1 = BEM.kernelM1(pontos, n_nos, func_radial_M1, normal)

    # HM = comp(KM, Xclt, Yclt_cont)
    # HM1 = comp(KM1, Xclt, Yclt_cont)
    HM = assemble_hmatrix(KM, Xclt, Yclt_cont; adm = adm_aux, comp)
    HM1 = assemble_hmatrix(KM1, Xclt, Yclt_cont; adm = adm_aux, comp)
    # Matrix(HMt) - Matrix(HM)
    # Matrix(HM1t) - Matrix(HM1)
    M = HM * intelems
    M1 = HM1 * intelems

    #matriz F
    KF = BEM.kernelF(pontos, n_nos, func_radial)

    HF = assemble_hmat(KF, Xclt, Yclt; adm = adm_aux, comp)
    #HF = comp(KF)
    # @show rank(HM), rank(HM1), rank(HF)

    #matriz S
    S, aux = gmres(HF, M, rtol = rtol, itmax = 100)

    #Matriz M final após multiplicação com a matriz D
    KD = BEM.kernelD(pontos, n_nos, dad.k, S)
    HD = assemble_hmat(KD, Xclt, Yclt; adm, comp)
    corrige_diagonal_M!(HD, M1)

    HD
end


function Monta_M_HddH2(
    dad::potencial,
    npg;
    tiporadial = "tps",
    atol = 1e-6,
    rtol = 1e-5,
    nmax = 10,
    eta = 3,
)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi, w = gausslegendre(npg)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]

    #vetor M e M1
    M = zeros(n_pontos)
    M1 = zeros(n_pontos)

    #intelems
    intelems = zeros(n_nos)
    for elem_j in dad.ELEM
        x = dad.NOS[elem_j.indices, :]
        intelem = integraelemd(x, qsi, w, elem_j)
        intelems[elem_j.indices] = intelem
    end

    #define variaiveis das hierarquicas
    splitter = BEM.PrincipalComponentSplitter(; nmax = nmax)
    Xclt, Yclt = ClusterTree(dad, splitter)

    V = [SVector{2,Float64}(row) for row in eachrow(dad.NOS)]
    Yclt_cont = ClusterTree(V, splitter)

    adm = StrongAdmissibilityStd(eta = eta)
    comp = PartialACA(; atol)

    pontos = [SVector{2,Float64}(row) for row in eachrow(nodes)]
    normal = [SVector{2,Float64}(row) for row in eachrow(dad.normal)]

    k = dad.k

    func_radial_M = int_interpolaρdρ_func(tiporadial)
    func_radial_M1 = R -> -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)

    #kernels e matrizes

    KM = BEM.kernelM(pontos, n_nos, func_radial_M, normal)
    KM1 = BEM.kernelM1(pontos, n_nos, func_radial_M1, normal)

    M_mat = HMatrix{typeof(Xclt)}(Xclt, Yclt_cont, adm)
    M1_mat = deepcopy(M_mat)

    KM = BEM.PermutedMatrix(KM, BEM.loc2glob(M_mat.rowtree), BEM.loc2glob(M_mat.coltree))
    KM1 = BEM.PermutedMatrix(KM1, BEM.loc2glob(M_mat.rowtree), BEM.loc2glob(M_mat.coltree))

    M = BEM.kernelvec!(M, KM, intelems[BEM.colperm(M_mat)], M_mat, BEM.ACABuffer(Float64))
    invpermute!(M, BEM.rowperm(M_mat))
    M1 = kernelvec!(M1, KM1, intelems[colperm(M_mat)], M1_mat, ACABuffer(Float64))
    invpermute!(M1, BEM.rowperm(M_mat))


    func_radial = radial_func(tiporadial)
    KF = BEM.kernelF(pontos, n_nos, func_radial)
    HF = assemble_hmat(KF, Xclt, Yclt; adm, comp)

    #@infiltrate
    # S = HF \ M
    S, _ = gmres(HF, M, rtol = rtol, itmax = 100) #ficou melhor com 10^-10

    # F, D = FeD(dad, nodes, tiporadial)
    # S = (M'/F)'
    #S = gmres(F, M, rtol = 1e-6, itmax = 100)[1] #ficou melhor com 10^-10

    KD = BEM.kernelD(pontos, n_nos, k, S)
    #KD = BEM.kernelD(pontos,n_nos,k,ones(size(S,1)))
    HD = assemble_hmat(KD, Xclt, Yclt; adm, comp)
    corrige_diagonal_M!(HD, M1)

    HD, M, M1, HF
end


function kernelvec!(
    C::AbstractVector,
    K::AbstractMatrix,
    B::AbstractVector,
    A::HMatrix,
    bufs,
    comp = PartialACA(; atol = 1e-6),
)
    #offset = pivot(A) .- 1
    if isleaf(A)
        # irange = rowrange(A) .- offset[1]
        # jrange = colrange(A) .- offset[2]

        # irange = rowrange(A) .- piv[1] .+ 1
        # jrange = colrange(A) .- piv[2] .+ 1

        irange = rowrange(A)
        jrange = colrange(A)

        if isadmissible(A)
            d = comp(K, A.rowtree, A.coltree, bufs)
            #d = _process_leaf!(A, K, comp, bufs)

        else

            out = Matrix{Float64}(undef, length(irange), length(jrange))
            getblock!(out, K, irange, jrange)
            d = out
            # mulblock!(view(C, irange), K, view(B, jrange), irange, jrange)
        end

        mul!(view(C, irange), d, view(B, jrange), 1, 1)

        # out = Matrix{Float64}(undef, length(irange), length(jrange))
        # getblock!(out, K, irange, jrange)
        # d = out

    else
        for block in children(A)
            kernelvec!(C, K, B, block, bufs)
        end
    end
    return C
end


function Monta_M_Hdd(
    dad::potencial,
    npg;
    tiporadial = "tps",
    atol = 1e-6,
    nmax = 10,
    eta = 3,
    rtol = 1e-5,
)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi, w = gausslegendre(npg)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(n_pontos)
    M1 = zeros(n_pontos)

    normal = dad.normal
    k = dad.k
    intelems = zeros(n_nos)

    teste = @elapsed begin
        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            intelem = integraelemd(x, qsi, w, elem_j)
            # @infiltrate
            intelems[elem_j.indices] = intelem
        end



        for i = 1:n_pontos #varia no ponto campo
            xi = nodes[i, 1]
            yi = nodes[i, 2]

            for j = 1:n_nos #varia nos pontos da integração de contorno
                xj = nodes[j, 1]
                yj = nodes[j, 2]

                r = [xj - xi, yj - yi]
                R = norm(r)

                if R == 0
                    continue
                else
                    m = int_interpolaρdρ(R, tipo = tiporadial)
                    m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)

                    M[i] = M[i] + dot(normal[j, :], r) / R^2 * m * intelems[j]
                    M1[i] = M1[i] + dot(normal[j, :], r) / R^2 * m1 * intelems[j]
                end
            end
        end
    end

    splitter = BEM.PrincipalComponentSplitter(; nmax = nmax)
    Xclt, Yclt = ClusterTree(dad, splitter)

    adm = StrongAdmissibilityStd(eta = eta)
    comp = PartialACA(; atol)

    pontos = [SVector{2,Float64}(row) for row in eachrow([dad.NOS; dad.pontos_internos])]

    func_radial = radial_func(tiporadial)

    KF = BEM.kernelF(pontos, n_nos, func_radial)

    HF = assemble_hmat(KF, Xclt, Yclt; adm, comp)

    #@infiltrate
    # S = HF \ M
    S, _ = gmres(HF, M, rtol = rtol, itmax = 100) #ficou melhor com 10^-10

    # F, D = FeD(dad, nodes, tiporadial)
    # S = (M'/F)'
    #S = gmres(F, M, rtol = 1e-6, itmax = 100)[1] #ficou melhor com 10^-10

    KD = BEM.kernelD(pontos, n_nos, k, S)
    #KD = BEM.kernelD(pontos,n_nos,k,ones(size(S,1)))
    HD = assemble_hmat(KD, Xclt, Yclt; adm, comp)
    corrige_diagonal_M!(HD, M1)

    HD, M, M1, HF
end



function corrige_CDC_M!(Hmat::HMatrix, index)

    for block in nodes(Hmat)
        hasdata(block) || continue
        !block.admissible || continue
        irange = rowrange(block) .- piv[1] .+ 1
        jrange = colrange(block) .- piv[2] .+ 1


        irangeg = rowperm(Hmat)[irange]
        jrangeg = colperm(Hmat)[jrange]

        inds = indexin(irange, jrange)
        for i = 1:size(irange, 1)
            inds[i] !== nothing || continue
            block.data[i, inds[i]] = -diag[irangeg[i]]
        end
    end
end

# function calcMs(dad::potencial, npg)
#   nodes = [dad.NOS; dad.pontos_internos]
#   n_pontos = size(nodes, 1)
#   M = zeros(n_pontos)
#   M1 = zeros(n_pontos)
#   qsi, w = gausslegendre(npg)

#   for i = 1:n_pontos #Laço dos pontos radiais
#     pf = nodes[i, :]
#     for elem_j in dad.ELEM  #Laço dos elementos
#       x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
#       m_el, m_el1 = calc_md(x, pf, dad.k, qsi, w, elem_j)
#       M[i] = M[i] + m_el
#       M1[i] = M1[i] + m_el1
#     end
#   end
#   M, M1
# end
# function calc_md(x, pf, k, qsi, w, elem)
#   npg = length(w)
#   m_el, m_el1 = 0, 0

#   for i = 1:npg
#     N, dN_geo = calc_fforma(qsi[i], elem)
#     pg = N' * x    # Ponto de gauss interpolador
#     r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
#     dxdqsi = dN_geo' * x   # dx/dξ & dy/dξ
#     dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
#     sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
#     sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

#     nx = sy # Componente x do vetor normal unit�rio
#     ny = -sx # Componente y do vetor normal unit�rio
#     # @infiltrate
#     R = norm(r)
#     m = int_interpolaρdρ(R)
#     m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)
#     # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
#     m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
#     m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
#   end
#   return m_el, m_el1
# end


function calc_Tid(dad::Union{potencial,helmholtz}, T, q, npg = 8)
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
            Ti[i] += h' * T[elem_j.indices] - g' * q[elem_j.indices]
        end
    end
    Ti
end

#=
function calc_Ti(dad::Union{potencial,helmholtz}, T, q, npg=8)
  nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.pontos_internos, 1)
  Ti = zeros(n)
  qsi, w = gausslegendre(npg)    # Quadratura de gauss

  for i = 1:n
    pf = dad.pontos_internos[i, :]   # Coordenada (x,y) dos pontos fonte
    for elem_j in dad.ELEM  #Laço dos elementos
      x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
      Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
      eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
      N_geo, ~ = calc_fforma(eet, elem_j)
      ps = N_geo' * x
      b = norm(ps' - pf) / norm(Δelem)
      eta, Jt = sinhtrans(qsi, eet, b)
      h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
      Ti[i] += h' * T[elem_j.indices] - g' * q[elem_j.indices]
    end
  end
  Ti
end

function calc_Ti(dad::potencial_iga, T, q, npg=8)
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

function calc_Aeb(dad::Union{potencial,helmholtz}, npg=8)
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
      eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
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


function calc_HeG(dad::Union{potencial,helmholtz}, b1, b2, npg=8)
  n1 = size(b1, 1)
  n2 = size(b2, 1)
  H = zeros(n1, 0)
  G = zeros(n1, 0)
  qsi, w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM[b2]  #Laço dos elementos
    x = [dad.NOS; dad.pontos_internos][elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos

    h1 = zeros(n1, size(elem_j))
    g1 = zeros(n1, size(elem_j))
    for i = 1:n1
      ind = b1[i]  #Laço dos pontos fontes
      pf = [dad.NOS; dad.pontos_internos][ind, :]   # Coordenada (x,y)  dos pontos fonte
      Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
      eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
      N_geo, ~ = calc_fforma(eet, elem_j)
      ps = N_geo' * x
      b = norm(ps' - pf) / norm(Δelem)
      eta, Jt = sinhtrans(qsi, eet, b)
      h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
      h[elem_j.indices.==ind] = h[elem_j.indices.==ind] .- 0.5

      h1[i, :] += h
      g1[i, :] += g
    end
    H = [H h1]
    G = [G g1]
  end
  H, G
end


function calc_HeG_interp(dad::Union{potencial,helmholtz}, b1, b2, npg=8, ninterp=3)
  collocCoord = [dad.NOS; dad.pontos_internos][b1, :]
  xmax = maximum(collocCoord, dims=1)
  xmin = minimum(collocCoord, dims=1)

  xs = criapontosinterp(ninterp)
  fontes, L, ninterp1, ninterp2 = gera_interpolação(ninterp, collocCoord, xmax, xmin, xs)

  H = zeros(ninterp1 * ninterp2, 0)
  G = zeros(ninterp1 * ninterp2, 0)
  n1, n2 = Nlinear(xs)
  xks = n1 * xmin + n2 * xmax

  qsi, w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM[b2]  #Laço dos elementos
    x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos

    h1 = zeros(ninterp1 * ninterp2, size(elem_j))
    g1 = zeros(ninterp1 * ninterp2, size(elem_j))
    ci = 0
    for i2 = 1:ninterp1
      for i1 = 1:ninterp2
        ci += 1

        pf = [xks[i1, 1], xks[i2, 2]]   # Coordenada (x,y)  dos pontos fonte
        Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
        eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
        N_geo, ~ = calc_fforma(eet, elem_j)
        ps = N_geo' * x
        b = norm(ps' - pf) / norm(Δelem)
        eta, Jt = sinhtrans(qsi, eet, b)
        h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)

        h1[ci, :] += h
        g1[ci, :] += g
      end
    end
    H = [H h1]
    G = [G g1]
  end
  L, H, G
end

function gera_interpolação(ninterp, NOS, xmax, xmin, xs, ϵ=1e-6)
  if (abs(xmax[1] - xmin[1]) < ϵ)
    fontes = (2.0 .* (NOS[:, 2] .- xmin[2]) ./ (xmax[2] - xmin[2]) .- 1)
    L = lagrange(fontes, xs, ninterp)
    ninterp2 = ninterp
    ninterp1 = 1
  elseif (abs(xmax[2] - xmin[2]) < ϵ)
    fontes = (2.0 .* (NOS[:, 1] .- xmin[1]) ./ (xmax[1] - xmin[1]) .- 1)
    L = lagrange(fontes, xs, ninterp)
    ninterp2 = 1
    ninterp1 = ninterp
  else
    fontes = [(2.0 .* (NOS[:, 1] .- xmin[1]) ./ (xmax[1] - xmin[1]) .- 1) (2.0 .* (NOS[:, 2] .- xmin[2]) ./ (xmax[2] - xmin[2]) .- 1)]
    L = lagrange(fontes, xs, ninterp, xs, ninterp)
    ninterp2 = ninterp
    ninterp1 = ninterp
  end
  fontes, L, ninterp1, ninterp2
end



function calc_HeG(dad::potencial_iga, npg=8)
  # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  nfonte = size(dad.NOS, 1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS, 1)
  H = zeros(n, n)
  G = zeros(n, n)
  qsi, w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM  #Laço dos elementos
    xf = elem_j.limites[:, 2]
    x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
    Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
    pc = dad.pontos_controle[:, elem_j.indices]
    # @infiltrate
    cf = pc[1:2, :] ./ pc[4, :]'
    for contafonte = 1:nfonte
      pf = dad.NOS[contafonte, :]   # Coordenada (x,y)  dos pontos fonte
      eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
      N, dN = calc_fforma(eet, elem_j, pc[4, :])
      ps = cf * N
      b = norm(ps - pf) / norm(Δelem)
      eta, Jt = sinhtrans(qsi, eet, b)
      h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)

      H[contafonte, elem_j.indices] += h
      G[contafonte, elem_j.indices] += g

    end
  end
  H - dad.E / 2, G
end
function integrabezier(pf, cf, we, eta, w, elem::bezier, prob::potencial_iga)
  h = zeros(Float64, size(elem))
  g = zeros(Float64, size(elem))
  for k = 1:size(w, 1)
    N, dN = calc_fforma(eta[k], elem, we)
    pg = cf * N    # Ponto de gauss interpolador
    r = pg - pf      # Distancia entre ponto de gauss e ponto fonte
    dxdqsi = cf * dN   # dx/dξ & dy/dξ
    dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
    sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
    sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
    Qast, Tast = calsolfund(r, [sy, -sx], prob)
    # h+=N*dgamadqsi*w[k]
    # g+=N*dgamadqsi*w[k]
    h += N * Qast * dgamadqsi * w[k] / 2
    g += N * Tast * dgamadqsi * w[k] / 2

  end
  h, g
end


function calc_HeG_interp(dad::potencial_iga, b1, b2, npg=8, ninterp=3)
  collocCoord = [dad.NOS; dad.pontos_internos][b1, :]
  xmax = maximum(collocCoord, dims=1)
  xmin = minimum(collocCoord, dims=1)

  xs = criapontosinterp(ninterp)
  fontes, L, ninterp1, ninterp2 = gera_interpolação(ninterp, collocCoord, xmax, xmin, xs)

  H = zeros(ninterp1 * ninterp2, 0)
  G = zeros(ninterp1 * ninterp2, 0)
  n1, n2 = Nlinear(xs)
  xks = n1 * xmin + n2 * xmax

  qsi, w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM[b2]  #Laço dos elementos
    xf = elem_j.limites[:, 2]
    x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
    Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
    pc = dad.pontos_controle[:, elem_j.indices]
    # @infiltrate
    cf = pc[1:2, :] ./ pc[4, :]'

    h1 = zeros(ninterp1 * ninterp2, size(elem_j))
    g1 = zeros(ninterp1 * ninterp2, size(elem_j))
    ci = 0
    for i2 = 1:ninterp1
      for i1 = 1:ninterp2
        ci += 1

        pf = [xks[i1, 1], xks[i2, 2]]   # Coordenada (x,y)  dos pontos fonte
        eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
        N, dN = calc_fforma(eet, elem_j, pc[4, :])
        ps = cf * N
        b = norm(ps - pf) / norm(Δelem)
        eta, Jt = sinhtrans(qsi, eet, b)
        h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)
        h1[ci, :] += h
        g1[ci, :] += g
      end
    end
    H = [H h1]
    G = [G g1]
  end
  L, H, G
end

function calc_HeG(dad::potencial_iga, b1, b2, npg=8)
  n1 = size(b1, 1)
  n2 = size(b2, 1)
  H = zeros(n1, 0)
  G = zeros(n1, 0)
  qsi, w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM[b2]  #Laço dos elementos
    xf = elem_j.limites[:, 2]
    x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
    Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
    pc = dad.pontos_controle[:, elem_j.indices]
    # @infiltrate
    cf = pc[1:2, :] ./ pc[4, :]'
    h1 = zeros(n1, size(elem_j))
    g1 = zeros(n1, size(elem_j))
    for i = 1:n1
      ind = b1[i]  #Laço dos pontos fontes
      pf = [dad.NOS; dad.pontos_internos][ind, :]   # Coordenada (x,y)  dos pontos fonte
      eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
      N, dN = calc_fforma(eet, elem_j, pc[4, :])
      ps = cf * N
      b = norm(ps - pf) / norm(Δelem)
      eta, Jt = sinhtrans(qsi, eet, b)
      h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)
      if ind in elem_j.sing
        h1[i, :] += h - 0.5 * dad.E[ind, elem_j.indices]
      else
        h1[i, :] += h
      end

      g1[i, :] += g
    end
    H = [H h1]
    G = [G g1]
  end
  H, G
end



function calsolfund(r, n, prob::Union{potencial,potencial_iga})
  R = norm(r)
  Qast = dot(r, n) / R^2 / (2 * π)       # Equação 4.36
  Tast = -log(R) / (2 * π * prob.k)
  Qast, Tast
end


=#
