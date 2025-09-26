function calc_HeGd(dad::elastico, npg = 3; interno = false)
    if interno
        n = size(dad.NOS, 1) + size(dad.pontos_internos, 1)
    else
        n = size(dad.NOS, 1)
    end
    intelems = zeros(n)
    H = zeros(2n, 2n)
    G = zeros(2n, 2n)
    qsi, w = BEM.gausslegendre(npg)    # Quadratura de gauss
    qsi2, w2 = BEM.gausslegendre(20)    # Quadratura de gauss

    @showprogress "Montando H e G" for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        # @infiltrate
        intelem = integraelemd(x, qsi, w, elem_j)
        intelems[elem_j.indices] = intelem

    end
    for i = 1:n
        # xi = dad.NOS[i, 1]
        # yi = dad.NOS[i, 2]
        if i > size(dad.NOS, 1)
            pf = dad.pontos_internos[i-nc(dad), :]
        else
            pf = dad.NOS[i, :]
        end
        for elem_j in dad.ELEM  #Laço dos elementos
            xj = dad.NOS[elem_j.indices[1], 1]
            yj = dad.NOS[elem_j.indices[1], 2]
            r = sqrt((pf[1] - xj)^2 + (pf[2] - yj)^2)
            if r > 2 * elem_j.comprimento
                for j in elem_j.indices
                    xj = dad.NOS[j, 1]
                    yj = dad.NOS[j, 2]
                    uast, tast = calsolfund([xj, yj], pf, dad.normal[j, :], dad)

                    H[2i-1:2i, 2j-1:2j] = tast * intelems[j]
                    G[2i-1:2i, 2j-1:2j] = uast * intelems[j]
                end
            else
                x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
                nosing = elem_j.indices .== i
                if sum(nosing) == 1
                    no_pf = findfirst(nosing)
                    eet = elem_j.ξs[no_pf]
                    b = 0.0
                else
                    Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
                    eet =
                        (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                        norm(Δelem)^2 + elem_j.ξs[1]
                    N_geo = calc_fforma(eet, elem_j, false)
                    ps = N_geo' * x
                    b = norm(ps' - pf) / norm(Δelem)
                end
                eta, Jt = sinhtrans(qsi2, eet, b)
                # eta,Jt=telles(qsi,eet)
                # @infiltrate
                h, g = integraelem(pf, x, eta, w2 .* Jt, elem_j, dad)
                # h, g = integraelem(pf, x, qsi, w , elem_j, dad)
                # nosing = elem_j.indices .== i

                # @infiltrate
                # @infiltrate isnan(sum(h))
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                H[2i-1:2i, cols] = h
                G[2i-1:2i, cols] = g
            end

        end
    end



    for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
        H[2i-1:2i, 2i-1:2i] .= 0
        H[2i-1:2i, 2i-1:2i] =
            -[sum(H[2i-1:2i, 1:2:end], dims = 2) sum(H[2i-1:2i, 2:2:end], dims = 2)]
        # H[2i-1, 2i-1] = 0.5
        # H[2i, 2i] = 0.5
    end
    # corrige_diagonais!(dad::elastico, H, G)
    # corrige_diagonais2!(dad::elastico, H, G)
    H, G
end
function calc_HeG_pre(dad::elastico, npg = 20; interno = false)
    if interno
        n = size(dad.NOS, 1) + size(dad.pontos_internos, 1)
    else
        n = size(dad.NOS, 1)
    end
    H = spzeros(2n, 2n)
    G = spzeros(2n, 2n)
    qsi, w = BEM.gausslegendre(npg)    # Quadratura de gauss

    @showprogress "Montando H e G" for i = 1:n
        if i > size(dad.NOS, 1)
            pf = dad.pontos_internos[i-nc(dad), :]
        else
            pf = dad.NOS[i, :]
        end
        for elem_j in dad.ELEM  #Laço dos elementos
            xj = dad.NOS[elem_j.indices[1], 1]
            yj = dad.NOS[elem_j.indices[1], 2]
            r = sqrt((pf[1] - xj)^2 + (pf[2] - yj)^2)
            if r <= 2 * elem_j.comprimento
                x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
                nosing = elem_j.indices .== i
                if sum(nosing) == 1
                    no_pf = findfirst(nosing)
                    eet = elem_j.ξs[no_pf]
                    
                    eta2, w2 = novelquad(1, eet, 1*npg)
                    # eta, Jt = sinhtrans(eta2, eet, 0.0)
                    h, g = integraelem(pf, x, eta2, w2, elem_j, dad)

                else
                    Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
                    eet =
                        (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                        norm(Δelem)^2 + elem_j.ξs[1]
                    N_geo = calc_fforma(eet, elem_j, false)
                    ps = N_geo' * x
                    b = norm(ps' - pf) / norm(Δelem)
                    eta, Jt = sinhtrans(qsi, eet, b)
                    # eta,Jt=telles(qsi,eet)
                    # @infiltrate
                    h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
                end
                # h, g = integraelem(pf, x, qsi, w , elem_j, dad)
                # nosing = elem_j.indices .== i

                # @infiltrate
                # @infiltrate isnan(sum(h))
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                H[2i-1:2i, cols] = h
                G[2i-1:2i, cols] = g
            end

        end
    end
 for i = 1:nc(dad)                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
        H[2i-1:2i, 2i-1:2i] +=
            [.5 0.0; 0.0 .5]
    end
 for i = 1:ni(dad)                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
        H[2*nc(dad)+2i-1:2*nc(dad)+2i, 2*nc(dad)+2i-1:2*nc(dad)+2i] +=
            [1.0 0.0; 0.0 1.0]
    end


    H, G
end

function calc_HeG_Hd(dad::elastico; npg = 3, atol = 1e-4)
    n = size(dad.NOS, 1)
    intelems = zeros(n)
    qsi, w = BEM.gausslegendre(npg)    # Quadratura de gauss

    for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        intelem = integraelemd(x, qsi, w, elem_j)
        # @infiltrate
        intelems[elem_j.indices] = intelem
    end

    # X = [[Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i in 1:nc(dad)]; [Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for i in 1:ni(dad)]]
    # # X = [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i in 1:nc(dad)]
    Hpre,Gpre=calc_HeG_pre(dad,interno = true)
    splitter = BEM.PrincipalComponentSplitter(; nmax = 10)
    # Xclt = Yclt = ClusterTree(X, splitter)
    # Xclt, Yclt = ClusterTree(dad, splitter)
#     elements = [
#     [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
#     [Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for i = 1:ni(dad)]
# ];
# Xclt = Yclt = ClusterTree(repeat(elements, inner = 2),splitter)
 Xclt, Yclt = ClusterTree(dad, splitter)
    # @infiltrate
    adm = StrongAdmissibilityStd(eta = 3)
    comp = PartialACA(; atol)

    KG = BEM.kernelGv(dad, intelems,Gpre)
    KH = BEM.kernelHv(dad, intelems,Hpre)
    HG = assemble_hmat(KG, Xclt, Yclt; adm, comp)
    HH = assemble_hmat(KH, Xclt, Yclt; adm, comp)
    # @infiltrate
    HH, HG
end
function corrige_diagonais!(Hmat::HMatrix)
    piv = pivot(Hmat)
    diag = Hmat * ones(size(Hmat, 2))

    for block in AbstractTrees.PreOrderDFS(Hmat)
        hasdata(block) || continue
        # !block.admissible || continue
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
end


function aplicaCDC(HH::HMatrix, HG::HMatrix, dad::DadosBEM)
    tipoCDC = BEM.tipoCDC(dad)
    valorCDC = [BEM.valorCDC(dad); zeros(ni(dad))]

    aux = Array
    piv = pivot(HH)
    nodesH = nodes(HH)
    nodesG = nodes(HG)
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
    HG * valorCDC
    # @infiltrate
end

# function calc_Tid(dad::potencial, T, q, npg = 3)
#     nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
#     nb = nc(dad)
#     npi = ni(dad)
#     intelems = zeros(nb)
#     qsi, w = gausslegendre(npg)    # Quadratura de gauss
#     normal_fonte = calc_normais(dad)
#     Ti = zeros(npi)
#     for elem_j in dad.ELEM  #Laço dos elementos
#         x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
#         intelem = integraelemd(x, qsi, w, elem_j)
#         # @infiltrate
#         intelems[elem_j.indices] = intelem

#     end
#     for i = 1:npi
#         xi = dad.pontos_internos[i, 1]
#         yi = dad.pontos_internos[i, 2]
#         for j = 1:nb
#             xj = dad.NOS[j, 1]
#             yj = dad.NOS[j, 2]
#             # r = sqrt((xi - xj)^2 + (yi - yj)^2)
#             Qast, Tast = calsolfund([xj - xi, yj - yi], normal_fonte[j, :], dad)
#             Ti[i] += Qast * intelems[j] * T[j] - Tast * intelems[j] * q[j]
#         end
#     end
#     Ti
# end




# #__________________________________________________________________________________________________________
# "Funcao para calcular a das funções de forma"
# function integraelemd(x, eta, w, elem)
#     h = zeros(Float64, size(elem))

#     for k in eachindex(w)
#         N, dN = calc_fforma(eta[k], elem)
#         # pg = N' * x    # Ponto de gauss interpolador
#         dxdqsi = dN' * x   # dx/dξ & dy/dξ
#         dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
#         # sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
#         # sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
#         # @infiltrate
#         # @show pg,pf,r,[sy,-sx],Qast
#         h += N * dgamadqsi * w[k]
#     end
#     h
# end
# function Monta_M_Hd(dad::potencial, npg; atol = 1e-4)
#     n_nos = size(dad.NOS, 1)
#     nelem = size(dad.ELEM, 1)
#     n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

#     qsi, w = gausslegendre(npg)
#     n_pontos = n_nos + n_noi
#     nodes = [dad.NOS; dad.pontos_internos]
#     M = zeros(n_pontos)
#     M1 = zeros(n_pontos)

#     splitter = BEM.PrincipalComponentSplitter(; nmax = 10)
#     Xclt, Yclt = ClusterTree(dad, splitter)

#     adm = StrongAdmissibilityStd()
#     comp = PartialACA(; atol)


#     KF = BEM.kernelF(dad)

#     HF = assemble_hmat(KF, Xclt, Yclt; adm, comp)
#     M, M1 = calcMs(dad, npg)
#     # @show size(M)
#     # @show length(M)
#     # @infiltrate
#     # S = HF \ M
#     S = gmres(HF, M)
#     KD = BEM.kernelD(dad, S)
#     HD = assemble_hmat(KD, Xclt, Yclt; adm, comp)
#     corrige_diagonal_M!(HD, M1)

#     # A = ones(length(M)) * M' / F .* D
#     # for i = 1:n_pontos #Laço dos pontos radiais
#     #   A[i, i] = 0
#     #   A[i, i] = -sum(A[i, :])
#     # end
#     # A + diagm(0 => M1)
#     HD
#     # M, M1, F, D
# end
# struct kernelF <: AbstractMatrix{Float64}
#     dad::DadosBEM
# end
# function Base.getindex(K::kernelF, i::Int, j::Int)
#     if i > nc(K.dad)
#         xi = K.dad.pontos_internos[i-nc(K.dad), 1]
#         yi = K.dad.pontos_internos[i-nc(K.dad), 2]
#     else
#         xi = K.dad.NOS[i, 1]
#         yi = K.dad.NOS[i, 2]
#     end

#     if i == j
#         return 0.0
#     end

#     if j > nc(K.dad)
#         xj = K.dad.pontos_internos[j-nc(K.dad), 1]
#         yj = K.dad.pontos_internos[j-nc(K.dad), 2]
#     else
#         xj = K.dad.NOS[j, 1]
#         yj = K.dad.NOS[j, 2]
#     end
#     r = sqrt((xi - xj)^2 + (yi - yj)^2)
#     interpola(r)
# end
# struct kernelD <: AbstractMatrix{Float64}
#     dad::DadosBEM
#     S::Vector{Float64}
# end
# function Base.getindex(K::kernelD, i::Int, j::Int)
#     if i > nc(K.dad)
#         xi = K.dad.pontos_internos[i-nc(K.dad), 1]
#         yi = K.dad.pontos_internos[i-nc(K.dad), 2]
#     else
#         xi = K.dad.NOS[i, 1]
#         yi = K.dad.NOS[i, 2]
#     end

#     if i == j
#         return 0.0
#     end

#     if j > nc(K.dad)
#         xj = K.dad.pontos_internos[j-nc(K.dad), 1]
#         yj = K.dad.pontos_internos[j-nc(K.dad), 2]
#     else
#         xj = K.dad.NOS[j, 1]
#         yj = K.dad.NOS[j, 2]
#     end
#     r = sqrt((xi - xj)^2 + (yi - yj)^2)
#     -log(r) / (2 * π * K.dad.k) * K.S[j]
# end
# function corrige_diagonal_M!(Hmat::HMatrix, M1::Vector{Float64})
#     piv = pivot(Hmat)
#     diag = Hmat * ones(size(Hmat, 2)) - M1

#     for block in AbstractTrees.PreOrderDFS(Hmat)
#         hasdata(block) || continue
#         !block.admissible || continue
#         irange = rowrange(block) .- piv[1] .+ 1
#         jrange = colrange(block) .- piv[2] .+ 1


#         irangeg = rowperm(Hmat)[irange]
#         jrangeg = colperm(Hmat)[jrange]

#         inds = indexin(irange, jrange)
#         for i = 1:size(irange, 1)
#             inds[i] !== nothing || continue
#             block.data[i, inds[i]] = -diag[irangeg[i]]
#         end
#     end

# end
# # function calcMs(dad::potencial, npg)
# #   nodes = [dad.NOS; dad.pontos_internos]
# #   n_pontos = size(nodes, 1)
# #   M = zeros(n_pontos)
# #   M1 = zeros(n_pontos)
# #   qsi, w = gausslegendre(npg)

# #   for i = 1:n_pontos #Laço dos pontos radiais
# #     pf = nodes[i, :]
# #     for elem_j in dad.ELEM  #Laço dos elementos
# #       x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
# #       m_el, m_el1 = calc_md(x, pf, dad.k, qsi, w, elem_j)
# #       M[i] = M[i] + m_el
# #       M1[i] = M1[i] + m_el1
# #     end
# #   end
# #   M, M1
# # end
# # function calc_md(x, pf, k, qsi, w, elem)
# #   npg = length(w)
# #   m_el, m_el1 = 0, 0

# #   for i = 1:npg
# #     N, dN_geo = calc_fforma(qsi[i], elem)
# #     pg = N' * x    # Ponto de gauss interpolador
# #     r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
# #     dxdqsi = dN_geo' * x   # dx/dξ & dy/dξ
# #     dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
# #     sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
# #     sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

# #     nx = sy # Componente x do vetor normal unit�rio
# #     ny = -sx # Componente y do vetor normal unit�rio
# #     # @infiltrate
# #     R = norm(r)
# #     m = int_interpolaρdρ(R)
# #     m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)
# #     # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
# #     m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
# #     m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
# #   end
# #   return m_el, m_el1
# # end


# #=
# function calc_Ti(dad::Union{potencial,helmholtz}, T, q, npg=8)
#   nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
#   n = size(dad.pontos_internos, 1)
#   Ti = zeros(n)
#   qsi, w = gausslegendre(npg)    # Quadratura de gauss

#   for i = 1:n
#     pf = dad.pontos_internos[i, :]   # Coordenada (x,y) dos pontos fonte
#     for elem_j in dad.ELEM  #Laço dos elementos
#       x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
#       Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#       eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
#       N_geo, ~ = calc_fforma(eet, elem_j)
#       ps = N_geo' * x
#       b = norm(ps' - pf) / norm(Δelem)
#       eta, Jt = sinhtrans(qsi, eet, b)
#       h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
#       Ti[i] += h' * T[elem_j.indices] - g' * q[elem_j.indices]
#     end
#   end
#   Ti
# end
# function calc_Ti(dad::potencial_iga, T, q, npg=8)
#   nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
#   n = size(dad.pontos_internos, 1)
#   Ti = zeros(n)
#   qsi, w = gausslegendre(npg)    # Quadratura de gauss

#   for i = 1:n
#     pf = dad.pontos_internos[i, :]   # Coordenada (x,y) dos pontos fonte
#     for elem_j in dad.ELEM  #Laço dos elementos
#       xf = elem_j.limites[:, 2]
#       x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#       Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
#       eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
#       pc = dad.pontos_controle[:, elem_j.indices]
#       # @infiltrate
#       cf = pc[1:2, :] ./ pc[4, :]'
#       N, dN = calc_fforma(eet, elem_j, pc[4, :])
#       ps = cf * N
#       b = norm(ps - pf) / norm(Δelem)

#       eta, Jt = sinhtrans(qsi, eet, b)
#       h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)

#       Ti[i] += h' * T[elem_j.indices] - g' * q[elem_j.indices]
#     end
#   end
#   Ti
# end

# function calc_Aeb(dad::Union{potencial,helmholtz}, npg=8)
#   nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
#   n = size(dad.NOS, 1)
#   ni = size(dad.pontos_internos, 1)
#   A = zeros(n + ni, n)
#   B = zeros(n + ni)
#   qsi, w = gausslegendre(npg)    # Quadratura de gauss

#   for i = 1:n+ni  #Laço dos pontos fontes
#     pf = [dad.NOS; dad.pontos_internos][i, :]   # Coordenada (x,y)  dos pontos fonte
#     for elem_j in dad.ELEM  #Laço dos elementos
#       x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
#       Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#       eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
#       N_geo, ~ = calc_fforma(eet, elem_j)
#       ps = N_geo' * x
#       b = norm(ps' - pf) / norm(Δelem)
#       eta, Jt = sinhtrans(qsi, eet, b)
#       h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
#       h[elem_j.indices.==i] = h[elem_j.indices.==i] .- 0.5
#       if elem_j.tipoCDC == 1
#         A[i, elem_j.indices] = h
#         B[i] += dot(g, elem_j.valorCDC)
#       else
#         A[i, elem_j.indices] = -g
#         B[i] += -dot(h, elem_j.valorCDC)
#       end

#     end
#   end
#   [A [zeros(n, ni); -diagm(ones(ni))]], B
# end


# function calc_HeG(dad::Union{potencial,helmholtz}, b1, b2, npg=8)
#   n1 = size(b1, 1)
#   n2 = size(b2, 1)
#   H = zeros(n1, 0)
#   G = zeros(n1, 0)
#   qsi, w = gausslegendre(npg)    # Quadratura de gauss
#   for elem_j in dad.ELEM[b2]  #Laço dos elementos
#     x = [dad.NOS; dad.pontos_internos][elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos

#     h1 = zeros(n1, size(elem_j))
#     g1 = zeros(n1, size(elem_j))
#     for i = 1:n1
#       ind = b1[i]  #Laço dos pontos fontes
#       pf = [dad.NOS; dad.pontos_internos][ind, :]   # Coordenada (x,y)  dos pontos fonte
#       Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#       eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
#       N_geo, ~ = calc_fforma(eet, elem_j)
#       ps = N_geo' * x
#       b = norm(ps' - pf) / norm(Δelem)
#       eta, Jt = sinhtrans(qsi, eet, b)
#       h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
#       h[elem_j.indices.==ind] = h[elem_j.indices.==ind] .- 0.5

#       h1[i, :] += h
#       g1[i, :] += g
#     end
#     H = [H h1]
#     G = [G g1]
#   end
#   H, G
# end


# function calc_HeG_interp(dad::Union{potencial,helmholtz}, b1, b2, npg=8, ninterp=3)
#   collocCoord = [dad.NOS; dad.pontos_internos][b1, :]
#   xmax = maximum(collocCoord, dims=1)
#   xmin = minimum(collocCoord, dims=1)

#   xs = criapontosinterp(ninterp)
#   fontes, L, ninterp1, ninterp2 = gera_interpolação(ninterp, collocCoord, xmax, xmin, xs)

#   H = zeros(ninterp1 * ninterp2, 0)
#   G = zeros(ninterp1 * ninterp2, 0)
#   n1, n2 = Nlinear(xs)
#   xks = n1 * xmin + n2 * xmax

#   qsi, w = gausslegendre(npg)    # Quadratura de gauss
#   for elem_j in dad.ELEM[b2]  #Laço dos elementos
#     x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos

#     h1 = zeros(ninterp1 * ninterp2, size(elem_j))
#     g1 = zeros(ninterp1 * ninterp2, size(elem_j))
#     ci = 0
#     for i2 = 1:ninterp1
#       for i1 = 1:ninterp2
#         ci += 1

#         pf = [xks[i1, 1], xks[i2, 2]]   # Coordenada (x,y)  dos pontos fonte
#         Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#         eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
#         N_geo, ~ = calc_fforma(eet, elem_j)
#         ps = N_geo' * x
#         b = norm(ps' - pf) / norm(Δelem)
#         eta, Jt = sinhtrans(qsi, eet, b)
#         h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)

#         h1[ci, :] += h
#         g1[ci, :] += g
#       end
#     end
#     H = [H h1]
#     G = [G g1]
#   end
#   L, H, G
# end

# function gera_interpolação(ninterp, NOS, xmax, xmin, xs, ϵ=1e-6)
#   if (abs(xmax[1] - xmin[1]) < ϵ)
#     fontes = (2.0 .* (NOS[:, 2] .- xmin[2]) ./ (xmax[2] - xmin[2]) .- 1)
#     L = lagrange(fontes, xs, ninterp)
#     ninterp2 = ninterp
#     ninterp1 = 1
#   elseif (abs(xmax[2] - xmin[2]) < ϵ)
#     fontes = (2.0 .* (NOS[:, 1] .- xmin[1]) ./ (xmax[1] - xmin[1]) .- 1)
#     L = lagrange(fontes, xs, ninterp)
#     ninterp2 = 1
#     ninterp1 = ninterp
#   else
#     fontes = [(2.0 .* (NOS[:, 1] .- xmin[1]) ./ (xmax[1] - xmin[1]) .- 1) (2.0 .* (NOS[:, 2] .- xmin[2]) ./ (xmax[2] - xmin[2]) .- 1)]
#     L = lagrange(fontes, xs, ninterp, xs, ninterp)
#     ninterp2 = ninterp
#     ninterp1 = ninterp
#   end
#   fontes, L, ninterp1, ninterp2
# end



# function calc_HeG(dad::potencial_iga, npg=8)
#   # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
#   nfonte = size(dad.NOS, 1)    # Quantidade de elementos discretizados no contorno
#   n = size(dad.NOS, 1)
#   H = zeros(n, n)
#   G = zeros(n, n)
#   qsi, w = gausslegendre(npg)    # Quadratura de gauss
#   for elem_j in dad.ELEM  #Laço dos elementos
#     xf = elem_j.limites[:, 2]
#     x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#     Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
#     pc = dad.pontos_controle[:, elem_j.indices]
#     # @infiltrate
#     cf = pc[1:2, :] ./ pc[4, :]'
#     for contafonte = 1:nfonte
#       pf = dad.NOS[contafonte, :]   # Coordenada (x,y)  dos pontos fonte
#       eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
#       N, dN = calc_fforma(eet, elem_j, pc[4, :])
#       ps = cf * N
#       b = norm(ps - pf) / norm(Δelem)
#       eta, Jt = sinhtrans(qsi, eet, b)
#       h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)

#       H[contafonte, elem_j.indices] += h
#       G[contafonte, elem_j.indices] += g

#     end
#   end
#   H - dad.E / 2, G
# end
# function integrabezier(pf, cf, we, eta, w, elem::bezier, prob::potencial_iga)
#   h = zeros(Float64, size(elem))
#   g = zeros(Float64, size(elem))
#   for k = 1:size(w, 1)
#     N, dN = calc_fforma(eta[k], elem, we)
#     pg = cf * N    # Ponto de gauss interpolador
#     r = pg - pf      # Distancia entre ponto de gauss e ponto fonte
#     dxdqsi = cf * dN   # dx/dξ & dy/dξ
#     dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
#     sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
#     sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
#     Qast, Tast = calsolfund(r, [sy, -sx], prob)
#     # h+=N*dgamadqsi*w[k]
#     # g+=N*dgamadqsi*w[k]
#     h += N * Qast * dgamadqsi * w[k] / 2
#     g += N * Tast * dgamadqsi * w[k] / 2

#   end
#   h, g
# end


# function calc_HeG_interp(dad::potencial_iga, b1, b2, npg=8, ninterp=3)
#   collocCoord = [dad.NOS; dad.pontos_internos][b1, :]
#   xmax = maximum(collocCoord, dims=1)
#   xmin = minimum(collocCoord, dims=1)

#   xs = criapontosinterp(ninterp)
#   fontes, L, ninterp1, ninterp2 = gera_interpolação(ninterp, collocCoord, xmax, xmin, xs)

#   H = zeros(ninterp1 * ninterp2, 0)
#   G = zeros(ninterp1 * ninterp2, 0)
#   n1, n2 = Nlinear(xs)
#   xks = n1 * xmin + n2 * xmax

#   qsi, w = gausslegendre(npg)    # Quadratura de gauss
#   for elem_j in dad.ELEM[b2]  #Laço dos elementos
#     xf = elem_j.limites[:, 2]
#     x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#     Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
#     pc = dad.pontos_controle[:, elem_j.indices]
#     # @infiltrate
#     cf = pc[1:2, :] ./ pc[4, :]'

#     h1 = zeros(ninterp1 * ninterp2, size(elem_j))
#     g1 = zeros(ninterp1 * ninterp2, size(elem_j))
#     ci = 0
#     for i2 = 1:ninterp1
#       for i1 = 1:ninterp2
#         ci += 1

#         pf = [xks[i1, 1], xks[i2, 2]]   # Coordenada (x,y)  dos pontos fonte
#         eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
#         N, dN = calc_fforma(eet, elem_j, pc[4, :])
#         ps = cf * N
#         b = norm(ps - pf) / norm(Δelem)
#         eta, Jt = sinhtrans(qsi, eet, b)
#         h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)
#         h1[ci, :] += h
#         g1[ci, :] += g
#       end
#     end
#     H = [H h1]
#     G = [G g1]
#   end
#   L, H, G
# end

# function calc_HeG(dad::potencial_iga, b1, b2, npg=8)
#   n1 = size(b1, 1)
#   n2 = size(b2, 1)
#   H = zeros(n1, 0)
#   G = zeros(n1, 0)
#   qsi, w = gausslegendre(npg)    # Quadratura de gauss
#   for elem_j in dad.ELEM[b2]  #Laço dos elementos
#     xf = elem_j.limites[:, 2]
#     x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
#     Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
#     pc = dad.pontos_controle[:, elem_j.indices]
#     # @infiltrate
#     cf = pc[1:2, :] ./ pc[4, :]'
#     h1 = zeros(n1, size(elem_j))
#     g1 = zeros(n1, size(elem_j))
#     for i = 1:n1
#       ind = b1[i]  #Laço dos pontos fontes
#       pf = [dad.NOS; dad.pontos_internos][ind, :]   # Coordenada (x,y)  dos pontos fonte
#       eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
#       N, dN = calc_fforma(eet, elem_j, pc[4, :])
#       ps = cf * N
#       b = norm(ps - pf) / norm(Δelem)
#       eta, Jt = sinhtrans(qsi, eet, b)
#       h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)
#       if ind in elem_j.sing
#         h1[i, :] += h - 0.5 * dad.E[ind, elem_j.indices]
#       else
#         h1[i, :] += h
#       end

#       g1[i, :] += g
#     end
#     H = [H h1]
#     G = [G g1]
#   end
#   H, G
# end



# function calsolfund(r, n, prob::Union{potencial,potencial_iga})
#   R = norm(r)
#   Qast = dot(r, n) / R^2 / (2 * π)       # Equação 4.36
#   Tast = -log(R) / (2 * π * prob.k)
#   Qast, Tast
# end


# =#