function define_SubRegioes(LINHA, DISCRE)
    nlinhas = size(LINHA, 1)
    j = 1
    k = 1
    kk = 1
    jj = 1
    SUBREGIAO = Vector{Vector{Int}}(undef, 1)
    SUBREGIAO[1] = []
    ponto_ini = LINHA[1, 2]

    for i = 1:nlinhas
        if k == 1
            ponto_ini = LINHA[i, 2]
        end
        p1 = LINHA[i, 2]
        p2 = LINHA[i, 3]
        push!(SUBREGIAO[j], i)
        k += 1
        if p2 == ponto_ini && i != nlinhas
            push!(SUBREGIAO, [])
            j += 1
            k = 1
        end
    end

    interfaces = Vector{Vector{Int}}(undef, 0)
    nreg = size(SUBREGIAO, 1)
    iint = 0
    for i = 1:nreg-1
        nlreg = size(SUBREGIAO[i], 1)
        for k = 1:nlreg
            ii = SUBREGIAO[i][k]
            if ii != 0
                p1 = LINHA[ii, 2]
                p2 = LINHA[ii, 3]
            elseif jj == 0
                p11 = 1000000
                p22 = 1000000
            end
            for j = i+1:nreg
                nlregj = size(SUBREGIAO[j], 1)
                for kk = 1:nlregj
                    jj = SUBREGIAO[j][kk]
                    if jj != 0
                        p11 = LINHA[jj, 2]
                        p22 = LINHA[jj, 3]
                    elseif jj == 0
                        p11 = 0
                        p22 = 0
                    end
                    if (p1 == p11 && p2 == p22) || (p1 == p22 && p2 == p11)
                        iint += 1
                        nelem = DISCRE[ii, 2]
                        nnos = DISCRE[ii, 3] * nelem
                        if p1 == p11 && p2 == p22
                            beta = 1
                        else
                            beta = -1
                        end
                        push!(interfaces, [iint, i, j, ii, jj, nnos, beta])

                    end
                end
            end
        end
    end
    nelem = [sum(DISCRE[s, 2]) for s in SUBREGIAO]
    celem = cumsum(nelem)
    cnnos = [0; cumsum(DISCRE[:, 2] .* DISCRE[:, 3])]
    n_i = size(interfaces, 1)
    equivale = zeros(Int64, 0, 5)
    for i = 1:n_i
        n1 = cnnos[interfaces[i][4]]+1:cnnos[interfaces[i][4]+1]
        n2 = cnnos[interfaces[i][5]]+1:cnnos[interfaces[i][5]+1]
        # @infiltrate
        if interfaces[i][7] == 1
            equivale = [equivale; [n1 n2 repeat([i interfaces[i][2] interfaces[i][3]], length(n1))]]
        else
            equivale = [equivale; [n1 n2[end:-1:1] repeat([i interfaces[i][2] interfaces[i][3]], length(n1))]]
        end
    end
    elem_reg = ones(celem[end])
    for i = 2:size(celem, 1)
        elem_reg[celem[i-1]+1:celem[i]] .= i
    end
    Hc = compatibilidade(cnnos[end], equivale)
    # @infiltrate
    return subregioes(elem_reg, equivale, Hc)
end
function compatibilidade(nc, equivale)
    neq = size(equivale, 1)
    Hc = zeros(4neq, 2nc + 4neq)
    for i = 1:neq
        Hc[4i-3, 2equivale[i, 1]-1] = 1
        Hc[4i-2, 2equivale[i, 1]] = 1
        Hc[4i-3, 2equivale[i, 2]-1] = -1
        Hc[4i-2, 2equivale[i, 2]] = -1
        # @infiltrate
        Hc[4i-1, 2nc+4i-3] = 1
        Hc[4i, 2nc+4i-2] = 1
        Hc[4i-1, 2nc+4i-1] = 1
        Hc[4i, 2nc+4i] = 1
    end
    Hc
end


temsubregioes(dad) = haskey(dad.k, :subregioes)

function fix_interface!(H, dad)
    ni = size(dad.k.interface, 1)
    for i = 1:ni
        indi = dad.k.subregioes.equivale[:, 3] .== i
        ind1 = dad.k.subregioes.equivale[indi, 1]
        ind2 = dad.k.subregioes.equivale[indi, 2]
        for elem in dad.ELEM[dad.k.subregioes.regiao.==dad.k.interface[2]]
            cols = [2elem.indices .- 1 2elem.indices]'[:]
            H[2ind1, cols] .== 0
            H[2ind1.-1, cols] .== 0
        end
        for elem in dad.ELEM[dad.k.subregioes.regiao.==dad.k.interface[3]]
            cols = [2elem.indices .- 1 2elem.indices]'[:]
            H[2ind2, cols] .== 0
            H[2ind2.-1, cols] .== 0
        end
    end
    H
end


function regiao_NO(dad)
    regs = ones(Int, size(dad.NOS, 1))
    for elem in dad.ELEM
        regs[elem.indices] .= elem.regiao
    end
    regs
end