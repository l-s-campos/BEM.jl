function divnode(X, t)
    # Realiza divisão binária dos nós da malha;
    # X = matrix que contêm as coordenadas dos nós;{t,2}
    # t = vetor com os números dos nós; [t]
    n = length(t)   # Quantidade de nós
    x = X[t, :] # Matrix com as coordenadas dos nós que pertencem ao bloco a ser dividido; {t,2}
    c = mean(x, dims = 1)   # Vetor com as coordenadas do centro geométrico do conjunto de nós;{1,2}
    #mean(x,1) = média ao longo da 1 dimensão da matrix (1dim = linhas).
    covx = cov(x)                           # Calcula matrix de covariância de x
    eig_valx, eig_vecx = eigen(covx)          # Calcula autovalores e autovetores da matriz de covariância de x
    ref = eig_vecx[:, argmax(eig_valx)]      # Define como referencia o autovetor relacionado ao maior autovalor
    # a direção desse autovetor e a direção e maior variabilidade dos dados

    attcond = [(x.-c)[i, :]' * ref for i = 1:n] # Condição que divide os nos em dois blocos diferentes.

    x1 = t[attcond.>=0]         # Bloco tal que a condição e >= 0
    x2 = t[attcond.<0]          # Bloco tal que a condição e < 0
    diam = 2 * maximum(sqrt.(((x.-c).*(x.-c))[:, 1] + ((x.-c).*(x.-c))[:, 2]))
    # Calcula diametro do conjunto de dados, centralizando eles; dia = 2*norma do ponto mais distante do centro
    return x1, x2, diam, c
end

function cluster(dad; max_elem = 4, η = 1.0)
    # X = Coordenadas (x,y) dos nós
    # max_elem = Define máximo de nós em cada folha, tal que: max_elem/2 =< nós em cada folha < max_elem
    centro = centros(dad)
    collocCoord = [dad.NOS; dad.pontos_internos]
    Tree1, child1, center_row1, diam1 = Tree(collocCoord, max_elem)
    if typeof(dad) == potencial_iga
        Tree2, child2, center_row2, diam2 =
            Tree(centro, max_elem, [dad.tipoCDC[i.indices[1]] for i in dad.ELEM])
    else
        Tree2, child2, center_row2, diam2 =
            Tree(centro, max_elem, [i.tipoCDC for i in dad.ELEM])
    end
    n1 = size(Tree1, 1)
    n2 = size(Tree2, 1)

    admiss = zeros(n1, n2)
    # Cria matriz para alocar o resultado da aplicação da condição de admissibilidade entre os blocos
    for i = 1:n1     # Para todos os nós da malha
        for j = 1:n2
            admiss[i, j] =
                η * norm(center_row1[i] - center_row2[j], 2) - max(diam1[i], diam2[j])
            # Condição de admissibilidade, para satisfazer deve ser > 0
        end
    end
    allow = admiss .>= 0 # Salva blocos onde a condição e satisfeita
    # @infiltrate
    block = blocks(Tree1, child1, Tree2, child2, allow) # Função que retorna os blocos admissiveis
    return Tree1, Tree2, block
end

function blocks(Tree1, child1, Tree2, child2, allow)
    fc1 = [2; 2; 3; 3]   # Primeiros Blocos a serem avaliados
    fc2 = [2; 3; 2; 3]   # Primeiros Blocos a serem avaliados
    # fc1(1) e fc(2) formam blocos a serem analisados -> (22, 23, 32, 33)
    block = zeros(Int, 0, 3)
    # Matrix que aloca os blocos admissiveis [:,1:2]
    # e se atende a condicao de admissiblidade [:,3]
    c1 = 0     # Contador
    # @infiltrate
    while c1 < length(fc1) / 2
        for i = 1:2
            # @infiltrate fc1[end]==0
            # @show fc1
            # @show fc2
            # @show c1,fc1[c1*2+i],fc2[c1*2+i]
            if allow[fc1[c1*2+i], fc2[c1*2+i]] == 1    # Se blocos são admissiveis
                block = vcat(block, [fc1[c1*2+i] fc2[c1*2+i] 1])
                # Adicionar blocos e identificador 1 (admissivel) a proxima linha matrix block
            else # Se blocos não são admissiveis
                if child1[fc1[c1*2+i], 1] == 0 && child2[fc2[c1*2+i], 1] == 0
                    # Se ambos os blocos não tem filhos, ou seja, se ambos sao folhas
                    block = vcat(block, [fc1[c1*2+i] fc2[c1*2+i] 0])
                    # Adicionar blocos e identificador 0 (não admissivel) a proxima linha matrix block
                else
                    if length(Tree1[fc1[c1*2+i]]) >= length(Tree2[fc2[c1*2+i]]) &&
                       child1[fc1[c1*2+i], 1] != 0
                        # Se a quantidade de elementos no bloco Tree[fc1[...]]] e >=  Tree[fc2[...]]
                        fc1 = [fc1; child1[fc1[c1*2+i], :]]        # Adiciona filhos a fc1[...]
                        fc2 = [fc2; fc2[c1*2+i]; fc2[c1*2+i]]    # Repete elemento de fc2[...]
                    else
                        fc1 = [fc1; fc1[c1*2+i]; fc1[c1*2+i]]    # Repete elemento de fc1[...]
                        fc2 = [fc2; child2[fc2[c1*2+i], :]]        # Adiciona filhos a fc2[...]
                    end
                end
            end
        end
        c1 = c1 + 1  # Atualiza contador
    end
    return block  #Matriz que contem os blocos analisados e sua condicao de admissibilidade
end

function matvec(hmat, b, block, Tree1, Tree2, dad)
    v = b * 0

    for i = 1:length(block[:, 3])
        b1 = Tree1[block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        cols = vcat([dad.ELEM[i].indices for i in b2]...)
        if block[i, 3] == 1
            v[b1] += hmat[i, 1] * (hmat[i, 2] * b[cols])
        else
            v[b1] += hmat[i, 1] * b[cols]
        end
    end
    v[nc(dad)+1:end] -= b[nc(dad)+1:end]
    v
end

function colsblocks(block, Tree2, dad)
    # cols= BEM.colsblocks(block,Tree2,dad)

    # v=b*0
    cols = Vector{Vector{Int64}}(undef, length(block[:, 3]))
    for i = 1:length(block[:, 3])
        # b1 =  @views Tree1[block[i,1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = @views Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        cols[i] = vcat([dad.ELEM[i].indices for i in b2]...)
        # if block[i,3]==1
        #     v[b1]+=hmat[i,1]*(hmat[i,2]*b[cols])
        # else
        #       v[b1]+=hmat[i,1]*b[cols]
        # end
    end
    cols
end
*(hmat::hmat, b) = matvec(hmat, b)
function matvec(hmat::hmat, b)
    v = b * 0
    for i = 1:length(hmat.block[:, 3])
        b1 = @views hmat.Tree1[hmat.block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)

        if hmat.block[i, 3] == 1
            v[b1] += @views hmat.A[i, 1] * (hmat.A[i, 2] * b[hmat.cols[i]])
            #  mul!(@view(v[b1]),hmat[i,1],hmat[i,2]*@view(b[cols[i]]),1,1)
        else
            #  v[b1]=v[b1]+mul!(v[b1],hmat[i,1],(b[cols[i]]))
            v[b1] += @views hmat.A[i, 1] * (b[hmat.cols[i]])
            # mul!(@view(v[b1]),hmat[i,1],@view(b[cols[i]]),1,1)

        end
    end
    v[nc(hmat.dad)+1:end] -= @views b[nc(hmat.dad)+1:end]
    v
end

function tamanho(hmat, block, Tree)
    A = 0
    for i = 1:length(block[:, 3])
        if block[i, 3] == 1
            A += length(hmat[i, 1]) + length(hmat[i, 2])
        else
            A += length(hmat[i, 1])
        end
    end
    A
end




function centros(dad)
    n = length(dad.ELEM)  # Number of curves
    # ind=Array{Array}(undef, n)
    centros = zeros(n, 2)
    for i = 1:n
        x = dad.NOS[dad.ELEM[i].indices, :]
        centros[i, :] = mean(x, dims = 1)
    end
    # centros,ind
    centros
end

function Tree(X, max_elem, tipoCDC = 1)
    m, n = size(X)                     # Tamanho da matriz contendo as coordernadas de cada nós {m,2}
    max_clt = ceil(Int, 2 * m / max_elem)  # Define limite superior para tamanho (nº de linhas) das matrizes e vetores utilizados
    child1 = zeros(1, 2 * max_clt)
    child2 = zeros(1, 2 * max_clt)
    t = collect(1:m)                # Nós
    inode = 1                       # Começa a contagem de nós da árvore
    ileaf = 1                       # Começa a contagem de folhas da árvore
    nodes = Array{Any}(undef, 2 * max_clt)   # Aloca um vetor de vetores para guardar os nós da malha pertencentes a cada nó da árvore
    # Nodes[i] = vetor com os nós da malha pertencentes ao nó i da árvore
    leaves = Array{Any}(undef, 2 * max_clt)  # Aloca um vetor para guardar as folhas
    child = zeros(Int, 2 * max_clt, 2)  # Aloca uma matriz para guardar os filhos de cada nó.
    # Child[i,:] = filhos do nó i
    if tipoCDC == 1
        nodes[1] = t                     # O 1º nó da árvore (raiz) contem todos os nós da malha.
        i = 1
    else
        nodes[1] = t
        nodes[2] = t[tipoCDC.==1]
        nodes[3] = t[tipoCDC.==0]
        inode = inode + 2
        child[1, 1] = 2
        child[1, 2] = 3
        i = 2
    end
    center_row = zeros(2 * max_clt, 2)  # Aloca um vetor para guardar o centro geometrico de cada bloco de nós da malha.
    diam = zeros(2 * max_clt)          # Aloca um vetor ppara guardar o diametro de cada bloco de nós da malha.
    while inode >= ileaf             # Enquanto o quantidade de nós for maior que a de folhas.
        # Observe que a condição só não vai ser satisfeita quando a árvore estiver completa.
        t1, t2, d, c = divnode(X, nodes[i])      # Executa a rotina que divide os nós da malha.
        center_row[i, :] = c                 # Salva centro geometrico do nó i da árvore
        diam[i] = d                         # Salva diametro do nó i da árvore
        if length(t1) > max_elem          # Se a quantidade de nós em t1 for maior que max_elem, definir como nó
            inode = inode + 1            # Chama proximo nó
            nodes[inode] = t1            # Define t1 como um nó
            child[i, 1] = inode           # Define t1 como filho do nó i
        else                             # Se a quantidade de nós for menor que max_elem, definir como folha
            leaves[ileaf] = t1           # Define t1 como uma folha
            ileaf = ileaf + 1            # Chama proxima folha
            child1[i] = ileaf            # Define t1 como folha do nó i
        end
        # Realiza o mesmo para t2---------------------------------------
        if length(t2) > max_elem
            inode = inode + 1
            nodes[inode] = t2
            child[i, 2] = inode
        else
            leaves[ileaf] = t2
            ileaf = ileaf + 1
            child2[i] = ileaf
        end
        # --------------------------------------------------------------
        i = i + 1
    end

    Tree = Array{Vector{Int64}}(undef, inode + ileaf - 1)    # Define tamanho da árvore
    for i = 1:inode                       # Para todos os nós
        Tree[i] = nodes[i]              # Insere os nós na árvore
        if child1[i] > 0                # Se aquele nó tem folhas
            child[i, 1] = child1[i] + inode - 1   # Adiciona as folhas pares na matriz child
        end
        if child2[i] > 0                         # Se aquele nó tem folhas
            child[i, 2] = child2[i] + inode - 1   # Adiciona as folhas impares na matriz child
        end
    end
    for i = 1:ileaf-1     # Para todos as folhas
        Tree[inode+i] = leaves[i]   # Insere as folhas na árvore
        x = X[leaves[i], :]
        c = mean(x, dims = 1)
        d = 2 * maximum(sqrt.(((x.-c).*(x.-c))[:, 1] + ((x.-c).*(x.-c))[:, 2]))        # Havia sido calculado somente os dos bloco que foram divididos, ou seja, os nós.
        center_row[inode+i, :] = c   # Adicona o c das folhas na matrixc center_row
        diam[inode+i] = d           # Adicona diam das folhas na matrix diam
        child[inode+i, :] = [0 0]    # Completa a matrix child com pares [0,0], pois folhas nao tem filhos
    end
    Tree, child, center_row, diam
end

function Hinterp(dad, Tree1, Tree2, block; ninterp = 3, compressão = true, ϵ = 1e-3)
    # arg = [NOS1,NOS_GEO1,tipoCDC,valorCDC,normal,ELEM1,k]
    #         1      2        3      4        5     6   7
    n = size(block, 1)               # Quantidade de Submatrizes
    HA = Array{Matrix{Float64}}(undef, n, 2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    HB = Array{Matrix{Float64}}(undef, n, 2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    for i = 1:n                       # Para cada Submatriz
        # @timeit to "Para cada Submatriz" begin
        b1 = Tree1[block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # Submatriz = Produto cartesiano I x J
        if block[i, 3] == 0                # Se esses blocos não são admissiveis
            HA[i, 1], HB[i, 1] = calc_HeG(dad, b1, b2, 16)
            # Salva na linha i da 1º coluna da matriz HA a matriz H e salva G na matriz B
        else                              # Caso contrario (Se blocos são admissiveis)
            HA[i, 1], HA[i, 2], HB[i, 1] = calc_HeG_interp(dad, b1, b2, 8, ninterp)
        end
        # end

    end
    return HA, HB
end

function Ainterp(dad, Tree1, Tree2, block, ninterp = 3; compressão = true, ϵ = 1e-3)
    # arg = [NOS1,NOS_GEO1,tipoCDC,valorCDC,normal,ELEM1,k]
    #         1      2        3      4        5     6   7
    nc = size(dad.NOS, 1)
    ni = size(dad.pontos_internos, 1)
    b = zeros(nc + ni)
    n = size(block, 1)               # Quantidade de Submatrizes
    HA = Array{Matrix{Float64}}(undef, n, 2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    # HB = Array{Any}(undef,n,2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]

    for i = 1:n                       # Para cada Submatriz
        # @timeit to "Para cada Submatriz" begin
        b1 = Tree1[block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # Submatriz = Produto cartesiano I x J

        if typeof(dad) == potencial_iga
            CDCs = vcat([dad.valorCDC[dad.ELEM[ii].indices] for ii in b2]...)
            tipo = dad.tipoCDC[dad.ELEM[b2[1]].indices[1]]
        else
            CDCs = vcat([dad.ELEM[ii].valorCDC for ii in b2]...)
            tipo = dad.ELEM[b2[1]].tipoCDC
        end

        if block[i, 3] == 0                # Se esses blocos não são admissiveis
            if tipo == 1
                HA[i, 1], temp = calc_HeG(dad, b1, b2, 16)
                b[b1] += temp * CDCs
            else
                temp, HA[i, 1] = calc_HeG(dad, b1, b2, 16)
                b[b1] -= temp * CDCs
                HA[i, 1] = -HA[i, 1]

            end
            # Salva na linha i da 1º coluna da matriz HA a matriz H e salva G na matriz B
        else
            # Caso contrario (Se blocos são admissiveis)
            if tipo == 1
                HA[i, 1], HA[i, 2], temp = calc_HeG_interp(dad, b1, b2, 8, ninterp)
                b[b1] += HA[i, 1] * (temp * CDCs)
            else
                HA[i, 1], temp, HA[i, 2] = calc_HeG_interp(dad, b1, b2, 8, ninterp)

                b[b1] -= HA[i, 1] * (temp * CDCs)
                HA[i, 2] = -HA[i, 2]
            end


        end
        # end

    end
    return HA, b
end
function montaAcheia(hmat, block, Tree1, Tree2, cols, dad::potencial_iga; full = true)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    nc = size(dad.NOS, 1)
    ni = size(dad.pontos_internos, 1)
    if full
        A = zeros(nc + ni, nc + ni)
    else
        A = spzeros(nc + ni, nc + ni)
    end
    for i = nc+1:nc+ni
        A[i, i] = -1
    end
    for i = 1:length(block[:, 3])
        b1 = Tree1[block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        @show i, cols[block[i, 2]]
        if block[i, 3] == 0
            # conta=1
            # cols=dad.ELEM[ii].indices
            # @show norm(A[b1,cols]-H[b1,cols])
            A[b1, cols[i]] += hmat[i, 1]
            # conta+=size(cols,1)

        elseif full
            # @infiltrate
            # conta=1
            # for ii in b2
            # cols=dad.ELEM[ii].indices
            # @show norm(A[b1,cols]-H[b1,cols])
            A[b1, cols[i]] += hmat[i, 1] * hmat[i, 2]
            # conta+=size(cols,1)
            # end
            # cols=vcat([dad.ELEM[i].indices for i in b2]...)
            # A[b1,cols]+=(hmat[i,1]*hmat[i,2])
            # A[b1,cols]+=H[b1,cols]
            # @show norm(A[b1,cols]-H[b1,cols])

        end
    end
    A
end

function montaAcheia(hmat, block, Tree1, Tree2, dad::potencial; full = true)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    nc = size(dad.NOS, 1)
    ni = size(dad.pontos_internos, 1)
    if full
        A = zeros(nc + ni, nc + ni)
    else
        A = spzeros(nc + ni, nc + ni)
    end
    for i = nc+1:nc+ni
        A[i, i] = -1
    end
    for i = 1:length(block[:, 3])
        b1 = Tree1[block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # @show i,block[i,3]
        cols = vcat([dad.ELEM[i].indices for i in b2]...)
        if block[i, 3] == 0
            A[b1, cols] += hmat[i, 1]
        elseif full
            A[b1, cols] += (hmat[i, 1] * hmat[i, 2])
        end
    end
    A
end



"indc=nosproximoskmeans(X,k=9)
qr1=pqrfact(M[:,indc],rtol=1e-8)
qr2=pqrfact(:c,qr1.Q,rtol=1e-8)
indl=qr2.p[1:size(qr2.Q,1)]
A1=M[:,indc]
A2=M[indl,indc]divide M[indl,:]
M=A1*A2"
function nosproximoskmeans(X, k = 9)
    res = kmeans(X', k; tol = 1e-2, max_iters = 30)
    I = zeros(Int, k)
    for i = 1:k
        inds = findall(res.assignments .== i)

        I[i] = inds[argmin(colwise(SqEuclidean(), res.centers[:, i], X[inds, :]'))]
    end
    I
end

function Akmeans(dad, Tree1, Tree2, block, nnucleos = 6; compressão = true, ϵ = 1e-3)
    # arg = [NOS1,NOS_GEO1,tipoCDC,valorCDC,normal,ELEM1,k]
    #         1      2        3      4        5     6   7

    nc = size(dad.NOS, 1)
    ni = size(dad.pontos_internos, 1)
    b = zeros(nc + ni)
    n = size(block, 1)               # Quantidade de Submatrizes
    HA = Array{Matrix{Float64}}(undef, n, 2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    # HB = Array{Any}(undef,n,2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    indices = Dict()
    for i in unique(block[:, 2])
        # push!(indices,nosproximoskmeans([dad.NOS;dad.pontos_internos][Tree2[i],:],nnucleos))
        indices[i] = nosproximoskmeans(centros(dad)[Tree2[i], :], nnucleos)
    end



    for i = 1:n                       # Para cada Submatriz
        # @timeit to "Para cada Submatriz" begin
        b1 = Tree1[block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # Submatriz = Produto cartesiano I x J
        if typeof(dad) == potencial_iga
            CDCs = vcat([dad.valorCDC[dad.ELEM[ii].indices] for ii in b2]...)
            tipo = dad.tipoCDC[dad.ELEM[b2[1]].indices[1]]
        else
            CDCs = vcat([dad.ELEM[ii].valorCDC for ii in b2]...)
            tipo = dad.ELEM[b2[1]].tipoCDC
        end


        if block[i, 3] == 0                # Se esses blocos não são admissiveis
            if tipo == 1
                HA[i, 1], temp = calc_HeG(dad, b1, b2, 16)
                b[b1] += temp * CDCs
            else
                temp, HA[i, 1] = calc_HeG(dad, b1, b2, 16)
                b[b1] -= temp * CDCs
                HA[i, 1] = -HA[i, 1]

            end
            # Salva na linha i da 1º coluna da matriz HA a matriz H e salva G na matriz B

        else
            # @infiltrate             # Caso contrario (Se blocos são admissiveis)
            H1, G1 = calc_HeG(dad, b1, b2[indices[block[i, 2]]], 16)
            qrH1 = pqrfact(H1, rtol = 1e-8)
            qrH2 = pqrfact(:c, qrH1.Q, rtol = 1e-8)
            indl = qrH2.p[1:size(qrH2.Q, 1)]
            H2, G2 = calc_HeG(dad, b1[indl], b2, 16)
            H2 = H1[indl, :] \ H2
            G2 = G1[indl, :] \ G2

            if tipo == 1
                HA[i, 1], HA[i, 2] = H1, H2
                b[b1] += G1 * (G2 * CDCs)
            else
                HA[i, 1], HA[i, 2] = -G1, G2
                b[b1] -= H1 * (H2 * CDCs)
            end


        end
        # end

    end
    return HA, b
end


function Akmeans2(dad, Tree1, Tree2, block, nnucleos = 6; compressão = true, ϵ = 1e-3)
    # arg = [NOS1,NOS_GEO1,tipoCDC,valorCDC,normal,ELEM1,k]
    #         1      2        3      4        5     6   7
    nc = size(dad.NOS, 1)
    ni = size(dad.pontos_internos, 1)
    b = zeros(nc + ni)
    n = size(block, 1)               # Quantidade de Submatrizes
    HA = Array{Matrix{Float64}}(undef, n, 2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    # HB = Array{Any}(undef,n,2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]

    indices2 = Dict()
    indices1 = Dict()
    for i in unique(block[:, 2])
        # push!(indices,nosproximoskmeans([dad.NOS;dad.pontos_internos][Tree2[i],:],nnucleos))
        indices2[i] = nosproximoskmeans(centros(dad)[Tree2[i], :], nnucleos)
    end
    for i in unique(block[:, 1])
        # push!(indices,nosproximoskmeans([dad.NOS;dad.pontos_internos][Tree2[i],:],nnucleos))
        indices1[i] =
            nosproximoskmeans([dad.NOS; dad.pontos_internos][Tree1[i], :], nnucleos)
    end


    for i = 1:n                       # Para cada Submatriz
        # @timeit to "Para cada Submatriz" begin
        b1 = Tree1[block[i, 1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i, 2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # Submatriz = Produto cartesiano I x J

        if typeof(dad) == potencial_iga
            CDCs = vcat([dad.valorCDC[dad.ELEM[ii].indices] for ii in b2]...)
            tipo = dad.tipoCDC[dad.ELEM[b2[1]].indices[1]]
        else
            CDCs = vcat([dad.ELEM[ii].valorCDC for ii in b2]...)
            tipo = dad.ELEM[b2[1]].tipoCDC
        end

        if block[i, 3] == 0                # Se esses blocos não são admissiveis
            if tipo == 1
                HA[i, 1], temp = calc_HeG(dad, b1, b2, 16)
                b[b1] += temp * CDCs
            else
                temp, HA[i, 1] = calc_HeG(dad, b1, b2, 16)
                b[b1] -= temp * CDCs
                HA[i, 1] = -HA[i, 1]

            end
            # Salva na linha i da 1º coluna da matriz HA a matriz H e salva G na matriz B
        else
            # @infiltrate             # Caso contrario (Se blocos são admissiveis)
            H1, G1 = calc_HeG(dad, b1, b2[indices2[block[i, 2]]], 16)
            H2, G2 = calc_HeG(dad, b1[indices1[block[i, 1]]], b2, 16)
            H2 = H1[indices1[block[i, 1]], :] \ H2
            G2 = G1[indices1[block[i, 1]], :] \ G2
            if tipo == 1
                HA[i, 1], HA[i, 2] = H1, H2
                b[b1] += G1 * (G2 * CDCs)
            else
                HA[i, 1], HA[i, 2] = -G1, G2
                b[b1] -= H1 * (H2 * CDCs)
            end


        end
        # end

    end
    return HA, b
end

function MatrizH(dad, f, p = 9)
    Tree1, Tree2, block = cluster(dad, η = 1, max_elem = 30)
    cols = colsblocks(block, Tree2, dad)
    Ai, bi = f(dad, Tree1, Tree2, block, p)
    hmat(Ai, block, Tree1, Tree2, dad, cols, size(bi, 1)), bi
end

size(A1::hmat) = (A1.n, A1.n)
