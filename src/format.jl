function calcula_arco(x1, y1, x2, y2, xc, yc, raio)
    # Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
    # horizontal direction 
    #  and 
    # the tet2 angle between the line from point (x2,y2) to (xc,yc)
    # and the horizontal direction
    dx1 = x1 - xc
    dy1 = y1 - yc
    dx2 = x2 - xc
    dy2 = y2 - yc

    if abs(dx1) < 1e-15         # <------- EmerEdited
        dx1 = 0
    elseif abs(dy1) < 1e-15
        dy1 = 0
    end
    # @show dy1,dx1           # <------- EmerEdited

    tet1 = atan(dy1, dx1)
    a = [dx1, dy1, 0]
    b = [dx2, dy2, 0]
    c = [0, 0, (dx1 * dy2 - dx2 * dy1)]

    # @show sqrt(sum(c.*c)), sum(a.*b)              # <------- EmerEdited
    aa = sqrt(sum(c .* c))
    bb = sum(a .* b)              # <------- EmerEdited
    if abs(aa) < 1e-15                  # <------- EmerEdited
        aa = 0
    elseif abs(bb) < 1e-15                 # <------- EmerEdited
        bb = 0
    end

    angle = atan(aa, bb)

    # angle = atan(norm(cross(a,b)),dot(a,b));  # <------- EmerEdited
    # @show tet1,angle           # <------- EmerEdited
    if (raio > 0)
        tet2 = tet1 + angle
        if abs(tet2) < 1e-15      # <------- EmerEdited
            tet2 = 0
        end
    else
        tet2 = tet1 - angle
        if abs(tet2) < 1e-15      # <------- EmerEdited
            tet2 = 0
        end
    end
    [tet1, tet2]
end

function calcula_centro(x1, y1, x2, y2, raio)
    # Compute the center of an arc given two points and the radius

    xm = (x1 + x2) / 2
    ym = (y1 + y2) / 2
    b = √((x2 - x1)^2 + (y2 - y1)^2)
    t1 = (x2 - x1) / b
    t2 = (y2 - y1) / b
    n1 = t2
    n2 = -t1
    h = √(abs(raio^2 - (b / 2)^2))
    if (raio > 0)
        if (n1 == 0)
            xc = xm
            yc = ym - n2 / abs(n2) * h
        else
            xc = -n1 / abs(n1) * √(h^2 * n1^2 / (n1^2 + n2^2)) + xm
            yc = n2 / n1 * (xc - xm) + ym
        end
    else
        if (n1 == 0)
            xc = xm
            yc = ym + n2 / abs(n2) * h
        else
            xc = n1 / abs(n1) * √(h^2 * n1^2 / (n1^2 + n2^2)) + xm
            yc = n2 / n1 * (xc - xm) + ym
        end
    end
    [xc, yc]
end
# function format_dad(entrada, NPX=2, NPY=2, afasta=1)

function format_dad(entrada, NPX=2, NPY=2, afasta=1; canto=false)
    # Programa para formatação dos dados de entrada
    # @infiltrate
    if length(entrada) == 6
        prob, PONTOS, SEGMENTOS, MALHA, CCSeg, k = entrada
        subregioes = 0
    elseif length(entrada) == 7
        prob, PONTOS, SEGMENTOS, MALHA, CCSeg, k, subregioes = entrada
        k = (; k..., subregioes)
    else
        println("Erro")
    end
    CORNERS = [PONTOS zeros(size(PONTOS, 1), 6)]
    num_elementos = sum(MALHA[:, 2])
    if prob == potencial
        ELEM = Vector{elemento}(undef, num_elementos)
    elseif prob == helmholtz
        ELEM = Vector{elemento}(undef, num_elementos)
    else
        # prob == elastico
        ELEM = Vector{elementov}(undef, num_elementos)
    end
    NOS = zeros(Float64, 0, 2)
    cont_nos = Int64(0)  # Counter to the physical nodes
    cont_el = Int64(0)# Counter to the elements (the number of physical and geometric elements is the same).
    num_lin = length(SEGMENTOS[:, 1])# Número de linhas/lados no contorno
    p_ini = round(Int64, SEGMENTOS[1, 2])
    pontosg, pesosg = gausslegendre(2)    # Quadratura de gauss

    # ______________________________________________________________________
    # Definition of the biggest dimension of the problem
    max_dl = 0
    for lin = 1:num_lin
        p1 = round(Int64, SEGMENTOS[lin, 2])
        p2 = round(Int64, SEGMENTOS[lin, 3])
        xp1 = PONTOS[p1, 2]
        yp1 = PONTOS[p1, 3]
        xp2 = PONTOS[p2, 2]
        yp2 = PONTOS[p2, 3]
        dl = √((xp1 - xp2)^2 + (yp1 - yp2)^2)
        if dl > max_dl
            max_dl = dl
        end
    end
    # _____________________________________________________________________

    no_ini = 1
    t = 1
    p2 = 0
    no1_prox = 0
    while (t <= num_lin)  # While over all lines
        tipo_elem = MALHA[t, 3]
        # while (p2 != p_ini)
        num_el_lin = MALHA[t, 2]# Number of the elements in the line t
        # Coordinates of the initial and final PONTOS of each line
        p1 = round(Int64, SEGMENTOS[t, 2])
        p2 = round(Int64, SEGMENTOS[t, 3])
        x1l = PONTOS[p1, 2]                 # [x1l,y1l,x2l,y2l]
        y1l = PONTOS[p1, 3]
        x2l = PONTOS[p2, 2]
        y2l = PONTOS[p2, 3]
        if (SEGMENTOS[t, 4] == 0)     # The segment is a straight line
            delta_x = x2l - x1l      # Increment in x and y direction
            delta_y = y2l - y1l

            length_ele = sqrt(delta_x^2 + delta_y^2)
            nor1l = delta_y / length_ele
            nor2l = -delta_x / length_ele

            CORNERS[p1, 6] = nor1l
            CORNERS[p1, 7] = nor2l
            CORNERS[p2, 4] = nor1l
            CORNERS[p2, 5] = nor2l
        else                       # The segment is an arc
            r = SEGMENTOS[t, 4]
            xc, yc = calcula_centro(x1l, y1l, x2l, y2l, r)    # Compute the center of the arc and its coordinates

            r1 = √((x1l - xc)^2 + (y1l - yc)^2) # Distance between p1 and center (r1)
            r2 = √((x2l - xc)^2 + (y2l - yc)^2) # and between p2 and center (r2)

            if abs(r1 - r2) < 0.00001 * max_dl
                # Compute the angle between the lines from point c to p1 [tet1) and c to p2 (tet2]
                tet1, tet2 = calcula_arco(x1l, y1l, x2l, y2l, xc, yc, r)
                if tet2 < tet1
                    tet2 = tet2 + 2 * π
                end

                # Angle of the sector defined by the arc
                if SEGMENTOS[t, 4] > 0
                    tet = abs(tet2 - tet1)
                    sig = 1
                else
                    tet = 2 * π - abs(tet2 - tet1)
                    sig = -1
                end

                # Angle between two nodes of the line
                divtet = tet / (num_el_lin)
            else
                error("Error in the data input file: Wrong central point")
            end
            CORNERS[p1, 6] = x1l - xc
            CORNERS[p1, 7] = y1l - yc
            CORNERS[p2, 4] = x1l - xc
            CORNERS[p2, 5] = y1l - yc
        end
        # Generation of elements and nodes
        CORNERS[p1, 9] = cont_el + 1

        qsi = range(0, stop=1, length=tipo_elem)
        for i = 1:num_el_lin
            if (SEGMENTOS[t, 4] == 0) # The segment is a straight line
                x_i = x1l + delta_x / num_el_lin * (i - 1)# initial x coordinate of the element
                y_i = y1l + delta_y / num_el_lin * (i - 1)# initial y coordinate of the element
                x_f = x1l + delta_x / num_el_lin * (i)     # final x coordinate of the element
                y_f = y1l + delta_y / num_el_lin * (i)     # final y coordinate of the element
                xs = x_i .+ (x_f - x_i) * qsi              # nos no elemento
                ys = y_i .+ (y_f - y_i) * qsi
            else  # The segment is an arc
                # Compute the node coordinates
                xs = xc .+ r1 * cos.(tet1 .+ (i - 1 .+ qsi) * sig * divtet)
                ys = yc .+ r1 * sin.(tet1 .+ (i - 1 .+ qsi) * sig * divtet)
            end
            NOS = [NOS; [xs ys]]

            cont_el = cont_el + 1
            nos = (cont_el-1)*tipo_elem+1:cont_el*tipo_elem
            # qsis = range(-1 + afasta / tipo_elem, stop=1 - afasta / tipo_elem, length=tipo_elem) # Parametrização de -1 a 1
            qsis, ~ = gausslegendre(tipo_elem) # Parametrização de -1 a 1
            # @infiltrate
            N1, dN1 = calc_fforma_gen(pontosg[1], qsis)
            N2, dN2 = calc_fforma_gen(pontosg[2], qsis)
            tamanho = norm(dN1' * [xs ys]) * pesosg[1] + norm(dN2' * [xs ys]) * pesosg[2]
            # tamanho=norm(N1'*[xs ys]+N2'*[xs ys])
            # @infiltrate
            # @infiltrate
            if prob == potencial || prob == helmholtz
                if subregioes == 0
                    ELEM[cont_el] = elemento(nos, CCSeg[t, 2], fill(CCSeg[t, 3], tipo_elem), qsis, tamanho, 0)
                else
                    ELEM[cont_el] = elemento(nos, CCSeg[t, 2], fill(CCSeg[t, 3], tipo_elem), qsis, tamanho, subregioes.regiao[cont_el])
                end

            elseif prob == elastico || prob == elastico_aniso || prob == placa_fina || prob == placa_fina_isotropica
                if subregioes == 0
                    ELEM[cont_el] = elementov(nos, CCSeg[t, [2, 4]], repeat(CCSeg[t, [3, 5]], 1, tipo_elem), qsis, tamanho, 0)
                else
                    ELEM[cont_el] = elementov(nos, CCSeg[t, [2, 4]], repeat(CCSeg[t, [3, 5]], 1, tipo_elem), qsis, tamanho, subregioes.regiao[cont_el])
                end
            else
                if subregioes == 0
                    ELEM[cont_el] = elementov(nos, CCSeg[t, [2, 4, 6]], repeat(CCSeg[t, [3, 5, 7]], 1, tipo_elem), qsis, tamanho, 0)
                else
                    ELEM[cont_el] = elementov(nos, CCSeg[t, [2, 4, 6]], repeat(CCSeg[t, [3, 5, 7]], 1, tipo_elem), qsis, tamanho, subregioes.regiao[cont_el])
                end
            end
        end
        CORNERS[p2, 8] = cont_el
        t = t + 1
        # end
    end
    if canto
        bc_canto = zeros(size(CORNERS, 1), 2)
        if CCSeg[end, 1] == 'f'
            bc_canto[:, 1] .= 1
        end
        # @infiltrate
    end
    #   vizinhos=false
    #   if vizinhos
    if afasta == 0
        NOSn = unique(NOS, dims=1)             # elimina pontos iguais
        ELEM1 = vcat([i.indices for i in ELEM]'...)
        ELEMn = zeros(size(ELEM1))

        for i = 1:size(NOSn, 1)
            equalind = findall(sum(isequal.(NOSn[i, :]', NOS), dims=2) .== 2)[:]   # true +true = 2
            for ind in equalind
                ELEMn[ELEM1.==ind[1]] .= i
            end
        end

        for i = 1:size(ELEM1, 1)
            if prob == potencial
                # @infiltrate
                ELEM[i] = elemento(ELEMn[i, :], ELEM[i].tipoCDC, ELEM[i].valorCDC, ELEM[i].ξs, 0)
            else
                prob == elastico
                ELEM[i] = elementov(ELEMn[i, :], ELEM[i].tipoCDC, ELEM[i].valorCDC, ELEM[i].ξs, 0)
            end
        end
        return prob(NOSn, gera_p_in(NPX, NPY, PONTOS, SEGMENTOS), ELEM, k)
    end
    # vizinho= zeros(Int64,size(ELEM,1),2) #matriz para armazenar o elemento inicial e o elemento final nos quais estao sendo calculado o salto


    #   for el=1:size(ELEM,1)
    #     # @infiltrate
    #     vizinho[el,:] = [ELEM[ELEMn[el,1].==ELEMn[:,end],1] ELEM[ELEMn[el,end].==ELEMn[:,1],1]]; #busca o no adjacente ao no el e monta a matriz SALTOSI[el no_adjacente]
    #   end
    # end

    for i = 1:num_elementos
        tipo = size(ELEM[i].indices, 1)
        ξs = range(-1, stop=1, length=tipo)

        N_geo = calc_fforma_gen.(ELEM[i].ξs, Ref(ξs)) # funções de forma generalizada
        xn = zeros(tipo, 2)
        for k = 1:tipo
            x = NOS[ELEM[i].indices, :] # coordenadas dos nos geometricos que compoem o elemento
            # @show x,N_geo
            # @infiltrate
            xn[k, :] = N_geo[k][1]' * x # coordenadas geometricas do no em qsi = 
        end
        NOS[ELEM[i].indices, :] = xn
    end
    # _____
    #   if vizinhos
    #     return dad,vizinho
    #     else
    #  return  ELEM
    #     end
    # @infiltrate

    if canto
        prob(NOS, gera_p_in(NPX, NPY, PONTOS, SEGMENTOS), ELEM, (; k..., cantos=CORNERS, bc_canto=bc_canto))
        # prob(NOS, gera_p_in(NPX, NPY, PONTOS, SEGMENTOS), ELEM, merge(k,(cantos=CORNERS),))
    else
        # @infiltrate
        prob(NOS, gera_p_in(NPX, NPY, PONTOS, SEGMENTOS), ELEM, k)
    end
end


function aplicaCDC(H, G, dad::Union{elastico,elastico_aniso})
    # nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    # n = size(dad.NOS, 1)
    A = deepcopy(H)
    b = zeros(size(H, 1))
    if typeof(dad) == elastico_aniso
        scale = dad.k.A3[3, 3]
    else
        scale = 1
    end
    for elem_i in dad.ELEM, i in 1:2  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[i] == 0
            A[:, 2ind_elem.+(i-2)] = -G[:, 2ind_elem.+(i-2)] * scale
            b += -H[:, 2ind_elem.+(i-2)] * elem_i.valorCDC[i, :]
        elseif elem_i.tipoCDC[i] == 1
            # A[:, 2ind_elem.+(i-2)] = H[:, 2ind_elem.+(i-2)]
            b += G[:, 2ind_elem.+(i-2)] * elem_i.valorCDC[i, :]
        end
    end
    if temsubregioes(dad)
        neq = size(dad.k.subregioes.equivale, 1)
        indeq = [2 * dad.k.subregioes.equivale[:, 1] .- 1 dad.k.subregioes.equivale[:, 1] * 2 2 * dad.k.subregioes.equivale[:, 2] .- 1 dad.k.subregioes.equivale[:, 2] * 2]'[:]
        @infiltrate
        A = [A -G[:, indeq]; dad.k.subregioes.Hc]
        b = [b; zeros(4neq)]
    end

    A, b
end
function aplicaCDC(H, G, dad::Union{elastico,elastico_aniso}, nosrestritos)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    A = zeros(2 * n, 2 * n)
    b = zeros(2 * n)
    if typeof(dad) == elastico_aniso
        scale = dad.k.A3[3, 3]
    else
        scale = 1
    end
    for elem_i in dad.ELEM, i in 1:2  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[i] == 0
            A[:, 2ind_elem.+(i-2)] = -G[:, 2ind_elem.+(i-2)] * scale
            b += -H[:, 2ind_elem.+(i-2)] * elem_i.valorCDC[i, :]
        elseif elem_i.tipoCDC[i] == 1
            A[:, 2ind_elem.+(i-2)] = H[:, 2ind_elem.+(i-2)]
            b += G[:, 2ind_elem.+(i-2)] * elem_i.valorCDC[i, :]
        end
    end
    for i = 1:size(nosrestritos, 1)
        no = nosrestritos[i, 1]
        tipo = nosrestritos[i, 2]
        # @infiltrate
        A[:, (no-1)*2+tipo] = -G[:, (no-1)*2+tipo]

    end
    A, b
end
function aplicaCDC(H::Matrix{Float64}, G::Matrix{Float64}, dad::Union{potencial,helmholtz})
    # nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    A = zeros(typeof(H[1]), n, n)
    b = zeros(typeof(H[1]), n)
    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        # @infiltrate
        if elem_i.tipoCDC == 0
            A[:, ind_elem] = -G[:, ind_elem]
            b += -H[:, ind_elem] * elem_i.valorCDC
        elseif elem_i.tipoCDC == 1
            A[:, ind_elem] = H[:, ind_elem]
            b += G[:, ind_elem] * elem_i.valorCDC
        end
    end
    A, b
end

function separa(dad::Union{elastico,elastico_aniso}, x)
    # Separa fluxo e temperatura

    # ncdc = número de linhas da matriz CDC
    # T = vetor que contêm as temperaturas nos nós
    # q = vetor que contêm o fluxo nos nós
    n = size(dad.NOS, 1)
    T = zeros(2n)
    q = zeros(2n)
    if typeof(dad) == elastico_aniso
        scale = dad.k.A3[3, 3]
    else
        scale = 1
    end

    for elem_i in dad.ELEM, i in 1:2   # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[i] == 0
            T[2ind_elem.+(i-2)] = elem_i.valorCDC[i, :] # A temperatura é a condição de contorno
            q[2ind_elem.+(i-2)] = x[2ind_elem.+(i-2)] * scale # O fluxo é o valor calculado
        elseif elem_i.tipoCDC[i] == 1
            T[2ind_elem.+(i-2)] = x[2ind_elem.+(i-2)] # 
            q[2ind_elem.+(i-2)] = elem_i.valorCDC[i, :]
        end
    end
    if temsubregioes(dad)
        neq = size(dad.k.subregioes.equivale, 1)
        indeq = [2 * dad.k.subregioes.equivale[:, 1] .- 1 dad.k.subregioes.equivale[:, 1] * 2 2 * dad.k.subregioes.equivale[:, 2] .- 1 dad.k.subregioes.equivale[:, 2] * 2]'[:]
        q[indeq] = x[2n+1:end]
        return [T[1:2:end] T[2:2:end]], [q[1:2:end] q[2:2:end]]

    end
    if length(x) > 2n
        uint = [x[2n+1:2:end] x[2n+2:2:end]]
        return [T[1:2:end] T[2:2:end]], [q[1:2:end] q[2:2:end]], uint
    else
        return [T[1:2:end] T[2:2:end]], [q[1:2:end] q[2:2:end]]
    end
end

function separa(dad::Union{elastico,elastico_aniso}, x, nosrestritos)
    # Separa fluxo e temperatura

    # ncdc = número de linhas da matriz CDC
    # T = vetor que contêm as temperaturas nos nós
    # q = vetor que contêm o fluxo nos nós
    n = size(dad.NOS, 1)
    T = zeros(2n)
    q = zeros(2n)
    if typeof(dad) == elastico_aniso
        scale = dad.k.A3[3, 3]
    else
        scale = 1
    end
    for elem_i in dad.ELEM, i in 1:2   # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[i] == 0
            T[2ind_elem.+(i-2)] = elem_i.valorCDC[i, :] # A temperatura é a condição de contorno
            q[2ind_elem.+(i-2)] = x[2ind_elem.+(i-2)] * scale # O fluxo é o valor calculado
        elseif elem_i.tipoCDC[i] == 1
            T[2ind_elem.+(i-2)] = x[2ind_elem.+(i-2)] # 
            q[2ind_elem.+(i-2)] = elem_i.valorCDC[i, :]
        end
    end
    for i = 1:size(nosrestritos, 1)
        no = nosrestritos[i, 1]
        tipo = nosrestritos[i, 2]
        ind = (no - 1) * 2 + tipo
        T[ind] = 0 # A temperatura é a condição de contorno
        q[ind] = x[ind] * scale # O fluxo é o valor calculado
    end
    [T[1:2:end] T[2:2:end]], [q[1:2:end] q[2:2:end]]
end

function separa(dad::Union{potencial,helmholtz}, x)
    # Separa fluxo e temperatura

    # ncdc = número de linhas da matriz CDC
    # T = vetor que contêm as temperaturas nos nós
    # q = vetor que contêm o fluxo nos nós
    n = size(dad.NOS, 1)
    T = zeros(typeof(x[1]), n)
    q = zeros(typeof(x[1]), n)
    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC == 0
            T[ind_elem] = elem_i.valorCDC # A temperatura é a condição de contorno
            q[ind_elem] = x[ind_elem] # O fluxo é o valor calculado
        elseif elem_i.tipoCDC == 1
            T[ind_elem] = x[ind_elem] # 
            q[ind_elem] = elem_i.valorCDC
        end
    end
    T, q
end


function inpoly(p, node, edge, reltol=1.0e-12)

    #   INPOLY: Point-in-polygon testing.
    # 
    #  Determine whether a series of points lie within the bounds of a polygon
    #  in the 2D plane. General non-convex, multiply-connected polygonal
    #  regions can be handled.
    # 
    #  SHORT SYNTAX:
    # 
    #    in = inpoly(p, node);
    # 
    #    p   : The points to be tested as an Nx2 array [x1 y1; x2 y2; etc].
    #    node: The vertices of the polygon as an Mx2 array [X1 Y1; X2 Y2; etc].
    #          The standard syntax assumes that the vertices are specified in
    #          consecutive order.
    # 
    #    in  : An Nx1 logical array with IN(i) = TRUE if P(i,:) lies within the
    #          region.
    # 
    #  LONG SYNTAX:
    # 
    #   [in, on] = inpoly(p, node, edge, tol);
    # 
    #   edge: An Mx2 array of polygon edges, specified as connections between
    #         the vertices in NODE: [n1 n2; n3 n4; etc]. The vertices in NODE
    #         do not need to be specified in connsecutive order when using the
    #         extended syntax.
    # 
    #   on  : An Nx1 logical array with ON(i) = TRUE if P(i,:) lies on a
    #         polygon edge. (A tolerance is used to deal with numerical
    #         precision, so that points within a distance of
    #         reltol*min(bbox(node)) from a polygon edge are considered "on" the 
    #         edge.
    # 
    #  EXAMPLE:
    # 
    #    polydemo;       #  Will run a few examples
    # 
    #  See also INPOLYGON

    #  The algorithm is based on the crossing number test, which counts the
    #  number of times a line that extends from each point past the right-most
    #  region of the polygon intersects with a polygon edge. Points with odd
    #  counts are inside. A simple implementation of this method requires each
    #  wall intersection be checked for each point, resulting in an O(N*M)
    #  operation count.
    # 
    #  This implementation does better in 2 ways:
    # 
    #    1. The test points are sorted by y-value and a binary search is used to
    #       find the first point in the list that has a chance of intersecting
    #       with a given wall. The sorted list is also used to determine when we
    #       have reached the last point in the list that has a chance of
    #       intersection. This means that in general only a small portion of
    #       points are checked for each wall, rather than the whole set.
    # 
    #    2. The intersection test is simplified by first checking against the
    #       bounding box for a given wall segment. Checking against the bbox is
    #       an inexpensive alternative to the full intersection test and allows
    #       us to take a number of shortcuts, minimising the number of times the
    #       full test needs to be done.
    # 
    #    Darren Engwirda: 2005-2009
    #    Email          : d_engwirda@hotmail.com
    #    Last updated   : 28/03/2009 with MATLAB 7.0
    # 
    #  Problems or suggestions? Email me.
    nnode = size(node, 1)

    # #  PRE-PROCESSING
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    n = size(p, 1)
    nc = size(edge, 1)

    #  Choose the direction with the biggest range as the "y-coordinate" for the
    #  test. This should ensure that the sorting is done along the best
    #  direction for long and skinny problems wrt either the x or y axes.
    dxy = maximum(p, 1) - minimum(p, 1)
    if dxy[1] > dxy[2]
        #  Flip co-ords if x range is bigger
        p = p[:, [2, 1]]
        node = node[:, [2, 1]]
    end

    #  Polygon bounding-box
    dxy = maximum(node, 1) - minimum(node, 1)
    tol = reltol * minimum(dxy)
    if tol == 0.0
        tol = reltol
    end

    #  Sort test points by y-value
    i = sortperm(p[:, 2])
    y = p[i, 2]
    x = p[i, 1]
    # #  MAIN LOOP
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    cn = falses(n)     #  Because we're dealing with mod(cn,2) we don't have
    #  to actually increment the crossing number, we can
    #  just flip a logical at each intersection (faster!)
    on = cn[:]
    for k = 1:nc         #  Loop through edges
        #  Nodes in current edge
        n1 = edge[k, 1]
        n2 = edge[k, 2]

        #  Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
        #            - also get xmin = min(x1,x2), xmax = max(x1,x2)
        y1 = node[n1, 2]
        y2 = node[n2, 2]
        if y1 < y2
            x1 = node[n1, 1]
            x2 = node[n2, 1]
        else
            yt = y1
            y1 = y2
            y2 = yt
            x1 = node[n2, 1]
            x2 = node[n1, 1]
        end
        if x1 > x2
            xmin = x2
            xmax = x1
        else
            xmin = x1
            xmax = x2
        end
        #  Binary search to find first point with y<=y1 for current edge
        if y[1] >= y1
            start = 1
        elseif y[n] < y1
            start = n + 1
        else
            lower = 1
            upper = n
            for j = 1:n
                start = convert(Int32, round(1 / 2 * (lower + upper)))
                if y[start] < y1
                    lower = start + 1
                elseif y[start-1] < y1
                    break
                else
                    upper = start - 1
                end
            end
        end
        #  Loop through points
        for j = start:n
            #  Check the bounding-box for the edge before doing the intersection
            #  test. Take shortcuts wherever possible!
            Y = y[j]   #  Do the array look-up once & make a temp scalar
            if Y <= y2
                X = x[j]   #  Do the array look-up once & make a temp scalar
                if X >= xmin
                    if X <= xmax

                        #  Check if we're "on" the edge
                        on[j] = on[j] || (abs((y2 - Y) * (x1 - X) - (y1 - Y) * (x2 - X)) <= tol)

                        #  Do the actual intersection test
                        if (Y < y2) && ((y2 - y1) * (X - x1) < (Y - y1) * (x2 - x1))
                            cn[j] = ~cn[j]
                        end

                    end
                elseif Y < y2   #  Deal with points exactly at vertices
                    #  Has to cross edge
                    cn[j] = ~cn[j]
                end
            else
                #  Due to the sorting, no points with >y
                #  value need to be checked
                break
            end
        end
    end
    #  Re-index to undo the sorting
    cn[i] = cn | on
    on[i] = on

    return cn

end      #  inpoly()






function aplicaCDC(H, G, dad::potencial_iga)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    A = deepcopy(H)
    b = zeros(n)
    valoresconhecidos = dad.E \ dad.valorCDC
    troca = dad.tipoCDC[:] .== 0
    # @infiltrate
    A[:, troca] = -G[:, troca]
    b += -H[:, troca] * valoresconhecidos[troca]
    b += G[:, .~troca] * valoresconhecidos[.~troca]

    A, b
end



function separa(dad::potencial_iga, x)
    # Separa fluxo e temperatura

    troca = dad.tipoCDC[:] .== 0
    # T = vetor que contêm as temperaturas nos nós
    # q = vetor que contêm o fluxo nos nós
    T = deepcopy(x[1:length(troca)])
    q = deepcopy(dad.valorCDC)
    #  @infiltrate
    T[troca] = dad.valorCDC[troca]
    q[troca] = x[1:length(troca)][troca]
    dad.E * T, dad.E * q
end

function aplicaCDC(H, G, dad::Union{elastico_iga,elastico_aniso_iga})
    # nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    A = deepcopy(H)
    b = zeros(2 * n)
    valoresconhecidos = dad.E \ dad.valorCDC
    troca = dad.tipoCDC'[:] .== 0
    # @infiltrate
    A[:, troca] = -G[:, troca]
    b += -H[:, troca] * valoresconhecidos'[troca]
    b += G[:, .~troca] * valoresconhecidos'[.~troca]

    A, b
end

function separa(dad::Union{elastico_iga,elastico_aniso_iga}, x)
    # nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno

    troca = dad.tipoCDC'[:] .== 0
    u = deepcopy(x[1:length(troca)])
    t = deepcopy(dad.valorCDC'[:])
    #  @infiltrate
    u[troca] = dad.valorCDC[troca]
    t[troca] = x[1:length(troca)][troca]
    [dad.E * u[1:2:end] dad.E * u[2:2:end]], [dad.E * t[1:2:end] dad.E * t[2:2:end]]
end

function calc_normais(dad; tangente=false)
    normais = deepcopy(dad.NOS)
    tangentes = deepcopy(dad.NOS)
    for elem_j in dad.ELEM#Laço dos elementos
        x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
        for i in 1:length(elem_j.indices) #
            ξ = elem_j.ξs[i]
            N_geo, dN = calc_fforma(ξ, elem_j)
            dxdqsi = dN' * x
            dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
            sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
            sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
            normais[elem_j.indices[i], :] = [sy, -sx]
            tangentes[elem_j.indices[i], :] = [sx, sy]
        end
    end
    if tangente
        return normais, tangentes
    end
    normais
end
function matrizrotação(dad)
    normais, tangentes = calc_normais(dad, tangente=true)
    R = spzeros(2nc(dad), 2nc(dad))
    for i = 1:nc(dad)
        R[2i-1:2i, 2i-1:2i] = [normais[i, :]; tangentes[i, :]]
    end
    R'
end
function aplicaCDC(H, G, dad::placa_fina)
    # nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    A = deepcopy(H)
    B = deepcopy(G)
    bcval = zeros(size(H, 2))
    ncanto = size(dad.k.bc_canto, 1)
    scale = dad.k.D22

    for elem_i in dad.ELEM, i in 1:2  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        # @show elem_i.tipoCDC[i]
        if elem_i.tipoCDC[i] == 0
            A[:, 2ind_elem.+(i-2)] = -G[:, 2ind_elem.+(i-2)] * scale
            B[:, 2ind_elem.+(i-2)] = -H[:, 2ind_elem.+(i-2)]
        end
        bcval[2ind_elem.+(i-2)] = elem_i.valorCDC[i, :]
    end
    for c = 1:ncanto
        if dad.k.bc_canto[c, 1] == 0
            A[:, end-ncanto+c] = -G[:, end-ncanto+c] * scale
            B[:, end-ncanto+c] = -H[:, end-ncanto+c]
        end
    end

    A, B, bcval
end


function aplicaCDC(H, G, dad::placa_fina_isotropica)
    # nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    A = deepcopy(H)
    B = deepcopy(G)
    bcval = zeros(size(H, 2))
    ncanto = size(dad.k.bc_canto, 1)
    scale = dad.k.D

    for elem_i in dad.ELEM, i in 1:2  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        # @show elem_i.tipoCDC[i]
        if elem_i.tipoCDC[i] == 0
            A[:, 2ind_elem.+(i-2)] = -G[:, 2ind_elem.+(i-2)] * scale
            B[:, 2ind_elem.+(i-2)] = -H[:, 2ind_elem.+(i-2)]
        end
        bcval[2ind_elem.+(i-2)] = elem_i.valorCDC[i, :]
    end
    for c = 1:ncanto
        if dad.k.bc_canto[c, 1] == 0
            A[:, end-ncanto+c] = -G[:, end-ncanto+c] * scale
            B[:, end-ncanto+c] = -H[:, end-ncanto+c]
        end
    end

    A, B, bcval
end


function aplicaCDC(H, G, dad::placa_espessa_isotropica)
    # nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    A = deepcopy(H)
    B = deepcopy(G)
    bcval = zeros(size(H, 2))
    scale = dad.k.D

    for elem_i in dad.ELEM, i in 1:3  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        # @show elem_i.tipoCDC[i]
        if elem_i.tipoCDC[i] == 0
            A[:, 3ind_elem.+(i-3)] = -G[:, 3ind_elem.+(i-3)] * scale
            B[:, 3ind_elem.+(i-3)] = -H[:, 3ind_elem.+(i-3)]
        end
        bcval[3ind_elem.+(i-3)] = elem_i.valorCDC[i, :]
    end

    A, B, bcval
end

function tipoCDC(dad)
    tipos = fill(true, ni(dad) + nc(dad))
    for i in dad.ELEM
        if i.tipoCDC == 0
            tipos[i.indices] .= false
        end
    end
    tipos
end
function valorCDC(dad)
    valores = fill(0.0,nc(dad))
    for i in dad.ELEM
        valores[i.indices] = i.valorCDC
    end
    valores
end
function valoresCDC(dad)
    valoresu = fill(0.0,ni(dad) +nc(dad))
    valoresq = fill(0.0,ni(dad) +nc(dad))
    for i in dad.ELEM
        if i.tipoCDC == 0
            valoresu[i.indices] = i.valorCDC
        else
            valoresq[i.indices] = i.valorCDC
        end
    end
    valoresu,valoresq
end