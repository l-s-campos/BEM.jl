#==========POTENCIAL==========#

function deriva_Ti_potencial(dad, Ti, dx, dy)

    n = length(Ti)

    ui = zeros(n)
    vi = zeros(n)

    global iterando = 0

    for i = 1:n

        global iterando
        iterando = i

        if dad.pontos_internos[i, 2] == minimum(dad.pontos_internos[:, 2])
            vi[i] = (Ti[indice_forward(dad, i, 'y', x_order, y_order)] - Ti[i]) / dy

        elseif dad.pontos_internos[i, 2] == maximum(dad.pontos_internos[:, 2])
            vi[i] = (Ti[i] - Ti[indice_back(dad, i, 'y', x_order, y_order)]) / dy

        else
            if indice_back(dad, i, 'y', x_order, y_order) !== nothing &&
               indice_forward(dad, i, 'y', x_order, y_order) !== nothing
                vi[i] =
                    (
                        Ti[indice_forward(dad, i, 'y', x_order, y_order)] -
                        Ti[indice_back(dad, i, 'y', x_order, y_order)]
                    ) / (2 * dy)

            elseif indice_forward(dad, i, 'y', x_order, y_order) !== nothing
                vi[i] = (Ti[indice_forward(dad, i, 'y', x_order, y_order)] - Ti[i]) / dy

            else
                vi[i] = (Ti[i] - Ti[indice_back(dad, i, 'y', x_order, y_order)]) / dy
            end
        end



        if dad.pontos_internos[i, 1] == minimum(dad.pontos_internos[:, 1])
            ui[i] = (Ti[i+1] - Ti[i]) / dx

        elseif dad.pontos_internos[i, 1] == maximum(dad.pontos_internos[:, 1])
            ui[i] = (Ti[i] - Ti[i-1]) / dx

        else
            if indice_back(dad, i, 'x', x_order, y_order) !== nothing &&
               indice_forward(dad, i, 'x', x_order, y_order) !== nothing
                ui[i] = (Ti[i+1] - Ti[i-1]) / (2 * dx)

            elseif indice_forward(dad, i, 'x', x_order, y_order) !== nothing
                ui[i] = (Ti[i+1] - Ti[i]) / dx

            else
                ui[i] = (Ti[i] - Ti[i-1]) / dx
            end
        end
    end
    return ui, vi
end

function deriva_T(dad, T, arg = 1)
    nT = length(T)
    dT = zeros(nT)

    if arg == 1
        dx = sqrt((dad.NOS[10, 1] - dad.NOS[9, 1])^2 + (dad.NOS[10, 2] - dad.NOS[9, 2])^2)
    else
        dx = sqrt(
            (dad.pontos_internos[10, 1] - dad.pontos_internos[9, 1])^2 +
            (dad.pontos_internos[10, 2] - dad.pontos_internos[9, 2])^2,
        )
    end

    for i = 1:nT
        if i == 1
            dT[i] = (T[i+1] - T[i]) / dx
        elseif i == nT
            dT[i] = (T[i] - T[i-1]) / dx
        else
            dT[i] = (T[i+1] - T[i-1]) / (2 * dx)
        end
    end
    return dT
end

#=========POISSON=============#

function monta_ϕ(dad)

    xj = [dad.NOS; dad.pontos_internos]
    xi = [dad.NOS; dad.pontos_internos]

    n = size(xi, 1)

    ϕ = zeros(2 * n, 2 * n)
    r = zeros(n, n)

    for i = 1:n
        for k = 1:n
            r[i, k] = norm(xj[i, :] - xi[k, :]) + 1e-9
        end
    end


    f = 1 .- r
    #f = r.^2 .* log.(r)
    #f = r.^3
    #f=r
    for i = 1:2*n
        if i % 2 == 1
            ϕ[i, 1:2:end] = f[Int((i + 1) / 2), :]
        else
            ϕ[i, 2:2:end] = f[Int(i / 2), :]
        end
    end

    return ϕ

end

function u_matrix(u_int, nx, ny)

    u_matrix = zeros(nx, ny)
    v_matrix = zeros(nx, ny)

    for i = 1:nx
        u_matrix[i, :] = (u_int[(i-1)*ny+1:i*ny, 1])'
        v_matrix[i, :] = (u_int[(i-1)*ny+1:i*ny, 2])'
    end
    u_matrix, v_matrix
end

function deriva_int(u, v, dx, dy, nx, ny)

    dudx = zeros(nx, ny)
    dudy = zeros(nx, ny)
    dvdx = zeros(nx, ny)
    dvdy = zeros(nx, ny)

    for i = 1:nx
        for j = 1:ny
            # Derivadas em relação a x
            if i == 1
                dudx[i, j] = (u[i+1, j] - u[i, j]) / dx  # Diferença frontal
                dvdx[i, j] = (v[i+1, j] - v[i, j]) / dx
            elseif i == nx
                dudx[i, j] = (u[i, j] - u[i-1, j]) / dx  # Diferença para trás
                dvdx[i, j] = (v[i, j] - v[i-1, j]) / dx
            else
                dudx[i, j] = (u[i+1, j] - u[i-1, j]) / (2dx)  # Diferença centrada
                dvdx[i, j] = (v[i+1, j] - v[i-1, j]) / (2dx)
            end

            # Derivadas em relação a y
            if j == 1
                dudy[i, j] = (u[i, j+1] - u[i, j]) / dy  # Diferença frontal
                dvdy[i, j] = (v[i, j+1] - v[i, j]) / dy
            elseif j == ny
                dudy[i, j] = (u[i, j] - u[i, j-1]) / dy  # Diferença para trás
                dvdy[i, j] = (v[i, j] - v[i, j-1]) / dy
            else
                dudy[i, j] = (u[i, j+1] - u[i, j-1]) / (2dy)  # Diferença centrada
                dvdy[i, j] = (v[i, j+1] - v[i, j-1]) / (2dy)
            end
        end
    end

    return dudx, dudy, dvdx, dvdy
end

function refina_global(dad, dx, dy)

    xi = minimum(dad.NOS[:, 1]) + dx
    yi = minimum(dad.NOS[:, 2]) + dy

    xf = maximum(dad.NOS[:, 1]) - dx
    yf = maximum(dad.NOS[:, 2]) - dy

    px = xi:dx:xf
    py = yi:dy:yf

    nx = length(px)
    ny = length(py)

    pontos_refinados = zeros(nx * ny, 2)

    for i = 1:ny
        for j = 1:nx
            pontos_refinados[(i-1)*nx+j, 1] = px[j]
            pontos_refinados[(i-1)*nx+j, 2] = py[i]
        end
    end
    return pontos_refinados
end

function refina_local(
    dad,
    dx,
    dy,
    retangulo = [1, 1, 1, 1],
    objeto = false,
    intervalo_objeto = 1:4*order*nelem,
)

    xi = retangulo[1]
    yi = retangulo[2]

    xf = retangulo[3]
    yf = retangulo[4]

    px = xi:dx:xf
    py = yi:dy:yf

    nx = length(px)
    ny = length(py)

    pontos_retirados = deepcopy(dad.pontos_internos)
    pontos_refinados = zeros(nx * ny, 2)

    # Remove pontos internos dessa região

    i_remover = []

    for i = 1:ni(dad)
        if pontos_retirados[i, 1] < xf &&
           pontos_retirados[i, 1] > xi &&
           pontos_retirados[i, 2] < yf &&
           pontos_retirados[i, 2] > yi
            push!(i_remover, i)
        end
    end

    pontos_retirados = pontos_retirados[setdiff(1:end, i_remover), :]

    # Cria a malha nesse local

    for i = 1:ny
        for j = 1:nx
            pontos_refinados[(i-1)*nx+j, 1] = px[j]
            pontos_refinados[(i-1)*nx+j, 2] = py[i]
        end
    end

    # Remove os pontos do objeto em torno no qual se refina

    i_interno = Int[]  # Vetor vazio de inteiros

    # Iterar sobre os pontos para identificar os internos

    tamanho = length(pontos_refinados[:, 1])
    if objeto
        for i = 1:tamanho
            dist =
                sqrt.(
                    (dad.NOS[intervalo_objeto, 1] .- pontos_refinados[i, 1]) .^ 2 +
                    (dad.NOS[intervalo_objeto, 2] .- pontos_refinados[i, 2]) .^ 2
                )
            indice = findfirst(x -> x == minimum(dist), dist)
            if dot(pontos_refinados[i, :], dad.normal[indice, :]) >
               dot(dad.NOS[indice, :], dad.normal[indice, :])
                push!(i_interno, i)
            end
        end
    end

    # Remover os pontos internos
    pontos_refinados = pontos_refinados[setdiff(1:end, i_interno), :]

    refinado_local = [pontos_retirados; pontos_refinados]
    #refinado_local = pontos_refinados

    return refinado_local
end


function poisson_Houbolt(dad, f, iter, u, u_dot, A, b, M)

    u_dot[:, iter] =
        (-18 * u[:, iter-1] + 9 * (u[:, iter-2]) - 2 * (u[:, iter-3])) / (6 * dt)
    x = A \ (b + M * (f + u_dot[:, iter]))

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    Ti = T[nc(dad)+1:end]
    Tc = T[1:nc(dad)]

    return Tc, Ti, q, u_dot
end

function poisson_Adams(dad, f, iter, u, u_dot, A1, A2, b1, b2, M)

    if iter < 7
        u_dot[:, iter] =
            (-18 * u[:, iter-1] + 9 * (u[:, iter-2]) - 2 * (u[:, iter-3])) / (6 * dt)
        x = A1 \ (b1 + M * (f + u_dot[:, iter]))
    else
        u_dot[:, iter] =
            -24 / 55 * u[:, iter-1] / dt + 59 / 55 * u_dot[:, iter-1] -
            37 / 55 * u_dot[:, iter-2] + 9 / 55 * u_dot[:, iter-3]
        x = A2 \ (b2 + M * (f + u_dot[:, iter]))
    end

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    Ti = T[nc(dad)+1:end]
    Tc = T[1:nc(dad)]

    return Tc, Ti, q, u_dot
end

function poisson_Euler(dad, f, iter, u, A, b, M)

    u_dot = -u[:, iter-1] / dt
    x = A \ (b + M * (f + u_dot))

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    Ti = T[nc(dad)+1:end]
    Tc = T[1:nc(dad)]

    return Tc, Ti, q
end

function poisson_Euler_unsteady(dad, f, u, A, b, M)

    u_dot = -u / tau
    x = A \ (b + M * (f + u_dot))

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    Ti = T[nc(dad)+1:end]
    Tc = T[1:nc(dad)]

    return Tc, Ti, q
end


function poisson_steady(dad, f, A, b, M)

    x = A \ (b + M * (f))

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    Ti = T[nc(dad)+1:end]
    Tc = T[1:nc(dad)]

    return Tc, Ti, q
end


function poisson_ψ(f, ψ_flux)

    # Não está eficiente pois calculada Ht,Gt e M a cada passo

    ## Formatação dos dados ________________________________________________

    dad = format_dad(cavity_ψ(nelem, order), NPX, NPY)

    Ht, Gt = BEM.calc_HeGt(dad)

    M = BEM.Monta_M_RIMd(dad, npg)
    M = M[nc(dad)+1:ni(dad)+nc(dad), :]

    A, b = BEM.aplicaCDC(Ht, Gt, dad)
    A = A[nc(dad)+1:ni(dad)+nc(dad), :]
    b = b[nc(dad)+1:ni(dad)+nc(dad)]

    for i = 1:ni(dad)
        A[i, 1:nc(dad)] = A[i, 1:nc(dad)] .* ψ_flux
    end

    x = A \ (b + M * (f))

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    Ti = T[nc(dad)+1:end]
    Tc = T[1:nc(dad)]

    return Tc, Ti, q
end

function poisson_ω_steady(dad, f, ψ_int, u_wall, j_int, j_normal)

    # Corrige CDC considerando dif finitas p derivada de 2 ordem nas paredes

    ψ_wall = 0

    for i = 1:nelem*4

        indices = dad.ELEM[i].indices

        for (index, j) in enumerate(indices)

            if dad.normal[j, 1] == 1 || dad.normal[j, 1] == -1
                dad.ELEM[i].valorCDC[index] =
                    2 * (ψ_wall - ψ_int[j_int[j]]) /
                    ((dad.NOS[j, 1] - dad.pontos_internos[j_int[j], 1])^2)
            elseif dad.normal[j, 2] == 1
                dad.ELEM[i].valorCDC[index] =
                    2 * (ψ_wall - ψ_int[j_int[j]]) /
                    ((dad.NOS[j, 2] - dad.pontos_internos[j_int[j], 2])) -
                    u_wall * 2 / (dad.NOS[j, 2] - dad.pontos_internos[j_int[j], 2])
            elseif dad.normal[j, 2] == -1
                dad.ELEM[i].valorCDC[index] =
                    2 * (ψ_wall - ψ_int[j_int[j]]) /
                    ((dad.NOS[j, 2] - dad.pontos_internos[j_int[j], 2])^2)
            else
                (ψ_wall - ψ_int[j_normal]) /
                (2 * (dad.NOS[j, 2] - dad.pontos_internos[j_normal, 2]))

            end
        end
    end

    #_________________________

    Ht, Gt = BEM.calc_HeGt(dad)
    A, b = BEM.aplicaCDC(Ht, Gt, dad)

    x = A \ (b + M * (f))

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    return T, q
end

function poisson_ω_transiente(dad, f, ψ_int, u, iter, u_wall)

    # Corrige CDC considerando dif finitas p derivada de 2 ordem nas paredes
    # Tentando fazer o caso von_karman, temos ψ = 0 no circulo, na face de entrada vai linearmente de 0 a 1

    ψ_wall = 0

    for i = 1:nelem*4
        indices = dad.ELEM[i].indices
        for (index, j) in enumerate(indices)

            dist =
                (dad.NOS[j, 1] .- dad.pontos_internos[:, 1]) .^ 2 .+
                (dad.NOS[j, 2] .- dad.pontos_internos[:, 2]) .^ 2
            j_int = findfirst(x -> x == minimum(dist), dist)

            if dad.normal[j, 1] == 1 || dad.normal[j, 1] == -1
                dad.ELEM[i].valorCDC[index] =
                    2 * (ψ_wall - ψ_int[j_int]) /
                    ((dad.NOS[j, 1] - dad.pontos_internos[j_int, 1])^2)
            elseif dad.normal[j, 2] == 1
                dad.ELEM[i].valorCDC[index] =
                    2 * (ψ_wall - ψ_int[j_int]) /
                    ((dad.NOS[j, 2] - dad.pontos_internos[j_int, 2])) -
                    u_wall * 2 / (dad.NOS[j, 2] - dad.pontos_internos[j_int, 2])
            elseif dad.normal[j, 2] == -1
                dad.ELEM[i].valorCDC[index] =
                    2 * (ψ_wall - ψ_int[j_int]) /
                    ((dad.NOS[j, 2] - dad.pontos_internos[j_int, 2])^2)
            else
                normal = dad.normal[i, :]
                dist_normal = zeros(ni(dad))
                for k = 1:ni(dad)
                    dist_normal[k] =
                        dot(dad.NOS[i, :], normal) - dot(dad.pontos_internos[k, :], normal)
                end
                j_normal = findfirst(x -> x == minimum(abs.(dist_normal)), dist_normal)

                (ψ_wall - ψ_int[j_normal]) /
                (2 * (dad.NOS[j, 2] - dad.pontos_internos[j_normal, 2]))

            end
        end
    end

    #_________________________

    Ht, Gt = BEM.calc_HeGt(dad)

    M = BEM.Monta_M_RIMd(dad, npg)# calc_HeG_potencial linha 310

    A, b = BEM.aplicaCDC(Ht - 11 * M / (6 * dt), Gt, dad) #Houbolt

    u_n = u[:, iter-1]
    u_n1 = u[:, iter-2]
    u_n2 = u[:, iter-3]

    u_dot = (-18 * u_n + 9 * (u_n1) - 2 * (u_n2)) / (6 * dt)

    x = A \ (b + M * (f + u_dot))

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    Ti = T[nc(dad)+1:end]
    Tc = T[1:nc(dad)]

    return Tc, Ti, q
end

function poisson_ω_transiente_von_karman(dad, f, ψ_int, u, iter)

    # Corrige CDC considerando dif finitas p derivada de 2 ordem nas paredes
    # Tentando fazer o caso von_karman, temos ψ = 0 no circulo, na face de entrada vai linearmente de 0 a 1

    #ψ_wall = 0

    ψ_wall = zeros(nc(dad))
    ψ_wall[1:5*order*nelem] .= 0
    ψ_wall[6*order*nelem+1:7*order*nelem] .= 0

    for i = 1:order*nelem
        ψ_wall[7*order*nelem+i:8*order*nelem] = i / (order * nelem)
        ψ_wall[5*order*nelem+i:6*order*nelem] = i / (order * nelem)
    end

    for i = 1:nelem*8
        indices = dad.ELEM[i].indices
        for (index, j) in enumerate(indices)

            dist =
                (dad.NOS[j, 1] .- dad.pontos_internos[:, 1]) .^ 2 .+
                (dad.NOS[j, 2] .- dad.pontos_internos[:, 2]) .^ 2
            j_int = findfirst(x -> x == minimum(dist), dist)

            if dad.normal[j, 1] == 1 || dad.normal[j, 1] == -1
                dad.ELEM[i].valorCDC[index] =
                    (ψ_wall[j] - ψ_int[j_int]) /
                    (2 * (dad.NOS[j, 1] - dad.pontos_internos[j_int, 1])^2)
            elseif dad.normal[j, 2] == 1
                dad.ELEM[i].valorCDC[index] =
                    (ψ_wall[j] - ψ_int[j_int]) /
                    (2 * (dad.NOS[j, 2] - dad.pontos_internos[j_int, 2]))
            elseif dad.normal[j, 2] == -1
                dad.ELEM[i].valorCDC[index] =
                    (ψ_wall[j] - ψ_int[j_int]) /
                    (2 * (dad.NOS[j, 2] - dad.pontos_internos[j_int, 2])^2)
            else
                normal = dad.normal[j, :]
                dist_normal =
                    dot(dad.NOS[i, :], normal) .- dot.(dad.pontos_internos[:, :], normal)
                j_normal = findfirst(x -> x == minimum(abs.(dist_normal)), dist_normal)

                (ψ_wall[j] - ψ_int[j_normal]) /
                (2 * (dad.NOS[j, 2] - dad.pontos_internos[j_normal, 2]))

            end
        end
    end

    #_________________________

    Ht, Gt = BEM.calc_HeGt(dad)

    M = BEM.Monta_M_RIMd(dad, npg)# calc_HeG_potencial linha 310

    A, b = BEM.aplicaCDC(Ht - 11 * M / (6 * dt), Gt, dad) #Houbolt

    u_n = u[:, iter-1]
    u_n1 = u[:, iter-2]
    u_n2 = u[:, iter-3]

    u_dot = (-18 * u_n + 9 * (u_n1) - 2 * (u_n2)) / (6 * dt)

    x = A \ (b + M * (f + u_dot))

    T = zeros(nc(dad) + ni(dad))
    T[1:nc(dad)], q = separa(dad, x)
    T[nc(dad)+1:end] = x[nc(dad)+1:end]

    Ti = T[nc(dad)+1:end]
    Tc = T[1:nc(dad)]

    return Tc, Ti, q
end


function deriva_contorno(dad, u, t, δ)

    # u = u[contorno]

    derivada_x = zeros(nc(dad))
    derivada_y = zeros(nc(dad))

    for i = 1:nc(dad)

        if dad.normal[i, 1] == 1
            if i != 1 && i != nc(dad)
                derivada_y[i] = (u[i+1] - u[i-1]) / (2 * δ[i])
            elseif i != 1
                derivada_y[i] = (u[i] - u[i-1]) / (δ[i])
            else
                derivada_y[i] = (u[i+1] - u[i]) / (δ[i])
            end

            derivada_x[i] = t[i]

        elseif dad.normal[i, 1] == -1

            if i != 1 && i != nc(dad)
                derivada_y[i] = -(u[i+1] - u[i-1]) / (2 * δ[i])
            elseif i != nc(dad)
                derivada_y[i] = -(u[i+1] - u[i]) / (δ[i])
            else
                derivada_y[i] = -(u[i] - u[i-1]) / (δ[i])
            end

            derivada_x[i] = -t[i]

        elseif dad.normal[i, 2] == 1
            if i != 1 && i != nc(dad)
                derivada_x[i] = -(u[i+1] - u[i-1]) / (2 * δ[i])
            elseif i != nc(dad)
                derivada_x[i] = -(u[i+1] - u[i]) / (δ[i])
            else
                derivada_x[i] = -(u[i] - u[i-1]) / (δ[i])
            end

            derivada_y[i] = t[i]

        elseif dad.normal[i, 2] == -1
            if i != 1 && i != nc(dad)
                derivada_x[i] = (u[i+1] - u[i-1]) / (2 * δ[i])
            elseif i != 1
                derivada_x[i] = (u[i] - u[i-1]) / (δ[i])
            else
                derivada_x[i] = (u[i+1] - u[i]) / (δ[i])
            end

            derivada_y[i] = -t[i]
        end
    end
    return derivada_x, derivada_y
end

function deriva_contorno_com_interno(
    dad,
    u,
    u_int,
    derivada_x_int,
    derivada_y_int,
    index_dist_cont_int,
    dist_cont_cont,
)

    # Tomar cuidado com derivada nas quinas

    derivada_x = zeros(nc(dad))
    derivada_y = zeros(nc(dad))

    for i = 1:nc(dad)

        if dad.normal[i, 1] == 1

            if i != 1 && i != nc(dad)
                derivada_y[i] = (u[i+1] - u[i-1]) / (2 * dist_cont_cont[i])
            elseif i != 1
                derivada_y[i] = (u[i] - u[i-1]) / (dist_cont_cont[i])
            else
                derivada_y[i] = (u[i+1] - u[i]) / (dist_cont_cont[i])
            end

            derivada_x[i] =
                (u[i] - u_int[index_dist_cont_int[i]]) /
                (dad.NOS[i, 1] - dad.pontos_internos[index_dist_cont_int[i], 1])

        elseif dad.normal[i, 1] == -1

            if i != 1 && i != nc(dad)
                derivada_y[i] = -(u[i+1] - u[i-1]) / (2 * dist_cont_cont[i])
            elseif i != nc(dad)
                derivada_y[i] = -(u[i+1] - u[i]) / (dist_cont_cont[i])
            else
                derivada_y[i] = -(u[i] - u[i-1]) / (dist_cont_cont[i])
            end

            derivada_x[i] =
                (u[i] - u_int[index_dist_cont_int[i]]) /
                (dad.NOS[i, 1] - dad.pontos_internos[index_dist_cont_int[i], 1])

        elseif dad.normal[i, 2] == 1
            if i != 1 && i != nc(dad)
                derivada_x[i] = -(u[i+1] - u[i-1]) / (2 * dist_cont_cont[i])
            elseif i != nc(dad)
                derivada_x[i] = -(u[i+1] - u[i]) / (dist_cont_cont[i])
            else
                derivada_x[i] = -(u[i] - u[i-1]) / (dist_cont_cont[i])
            end

            derivada_y[i] =
                (u[i] - u_int[index_dist_cont_int[i]]) /
                (dad.NOS[i, 2] - dad.pontos_internos[index_dist_cont_int[i], 2])

        elseif dad.normal[i, 2] == -1
            if i != 1 && i != nc(dad)
                derivada_x[i] = (u[i+1] - u[i-1]) / (2 * dist_cont_cont[i])
            elseif i != 1
                derivada_x[i] = (u[i] - u[i-1]) / (dist_cont_cont[i])
            else
                derivada_x[i] = (u[i+1] - u[i]) / (dist_cont_cont[i])
            end

            derivada_y[i] =
                (u[i] - u_int[index_dist_cont_int[i]]) /
                (dad.NOS[i, 2] - dad.pontos_internos[index_dist_cont_int[i], 2])
        else
            derivada_y[i] = derivada_y_int[index_dist_cont_int[i]]
            derivada_x[i] = derivada_x_int[index_dist_cont_int[i]]
        end
    end
    return derivada_x, derivada_y
end

function deriva_interno_generico(
    dad,
    u,
    dx,
    dy,
    index_backward,
    index_forward,
    index_dist_int_int,
)

    n_int = length(dad.pontos_internos[:, 1])

    derivada_x = zeros(n_int)
    derivada_y = zeros(n_int)

    for i = 1:n_int

        if index_backward[i, 1] !== nothing && index_forward[i, 1] !== nothing
            derivada_x[i] =
                (u[Int(index_forward[i, 1])] - u[Int(index_backward[i, 1])]) / (2 * dx)

        elseif index_backward[i, 1] !== nothing
            derivada_x[i] = (u[i] - u[Int(index_backward[i, 1])]) / (dx)

        elseif index_forward[i, 1] !== nothing
            derivada_x[i] = (u[Int(index_forward[i, 1])] - u[i]) / (dx)

        else
            derivada_x[i] =
                (u[Int(index_dist_int_int[i])] - u[i]) / (
                    dad.pontos_internos[Int(index_dist_int_int[i]), 1] -
                    dad.pontos_internos[i, 1]
                )
        end

        if index_backward[i, 2] !== nothing && index_forward[i, 2] !== nothing
            derivada_y[i] =
                (u[Int(index_forward[i, 2])] - u[Int(index_backward[i, 2])]) / (2 * dy)

        elseif index_backward[i, 2] !== nothing
            derivada_y[i] = (u[i] - u[Int(index_backward[i, 2])]) / (dy)

        elseif index_forward[i, 2] !== nothing
            derivada_y[i] = (u[Int(index_forward[i, 2])] - u[i]) / (dy)

        else
            derivada_y[i] =
                (u[Int(index_dist_int_int[i])] - u[i]) / (
                    dad.pontos_internos[Int(index_dist_int_int[i]), 2] -
                    dad.pontos_internos[i, 2]
                )
        end

    end
    return derivada_x, derivada_y
end

function deriva_interno_generico_refinado(
    dad,
    u,
    index_backward,
    index_forward,
    index_dist_int_int,
)

    n_int = length(dad.pontos_internos[:, 1])

    derivada_x = zeros(n_int)
    derivada_y = zeros(n_int)

    #Tomar cuidado na transição da zona n refinada para a refinada

    for i = 1:n_int

        if index_backward[i, 1] !== nothing && index_forward[i, 1] !== nothing
            derivada_x[i] =
                (u[Int(index_forward[i, 1])] - u[Int(index_backward[i, 1])]) / (
                    2 * (
                        dad.pontos_internos[Int(index_forward[i, 1]), 1] -
                        dad.pontos_internos[Int(index_backward[i, 1]), 1]
                    )
                )

        elseif index_backward[i, 1] !== nothing
            derivada_x[i] =
                (u[i] - u[Int(index_backward[i, 1])]) / ((
                    dad.pontos_internos[i, 1] -
                    dad.pontos_internos[Int(index_backward[i, 1]), 1]
                ))

        elseif index_forward[i, 1] !== nothing
            derivada_x[i] =
                (u[Int(index_forward[i, 1])] - u[i]) / ((
                    dad.pontos_internos[Int(index_forward[i, 1]), 1] -
                    dad.pontos_internos[i, 1]
                ))

        else
            derivada_x[i] =
                (u[Int(index_dist_int_int[i])] - u[i]) / (
                    dad.pontos_internos[Int(index_dist_int_int[i]), 1] -
                    dad.pontos_internos[i, 1]
                )
        end

        if index_backward[i, 2] !== nothing && index_forward[i, 2] !== nothing
            derivada_y[i] =
                (u[Int(index_forward[i, 2])] - u[Int(index_backward[i, 2])]) / (
                    2 * (
                        dad.pontos_internos[Int(index_forward[i, 2]), 2] -
                        dad.pontos_internos[Int(index_backward[i, 2]), 2]
                    )
                )

        elseif index_backward[i, 2] !== nothing
            derivada_y[i] =
                (u[i] - u[Int(index_backward[i, 2])]) / ((
                    dad.pontos_internos[i, 2] -
                    dad.pontos_internos[Int(index_backward[i, 2]), 2]
                ))

        elseif index_forward[i, 2] !== nothing
            derivada_y[i] =
                (u[Int(index_forward[i, 2])] - u[i]) / ((
                    dad.pontos_internos[Int(index_forward[i, 2]), 2] -
                    dad.pontos_internos[i, 2]
                ))

        else
            derivada_y[i] =
                (u[Int(index_dist_int_int[i])] - u[i]) / (
                    dad.pontos_internos[Int(index_dist_int_int[i]), 2] -
                    dad.pontos_internos[i, 2]
                )
        end

    end
    return derivada_x, derivada_y
end

function deriva_contorno_generico(dad, u, t, δ, transição)

    # Derivada nas quinas incorreta pois o delta é calculado errado e 

    derivada_x = zeros(nc(dad))
    derivada_y = zeros(nc(dad))

    for i = 1:nc(dad)

        if dad.normal[i, 1] == 1

            if i != 1 && i != nc(dad) && i != transição
                derivada_y[i] = (u[i+1] - u[i-1]) / (2 * δ[i])
            elseif i != 1
                derivada_y[i] = (u[i] - u[i-1]) / (δ[i])
            else
                derivada_y[i] = (u[i+1] - u[i]) / (δ[i])
            end

            derivada_x[i] = t[i]

        elseif dad.normal[i, 1] == -1

            if i != 1 && i != nc(dad) && i != transição
                derivada_y[i] = -(u[i+1] - u[i-1]) / (2 * δ[i])
            elseif i != 1
                derivada_y[i] = -(u[i] - u[i-1]) / (δ[i])
            else
                derivada_y[i] = -(u[i+1] - u[i]) / (δ[i])
            end

            derivada_x[i] = -t[i]

        elseif dad.normal[i, 2] == 1
            if i != 1 && i != nc(dad) && i != transição
                derivada_x[i] = -(u[i+1] - u[i-1]) / (2 * δ[i])
            elseif i != 1
                derivada_x[i] = -(u[i] - u[i-1]) / (δ[i])
            else
                derivada_x[i] = -(u[i+1] - u[i]) / (δ[i])
            end

            derivada_y[i] = t[i]

        elseif dad.normal[i, 2] == -1
            if i != 1 && i != nc(dad) && i != transição
                derivada_x[i] = (u[i+1] - u[i-1]) / (2 * δ[i])
            elseif i != 1
                derivada_x[i] = (u[i] - u[i-1]) / (δ[i])
            else
                derivada_x[i] = (u[i+1] - u[i]) / (δ[i])
            end

            derivada_y[i] = -t[i]
        else
            derivada_normal = t[i]

            normal = dad.normal[i, :]
            tangencial = [-normal[2], normal[1]]

            if i != 1 && i != nc(dad) && i != transição
                derivada_tangencial = (u[i+1] - u[i-1]) / (2 * δ[i])

            elseif i != 1
                derivada_tangencial = (u[i] - u[i-1]) / (δ[i])
            else
                derivada_tangencial = (u[i+1] - u[i]) / (δ[i])
            end

            #Projeção da derivada
            derivada_x[i] =
                dot([1, 0], normal) * derivada_normal +
                dot([1, 0], tangencial) * derivada_tangencial
            derivada_y[i] =
                dot([0, 1], normal) * derivada_normal +
                dot([0, 1], tangencial) * derivada_tangencial
        end
    end
    return derivada_x, derivada_y
end

function indice_back(dad, i, arg, x_order, y_order)

    if arg == "x"

        idx_order = findfirst(x -> x == dad.pontos_internos[i, 1], x_order)

        if idx_order == 1
            indice = nothing

        else
            idx_dad = findfirst(x -> x == x_order[idx_order-1], dad.pontos_internos[:, 1])
            indice = findfirst(
                row ->
                    row[1] == dad.pontos_internos[idx_dad, 1] &&
                        row[2] == dad.pontos_internos[i, 2],
                eachrow(dad.pontos_internos),
            )
        end

    else

        idx_order = findfirst(x -> x == dad.pontos_internos[i, 2], y_order)

        if idx_order == 1
            indice = nothing
        else
            idx_dad = findfirst(x -> x == y_order[idx_order-1], dad.pontos_internos[:, 2])
            indice = findfirst(
                row ->
                    row[1] == dad.pontos_internos[i, 1] &&
                        row[2] == dad.pontos_internos[idx_dad, 2],
                eachrow(dad.pontos_internos),
            )
        end
    end

    return indice
end

function indice_forward(dad, i, arg, x_order, y_order)

    if arg == "x"

        idx_order = findfirst(x -> x == dad.pontos_internos[i, 1], x_order)

        if idx_order == length(x_order)

            indice = nothing
        else
            idx_dad = findfirst(x -> x == x_order[idx_order+1], dad.pontos_internos[:, 1])

            indice = findfirst(
                row ->
                    row[1] == dad.pontos_internos[idx_dad, 1] &&
                        row[2] == dad.pontos_internos[i, 2],
                eachrow(dad.pontos_internos),
            )
        end
    else
        idx_order = findfirst(x -> x == dad.pontos_internos[i, 2], y_order)

        if idx_order == length(y_order)
            indice = nothing
        else
            idx_dad = findfirst(x -> x == y_order[idx_order+1], dad.pontos_internos[:, 2])

            indice = findfirst(
                row ->
                    row[1] == dad.pontos_internos[i, 1] &&
                        row[2] == dad.pontos_internos[idx_dad, 2],
                eachrow(dad.pontos_internos),
            )
        end
    end

    return indice
end

# Derivadas por diferenças finitas

function Dₓ_1D(n::Int, dx::Float64)
    # Inicializar a matriz como uma matriz esparsa para eficiência
    D = spzeros(n, n)

    # Preencher as fórmulas para as bordas
    D[1, 1] = -3 / 2
    D[1, 2] = 2
    D[1, 3] = -1 / 2
    D[n, n-2] = 1 / 2
    D[n, n-1] = -2
    D[n, n] = 3 / 2

    # Preencher diferenças centradas para o interior
    for i = 2:n-1
        D[i, i-1] = -1 / 2
        D[i, i+1] = 1 / 2
    end

    # Dividir por dx
    D /= dx
    return D
end

function Dₓ_1D_backwind(n::Int, dx::Float64)
    # Inicializar a matriz como uma matriz esparsa para eficiência
    D = spzeros(n, n)

    # Preencher as fórmulas para backwind (diferenciação para trás)
    for i = 2:n
        D[i, i] = 1
        D[i, i-1] = -1
    end

    # Ajustar a borda
    D[1, 1] = 1
    D[1, 2] = -1

    # Dividir por dx
    D /= dx
    return D
end

function Dₓ_1D_upwind(n::Int, dx::Float64)
    # Inicializar a matriz como uma matriz esparsa para eficiência
    D = spzeros(n, n)

    # Preencher as fórmulas para upwind (diferenciação para frente)
    for i = 1:n-1
        D[i, i] = -1
        D[i, i+1] = 1
    end

    # Ajustar a borda
    D[n, n] = -1
    D[n, n-1] = 1

    # Dividir por dx
    D /= dx
    return D
end

function Dₓ(nx, ny, h)
    N = nx * ny  # Número total de nós na malha 2D
    Dx = spzeros(N, N)  # Matriz esparsa para armazenar diferenças finitas

    for j = 1:ny  # Itera ao longo de cada linha da malha (direção y)
        for i = 1:nx  # Itera ao longo de cada coluna da malha (direção x)
            idx = (j - 1) * nx + i  # Índice linear baseado em i (coluna) e j (linha)

            if i == 1
                # Primeira coluna (condição de contorno à esquerda)
                Dx[idx, idx] = -3 / (2 * h)
                Dx[idx, idx+1] = 2 / h
                if i + 2 <= nx
                    Dx[idx, idx+2] = -1 / (2 * h)
                end
            elseif i == nx
                # Última coluna (condição de contorno à direita)
                if i - 2 > 0
                    Dx[idx, idx-2] = 1 / (2 * h)
                end
                Dx[idx, idx-1] = -2 / h
                Dx[idx, idx] = 3 / (2 * h)
            else
                # Ponto interno
                Dx[idx, idx-1] = -1 / (2 * h)
                Dx[idx, idx+1] = 1 / (2 * h)
            end
        end
    end

    return Dx
end

function Dₓ_forward(nx, ny, h)
    N = nx * ny  # Número total de nós na malha 2D
    Dx = spzeros(N, N)  # Matriz esparsa para armazenar diferenças finitas

    for j = 1:ny  # Itera ao longo de cada linha da malha (direção y)
        for i = 1:nx  # Itera ao longo de cada coluna da malha (direção x)
            idx = (j - 1) * nx + i  # Índice linear baseado em i (coluna) e j (linha)

            if i == nx
                # Última coluna: diferenças regressivas (backward)
                if i - 2 > 0
                    Dx[idx, idx-2] = -1 / (2 * h)
                end
                Dx[idx, idx-1] = 1 / h
                Dx[idx, idx] = -1 / h
            elseif i == nx - 1
                # Penúltima coluna: diferenças centrais, pois forward não é possível para todos os pontos
                Dx[idx, idx-1] = -1 / (2 * h)
                Dx[idx, idx+1] = 1 / (2 * h)
            else
                # Ponto interno: diferenças forward
                Dx[idx, idx] = -1 / h
                Dx[idx, idx+1] = 1 / h
            end
        end
    end

    return Dx
end

function Dy(nx, ny, h)
    N = nx * ny  # Número total de nós na malha 2D
    Dy = spzeros(N, N)  # Matriz esparsa para armazenar diferenças finitas

    for j = 1:ny  # Itera ao longo de cada linha da malha (direção y)
        for i = 1:nx  # Itera ao longo de cada coluna da malha (direção x)
            idx = (j - 1) * nx + i  # Índice linear baseado em i (coluna) e j (linha)

            if j == 1
                # Primeira linha (condição de contorno inferior)
                Dy[idx, idx] = -3 / (2 * h)
                Dy[idx, idx+nx] = 2 / h
                if j + 2 <= ny
                    Dy[idx, idx+2*nx] = -1 / (2 * h)
                end
            elseif j == ny
                # Última linha (condição de contorno superior)
                if j - 2 > 0
                    Dy[idx, idx-2*nx] = 1 / (2 * h)
                end
                Dy[idx, idx-nx] = -2 / h
                Dy[idx, idx] = 3 / (2 * h)
            else
                # Ponto interno
                Dy[idx, idx-nx] = -1 / (2 * h)
                Dy[idx, idx+nx] = 1 / (2 * h)
            end
        end
    end

    return Dy
end

function Dy_forward(nx, ny, h)
    N = nx * ny  # Número total de nós na malha 2D
    Dy = spzeros(N, N)  # Matriz esparsa para armazenar diferenças finitas

    for j = 1:ny  # Itera ao longo de cada linha da malha (direção y)
        for i = 1:nx  # Itera ao longo de cada coluna da malha (direção x)
            idx = (j - 1) * nx + i  # Índice linear baseado em i (coluna) e j (linha)

            if j == ny
                # Última linha: diferenças regressivas (backward)
                if j - 2 > 0
                    Dy[idx, idx-2*nx] = -1 / (2 * h)
                end
                Dy[idx, idx-nx] = 1 / h
                Dy[idx, idx] = -1 / h
            elseif j == ny - 1
                # Penúltima linha: diferenças centrais, pois forward não é possível
                Dy[idx, idx-nx] = -1 / (2 * h)
                Dy[idx, idx+nx] = 1 / (2 * h)
            else
                # Ponto interno: diferenças forward
                Dy[idx, idx] = -1 / h
                Dy[idx, idx+nx] = 1 / h
            end
        end
    end

    return Dy
end

function Dₓ_backward(nx, ny, h)
    N = nx * ny  # Número total de nós na malha 2D
    Dx = spzeros(N, N)  # Matriz esparsa para armazenar diferenças finitas

    for j = 1:ny  # Itera ao longo de cada linha da malha (direção y)
        for i = 1:nx  # Itera ao longo de cada coluna da malha (direção x)
            idx = (j - 1) * nx + i  # Índice linear baseado em i (coluna) e j (linha)

            if i == 1
                # Primeira coluna: diferenças progressivas (forward)
                Dx[idx, idx] = -1 / h
                Dx[idx, idx+1] = 1 / h
            elseif i == 2
                # Segunda coluna: diferenças centrais
                Dx[idx, idx-1] = -1 / (2 * h)
                Dx[idx, idx+1] = 1 / (2 * h)
            else
                # Ponto interno: diferenças regressivas (backward)
                Dx[idx, idx-1] = 1 / h
                Dx[idx, idx] = -1 / h
            end
        end
    end

    return Dx
end

function Dy_backward(nx, ny, h)
    N = nx * ny  # Número total de nós na malha 2D
    Dy = spzeros(N, N)  # Matriz esparsa para armazenar diferenças finitas

    for j = 1:ny  # Itera ao longo de cada linha da malha (direção y)
        for i = 1:nx  # Itera ao longo de cada coluna da malha (direção x)
            idx = (j - 1) * nx + i  # Índice linear baseado em i (coluna) e j (linha)

            if j == 1
                # Primeira linha: diferenças centrais (não é possível usar backward aqui)
                if j + 2 <= ny
                    Dy[idx, idx+2*nx] = 1 / (2 * h)
                end
                Dy[idx, idx+nx] = -2 / h
                Dy[idx, idx] = 3 / (2 * h)
            elseif j == 2
                # Segunda linha: diferenças centrais
                Dy[idx, idx-nx] = -1 / (2 * h)
                Dy[idx, idx+nx] = 1 / (2 * h)
            else
                # Ponto interno: diferenças regressivas (backward)
                Dy[idx, idx] = 1 / h
                Dy[idx, idx-nx] = -1 / h
            end
        end
    end

    return Dy
end

function Dₓₓ(nx, ny, h)
    n = nx * ny  # Total de pontos no grid 2D
    Dxx = zeros(Float64, n, n)

    for j = 1:ny
        for i = 1:nx
            idx = (j - 1) * nx + i  # Índice global do ponto (i, j)

            if i == 1
                # Forward difference para o limite esquerdo
                Dxx[idx, idx] = 2
                Dxx[idx, idx+1] = -5
                Dxx[idx, idx+2] = 4
                Dxx[idx, idx+3] = -1
            elseif i == 2 || i == nx - 1
                # Diferença central para os pontos próximos aos limites
                Dxx[idx, idx-1] = 1
                Dxx[idx, idx] = -2
                Dxx[idx, idx+1] = 1
            elseif i == nx
                # Backward difference para o limite direito
                Dxx[idx, idx] = 2
                Dxx[idx, idx-1] = -5
                Dxx[idx, idx-2] = 4
                Dxx[idx, idx-3] = -1
            else
                # Diferença central para os pontos internos
                Dxx[idx, idx-1] = 1
                Dxx[idx, idx] = -2
                Dxx[idx, idx+1] = 1
            end
        end
    end

    return Dxx / h^2
end

function Dyy(nx, ny, h)
    n = nx * ny  # Total de pontos no grid 2D
    Dyy = zeros(Float64, n, n)

    for j = 1:ny
        for i = 1:nx
            idx = (j - 1) * nx + i  # Índice global do ponto (i, j)

            if j == 1
                # Forward difference para o limite inferior
                Dyy[idx, idx] = 2
                Dyy[idx, idx+nx] = -5
                Dyy[idx, idx+2*nx] = 4
                Dyy[idx, idx+3*nx] = -1
            elseif j == 2 || j == ny - 1
                # Diferença central para os pontos próximos aos limites
                Dyy[idx, idx-nx] = 1
                Dyy[idx, idx] = -2
                Dyy[idx, idx+nx] = 1
            elseif j == ny
                # Backward difference para o limite superior
                Dyy[idx, idx] = 2
                Dyy[idx, idx-nx] = -5
                Dyy[idx, idx-2*nx] = 4
                Dyy[idx, idx-3*nx] = -1
            else
                # Diferença central para os pontos internos
                Dyy[idx, idx-nx] = 1
                Dyy[idx, idx] = -2
                Dyy[idx, idx+nx] = 1
            end
        end
    end

    return Dyy / h^2
end

function Dy_ordem_4(nx, ny, h)
    N = nx * ny  # Número total de nós na malha 2D
    Dy = spzeros(N, N)  # Matriz esparsa para armazenar diferenças finitas

    for j = 1:ny  # Itera ao longo de cada linha da malha (direção y)
        for i = 1:nx  # Itera ao longo de cada coluna da malha (direção x)
            idx = (j - 1) * nx + i  # Índice linear baseado em i (coluna) e j (linha)

            if j == 1
                # Primeira linha (condição de contorno inferior)
                Dy[idx, idx] = -25 / (12 * h)
                Dy[idx, idx+nx] = 4 / h
                if j + 2 <= ny
                    Dy[idx, idx+2*nx] = -3 / (2 * h)
                end
                if j + 3 <= ny
                    Dy[idx, idx+3*nx] = 4 / (3 * h)
                end
                if j + 4 <= ny
                    Dy[idx, idx+4*nx] = -1 / (4 * h)
                end
            elseif j == 2
                # Segunda linha (uma linha acima da borda inferior)
                Dy[idx, idx-nx] = -1 / (2 * h)
                Dy[idx, idx] = 0
                Dy[idx, idx+nx] = 1 / (2 * h)
            elseif j == ny
                # Última linha (condição de contorno superior)
                Dy[idx, idx-nx]
                Dy[idx, idx] = -4 / h
                if j - 2 > 0
                    Dy[idx, idx-2*nx] = 3 / (2 * h)
                end
                if j - 3 > 0
                    Dy[idx, idx-3*nx] = -4 / (3 * h)
                end
                if j - 4 > 0
                    Dy[idx, idx-4*nx] = 1 / (4 * h)
                end
            elseif j == ny - 1
                # Penúltima linha (uma linha abaixo da borda superior)
                Dy[idx, idx-nx] = -1 / (2 * h)
                Dy[idx, idx] = 0
                Dy[idx, idx+nx] = 1 / (2 * h)
            else
                # Ponto interno (diferenças centrais de 4ª ordem)
                Dy[idx, idx-2*nx] = 1 / (12 * h)
                Dy[idx, idx-nx] = -8 / (12 * h)
                Dy[idx, idx+nx] = 8 / (12 * h)
                Dy[idx, idx+2*nx] = -1 / (12 * h)
            end
        end
    end

    return Dy
end

function Dₓ_ordem_4(nx, ny, h)
    N = nx * ny  # Número total de nós na malha 2D
    Dx = spzeros(N, N)  # Matriz esparsa para armazenar diferenças finitas

    for j = 1:ny  # Itera ao longo de cada linha da malha (direção y)
        for i = 1:nx  # Itera ao longo de cada coluna da malha (direção x)
            idx = (j - 1) * nx + i  # Índice linear baseado em i (coluna) e j (linha)

            if i == 1
                # Primeira coluna (condição de contorno à esquerda)
                Dx[idx, idx] = -25 / (12 * h)
                Dx[idx, idx+1] = 4 / h
                if i + 2 <= nx
                    Dx[idx, idx+2] = -3 / (2 * h)
                end
                if i + 3 <= nx
                    Dx[idx, idx+3] = 4 / (3 * h)
                end
                if i + 4 <= nx
                    Dx[idx, idx+4] = -1 / (4 * h)
                end
            elseif i == 2
                # Segunda coluna (uma coluna à direita da borda esquerda)
                Dx[idx, idx-1] = -1 / (2 * h)
                Dx[idx, idx] = 0
                Dx[idx, idx+1] = 1 / (2 * h)
            elseif i == nx
                # Última coluna (condição de contorno à direita)
                Dx[idx, idx] = 25 / (12 * h)
                Dx[idx, idx-1] = -4 / h
                if i - 2 > 0
                    Dx[idx, idx-2] = 3 / (2 * h)
                end
                if i - 3 > 0
                    Dx[idx, idx-3] = -4 / (3 * h)
                end
                if i - 4 > 0
                    Dx[idx, idx-4] = 1 / (4 * h)
                end
            elseif i == nx - 1
                # Penúltima coluna (uma coluna à esquerda da borda direita)
                Dx[idx, idx-1] = -1 / (2 * h)
                Dx[idx, idx] = 0
                Dx[idx, idx+1] = 1 / (2 * h)
            else
                # Ponto interno (diferenças centrais de 4ª ordem)
                Dx[idx, idx-2] = 1 / (12 * h)
                Dx[idx, idx-1] = -8 / (12 * h)
                Dx[idx, idx+1] = 8 / (12 * h)
                Dx[idx, idx+2] = -1 / (12 * h)
            end
        end
    end

    return Dx
end

function ψ_dif_finitas(ψ, ω, nx, ny, Δx, Δy)
    # Precomputação de fatores constantes
    c1 = Δx^2 * Δy^2 / (2 * (Δx^2 + Δy^2))
    c2 = Δy^2 / (2 * (Δx^2 + Δy^2))
    c3 = Δx^2 / (2 * (Δx^2 + Δy^2))

    # Cria uma cópia de ψ para atualização
    ψ_new = copy(ψ)

    # Itera sobre a malha excluindo bordas
    for j = 2:(ny-1)
        for i = 2:(nx-1)
            idx = (j - 1) * nx + i  # Índice linear para (i, j)

            # Índices dos vizinhos
            idx_ip1 = idx + 1       # (i+1, j)
            idx_im1 = idx - 1       # (i-1, j)
            idx_jp1 = idx + nx      # (i, j+1)
            idx_jm1 = idx - nx      # (i, j-1)

            # Atualiza ψ usando a equação fornecida
            ψ_new[idx] =
                c1 * ω[idx] +
                c2 * (ψ_new[idx_ip1] + ψ_new[idx_im1]) +
                c3 * (ψ_new[idx_jp1] + ψ_new[idx_jm1])
        end
    end

    return ψ_new
end

function solve_poisson(nx, ny, dx, dy, q, psi)
    # nx: número de pontos em x
    # ny: número de pontos em y
    # dx: espaçamento da malha em x
    # dy: espaçamento da malha em y
    # q: vetor linearizado correspondente a -ω
    # psi: vetor linearizado inicial de ψ


    factor = 1 / (2 * (1 / dx^2 + 1 / dy^2))

    # Iterar sobre os pontos internos da malha
    for j = 2:(ny-1)
        for i = 2:(nx-1)
            idx = (j - 1) * nx + i  # Índice linear para (i, j)

            # Índices dos vizinhos
            idx_ip1 = idx + 1       # (i+1, j)
            idx_im1 = idx - 1       # (i-1, j)
            idx_jp1 = idx + nx      # (i, j+1)
            idx_jm1 = idx - nx      # (i, j-1)

            # Atualizar psi usando diferenças finitas centradas
            psi[idx] =
                factor * (
                    (psi[idx_ip1] + psi[idx_im1]) / dx^2 +
                    (psi[idx_jp1] + psi[idx_jm1]) / dy^2 +
                    q[idx]
                )
        end
    end

    return psi
end


function solve_poisson_with_edges(nx, ny, dx, dy, q, psi)
    # nx: número de pontos em x
    # ny: número de pontos em y
    # dx: espaçamento da malha em x
    # dy: espaçamento da malha em y
    # q: vetor linearizado correspondente a -omega
    # psi: vetor linearizado inicial de psi

    # Pré-cálculo de fatores devido a dx e dy diferentes
    dx2 = dx^2
    dy2 = dy^2
    factor = 1 / (2 * (1 / dx2 + 1 / dy2))

    # Iterar sobre todos os pontos da malha
    for j = 1:ny
        for i = 1:nx
            idx = (j - 1) * nx + i  # Índice linear para (i, j)

            # Atualizar psi usando diferenças finitas, com tratamento especial nas bordas
            if i == 1  # Borda esquerda
                if j == 1  # Canto inferior esquerdo
                    psi[idx] =
                        factor * (
                            (psi[idx+1] - psi[idx]) / dx2 +  # Vizinho à direita
                            (psi[idx+nx] - psi[idx]) / dy2 +  # Vizinho acima
                            q[idx]
                        )
                elseif j == ny  # Canto superior esquerdo
                    psi[idx] =
                        factor * (
                            (psi[idx+1] - psi[idx]) / dx2 +  # Vizinho à direita
                            (psi[idx] - psi[idx-nx]) / dy2 +  # Vizinho abaixo
                            q[idx]
                        )
                else  # Borda esquerda (sem ser canto)
                    psi[idx] =
                        factor * (
                            (psi[idx+1] - psi[idx]) / dx2 +  # Vizinho à direita
                            (psi[idx+nx] + psi[idx-nx] - 2 * psi[idx]) / dy2 +  # Vizinhos acima e abaixo
                            q[idx]
                        )
                end
            elseif i == nx  # Borda direita
                if j == 1  # Canto inferior direito
                    psi[idx] =
                        factor * (
                            (psi[idx] - psi[idx-1]) / dx2 +  # Vizinho à esquerda
                            (psi[idx+nx] - psi[idx]) / dy2 +  # Vizinho acima
                            q[idx]
                        )
                elseif j == ny  # Canto superior direito
                    psi[idx] =
                        factor * (
                            (psi[idx] - psi[idx-1]) / dx2 +  # Vizinho à esquerda
                            (psi[idx] - psi[idx-nx]) / dy2 +  # Vizinho abaixo
                            q[idx]
                        )
                else  # Borda direita (sem ser canto)
                    psi[idx] =
                        factor * (
                            (psi[idx] - psi[idx-1]) / dx2 +  # Vizinho à esquerda
                            (psi[idx+nx] + psi[idx-nx] - 2 * psi[idx]) / dy2 +  # Vizinhos acima e abaixo
                            q[idx]
                        )
                end
            elseif j == 1  # Borda inferior (sem cantos)
                psi[idx] =
                    factor * (
                        (psi[idx+1] + psi[idx-1] - 2 * psi[idx]) / dx2 +  # Vizinhos à direita e esquerda
                        (psi[idx+nx] - psi[idx]) / dy2 +  # Vizinho acima
                        q[idx]
                    )
            elseif j == ny  # Borda superior (sem cantos)
                psi[idx] =
                    factor * (
                        (psi[idx+1] + psi[idx-1] - 2 * psi[idx]) / dx2 +  # Vizinhos à direita e esquerda
                        (psi[idx] - psi[idx-nx]) / dy2 +  # Vizinho abaixo
                        q[idx]
                    )
            else  # Região interna
                psi[idx] =
                    factor * (
                        (psi[idx+1] + psi[idx-1] - 2 * psi[idx]) / dx2 +  # Vizinhos à direita e esquerda
                        (psi[idx+nx] + psi[idx-nx] - 2 * psi[idx]) / dy2 +  # Vizinhos acima e abaixo
                        q[idx]
                    )
            end
        end
    end

    return psi
end



#________________________

function passo_não_linear(
    x,
    dad,
    H,
    G,
    M,
    dNx,
    dNy;
    Re = 1,
    relaxation = 0.01,
    erro_vel_min = 1e-6,
    maxiter = 1000,
)
    iter = 0
    erro_vel = 1
    nolinear = zeros(typeof(x[1]), 2 * (ni(dad) + nc(dad)))

    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)
    u = zeros(nc(dad) + ni(dad), 2)
    u[contorno, :], t = separa(dad, x)
    u[interno, 1] = x[2*nc(dad)+1:2:end]
    u[interno, 2] = x[2*nc(dad)+2:2:end]

    A, b = BEM.aplicaCDC(H, Re * G, dad)
    u_before = deepcopy(u)

    prog = BEM.ProgressThresh(erro_vel_min; desc = "Re: $Re; erro_vel =")

    while erro_vel > erro_vel_min && iter < maxiter

        dudx = dNx * u[:, 1]
        dudy = dNy * u[:, 1]

        dvdx = dNx * u[:, 2]
        dvdy = dNy * u[:, 2]

        nolinear[1:2:2*n] = u[:, 1] .* dudx + u[:, 2] .* dudy
        nolinear[2:2:2*n] = u[:, 1] .* dvdx + u[:, 2] .* dvdy

        x = A \ (b - Re * M * nolinear)


        u[contorno, :] = u[contorno, :] * (1 - relaxation) + relaxation * separa(dad, x)[1]
        t = t * (1 - relaxation) + relaxation * separa(dad, x)[2]

        u[interno, 1] = u[interno, 1] * (1 - relaxation) + relaxation * x[2*nc(dad)+1:2:end]
        u[interno, 2] = u[interno, 2] * (1 - relaxation) + relaxation * x[2*nc(dad)+2:2:end]

        erro_vel = nrmse(u_before, u)

        BEM.update!(prog, erro_vel)

        # println("Re: $Re; Iteração: $iter; erro_vel = $erro_vel")
        # @show norm(x)
        u_before = deepcopy(u)
        iter = iter + 1
    end
    x
end



function passo_não_linear2(
    x,
    dad,
    H,
    G,
    M,
    dNx,
    dNy;
    Re = 1,
    relaxation = 1,
    erro_vel_min = 1e-6,
    maxiter = 1000,
)
    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)
    u = zeros(typeof(x[1]), nc(dad) + ni(dad), 2)
    nolinear = zeros(typeof(x[1]), 2 * (ni(dad) + nc(dad)))
    dnl = zeros(typeof(x[1]), 2 * (ni(dad) + nc(dad)), 2 * (ni(dad) + nc(dad)))
    dnl2 = zeros(typeof(x[1]), 2 * (ni(dad) + nc(dad)))

    u[contorno, :], t = separa(dad, x)
    u[interno, 1] = x[2*nc(dad)+1:2:end]
    u[interno, 2] = x[2*nc(dad)+2:2:end]

    A, b = BEM.aplicaCDC(H, Re * G, dad)
    u_before = deepcopy(u)

    iter = 0
    erro_vel = 1
    prog = BEM.ProgressThresh(erro_vel_min; desc = "Re: $Re; erro_vel =")
    while erro_vel > erro_vel_min && iter < maxiter

        dudx = dNx * u[:, 1]
        dudy = dNy * u[:, 1]

        dvdx = dNx * u[:, 2]
        dvdy = dNy * u[:, 2]

        nolinear[1:2:2*n] = u[:, 1] .* dudx + u[:, 2] .* dudy
        nolinear[2:2:2*n] = u[:, 1] .* dvdx + u[:, 2] .* dvdy
        for i = 1:n
            if BEM.tipoCDC(dad)[i] .== 0
                continue
            end
            dnl[2i-1, 2i-1] = dudx[i]
            dnl[2i-1, 2i] = dvdx[i]
            dnl[2i, 2i-1] = dudy[i]
            dnl[2i, 2i] = dvdy[i]
        end

        Residuo = A * x - (b - Re * M * nolinear)
        jacobiano = A - Re * M * dnl
        x = x - jacobiano \ Residuo


        u[contorno, :] = u[contorno, :] * (1 - relaxation) + relaxation * separa(dad, x)[1]
        t = t * (1 - relaxation) + relaxation * separa(dad, x)[2]

        u[interno, 1] = u[interno, 1] * (1 - relaxation) + relaxation * x[2*nc(dad)+1:2:end]
        u[interno, 2] = u[interno, 2] * (1 - relaxation) + relaxation * x[2*nc(dad)+2:2:end]

        erro_vel = nrmse(u_before, u)

        BEM.update!(prog, erro_vel)

        # println("Re: $Re; Iteração: $iter; erro_vel = $erro_vel")
        # @show norm(x)
        u_before .= u
        iter = iter + 1
    end
    # Residuo = x - A \ (b - Re * M * nolinear)
    # jacobiano = I - A \ (Re * M * dnl)
    x
end