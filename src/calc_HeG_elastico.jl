function calc_HeG(dad::elastico, npg = 8; interno = false, subcheia = false)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    H = zeros(2 * n, 2 * n)
    G = zeros(2 * n, 2 * n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    if temsubregioes(dad)
        regs = regiao_NO(dad)
    end
    @showprogress "Montando H e G" for i = 1:n
        pf = dad.NOS[i, :]   # Coordenada (x,y)  dos pontos fonte
        for elem_j in dad.ELEM  #Laço dos elementos

            if temsubregioes(dad)# matriz esparsa
                if subcheia
                    if i in dad.k.subregioes.equivale[:, 1:2]
                        if regs[i] != elem_j.regiao
                            continue
                        end
                    end
                else
                    if regs[i] != elem_j.regiao
                        continue
                    end
                end
            end

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
            eta, Jt = sinhtrans(qsi, eet, b)
            # eta,Jt=telles(qsi,eet)
            # @infiltrate
            h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
            # h, g = integraelem(pf, x, qsi, w , elem_j, dad)
            # nosing = elem_j.indices .== i

            # @infiltrate
            # @infiltrate isnan(sum(h))
            cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
            H[2i-1:2i, cols] = h
            G[2i-1:2i, cols] = g
        end
    end

    for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
        H[2i-1:2i, 2i-1:2i] .= 0
        H[2i-1:2i, 2i-1:2i] =
            -[sum(H[2i-1:2i, 1:2:end], dims = 2) sum(H[2i-1:2i, 2:2:end], dims = 2)]
        # H[2i-1, 2i] = -sum(H[2i-1, 2:2:end])
        # H[2i, 2i-1] = -sum(H[2i, 1:2:end])
        # H[2i-1, 2i-1] += 0.5
        # H[2i, 2i] += 0.5
    end
    if interno
        ni = size(dad.pontos_internos, 1)
        Hi = zeros(2 * ni, 2 * n)
        Gi = zeros(2 * ni, 2 * n)
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
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                # @infiltrate
                Hi[2i-1:2i, cols] = h
                Gi[2i-1:2i, cols] = g
            end
        end
        return [H zeros(2n, 2ni); Hi I], [G; Gi]
    else
        return H, G
    end
end

function calc_HeG(dad::elastico_aniso, npg = 8; interno = false)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    H = zeros(2 * n, 2 * n)
    G = zeros(2 * n, 2 * n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    for i = 1:n
        pf = dad.NOS[i, :]   # Coordenada (x,y)  dos pontos fonte
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
            nosing = elem_j.indices .== i
            if sum(nosing) == 1
                no_pf = findfirst(nosing)
                xi0 = elem_j.ξs[no_pf]
                Hesing = calc_Hsing(xi0, x, dad, elem_j, qsi, w)
                h[:, 2*no_pf-1:2*no_pf] = Hesing
            end
            cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
            H[2i-1:2i, cols] = h
            G[2i-1:2i, cols] = g
        end
    end
    if interno
        ni = size(dad.pontos_internos, 1)
        Hi = zeros(2 * ni, 2 * n)
        Gi = zeros(2 * ni, 2 * n)
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
                cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
                # @infiltrate
                Hi[2i-1:2i, cols] = h
                Gi[2i-1:2i, cols] = g
            end
        end
        return [H zeros(2n, 2ni); Hi I], [G; Gi]
    else
        return H, G
    end
end


function calc_Hsing(xi0, x, dad, elem, xi, w)
    mi = dad.k.mi
    A = dad.k.A
    q = dad.k.q
    g = dad.k.g

    N, dN = calc_fforma(xi0, elem)
    x0 = N' * x    # Ponto fonte
    npontos = length(xi)
    integral = log(1 - xi0) - log(1 + xi0)
    Hsing = zeros(Complex, 2, 2)
    Hesing = zeros(Float64, 2, 2)
    for i = 1:npontos
        delta = xi[i] - xi0
        kernelsing = (1 / delta)
        N, dN = calc_fforma(xi[i], elem)
        # @infiltrate
        if (xi0 ≈ -2 / 3)
            Nno = N[1]
            no = 1
        elseif (xi0 ≈ 0)
            Nno = N[2]
            no = 2
        else
            Nno = N[3]
            no = 3
        end
        pg = N' * x    # Ponto de gauss interpolador
        dxdqsi = dN' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        nx = sy
        ny = -sx

        # Compute the distance from the source point (xf,yf) to the field point (xcampo, ycampo)
        z1 = pg[1] - x0[1] + mi[1] * (pg[2] - x0[2])
        mi_n_z1 = (mi[1] * nx - ny) / z1
        z2 = pg[1] - x0[1] + mi[2] * (pg[2] - x0[2])
        mi_n_z2 = (mi[2] * nx - ny) / z2
        kernel1 = mi_n_z1 * dgamadqsi * Nno
        kernel2 = mi_n_z2 * dgamadqsi * Nno
        kernelreg1 = kernel1 - kernelsing
        kernelreg2 = kernel2 - kernelsing
        Hsing = Hsing + [kernelreg1 0; 0 kernelreg2] * w[i]
    end
    Hsing = Hsing + [integral 0; 0 integral]
    Hesing = 2 * real(A * Hsing * conj(g)') + [1/2 0; 0 1/2]
    return Hesing
end






#__________________________________________________________________________________________________________
"Funcao para calcular fazer integracao no contorno "
function integraelem(pf, x, eta, w, elem, dad::Union{elastico,elastico_aniso})
    h = zeros(Float64, 2, 2 * size(elem))
    g = zeros(Float64, 2, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x    # Ponto de gauss interpolador
        dxdqsi = dN' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        uast, tast = calsolfund(pg', pf, [sy, -sx], dad, elem.regiao)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        h += tast * Nm * dgamadqsi * w[k]
        g += uast * Nm * dgamadqsi * w[k]
        # @infiltrate

    end
    h, g
end



function calc_Aeb(dad::Union{elastico,elastico_aniso}, npg = 8)
    nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    A = zeros(2 * n, 2 * n)
    B = zeros(2 * n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    for i = 1:n
        pf = [dad.NOS; dad.pontos_internos][i, :]   # Coordenada (x,y)  dos pontos fonte
        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
            eet =
                (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) /
                norm(Δelem)^2 + elem_j.ξs[1]
            N = calc_fforma(eet, elem_j, false)
            ps = N_geo' * x
            b = norm(ps' - pf) / norm(Δelem)
            eta, Jt = sinhtrans(qsi, eet, b)
            # eta,Jt=telles(qsi,eet)
            h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
            # falta a integral singular
            for ii = 1:2
                if elem_j.tipoCDC[ii] == 1
                    A[2i-1:2i, 2elem_j.indices.-(2-ii)] = h[:, ii:2:end]
                    B[2i-1:2i] = g[:, ii:2:end] * elem_j.valorCDC[ii, :]
                else
                    A[2i-1:2i, 2elem_j.indices] = -g[:, ii:2:end]
                    B[2i-1:2i] += -h[:, ii:2:end] * elem_j.valorCDC[ii, :]
                end
            end
        end
    end
    A, B
end

function calc_HeG(dad::Union{elastico_iga,elastico_aniso_iga}, npg = 8)
    # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    H = zeros(2 * n, 2 * n)
    G = zeros(2 * n, 2 * n)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    for elem_j in dad.ELEM  #Laço dos elementos
        xf = elem_j.limites[:, 2]
        x0 = elem_j.limites[:, 1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
        Δelem = xf - x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
        pc = dad.pontos_controle[:, elem_j.indices]
        # @infiltrate
        cf = pc[1:2, :] ./ pc[4, :]'
        # for i =1 : 1
        for i = 1:n
            pf = dad.NOS[i, :]   # Coordenada (x,y)  dos pontos fonte
            eet = 2 * dot(Δelem, pf .- x0) / norm(Δelem)^2 - 1
            N, dN = calc_fforma(eet, elem_j, pc[4, :])
            ps = cf * N
            b = norm(ps - pf) / norm(Δelem)
            eta, Jt = sinhtrans(qsi, eet, b)
            h, g = integrabezier(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad)
            if i in elem_j.sing
                # @infiltrate
                ind = findfirst(i .== elem_j.sing)
                eet = elem_j.eets[ind]
                eta, Jt = sinhtrans(qsi, eet, b)
                h = integrabeziersing(pf, cf, pc[4, :], eta, w .* Jt, elem_j, dad, eet)
            end
            cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
            # @infiltrate
            H[2i-1:2i, cols] = h
            G[2i-1:2i, cols] = g
        end
    end


    H[1:2:end, 1:2:end] += dad.E / 2
    H[2:2:end, 2:2:end] += dad.E / 2
    H, G
end


function integrabeziersing(pf, cf, we, eta, w, elem::bezier, dad::elastico_iga, eet)
    h = zeros(Float64, 2, 2 * size(elem))
    basisrc, dN = calc_fforma(eet, elem, we)
    dxdqsi = cf * dN   # dx/dξ & dy/dξ
    dgamadqsif = norm(dxdqsi) / 2  # dΓ/dξ = J(ξ) Jacobiano
    matbasisrc = zeros(2, 2 * size(elem))
    matbasisrc[1, 1:2:end] = basisrc
    matbasisrc[2, 2:2:end] = basisrc

    E, ν = dad.Ev[1], dad.Ev[2]


    hterm = [
        0 -(1 - 2 * ν)/(4*pi*(1-ν))
        (1-2*ν)/(4*pi*(1-ν)) 0
    ]
    htermMatrix = hterm * matbasisrc

    Nm = zeros(Float64, 2, 2 * size(elem))


    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem, we)
        pg = cf * N    # Ponto de gauss interpolador
        dxdqsi = cf * dN   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        uast, tast = calsolfund(pg, pf, [sy, -sx], dad)
        # h+=N*dgamadqsi*w[k]
        # g+=N*dgamadqsi*w[k]
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        # @infiltrate
        h += (tast * Nm * dgamadqsi / 2 - htermMatrix / (eta[k] - eet)) * w[k]
    end
    # @show h
    if abs(eet) == 1
        beta_m = 1 / dgamadqsif
        h += htermMatrix * log(abs(2 / beta_m)) * sign(-eet)
        #        println("h = $(htermMatrix*log(abs(2/beta_m))*sign(-eet))")
    else
        # @infiltrate
        h += htermMatrix * log(abs((1 - eet) / (1 + eet)))
        # @show htermMatrix * log(abs((1 - eet) / (1 + eet)))
    end
    h
end




function integrabezier(
    pf,
    cf,
    we,
    eta,
    w,
    elem::bezier,
    dad::Union{elastico_iga,elastico_aniso_iga},
)
    h = zeros(Float64, 2, 2 * size(elem))
    g = zeros(Float64, 2, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))

    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem, we)
        pg = cf * N    # Ponto de gauss interpolador
        #   r = pg-pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = cf * dN   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        uast, tast = calsolfund(pg, pf, [sy, -sx], dad)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        h += tast * Nm * dgamadqsi * w[k] / 2
        g += uast * Nm * dgamadqsi * w[k] / 2

    end
    h, g
end


function calsolfund(pg, pf, n, dad::Union{elastico,elastico_iga}, regiao = 0)
    # @infiltrate
    if regiao == 0
        E, v = dad.k.E, dad.k.nu
    else
        E, v = dad.k.E[regiao], dad.k.nu[regiao]
    end
    r = pg - pf      # Distancia entre ponto de gauss e ponto fonte
    # @infiltrate
    GE = E / (2 * (1 + v))
    # Distance of source and field points
    R = norm(r)

    # Components of the unity vector in the radial direction
    r1 = r / R

    # Plane elasticity fundamental solutions
    prod1 = 4 * pi * (1 - v)
    prod2 = (3 - 4 * v) * log(1 / R)
    prod3 = dot(r1, n)

    u11 = (prod2 + r1[1]^2) / (2 * prod1 * GE)
    u22 = (prod2 + r1[2]^2) / (2 * prod1 * GE)
    u12 = (r1[1] * r1[2]) / (2 * prod1 * GE)
    u21 = u12

    t11 = -(prod3 * ((1 - 2 * v) + 2 * r1[1]^2)) / (prod1 * R)
    t22 = -(prod3 * ((1 - 2 * v) + 2 * r1[2]^2)) / (prod1 * R)
    t12 =
        -((prod3 * 2 * r1[1] * r1[2]) - (1 - 2 * v) * (r1[1] * n[2] - r1[2] * n[1])) /
        (prod1 * R)
    t21 =
        -((prod3 * 2 * r1[1] * r1[2]) - (1 - 2 * v) * (r1[2] * n[1] - r1[1] * n[2])) /
        (prod1 * R)
    uast = [
        u11 u12
        u21 u22
    ]

    tast = [
        t11 t12
        t21 t22
    ]

    uast, tast
end

function calc_tens_cont(dad::elastico_iga, u, t)
    # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    tens_cont = zeros(n, 3)
    tens_nt = zeros(n, 3)
    uc = dad.E \ u
    for elem_j in dad.ELEM  #Laço dos elementos

        pc = dad.pontos_controle[:, elem_j.indices]
        # @infiltrate
        cf = pc[1:2, :] ./ pc[4, :]'
        ni = dad.Ev[2]
        E = dad.Ev[1]
        for i in elem_j.sing

            # pf = dad.NOS[i,:]   # Coordenada (x,y)  dos pontos fonte
            normal = dad.normal[i, :]
            lij = [
                normal
                -normal[2] normal[1]
            ]
            tr_loc = lij * t[i, :]
            s11 = tr_loc[1]
            s12 = tr_loc[2]
            du = uc[elem_j.indices, :] * dad.base[i, 2]
            dudqsi = lij * du
            e22 = dudqsi ./ dad.jacobiano[1]
            s22 = E * e22 / (1 - ni^2) + ni * (1 + ni) / (1 - ni^2) * s11
            s = [s11 s12; s12 s22]
            sigma = lij * s * lij'
            tens_cont[i, :] = [sigma[1, 1], sigma[2, 2], sigma[1, 2]]
            tens_nt[i, :] = [s11, s22, s12]  #normal,  tangente, normal tangente
        end
    end
    tens_cont
end

function calc_tens_cont(dad::elastico_aniso_iga, u, t)
    # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
    rigidez = dad.k.A3


    n = size(dad.NOS, 1)
    tens_nt = zeros(n, 3)
    tens_cont = zeros(n, 3)
    eps_cont = zeros(n, 3)
    uc = dad.E \ u
    for elem_j in dad.ELEM  #Laço dos elementos

        pc = dad.pontos_controle[:, elem_j.indices]
        # @infiltrate
        cf = pc[1:2, :] ./ pc[4, :]'
        for i in elem_j.sing

            # pf = dad.NOS[i,:]   # Coordenada (x,y)  dos pontos fonte
            normal = dad.normal[i, :]
            lij = [
                normal[1] normal[2]
                -normal[2] normal[1]
            ]
            cs = normal[1]
            si = normal[2]

            tr_loc = lij * t[i, :]
            s11 = tr_loc[1]
            s12 = tr_loc[2]

            Transf = [
                cs^2 si^2 2*cs*si
                si^2 cs^2 -2*cs*si
                -cs*si cs*si cs^2-si^2
            ]

            rigidez_loc = inv(Transf) * rigidez * inv(Transf')

            matrix = deepcopy(rigidez_loc)
            du = uc[elem_j.indices, :]' * dad.base[i, 2]
            dudqsi = lij * du
            e22 = dudqsi[2] ./ dad.base[i, 3][1]
            vector = [
                s11 - rigidez_loc[1, 2] * e22
                -rigidez_loc[2, 2] * e22
                s12 - rigidez_loc[3, 2] * e22
            ]
            matrix[:, 2] = [0 -1 0]'
            vector2 = inv(matrix) * vector
            e11 = vector2[1]
            s22 = vector2[2]
            e12 = vector2[3]

            s = [s11 s12; s12 s22]
            ep = [e11 e12; e12 e22]
            sigma = lij * s * lij'
            epsilon = lij * ep * lij'

            # @infiltrate
            tens_cont[i, :] = [sigma[1, 1], sigma[2, 2], sigma[1, 2]]
            eps_cont[i, :] = [epsilon[1, 1], epsilon[2, 2], epsilon[1, 2]]
            tens_nt[i, :] = [s11, s22, s12]  #normal,  tangente, normal tangente
        end
    end
    tens_cont, eps_cont, tens_nt
end

function calc_tens_cont(dad::elastico, uc, t)
    # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS, 1)
    tens_cont = zeros(n, 3)
    tens_nt = zeros(n, 3)
    for elem_j in dad.ELEM  #Laço dos elementos

        x = dad.NOS[elem_j.indices, :]
        ni = dad.k.nu
        E = dad.k.E
        for i = 1:size(elem_j.indices, 1)
            N, dN = calc_fforma(elem_j.ξs[i], elem_j)
            dxdqsi = dN' * x
            dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
            sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
            sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

            normal = [sy -sx]

            lij = [
                normal
                -normal[2] normal[1]
            ]
            tr_loc = lij * t[elem_j.indices[i], :]
            s11 = tr_loc[1]
            s12 = tr_loc[2]
            du = uc[elem_j.indices, :]' * dN

            dudqsi = lij * du
            e22 = dudqsi[2] ./ dgamadqsi
            # @infiltrate
            s22 = E * e22 + ni * s11
            # s22 = E * e22 / (1 - ni^2) + ni * (1 + ni) / (1 - ni^2) * s11
            s = [s11 s12; s12 s22]
            sigma = lij * s * lij'
            tens_cont[elem_j.indices[i], :] = [sigma[1, 1], sigma[2, 2], sigma[1, 2]]
            tens_nt[elem_j.indices[i], :] = [s11, s22, s12]  #normal,  tangente, normal tangente
        end
    end
    tens_cont, tens_nt
end

function calc_tens_cont(dad::elastico_aniso, uc, t)
    # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
    rigidez = dad.k.A3
    # rigidez = dad.k.rigidez

    n = size(dad.NOS, 1)
    tens_nt = zeros(n, 3)
    tens_cont = zeros(n, 3)
    eps_cont = zeros(n, 3)
    for elem_j in dad.ELEM  #Laço dos elementos

        x = dad.NOS[elem_j.indices, :]

        for i = 1:length(elem_j.indices)
            N, dN = calc_fforma(elem_j.ξs[i], elem_j)
            dxdqsi = dN' * x
            dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
            sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
            sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

            normal = [sy, -sx]

            lij = [
                normal[1] normal[2]
                -normal[2] normal[1]
            ]
            cs = normal[1]
            si = normal[2]

            tr_loc = lij * t[elem_j.indices[i], :]
            s11 = tr_loc[1]
            s12 = tr_loc[2]

            Transf = [
                cs^2 si^2 2*cs*si
                si^2 cs^2 -2*cs*si
                -cs*si cs*si cs^2-si^2
            ]

            rigidez_loc = inv(Transf) * rigidez * inv(Transf')

            matrix = deepcopy(rigidez_loc)
            du = uc[elem_j.indices, :]' * dN
            dudqsi = lij * du
            e22 = dudqsi[2] ./ dgamadqsi
            vector = [
                s11 - rigidez_loc[1, 2] * e22
                -rigidez_loc[2, 2] * e22
                s12 - rigidez_loc[3, 2] * e22
            ]
            matrix[:, 2] = [0 -1 0]'
            vector2 = inv(matrix) * vector
            e11 = vector2[1]
            s22 = vector2[2]
            e12 = vector2[3]

            s = [s11 s12; s12 s22]
            ep = [e11 e12; e12 e22]
            sigma = lij * s * lij'
            epsilon = lij * ep * lij'

            tens_cont[elem_j.indices[i], :] = [sigma[1, 1], sigma[2, 2], sigma[1, 2]]
            # @infiltrate
            eps_cont[elem_j.indices[i], :] = [epsilon[1, 1], epsilon[2, 2], epsilon[1, 2]]
            tens_nt[elem_j.indices[i], :] = [s11, s22, s12]  #normal,  tangente, normal tangente
        end
    end
    tens_cont, tens_nt
end


function calc_tens_int(dad::Union{elastico,elastico_aniso}, u, t, regiao = 0, npg = 8)
    n = size(dad.pontos_internos, 1)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    stress = zeros(n, 3)
    for i = 1:n
        pf = dad.pontos_internos[i, :]   # Coordenada (x,y)  dos pontos fonte
        for elem_j in dad.ELEM  #Laço dos elementos
            if temsubregioes(dad)
                if elem_j.regiao != regiao[i]
                    continue
                end
            end
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
            d, s = integradelem(pf, x, eta, w .* Jt, elem_j, dad)
            # nosing = elem_j.indices .== i
            cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
            # @infiltrate
            stress[i, :] += d * t'[:][cols] - s * u'[:][cols]
        end
    end

    # for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    #   H[2i-1:2i,2i-1:2i].=0
    #   H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
    # end
    stress
end

function calc_SeD(dad::Union{elastico,elastico_aniso}, npg = 8)
    np = nc(dad)    # Quantidade de elementos discretizados no contorno
    npi = ni(dad)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    n = npi + np
    S = zeros(3n, 2np)
    D = zeros(3n, 2np)
    for i = 1:n
        if i <= np
            pf = dad.NOS[i, :]   # Coordenada (x,y)  dos pontos fonte
        else
            pf = dad.pontos_internos[i-np, :]   # Coordenada (x,y)  dos pontos fonte
        end
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
            nosing = elem_j.indices .== i
            if sum(nosing) == 1
                no_pf = findfirst(nosing)
                xi0 = elem_j.ξs[no_pf]
                d, s = integradelemsing_num(pf, x, elem_j, dad, xi0, 100)
                # @infiltrate
            else
                d, s = integradelem(pf, x, eta, w .* Jt, elem_j, dad)
            end
            # nosing = elem_j.indices .== i
            cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
            # @infiltrate
            S[3i-2:3i, cols] = s
            D[3i-2:3i, cols] = d
        end
    end

    # for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    #   H[2i-1:2i,2i-1:2i].=0
    #   H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
    # end
    S, D
end

function calc_desl_int(dad::Union{elastico}, u, t, npg = 8)
    n = size(dad.pontos_internos, 1)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    desl = zeros(n, 2)
    for i = 1:n
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
            cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
            # @infiltrate
            desl[i, :] += g * t'[:][cols] - h * u'[:][cols]
        end
    end

    # for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    #   H[2i-1:2i,2i-1:2i].=0
    #   H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
    # end
    desl
end
function calc_desl_int(dad::Union{elastico_aniso}, u, t, npg = 8)
    n = size(dad.pontos_internos, 1)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    desl = zeros(n, 2)
    scale = dad.k.A3[3, 3]

    for i = 1:n
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
            cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
            # @infiltrate
            desl[i, :] += g * t'[:][cols] - h * u'[:][cols]
        end
    end
    # for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    #   H[2i-1:2i,2i-1:2i].=0
    #   H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
    # end
    desl
end
function integradelem(pf, x, eta, w, elem, dad::Union{elastico,elastico_aniso})
    d = zeros(Float64, 3, 2 * size(elem))
    s = zeros(Float64, 3, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x    # Ponto de gauss interpolador
        dxdqsi = dN' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        D, S = caldsolfund(pg', pf, [sy, -sx], dad, elem.regiao)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        D_2 = [
            D[1, 1, 1] D[2, 1, 1]
            D[1, 2, 2] D[2, 2, 2]
            D[1, 1, 2] D[2, 1, 2]
        ]

        S_2 = [
            S[1, 1, 1] S[2, 1, 1]
            S[1, 2, 2] S[2, 2, 2]
            S[1, 1, 2] S[2, 1, 2]
        ]
        # @infiltrate
        d += D_2 * Nm * dgamadqsi * w[k]
        s += S_2 * Nm * dgamadqsi * w[k]

    end
    d, s
end

function caldsolfund(pg, pf, n, dad::Union{elastico,elastico_iga}, regiao = 0)
    # @infiltrate
    if regiao == 0
        E, v = dad.k.E, dad.k.nu
    else
        E, v = dad.k.E[regiao], dad.k.nu[regiao]
    end
    R = pg - pf      # Distancia entre ponto de gauss e ponto fonte
    # @infiltrate
    GEL = E / (2 * (1 + v))
    # Distance of source and field points
    r = norm(R)

    # Components of the unity vector in the radial direction
    dr = R / r
    fat1 = 4 * pi * (1 - v)
    fat2 = 1 - 2 * v
    drdn = dot(dr, n)
    kro = [1 0; 0 1]  # Kronecker delta

    # Initialization of D and S tensors
    D = zeros(2, 2, 2)
    S = zeros(2, 2, 2)

    # Calculation of D and S tensors
    for k = 1:2
        for i = 1:2
            for j = 1:2
                # Calculation of D tensor
                d1 = fat2 * (kro[k, i] * dr[j] + kro[k, j] * dr[i] - kro[i, j] * dr[k])
                d2 = 2 * dr[i] * dr[j] * dr[k]
                D[k, i, j] = (d1 + d2) / (fat1 * r)

                # Calculation of S tensor
                t1 =
                    2 *
                    drdn *
                    (
                        fat2 * kro[i, j] * dr[k] +
                        v * (kro[i, k] * dr[j] + kro[j, k] * dr[i]) -
                        4 * dr[i] * dr[j] * dr[k]
                    )
                t2 = 2 * v * (n[i] * dr[j] * dr[k] + n[j] * dr[i] * dr[k])
                t3 =
                    fat2 *
                    (2 * n[k] * dr[i] * dr[j] + n[j] * kro[i, k] + n[i] * kro[j, k]) -
                    (1 - 4 * v) * n[k] * kro[i, j]
                S[k, i, j] = (t1 + t2 + t3) * 2 * GEL / (fat1 * r^2)
                # @infiltrate
            end
        end
        D[k, 2, 1] = D[k, 1, 2]
        S[k, 2, 1] = S[k, 1, 2]
    end
    return D, S
end


vonmises(tens) =
    sqrt.(tens[:, 1] .^ 2 - tens[:, 1] * tens[:, 2] + tens[:, 2] .^ 2 + 2 * tens[:, 3] .^ 2)


function Monta_M_RIM(dad::Union{elastico,elastico_aniso}, npg1 = 10)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi1, w1 = gausslegendre(npg1)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(2n_pontos, 2n_pontos)
    F = zeros(n_pontos, n_pontos)
    normal_fonte = calc_normais(dad)

    # Cálculo da matriz [F]
    @showprogress "Montando F" for i = 1:n_pontos
        xi = nodes[i, 1]
        yi = nodes[i, 2]
        for j = 1:n_pontos
            xj = nodes[j, 1]
            yj = nodes[j, 2]
            r = sqrt((xi - xj)^2 + (yi - yj)^2)
            F[i, j] = interpola(r)
        end
    end
    @showprogress "Montando M" for i = 1:n_pontos
        if i <= n_nos
            pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
            nf = normal_fonte[i, :]
            caso = "contorno"
        elseif i <= n_nos + n_noi
            pf = dad.pontos_internos[i-n_nos, :] # Coordenada (x,y)
            nf = zeros(2)
            caso = "interno"
        end
        pr = nodes[i, :]
        for j = 1:n_pontos #Laço dos pontos fontes
            pr = nodes[j, :]
            for el = 1:nelem
                elem_j = dad.ELEM[el]#Laço dos elementos
                x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
                m_el = calc_m(x, pf, pr, qsi1, w1, elem_j, dad)
                M[2i-1:2i, 2j-1:2j] += m_el
            end
        end
    end
    F2 = zeros(2n_pontos, 2n_pontos)
    for i = 1:n_pontos
        for j = 1:n_pontos
            F2[2i-1, 2j-1] = F[i, j]
            F2[2i, 2j] = F[i, j]
        end
    end
    # @infiltrate
    M = M / F2
    M
end

function calc_m(x, pf, pr, qsi1, w1, elem, dad::Union{elastico,elastico_aniso})

    npg = length(w1)
    m_el = zeros(2, 2)

    for i = 1:npg
        N, dN = calc_fforma(qsi1[i], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        ny = -dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        nx = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        r = pg' - pf
        m = calcula_F(pr, pf, pg, [nx, ny], qsi1, w1, dad)
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w1[i]
    end
    return m_el
end

function calcula_F(pr, pf, pg, n, qsi, w, dad::Union{elastico,elastico_aniso}) #
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(2, 2) # Inicializa a integral de F_area
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
        rline = norm(xc - pr)
        ro = r / 2 * (qsi[i] + 1)
        term, ~ = calsolfund(xc, pf, n, dad)
        # term = ones(2, 2)
        f = interpola(rline)
        # @infiltrate
        F = F + term * f * ro * drodqsi * w[i]# Integral de F_area
    end
    return F
end

function Monta_M_RIMd(dad::Union{elastico,elastico_aniso}, npg)
    n_nos = size(dad.NOS, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    # n_canto = size(dad.k.cantos, 1)


    n_pontos = n_nos + n_noi
    if haskey(dad.k, :cantos)
        nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    else
        nodes = [dad.NOS; dad.pontos_internos]
    end

    F = zeros(n_pontos, n_pontos)
    # D11 = zeros(n_pontos, n_pontos)
    D = zeros(2n_pontos, 2n_pontos)
    # Dx = zeros(2n_pontos, 2n_pontos)
    # Dy = zeros(2n_pontos, 2n_pontos)
    # dA = zeros(2n_pontos, 4n_pontos)
    M1 = zeros(n_pontos)
    # M11 = zeros(n_pontos)
    M2 = zeros(2n_pontos, 2)
    # M2x = zeros(2n_pontos, 2)
    # M2y = zeros(2n_pontos, 2)
    # normal_fonte = calc_normais(dad)
    # Cálculo da matriz [F]
    @showprogress "Montando F e D" for i = 1:n_pontos
        if i <= n_nos
            pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
            # nf = normal_fonte[i, :]
            caso = "contorno"
        elseif i <= n_nos + n_noi
            pf = dad.pontos_internos[i-n_nos, :] # Coordenada (x,y)
            caso = "interno"
        end
        for j = 1:n_pontos
            # @show i, j
            pr = nodes[j, :]
            r = pr - pf
            R = norm(r)
            # @infiltrate
            F[i, j] = interpola(R)
            if i == j
                continue
            end
            # D[2i-1:2i, 2j-1:2j] = [1 2; 3 4]
            D[2i-1:2i, 2j-1:2j], ~ = calsolfund(pr, pf, [0, 0], dad)
            # D[2i-1:2i, 2j-1:2j], Dx[2i-1:2i, 2j-1:2j], Dy[2i-1:2i, 2j-1:2j] = caldusolfund(pr, pf, nf, dad)
            # D11[i, j] = D[2i-1, 2j-1]
        end
        qsi, w = gausslegendre(npg)

        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            m_el, m_el1 = calc_md(x, pf, qsi, w, elem_j, dad)
            # @infiltrate
            # m_el, m_el1, m_el1x, m_el1y = calc_md(x, pf, nf, qsi, w, elem_j, dad, pre)
            M1[i] = M1[i] + m_el
            M2[2i-1:2i, :] = M2[2i-1:2i, :] + m_el1
            # M2[i] += m_el1[1]
            # M2x[2i-1:2i, :] = M2x[2i-1:2i, :] + m_el1x
            # M2y[2i-1:2i, :] = M2y[2i-1:2i, :] + m_el1y

        end
    end
    # @show size(M)
    # @show length(M)

    # F2 = zeros(2n_pontos, 2n_pontos)
    # M1d = zeros(2, 2n_pontos)
    # for i = 1:n_pontos
    #   for j = 1:n_pontos
    #     F2[2i-1, 2j-1] = F[i, j]
    #     F2[2i, 2j] = F[i, j]
    #   end
    #   M1d[1, 2i-1] = M1[i]
    #   M1d[2, 2i] = M1[i]

    # end

    aux = M1' / F
    aux2 = [aux; aux][:]'
    A = aux2 .* D
    # aux1 = M1d / F2
    # # aux = [aux aux]'
    # @infiltrate
    # A2 = zeros(2n_pontos, 2n_pontos)
    # A2[:, 1:2:end] = aux1[1, :] .* D[:, 1:2:end]
    # A2[:, 2:2:end] = aux1[2, :] .* D[:, 2:2:end]

    # A11 = aux .* D[1:2:end, 1:2:end]
    # A12 = aux .* D[1:2:end, 2:2:end]
    # A21 = aux .* D[2:2:end, 1:2:end]
    # A22 = aux .* D[2:2:end, 2:2:end]
    # Ax = aux .* Dx
    # Ay = aux .* Dy
    for i = 1:n_pontos#Laço dos pontos radiais
        A[2i-1:2i, 2i-1:2i] .= 0
        A[2i-1:2i, 2i-1:2i] =
            -[sum(A[2i-1:2i, 1:2:end], dims = 2) sum(A[2i-1:2i, 2:2:end], dims = 2)] +
            M2[2i-1:2i, :]
        # Ax[2i-1:2i, 2i-1:2i] .= 0
        # Ax[2i-1:2i, 2i-1:2i] = -[sum(Ax[2i-1:2i, 1:2:end], dims=2) sum(Ax[2i-1:2i, 2:2:end], dims=2)] + M2x[2i-1:2i, :]
        # Ay[2i-1:2i, 2i-1:2i] .= 0
        # Ay[2i-1:2i, 2i-1:2i] = -[sum(Ay[2i-1:2i, 1:2:end], dims=2) sum(Ay[2i-1:2i, 2:2:end], dims=2)] + M2y[2i-1:2i, :]
    end
    A
end
function calc_md(x, pf, qsi, w, elem, dad::Union{elastico,elastico_aniso})
    npg = length(w)
    m_el, m_el1 = 0, zeros(2, 2)
    for i = 1:npg
        N, dN_geo = calc_fforma(qsi[i], elem)
        pg = N' * x    # Ponto de gauss interpolador
        r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = dN_geo' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

        nx = sy # Componente x do vetor normal unit�rio
        ny = -sx # Componente y do vetor normal unit�rio
        # @infiltrate
        R = norm(r)
        m = int_interpolaρdρ(R)
        m1 = calcula_F1(pf, pg, qsi, w, dad)
        # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
        m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]

    end
    return m_el, m_el1
end

function calcula_F1(pf, pg, qsi, w, dad::Union{elastico,elastico_aniso}) #
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(2, 2) # Inicializa a integral de F_area
    # Fx = zeros(2, 2) # Inicializa a integral de F_area
    # Fy = zeros(2, 2) # Inicializa a integral de F_area
    # theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        ro = r / 2 * (qsi[i] + 1)
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)

        # term = [1 2; 3 4]
        term, ~ = calsolfund(xc, pf, [0, 0], dad)
        F = F + term * ro * drodqsi * w[i]# Integral de F_area
        # Fx = F + termx * ro * drodqsi * w[i]# Integral de F_area
        # Fy = F + termy * ro * drodqsi * w[i]# Integral de F_area
    end
    return F
end
function Monta_dM_RIM(dad::Union{elastico,elastico_aniso}, npg1 = 10)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi1, w1 = gausslegendre(npg1)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    M = zeros(3n_pontos, 2n_pontos)
    F = zeros(n_pontos, n_pontos)
    normal_fonte = calc_normais(dad)

    # Cálculo da matriz [F]
    @showprogress "Montando F" for i = 1:n_pontos
        xi = nodes[i, 1]
        yi = nodes[i, 2]
        for j = 1:n_pontos
            xj = nodes[j, 1]
            yj = nodes[j, 2]
            r = sqrt((xi - xj)^2 + (yi - yj)^2)
            F[i, j] = interpola(r)
        end
    end
    @showprogress "Montando M" for i = 1:n_pontos
        if i <= n_nos
            pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
            nf = normal_fonte[i, :]
            caso = "contorno"
        elseif i <= n_nos + n_noi
            pf = dad.pontos_internos[i-n_nos, :] # Coordenada (x,y)
            nf = zeros(2)
            caso = "interno"
        end
        pr = nodes[i, :]
        for j = 1:n_pontos #Laço dos pontos fontes
            pr = nodes[j, :]
            for el = 1:nelem
                elem_j = dad.ELEM[el]#Laço dos elementos
                x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
                m_el = calc_dm(x, pf, pr, qsi1, w1, elem_j, dad)
                M[3i-2:3i, 2j-1:2j] += m_el
            end
        end
    end
    F2 = zeros(2n_pontos, 2n_pontos)
    for i = 1:n_pontos
        for j = 1:n_pontos
            F2[2i-1, 2j-1] = F[i, j]
            F2[2i, 2j] = F[i, j]
        end
    end
    # @infiltrate
    M = M / F2
    M
end


function Monta_dM_RIMd(dad::Union{elastico,elastico_aniso}, npg = 10)
    n_nos = size(dad.NOS, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    # n_canto = size(dad.k.cantos, 1)

    n_pontos = n_nos + n_noi
    if haskey(dad.k, :cantos)
        nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    else
        nodes = [dad.NOS; dad.pontos_internos]
    end

    F = zeros(n_pontos, n_pontos)
    D = zeros(3n_pontos, 2n_pontos)

    M1 = zeros(n_pontos)
    M2 = zeros(3n_pontos, 2)

    normal_fonte = calc_normais(dad)
    # Cálculo da matriz [F]
    @showprogress "Montando F e D" for i = 1:n_pontos
        if i <= n_nos
            pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
            nf = normal_fonte[i, :]
            caso = "contorno"
        elseif i <= n_nos + n_noi
            pf = dad.pontos_internos[i-n_nos, :] # Coordenada (x,y)
            nf = zeros(2)
            caso = "interno"
        end
        for j = 1:n_pontos
            # @show i, j
            pr = nodes[j, :]
            r = pr - pf
            R = norm(r)
            # @infiltrate
            F[i, j] = interpola(R)
            if i == j
                continue
            end
            Dt, ~ = caldsolfund(pr, pf, [0, 0], dad)
            term = [
                Dt[1, 1, 1] Dt[2, 1, 1]
                Dt[1, 2, 2] Dt[2, 2, 2]
                Dt[1, 1, 2] Dt[2, 1, 2]
            ]
            D[3i-2:3i, 2j-1:2j] = term
        end
        qsi, w = gausslegendre(npg)

        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            # @infiltrate
            m_el, m_el1 = calc_dmd(x, pf, qsi, w, elem_j, dad)
            M1[i] = M1[i] + m_el
            M2[3i-2:3i, :] = M2[3i-2:3i, :] + m_el1
        end
    end

    aux = M1' / F
    aux = [aux; aux][:]'
    A = aux .* D

    for i = 1:n_pontos #Laço dos pontos radiais
        A[3i-2:3i, 2i-1:2i] =
            -[sum(A[3i-2:3i, 1:2:end], dims = 2) sum(A[3i-2:3i, 2:2:end], dims = 2)] +
            M2[3i-2:3i, :]
    end

    A
end

function calc_dmd(x, pf, qsi, w, elem, dad::elastico)
    npg = length(w)
    m_el, m_el1, m_el1rx, m_el1ry = 0, zeros(3, 2), zeros(3, 2), zeros(3, 2)

    for i = 1:npg
        N, dN_geo = calc_fforma(qsi[i], elem)
        pg = N' * x    # Ponto de gauss interpolador
        r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = dN_geo' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

        nx = sy # Componente x do vetor normal unit�rio
        ny = -sx # Componente y do vetor normal unit�rio
        # @infiltrate
        R = norm(r)
        m = int_interpolaρdρ(R)
        m1, m1x, m1y = calcula_dF1(pf, pg, qsi, w, dad)
        # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
        m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
        m_el1rx += dot([nx, ny], r) / norm(r)^2 * m1x * dgamadqsi * w[i]
        m_el1ry += dot([nx, ny], r) / norm(r)^2 * m1y * dgamadqsi * w[i]
    end
    return m_el, m_el1#, m_el1rx, m_el1ry
end
function calc_dm(x, pf, pr, qsi, w, elem, dad::elastico)
    npg = length(w)
    m_el1 = zeros(3, 2)

    for i = 1:npg
        N, dN_geo = calc_fforma(qsi[i], elem)
        pg = N' * x    # Ponto de gauss interpolador
        r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = dN_geo' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

        nx = sy # Componente x do vetor normal unit�rio
        ny = -sx # Componente y do vetor normal unit�rio
        # @infiltrate
        r = pg' - pf
        m1 = calcula_dF(pr, pf, pg, qsi, w, dad)
        # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
        # @infiltrate
        m_el1 += dot([nx, ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
    end
    return m_el1
end
function calcula_dF1(pf, pg, qsi, w, dad::elastico) #
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(3, 2) # Inicializa a integral de F_area
    Fx = zeros(3, 2) # Inicializa a integral de F_area
    Fy = zeros(3, 2) # Inicializa a integral de F_area
    # theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        ro = r / 2 * (qsi[i] + 1)
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)

        D, ~ = caldsolfund(xc, pf, [0, 0], dad)
        term = [
            D[1, 1, 1] D[2, 1, 1]
            D[1, 2, 2] D[2, 2, 2]
            D[1, 1, 2] D[2, 1, 2]
        ]
        # @infiltrate
        F = F + term * ro * drodqsi * w[i]# Integral de F_area
        Fx = F + term * ro * drodqsi * w[i] * R[1] / 2 * (qsi[i] + 1)# Integral de F_area
        Fy = F + term * ro * drodqsi * w[i] * R[2] / 2 * (qsi[i] + 1)# Integral de F_area
        # F3 = F3 + term3 * ro * drodqsi * w[i]# Integral de F_area
    end
    return F, Fx, Fy
end
function calcula_dF(pr, pf, pg, qsi, w, dad::Union{elastico,elastico_aniso}) #

    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    F = zeros(3, 2) # Inicializa a integral de F_area
    # theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
        ro = r / 2 * (qsi[i] + 1)

        rline = norm(xc - pr)
        f = interpola(rline)

        D, ~ = caldsolfund(pg', pf, [0, 1], dad)
        term = [
            D[1, 1, 1] D[2, 1, 1]
            D[1, 2, 2] D[2, 2, 2]
            D[1, 1, 2] D[2, 1, 2]
        ]

        # @infiltrate
        F = F + f * term * ro * drodqsi * w[i]# Integral de F_area
        # F3 = F3 + term3 * ro * drodqsi * w[i]# Integral de F_area
    end
    return F
end



function calc_forçacorpo(dad::Union{elastico,elastico_aniso}, f, npg1 = 10)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes

    qsi1, w1 = gausslegendre(npg1)
    n_pontos = n_nos + n_noi
    nodes = [dad.NOS; dad.pontos_internos]
    q = zeros(2n_pontos)
    dq = zeros(3n_pontos)
    # Cálculo da matriz [F]
    @showprogress "Montando q" for i = 1:n_pontos #Laço dos pontos fontes
        pf = nodes[i, :]
        for el = 1:nelem
            elem_j = dad.ELEM[el]#Laço dos elementos
            x = dad.NOS[elem_j.indices, :] # Coordenada (x,y) dos nós geométricos
            m_el, dm_el = calc_corpo(x, pf, qsi1, w1, elem_j, dad, f)
            q[2i-1:2i] += m_el
            dq[3i-2:3i] += dm_el
        end
    end
    q, dq
end
function calc_corpo(x, pf, qsi, w, elem, dad::Union{elastico,elastico_aniso}, f)

    npg = length(w)
    m_el = zeros(2)
    dm_el = zeros(3)

    for i = 1:npg
        N, dN = calc_fforma(qsi[i], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        ny = -dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        nx = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        R = pg' - pf
        r = norm(R)
        drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
        #    para qsi (dr/dqsi)
        m = zeros(2)
        dm = zeros(3)
        for i = 1:npg # Percorre os pontos de integra��o
            xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
            ro = r / 2 * (qsi[i] + 1)
            term, ~ = calsolfund(xc, pf, [0, 0], dad)
            D, ~ = caldsolfund(xc, pf, [0, 0], dad)
            dterm = [
                D[1, 1, 1] D[2, 1, 1]
                D[1, 2, 2] D[2, 2, 2]
                D[1, 1, 2] D[2, 1, 2]
            ]
            # term = ones(2, 2)
            # @infiltrate
            m = m + term * f(xc[1], xc[2]) * ro * drodqsi * w[i]# Integral de F_area
            dm = dm + dterm * f(xc[1], xc[2]) * ro * drodqsi * w[i]# Integral de F_area
        end
        # @infiltrate
        m_el += dot([nx, ny], R) / r^2 * m * dgamadqsi * w[i]
        dm_el += dot([nx, ny], R) / r^2 * dm * dgamadqsi * w[i]
    end
    return m_el, dm_el
end



function muda_nt!(dad, H, G)
    # matriz de mudança de coordenadas
    R = spzeros(2 * nc(dad), 2 * nc(dad))
    for i = 1:nc(dad)
        R[2i-1:2i, 2i-1:2i] =
            [dad.normal[i, 1] dad.normal[i, 2]; -dad.normal[i, 2] dad.normal[i, 1]]
    end
    H .= H * R
    G .= G * R
    return nothing
end

function muda_nt(dad, H, G)
    # matriz de mudança de coordenadas
    R = spzeros(2 * nc(dad), 2 * nc(dad))
    for i = 1:nc(dad)
        R[2i-1:2i, 2i-1:2i] =
            [dad.normal[i, 1] dad.normal[i, 2]; -dad.normal[i, 2] dad.normal[i, 1]]
    end
    return H * R, G * R
end
"""
verifica contato com superfície rígida
"""
function verifica_contato(x, h, dad)
    # Verifica a condição de contato de cada nó
    nnoscontato = size(h[1], 1)
    contato = zeros(Int, nnoscontato)
    nlinhas = size(x, 1) - size(h[1], 1) * 2

    for k = 1:nnoscontato  # Percorre todos os elementos que podem entrar
        no_contato = h[1][k]
        posux = 2 * no_contato - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy = 2 * no_contato        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx = nlinhas + 2 * k - 1   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty = nlinhas + 2 * k       # Posição da coluna da matriz A2 referente à força de superfície na direção t
        tn = -x[posty]
        tt = x[postx]
        un = -x[posuy]
        ut = x[posux]

        # tn = x[postx] * dad.normal[no_contato, 1] + x[posty] * dad.normal[no_contato, 2]
        # un = x[posux] * dad.normal[no_contato, 1] + x[posuy] * dad.normal[no_contato, 2]
        # Verifica se un é maior que o gap e se tn é de tração. Caso seja
        # verdade, o nó não se encontra em contato
        # @infiltrate
        # @show un, h[2][k], un - h[2][k]

        if un < h[2][k] - 1e-10 || tn > 0
            # Se entrou aqui o nó não encontra-se em contato
            contato[k] = 1         # O nó está livre
        elseif abs(tt / tn) + 1e-10 > dad.k.μ || ut * tt < 0
            # elseif dad.k.μ == 0 || (abs(tt / tn) + 1e-10 > dad.k.μ && ut * tt < 0)
            # Se entrou aqui o nó encontra-se em contato  e escorregamento
            contato[k] = sign(tt) * 2         # O nó está em contato
        else
            contato[k] = 3         # O nó está em contato e adesão
            # @show abs(tt / tn) / dad.k.μ
        end
    end
    # @show contato
    return contato
end
"""
aplica contato com superfície rígida
"""
function aplica_contato!(h, contato, A, b, dad)
    nnoscontato = size(contato, 1)  # Número de nós na região de possível contato
    nlinhas = size(b, 1) - size(h[1], 1) * 2

    # Imposição das condições de contato na matriz A2 e b2
    for k = 1:nnoscontato  # Percorre todos os elementos que podem entrar
        tipocontato = contato[k]  # tipo da condição de contato
        no_contato = h[1][k]         # número do nó em contato
        posux = 2 * no_contato - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy = 2 * no_contato        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx = nlinhas + 2 * k - 1   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty = nlinhas + 2 * k       # Posição da coluna da matriz A2 referente à força de superfície na direção t
        # @show posux, posuy, postx, posty

        A[[postx, posty], :] .= 0
        b[[postx, posty]] .= 0

        if tipocontato == 1  # Zona livre de contato
            A[posty, posty] = 1
            A[postx, postx] = 1

        elseif tipocontato == 3  # Zona em contato e adesão ut=0 un=gap
            A[postx, posuy] = 1#
            A[posty, posux] = 1
            b[postx] = -h[2][k]  # uy = -gap
        else #escorregamento
            A[posty, posty] = -dad.k.μ * sign(tipocontato)
            A[posty, postx] = 1
            A[postx, posuy] = -1
            b[postx] = h[2][k]
        end
    end
end

"""
verifica contato 2corpos
"""
function verifica_contato2(x, h, dad)
    # Verifica a condição de contato de cada nó
    nnoscontato = size(h[1], 1)
    contato = zeros(Int, nnoscontato)
    nlinhas = nc(dad) * 2

    for k = 1:nnoscontato  # Percorre todos os elementos que podem entrar
        no_contato = h[1][k]
        no_contato2 = h[2][k]

        posux = 2 * no_contato - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy = 2 * no_contato        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx = nlinhas + 4 * k - 3   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty = nlinhas + 4 * k - 2      # Posição da coluna da matriz A2 referente à força de superfície na direção t

        posux2 = 2 * no_contato2 - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy2 = 2 * no_contato2        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        # postx2 = nlinhas + 4 * k - 1   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        # posty2 = nlinhas + 4 * k
        # tx = x[postx]
        # ty = x[posty]
        # ux = x[posux]
        # uy = x[posuy]
        un = (x[posux] - x[posux2]) * h[4][k, 1] + (x[posuy] - x[posuy2]) * h[4][k, 2]
        ut = -(x[posux] - x[posux2]) * h[4][k, 2] + (x[posuy] - x[posuy2]) * h[4][k, 1]
        tn = (x[postx]) * h[4][k, 1] + (x[posty]) * h[4][k, 2]
        tt = -(x[postx]) * h[4][k, 2] + (x[posty]) * h[4][k, 1]

        # Verifica se un é maior que o gap e se tn é de tração. Caso seja
        # verdade, o nó não se encontra em contato
        if un < h[3][k] || tn > 0
            # @infiltrate
            # Se entrou aqui o nó não encontra-se em contato
            contato[k] = 1         # O nó está livre
        elseif dad.k.μ == 0 || (abs(tt / tn) - 1e-10 < dad.k.μ && ut * tt > 0)
            # Se entrou aqui o nó encontra-se em contato  e escorregamento
            contato[k] = sign(tt) * 2         # O nó está em contato
        else
            contato[k] = 3         # O nó está em contato e adesão

        end
    end
    # @show contato
    return contato
end
"""
aplica contato 2corpos
"""
function aplica_contato2!(h, contato, A, b, dad)
    nnoscontato = size(contato, 1)  # Número de nós na região de possível contato
    nlinhas = nc(dad) * 2

    # Imposição das condições de contato na matriz A2 e b2
    for k = 1:nnoscontato  # Percorre todos os elementos que podem entrar
        tipocontato = contato[k]  # tipo da condição de contato
        no_contato = h[1][k]         # número do nó em contato
        posux = 2 * no_contato - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy = 2 * no_contato        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx = nlinhas + 4 * k - 3   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty = nlinhas + 4 * k - 2      # Posição da coluna da matriz A2 referente à força de superfície na direção t
        # @show posux, posuy, postx, posty

        no_contato2 = h[2][k]         # número do nó em contato
        posux2 = 2 * no_contato2 - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy2 = 2 * no_contato2        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx2 = nlinhas + 4 * k - 1  # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty2 = nlinhas + 4 * k       # Posição da coluna da matriz A2 referente à força de superfície na direção t
        # @show [postx, posty, postx2, posty2]
        A[[postx, posty, postx2, posty2], :] .= 0
        b[[postx, posty, postx2, posty2]] .= 0
        if tipocontato == 1  # Zona livre de contato tx e ty=0
            A[postx, postx] = 1

            A[posty, posty] = 1

            A[postx2, postx2] = 1

            A[posty2, posty2] = 1

        elseif tipocontato == 3  # Zona em contato e adesão un = gap t1=t2 ut=0
            A[postx, postx] = 1
            A[postx, postx2] = 1


            A[posty, posty] = 1
            A[posty, posty2] = 1


            A[postx2, posux] = h[4][k, 1]
            A[postx2, posuy] = h[4][k, 2]
            A[postx2, posux2] = -h[4][k, 1]
            A[postx2, posuy2] = -h[4][k, 2]
            b[postx2] = h[3][k]  # un = gap

            A[posty2, posux] = -h[4][k, 2]#m
            A[posty2, posuy] = h[4][k, 1]
            A[posty2, posux2] = h[4][k, 2]
            A[posty2, posuy2] = -h[4][k, 1]

        else #escorregamento un = gap t1=t2 tt=μtn
            A[postx, postx] = 1
            A[postx, postx2] = 1

            A[posty, posty] = 1
            A[posty, posty2] = 1

            A[postx2, posux] = h[4][k, 1]
            A[postx2, posuy] = h[4][k, 2]
            A[postx2, posux2] = -h[4][k, 1]
            A[postx2, posuy2] = -h[4][k, 2]
            b[postx2] = h[3][k]  # un = gap

            # A[posty2, postx] = -h[4][k, 2]
            # A[posty2, posty] = h[4][k, 2]
            # A[posty2, postx] = -dad.k.μ * sign(tipocontato) * h[4][k, 1]
            # A[posty2, posty] = -dad.k.μ * sign(tipocontato) * h[4][k, 2]

            A[posty2, postx] = -h[4][k, 2] - dad.k.μ * sign(tipocontato) * h[4][k, 1]
            A[posty2, posty] = h[4][k, 1] + dad.k.μ * sign(tipocontato) * h[4][k, 2]

        end
    end
end

function Contato_NL(x0, p)
    A2, b2, h, dad = p
    # Verifica a condição de contato de cada nó
    contato = verifica_contato(x0, h, dad)
    # Aplica a condição de contato em cada nó
    aplica_contato!(h, contato, A2, b2, dad)
    A2 * x0 - b2
end

function Contato_NL_newton(dad, x0, A2, b2, h; maxiter = 10, tol = 1e-8)
    #Tolerância para o erro para parar do método de Newton
    x = deepcopy(x0)
    for i = 1:maxiter
        # Verifica a condição de contato de cada nó
        contato = verifica_contato(x0, h, dad)
        # @show contato
        # @infiltrate
        # Aplica a condição de contato em cada nó
        aplica_contato!(h, contato, A2, b2, dad)
        # y0 = A2 * x0 - b2
        # x = x0 - A2 \ y0
        x = A2 \ b2
        e = norm(x - x0) / norm(x)
        x0 = x
        if e < tol
            return x0
        end
    end
    x0
end
"""
Resolve numericamente um problema não-linear com condições de contato utilizando o método de Newton.
# Argumentos:
    dad: Estrutura de dados
    x0: Chute inicial para a solução.
    A2: Matriz do BEM
    b2: Vetor do BEM
    h: variável com as informações sobre o contato
    maxiter (opcional): Número máximo de iterações do método de Newton.
    tol (opcional): Tolerância para o erro, utilizada como critério de parada.

# Retorno:
    x: Solução aproximada do problema, ou o último chute se o critério de parada não for satisfeito.
"""
function Contato_NL_newton2(dad, x0, A2, b2, h; maxiter = 10, tol = 1e-8)
    #Tolerância para o erro para parar do método de Newton
    x = deepcopy(x0)
    for i = 1:maxiter
        # Verifica a condição de contato de cada nó
        contato = verifica_contato2(x0, h, dad)
        # @infiltrate
        # Aplica a condição de contato em cada nó
        aplica_contato2!(h, contato, A2, b2, dad)
        # y0 = A2 * x0 - b2
        # x = x0 - A2 \ y0
        x = A2 \ b2
        @show e = norm(x - x0) / norm(x)
        x0 = x
        if e < tol
            return x0
        end
    end
    x0
end





function Contato_NL_newton_incremental(
    dad,
    x0,
    A2,
    b2,
    h,
    maxiter,
    tol,
    npassos,
    nosrestritos,
    contato_certo,
)
    #Tolerância para o erro para parar do método de Newton
    nnos = size(dad.NOS)[1]
    nelem = size(dad.ELEM)[1]
    u = zeros(nnos, 2)
    t = zeros(nnos, 2)
    bdiv = b2 / npassos
    for i = 1:npassos

        # Resolve a equação não linear usando o método de Newton


        x = newton_passo_de_carga(
            dad,
            x0,
            A2,
            bdiv,
            h,
            u,
            t,
            maxiter,
            tol,
            i,
            contato_certo,
        )

        # Reordena deslocamentos e forças de superfície de acordo com as condições
        # de contorno
        u_dt, t_dt = separa(dad, x, nosrestritos, h)
        u = u + u_dt
        t = t + t_dt

        x0 = deepcopy(x)

    end
    u, t
end
function newton_passo_de_carga(dad, x0, A2, b2, h, u, t, maxiter, tol, ii, contato_certo)
    #Tolerância para o erro para parar do método de Newton
    x = deepcopy(x0)
    # print("\n newton_passo_de_carga\n")
    for i = 1:maxiter
        # Verifica a condição de contato de cada nó
        contato = verifica_contato_incremental(x0, h, dad, u, t, i)
        if ii == 2
            contato = contato_certo
        end        # print("\n Contato = ", contato)
        # Aplica a condição de contato em cada nó
        aplica_contato_incremental!(h, contato, A2, b2, dad, u, t, i)
        y0 = A2 * x0 - b2
        x = x0 - A2 \ y0
        e = norm(x - x0) / norm(x)
        # x = A2 \ b2
        x0 = x
        if e < tol
            return x0
        end
    end
    x0
end

"""
aplica contato com superfície rígida
"""
function aplica_contato_incremental!(h, contato, A, b, dad, u, t, i)
    nnoscontato = size(contato, 1)  # Número de nós na região de possível contato
    nlinhas = size(b, 1) - size(h[1], 1) * 2

    # Imposição das condições de contato na matriz A2 e b2
    for k = 1:nnoscontato  # Percorre todos os elementos que podem entrar
        tipocontato = contato[k]  # tipo da condição de contato
        no_contato = h[1][k]         # número do nó em contato
        posux = 2 * no_contato - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy = 2 * no_contato        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx = nlinhas + 2 * k - 1   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty = nlinhas + 2 * k       # Posição da coluna da matriz A2 referente à força de superfície na direção t

        if (i == 2)
            # @show posux, posuy, postx, posty
            # @infiltrate
        end


        A[[postx, posty], :] .= 0
        b[[postx, posty]] .= 0

        if tipocontato == 1  # Zona livre de contato
            A[posty, posty] = 1
            A[postx, postx] = 1

        elseif tipocontato == 3  # Zona em contato e adesão ut=0 un=gap
            A[postx, posuy] = 1#
            A[posty, posux] = 1
            b[postx] = -h[2][k] - u[no_contato, 2] # un = gap
        # b[postx] = -h[2][k]  # uy = -gap

        # b[posty] = u[no_contato, 1]  # un = gap
        else #escorregamento
            A[posty, posty] = -dad.k.μ * sign(tipocontato)
            A[posty, postx] = 1
            A[postx, posuy] = 1
            b[postx] = -h[2][k] - u[no_contato, 2]
            b[posty] = t[no_contato, 2] * dad.k.μ * sign(tipocontato) - t[no_contato, 1]
            # b[postx] = h[2][k]


        end
    end
end


"""
verifica contato com superfície rígida
"""
function verifica_contato_incremental(x, h, dad, u, t, i)
    # Verifica a condição de contato de cada nó
    nnoscontato = size(h[1], 1)
    contato = zeros(Int, nnoscontato)
    nlinhas = size(x, 1) - size(h[1], 1) * 2#nlinhas = nc(dad) * 2

    for k = 1:nnoscontato  # Percorre todos os elementos que podem entrar
        no_contato = h[1][k]
        posux = 2 * no_contato - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy = 2 * no_contato        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx = nlinhas + 2 * k - 1   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty = nlinhas + 2 * k       # Posição da coluna da matriz A2 referente à força de superfície na direção t
        deltatn = -x[posty]
        deltatt = x[postx]
        deltaun = -x[posuy]
        deltaut = x[posux]

        tn = -t[no_contato, 2] + deltatn
        tt = t[no_contato, 1] + deltatt
        un = -u[no_contato, 2] + deltaun
        ut = u[no_contato, 1] + deltaut



        # tn = x[postx] * dad.normal[no_contato, 1] + x[posty] * dad.normal[no_contato, 2]
        # un = x[posux] * dad.normal[no_contato, 1] + x[posuy] * dad.normal[no_contato, 2]
        # @infiltrate
        # Verifica se un é maior que o gap e se tn é de tração. Caso seja
        # verdade, o nó não se encontra em contato
        if un < h[2][k] || tn > 0
            # Se entrou aqui o nó não encontra-se em contato
            contato[k] = 1         # O nó está livre
        elseif dad.k.μ == 0 || abs(tt / tn) + 1e-10 > dad.k.μ || ut * tt < 0
            # elseif abs(tt / tn) + 1e-10 > dad.k.μ
            # Se entrou aqui o nó encontra-se em contato  e escorregamento
            contato[k] = sign(tt) * 2         # O nó está em contato
        else
            contato[k] = 3         # O nó está em contato e adesão

        end
    end

    return contato
end



function Contato_NL_newton_incremental2(
    dad,
    x0,
    A2,
    b2,
    h,
    maxiter,
    tol,
    npassos,
    nosrestritos,
)
    #Tolerância para o erro para parar do método de Newton
    nnos = size(dad.NOS)[1]
    nelem = size(dad.ELEM)[1]
    u = zeros(nnos, 2)
    t = zeros(nnos, 2)
    for i = 1:npassos

        # Resolve a equação não linear usando o método de Newton


        x = newton_passo_de_carga2(dad, x0, A2, b2, h, u, t, maxiter, tol, i)

        # Reordena deslocamentos e forças de superfície de acordo com as condições
        # de contorno
        u_dt, t_dt = separa(dad, x, nosrestritos, h)
        u = u + u_dt
        t = t + t_dt

        x0 = deepcopy(x)

    end
    u, t
end
function newton_passo_de_carga2(dad, x0, A2, b2, h, u, t, maxiter, tol, i)
    #Tolerância para o erro para parar do método de Newton
    x = deepcopy(x0)
    print("\n newton_passo_de_carga\n")
    for i = 1:maxiter
        # Verifica a condição de contato de cada nó
        contato = verifica_contato_incremental2(x0, h, dad, u, t)

        print("\n Contato = ", contato)

        # Aplica a condição de contato em cada nó
        aplica_contato_incremental2!(h, contato, A2, b2, dad, u, t)
        y0 = A2 * x0 - b2
        x = x0 - A2 \ y0
        e = norm(x - x0) / norm(x)
        # x = A2 \ b2
        x0 = x
        if e < tol
            return x0
        end
    end
    x0
end




"""
aplica contato com superfície rígida
"""
function aplica_contato_incremental2!(h, contato, A, b, dad, u, t)
    nnoscontato = size(contato, 1)  # Número de nós na região de possível contato
    nlinhas = nc(dad) * 2
    # Imposição das condições de contato na matriz A2 e b2
    for k = 1:nnoscontato  # Percorre todos os elementos que podem entrar
        tipocontato = contato[k]  # tipo da condição de contato
        no_contato = h[1][k]         # número do nó em contato
        posux = 2 * no_contato - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy = 2 * no_contato        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx = nlinhas + 4 * k - 3   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty = nlinhas + 4 * k - 2      # Posição da coluna da matriz A2 referente à força de superfície na direção t
        # @show posux, posuy, postx, posty

        no_contato2 = h[2][k]         # número do nó em contato
        posux2 = 2 * no_contato2 - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy2 = 2 * no_contato2        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx2 = nlinhas + 4 * k - 1  # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty2 = nlinhas + 4 * k       # Posição da coluna da matriz A2 referente à força de superfície na direção t


        A[[postx, posty, postx2, posty2], :] .= 0
        b[[postx, posty, postx2, posty2]] .= 0


        if tipocontato == 1  # Zona livre de contato
            A[postx, postx] = 1

            A[posty, posty] = 1

            A[postx2, postx2] = 1

            A[posty2, posty2] = 1

        elseif tipocontato == 3  # Zona em contato e adesão ut=0 un=gap
            A[postx, postx] = 1
            A[postx, postx2] = 1

            A[posty, posty] = 1
            A[posty, posty2] = 1

            A[postx2, posux] = h[4][k, 1]
            A[postx2, posuy] = h[4][k, 2]
            A[postx2, posux2] = -h[4][k, 1]
            A[postx2, posuy2] = -h[4][k, 2]
            b[postx2] = h[3][k] + u[no_contato, 2] # un = gap

            A[posty2, posux] = -h[4][k, 2]#m
            A[posty2, posuy] = h[4][k, 1]
            A[posty2, posux2] = h[4][k, 2]
            A[posty2, posuy2] = -h[4][k, 1]
            b[posty2] = u[no_contato, 1]  # un = gap
        else #escorregamento
            A[postx, postx] = 1
            A[postx, postx2] = 1

            A[posty, posty] = 1
            A[posty, posty2] = 1

            A[postx2, posux] = h[4][k, 1]
            A[postx2, posuy] = h[4][k, 2]
            A[postx2, posux2] = -h[4][k, 1]
            A[postx2, posuy2] = -h[4][k, 2]
            b[postx2] = h[3][k] + u[no_contato, 2] # un = gap

            A[posty2, postx] = -h[4][k, 2] - dad.k.μ * sign(tipocontato) * h[4][k, 1]
            A[posty2, posty] = h[4][k, 1] + dad.k.μ * sign(tipocontato) * h[4][k, 2]

            b[posty] = -t[no_contato, 2] * dad.k.μ * sign(tipocontato) - t[no_contato, 1]

        end
    end
end


"""
verifica contato com superfície rígida
"""
function verifica_contato_incremental2(x, h, dad, u, t)
    # Verifica a condição de contato de cada nó
    nnoscontato = size(h[1], 1)
    contato = zeros(Int, nnoscontato)
    nlinhas = nc(dad) * 2

    for k = 1:nnoscontato  # Percorre todos os elementos que podem entrar
        no_contato = h[1][k]
        no_contato2 = h[2][k]

        posux = 2 * no_contato - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy = 2 * no_contato        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        postx = nlinhas + 4 * k - 3   # Posição da coluna da matriz A2 referente à força de superfície na direção n
        posty = nlinhas + 4 * k - 2      # Posição da coluna da matriz A2 referente à força de superfície na direção t

        posux2 = 2 * no_contato2 - 1    # Posição da coluxa da matriz A2 referente ao deslocamento na direção n
        posuy2 = 2 * no_contato2        # Posição da coluna da matriz A2 referente ao deslocamento na direção t
        # @infiltrate
        # @show posux, posux2, posuy, posuy2, postx, posty

        deltaun = (x[posux] - x[posux2]) * h[4][k, 1] + (x[posuy] - x[posuy2]) * h[4][k, 2]
        deltaut = -(x[posux] - x[posux2]) * h[4][k, 2] + (x[posuy] - x[posuy2]) * h[4][k, 1]
        deltatn = (x[postx]) * h[4][k, 1] + (x[posty]) * h[4][k, 2]
        deltatt = -(x[postx]) * h[4][k, 2] + (x[posty]) * h[4][k, 1]

        un0 =
            (u[no_contato, 1] - u[no_contato2, 1]) * h[4][k, 1] +
            (u[no_contato, 2] - u[no_contato2, 2]) * h[4][k, 2]
        ut0 =
            -(u[no_contato, 1] - u[no_contato2, 1]) * h[4][k, 2] +
            (u[no_contato, 2] - u[no_contato2, 2]) * h[4][k, 1]
        tn0 = (t[no_contato, 1]) * h[4][k, 1] + (t[no_contato, 2]) * h[4][k, 2]
        tt0 = -(t[no_contato, 1]) * h[4][k, 2] + (t[no_contato, 2]) * h[4][k, 1]

        tn = tn0 + deltatn
        tt = tt0 + deltatt
        un = un0 + deltaun
        ut = ut0 + deltaut

        # @infiltrate
        # Verifica se un é maior que o gap e se tn é de tração. Caso seja
        # verdade, o nó não se encontra em contato
        if un < h[3][k] || tn > 0
            # Se entrou aqui o nó não encontra-se em contato
            contato[k] = 1         # O nó está livre
        elseif dad.k.μ == 0 || (abs(tt / tn) - 1e-10 < dad.k.μ && ut * tt > 0)
            # Se entrou aqui o nó encontra-se em contato  e escorregamento
            contato[k] = sign(tt) * 2         # O nó está em contato
        else
            contato[k] = 3         # O nó está em contato e adesão

        end
    end

    return contato
end


################# Função calc_derivX_solfund ###############################
function calc_deriv_solfund(pg, pf, n, dad::Union{elastico,elastico_iga}, regiao = 0)
    # @infiltrate
    if regiao == 0
        E, ni = dad.k.E, dad.k.nu
    else
        E, ni = dad.k.E[regiao], dad.k.nu[regiao]
    end
    r = pg - pf      # Distancia entre ponto de gauss e ponto fonte

    #Calcula as derivadas em x das soluções fundamentais

    mi = E / (2 * (1 + ni))

    # Distance of source and field points
    r1 = r[1]
    r2 = r[2]
    R = norm(r)
    nx = n[1]
    ny = n[2]

    # Derivada em x das soluções fundamentais

    C1 = 1 / (8 * pi * (1 - ni) * mi)

    u11x = (C1 * (2 * r1 - (2 * r1^3) / R^2 + r1 * (4 * ni - 3))) / R^2
    u12x = (C1 * (r2 - (2 * r1^2 * r2) / R^2)) / R^2
    u21x = (C1 * (r2 - (2 * r1^2 * r2) / R^2)) / R^2
    u22x = (C1 * (r1 * (4 * ni - 3) - (2 * r1 * r2^2) / R^2)) / R^2

    C2 = -1 / (4 * pi * (1 - ni))
    C3 = (1 / R) * (r1 * nx + r2 * ny)

    t11x =
        (
            C2 * (
                (2 * nx * r1^2) / R^2 - nx * (2 * ni - 1) +
                (2 * C3 * (2 * r1 - (4 * r1^3) / R^2 + r1 * (2 * ni - 1))) / R
            )
        ) / R^2
    t12x =
        (
            C2 * (
                ny * (2 * ni - 1) -
                (2 * ((2 * ni - 1) * (ny * r1^2 - nx * r2 * r1) - nx * r1 * r2)) / R^2 +
                (2 * C3 * (r2 - (4 * r1^2 * r2) / R^2)) / R
            )
        ) / R^2
    t21x =
        (
            C2 * (
                (2 * ((2 * ni - 1) * (ny * r1^2 - nx * r2 * r1) + nx * r1 * r2)) / R^2 -
                ny * (2 * ni - 1) + (2 * C3 * (r2 - (4 * r1^2 * r2) / R^2)) / R
            )
        ) / R^2
    t22x =
        (
            C2 * (
                (2 * C3 * (r1 * (2 * ni - 1) - (4 * r1 * r2^2) / R^2)) / R -
                nx * (2 * ni - 1) + (2 * nx * r2^2) / R^2
            )
        ) / R^2

    C1 = 1 / (8 * pi * (1 - ni) * mi)

    u11y = (C1 * (r2 * (4 * ni - 3) - (2 * r1^2 * r2) / R^2)) / R^2
    u12y = (C1 * (r1 - (2 * r1 * r2^2) / R^2)) / R^2
    u21y = (C1 * (r1 - (2 * r1 * r2^2) / R^2)) / R^2
    u22y = (C1 * (2 * r2 - (2 * r2^3) / R^2 + r2 * (4 * ni - 3))) / R^2

    C2 = -1 / (4 * pi * (1 - ni))
    C3 = (1 / R) * (r1 * nx + r2 * ny)

    t11y =
        (
            C2 * (
                (2 * C3 * (r2 * (2 * ni - 1) - (4 * r1^2 * r2) / R^2)) / R -
                ny * (2 * ni - 1) + (2 * ny * r1^2) / R^2
            )
        ) / R^2
    t12y =
        (
            C2 * (
                (2 * ((2 * ni - 1) * (nx * r2^2 - ny * r1 * r2) + ny * r1 * r2)) / R^2 -
                nx * (2 * ni - 1) + (2 * C3 * (r1 - (4 * r1 * r2^2) / R^2)) / R
            )
        ) / R^2
    t21y =
        (
            C2 * (
                nx * (2 * ni - 1) -
                (2 * ((2 * ni - 1) * (nx * r2^2 - ny * r1 * r2) - ny * r1 * r2)) / R^2 +
                (2 * C3 * (r1 - (4 * r1 * r2^2) / R^2)) / R
            )
        ) / R^2
    t22y =
        (
            C2 * (
                (2 * ny * r2^2) / R^2 - ny * (2 * ni - 1) +
                (2 * C3 * (2 * r2 - (4 * r2^3) / R^2 + r2 * (2 * ni - 1))) / R
            )
        ) / R^2

    # Assembly of matrices that contain fundamental solutions.
    uasty = [
        u11y u12y
        u21y u22y
    ]

    tasty = [
        t11y t12y
        t21y t22y
    ]

    # Assembly of matrices that contain fundamental solutions.
    uastx = [
        u11x u12x
        u21x u22x
    ]

    tastx = [
        t11x t12x
        t21x t22x
    ]

    return uastx, tastx, uasty, tasty
end


function calc_dHedG(dad::Union{elastico,elastico_aniso}, npg = 8)
    np = nc(dad)    # Quantidade de elementos discretizados no contorno
    npi = ni(dad)
    qsi, w = gausslegendre(npg)    # Quadratura de gauss
    n = npi + np
    Hx = zeros(2n, 2np)
    Hy = zeros(2n, 2np)
    Gx = zeros(2n, 2np)
    Gy = zeros(2n, 2np)
    for i = 1:n
        if i <= np
            pf = dad.NOS[i, :]   # Coordenada (x,y)  dos pontos fonte
        else
            pf = dad.pontos_internos[i-np, :]   # Coordenada (x,y)  dos pontos fonte
        end
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
            nosing = elem_j.indices .== i
            if sum(nosing) == 1
                no_pf = findfirst(nosing)
                xi0 = elem_j.ξs[no_pf]
                hx, gx, hy, gy = integra_deriv_elemsing_num(pf, x, elem_j, dad, xi0, 100)
                # @infiltrate
            else
                hx, gx, hy, gy = integra_deriv_elem(pf, x, eta, w .* Jt, elem_j, dad)
            end
            # nosing = elem_j.indices .== i
            cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
            # @infiltrate
            if i <= np
                Hx[2i-1:2i, cols] = hx * 2
                Gx[2i-1:2i, cols] = gx * 2
                Hy[2i-1:2i, cols] = hy * 2
                Gy[2i-1:2i, cols] = gy * 2
            else
                Hx[2i-1:2i, cols] = hx
                Gx[2i-1:2i, cols] = gx
                Hy[2i-1:2i, cols] = hy
                Gy[2i-1:2i, cols] = gy

            end
        end
    end

    # for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    #   H[2i-1:2i,2i-1:2i].=0
    #   H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
    # end
    Hx, Gx, Hy, Gy
end


function integra_deriv_elem(pf, x, eta, w, elem, dad::Union{elastico,elastico_aniso})
    hx = zeros(Float64, 2, 2 * size(elem))
    gx = zeros(Float64, 2, 2 * size(elem))

    hy = zeros(Float64, 2, 2 * size(elem))
    gy = zeros(Float64, 2, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    for k = 1:size(w, 1)
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x    # Ponto de gauss interpolador
        dxdqsi = dN' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        Ux, Tx, Uy, Ty = calc_deriv_solfund(pg', pf, [sy, -sx], dad, elem.regiao)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N

        # @infiltrate
        hx += Tx * Nm * dgamadqsi * w[k]
        gx += Ux * Nm * dgamadqsi * w[k]
        hy += Ty * Nm * dgamadqsi * w[k]
        gy += Uy * Nm * dgamadqsi * w[k]

    end
    hx, gx, hy, gy
end



function Monta_deriv_M_RIMd(dad::Union{elastico,elastico_aniso}, npg = 10)
    n_nos = size(dad.NOS, 1)
    n_noi = size(dad.pontos_internos, 1) #Number of internal nodes
    # n_canto = size(dad.k.cantos, 1)

    n_pontos = n_nos + n_noi
    if haskey(dad.k, :cantos)
        nodes = [dad.NOS; dad.pontos_internos; dad.k.cantos[:, 2:3]]
    else
        nodes = [dad.NOS; dad.pontos_internos]
    end

    F = zeros(n_pontos, n_pontos)
    Dx = zeros(2n_pontos, 2n_pontos)
    Dy = zeros(2n_pontos, 2n_pontos)

    M1 = zeros(n_pontos)
    M2x = zeros(2n_pontos, 2)
    M2y = zeros(2n_pontos, 2)

    # Cálculo da matriz [F]
    @showprogress "Montando F e D" for i = 1:n_pontos
        if i <= n_nos
            pf = dad.NOS[i, :] # Coordenada (x,y)dos pontos fonte
            nf = dad.normal[i, :]
            caso = "contorno"
        elseif i <= n_nos + n_noi
            pf = dad.pontos_internos[i-n_nos, :] # Coordenada (x,y)
            nf = zeros(2)
            caso = "interno"
        end
        for j = 1:n_pontos
            # @show i, j
            pr = nodes[j, :]
            r = pr - pf
            R = norm(r)
            # @infiltrate
            F[i, j] = interpola(R)
            if i == j
                continue
            end
            Ux, Uy = calc_deriv_solfund(pr, pf, [0, 0], dad)[[1, 3]]

            Dx[2i-1:2i, 2j-1:2j] = Ux
            Dy[2i-1:2i, 2j-1:2j] = Uy
        end
        qsi, w = gausslegendre(npg)

        for elem_j in dad.ELEM  #Laço dos elementos
            x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
            # @infiltrate
            m_el, m_el1x, m_el1y = calc_deriv_md(x, pf, qsi, w, elem_j, dad)
            M1[i] = M1[i] + m_el
            M2x[2i-1:2i, :] = M2x[2i-1:2i, :] + m_el1x
            M2y[2i-1:2i, :] = M2y[2i-1:2i, :] + m_el1y
        end
    end

    aux = M1' / F
    aux = [aux; aux][:]'
    Ax = aux .* Dx
    Ay = aux .* Dy

    for i = 1:n_pontos #Laço dos pontos radiais
        Ax[2i-1:2i, 2i-1:2i] =
            -[sum(Ax[2i-1:2i, 1:2:end], dims = 2) sum(Ax[2i-1:2i, 2:2:end], dims = 2)] +
            M2x[2i-1:2i, :]
        Ay[2i-1:2i, 2i-1:2i] =
            -[sum(Ay[2i-1:2i, 1:2:end], dims = 2) sum(Ay[2i-1:2i, 2:2:end], dims = 2)] +
            M2y[2i-1:2i, :]
    end

    Ax, Ay
end


function calc_deriv_md(x, pf, qsi, w, elem, dad::elastico)
    npg = length(w)
    m_el, m_el1x, m_el1y = 0, zeros(2, 2), zeros(2, 2)

    for i = 1:npg
        N, dN_geo = calc_fforma(qsi[i], elem)
        pg = N' * x    # Ponto de gauss interpolador
        r = pg' - pf      # Distancia entre ponto de gauss e ponto fonte
        dxdqsi = dN_geo' * x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ

        nx = sy # Componente x do vetor normal unit�rio
        ny = -sx # Componente y do vetor normal unit�rio
        # @infiltrate
        R = norm(r)
        m = int_interpolaρdρ(R)
        # m1x, m1y = calcula_deriv_F(pf, pg, qsi, w, dad)
        m1x, m1y = intradial_deriv_solfund(pg', pf, dad)
        # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
        # @infiltrate
        m_el += dot([nx, ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
        m_el1x += dot([nx, ny], r) / norm(r)^2 * m1x * dgamadqsi * w[i]
        m_el1y += dot([nx, ny], r) / norm(r)^2 * m1y * dgamadqsi * w[i]

    end
    return m_el, m_el1x, m_el1y
end

function calcula_deriv_F(pf, pg, qsi, w, dad::elastico) #
    npg = length(w)
    R = (pg' - pf)
    r = norm(R)
    drodqsi = r / 2 # Jacobiano da transforma��o de vari�vel de r
    #    para qsi (dr/dqsi)
    Fx = zeros(2, 2) # Inicializa a integral de F_area
    Fy = zeros(2, 2) # Inicializa a integral de F_area
    # theta = atan(R[2], R[1])
    for i = 1:npg # Percorre os pontos de integra��o
        ro = r / 2 * (qsi[i] + 1)
        xc = pf + R / 2 * (qsi[i] + 1) # ro=ro(qsi)
        Ux, Uy = calc_deriv_solfund(xc, pf, [0, 0], dad)[[1, 3]]

        # @infiltrate
        Fx = Fx + Ux * ro * drodqsi * w[i]# Integral de F_area
        Fy = Fy + Uy * ro * drodqsi * w[i]# Integral de F_area
    end
    return Fx, Fy
end



################# Função calc_derivX_solfund ###############################
function intradial_deriv_solfund(pg, pf, dad::Union{elastico,elastico_iga}, regiao = 0)
    # @infiltrate
    if regiao == 0
        E, ni = dad.k.E, dad.k.nu
    else
        E, ni = dad.k.E[regiao], dad.k.nu[regiao]
    end
    r = pg - pf      # Distancia entre ponto de gauss e ponto fonte

    #Calcula as derivadas em x das soluções fundamentais

    mi = E / (2 * (1 + ni))

    # Distance of source and field points
    r1 = r[1]
    r2 = r[2]
    R = norm(r)

    # Derivada em x das soluções fundamentais

    C1 = 1 / (8 * pi * (1 - ni) * mi)

    u11x = (C1 * (2 * r1 - (2 * r1^3) / R^2 + r1 * (4 * ni - 3)))
    u12x = (C1 * (r2 - (2 * r1^2 * r2) / R^2))
    u21x = (C1 * (r2 - (2 * r1^2 * r2) / R^2))
    u22x = (C1 * (r1 * (4 * ni - 3) - (2 * r1 * r2^2) / R^2))

    u11y = (C1 * (r2 * (4 * ni - 3) - (2 * r1^2 * r2) / R^2))
    u12y = (C1 * (r1 - (2 * r1 * r2^2) / R^2))
    u21y = (C1 * (r1 - (2 * r1 * r2^2) / R^2))
    u22y = (C1 * (2 * r2 - (2 * r2^3) / R^2 + r2 * (4 * ni - 3)))

    # Assembly of matrices that contain fundamental solutions.
    uasty = [
        u11y u12y
        u21y u22y
    ]

    # Assembly of matrices that contain fundamental solutions.
    uastx = [
        u11x u12x
        u21x u22x
    ]
    return uastx, uasty
end
