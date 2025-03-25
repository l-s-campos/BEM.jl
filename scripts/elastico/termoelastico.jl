## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
function testagao(i, metodo)
    nelem = i  #Numero de elementos
    NPY = 2i#pontos internos na direção y
    NPX = 2NPY#pontos internos na direção x
    npg = 10    #apenas números pares
    ## Formatação dos dados ________________________________________________
    # println("1. Formatando os dados")
    # dad = format_dad(telastico1d(nelem, 3), NPX, NPY) # dados
    dad = format_dad(telasticogao(nelem, 3), NPX, NPY) # dados
    # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
    # u, ss = termoelasticidade(dad, npg, k=2e-6, θ=1)#temperatura constante

    # u, ss = termoelasticidade(dad, npg, k=2e-6, θ=ones(nc(dad) + ni(dad)))#valores nos pontos
    c0 = 0
    c1 = -60
    c2 = 40
    f(x, y) = c2 * y^2 + c1 * y + c0

    EE = dad.k.E
    v = dad.k.nu
    k = dad.k.k
    # fc(x, y) = [x, y]
    # df(x, y) = kchap * [0, 2c2 * y^1 + c1]
    # q2, dq2 = BEM.calc_forçacorpo(dad, fc, 10)
    # dMpe = BEM.Monta_dM_RIMd(dad, npg)
    # dq2t = dMpe * vcat(fc.(x, y)...)
    # [dq2 dq2t]

    # F, Fx, Fy = BEM.montaFs([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])
    # dF = [zeros(size(Fx)); zeros(size(Fy))]
    # for i = 1:ni(dad)+nc(dad)
    #     dF[2i-1, :] = Fx[i, :]
    #     dF[2i, :] = Fy[i, :]
    # end
    uana(y) =
        (c1 / 2 * (y^2 - 0.5^2) + c2 / 3 * (y^3 + 0.5^3) + c0 * (y + 0.5)) * 1.3 / 0.7 * k
    tdad = @timed u, ss = termoelasticidade(dad, npg, θ = f, metodo = metodo)#função

    # x = [dad.NOS[:, 1]; dad.pontos_internos[:, 1]];
    y = [dad.NOS[:, 2]; dad.pontos_internos[:, 2]]
    ssa = -k / 0.7 * f.(0, y) * dad.k.E
    ua = uana.(y)
    eu = nrmse(ua, u[:, 2])
    es = nrmse(ssa, ss[:, 1])
    [nc(dad) ni(dad) eu es tdad.time]
end
function testacil(i, metodo)
    nelem = i  #Numero de elementos
    NPY = 5i#pontos internos na direção y
    NPX = NPY#pontos internos na direção x
    npg = 12    #apenas números pares
    ## Formatação dos dados ________________________________________________
    # println("1. Formatando os dados")
    # dad = format_dad(telastico1d(nelem, 3), NPX, NPY) # dados
    entradaele, entradapot = telasticocilindro(nelem, 3)
    dadel = format_dad(entradaele, NPX, NPY) # dados
    dadpot = format_dad(entradapot, NPX, NPY) # dados
    # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
    # u, ss = termoelasticidade(dad, npg, k=2e-6, θ=1)#temperatura constante
    # H, G = calc_HeG(dadpot, npg)  #importante
    # H1,G1 = calc_HeG(dad,10*npg)  #importante
    # @infiltrate
    A, b = calc_Aeb(dadpot) # Calcula a matriz A e o vetor b
    println("3. Resolvendo o sistema linear")
    x = A \ b
    println("4. Separando fluxo e temperatura")
    T, q = separa(dadpot, x) #importante
    Ti = calc_Ti(dadpot, T, q, npg)
    # u, ss = termoelasticidade(dad, npg, k=2e-6, θ=ones(nc(dad) + ni(dad)))#valores nos pontos
    temp = [T; Ti]
    # tdad = @timed u, ss = termoelasticidade(dadel, npg, θ=1, metodo=metodo)#função
    r(x, y) = sqrt(x^2 + y^2)
    Ti = 100
    Te = 0
    re = 80
    ri = 30

    EE = dadel.k.E
    v = dadel.k.nu
    ka = dadel.k.k
    k = 60

    # c = 0
    g(x, y) =
        1 / 2 * (Te * r(x, y)^2 - Ti * ri^2) +
        (Ti - Te) / (2 * log(ri / re)) *
        ((ri^2 - r(x, y)^2) / 2 + r(x, y)^2 * log(r(x, y) / re))

    ua(x, y) =
        ka * (1 + v) / (1 - v) / r(x, y) *
        (g(x, y) + ((1 - 2v) * r(x, y)^2 + ri^2) / (re^2 - ri^2) * g(0, re))
    ta(x, y) = Te + (Ti - Te) / log(ri / re) * log(r(x, y) / re)
    function dta(x, y)
        r_xy = r(x, y)

        dtadr = (Ti - Te) / log(ri / re) / r_xy

        drdx = x / r_xy
        drdy = y / r_xy

        dta_dy = dtadr * drdy
        dta_dx = dtadr * drdx

        return [dta_dx, dta_dy]
    end


    # ssa(x, y) = 2ka * EE * v / (1 - v) / (re^2 - ri^2) * g(0, re) - ka * EE * v / (1 - v) * ta(x, y)
    ssa(x, y) =
        ka * EE * (Ti - Te) / (1 - v) / 2 / log(re / ri) *
        (v - 2 * ri^2 * v / (re^2 - ri^2) * log(re / ri) - 2 * log(re / r(x, y)))
    tdad = @timed u, ss = termoelasticidade(dadel, npg, θ = ta, metodo = metodo)#função

    x = [dadpot.NOS[:, 1]; dadpot.pontos_internos[:, 1]]
    y = [dadpot.NOS[:, 2]; dadpot.pontos_internos[:, 2]]
    tensa = ssa.(x, y)
    dispa = ua.(x, y)
    tempa = ta.(x, y)
    # dtempa = reduce(hcat, dta.(x, y))'
    un = r.(u[:, 1], u[:, 2])
    ssn = v * (ss[:, 1] + ss[:, 2]) - ka * EE * tempa
    rs = r.(x, y)
    # @infiltrate
    # lines(rs, dispa)
    # etemp = nrmse(tempa, temp)
    es = nrmse(tensa, ssn)
    ed = nrmse(dispa, un)
    [nc(dadpot) ni(dadpot) ed es tdad.time]
end
function testapatch(i, metodo)
    nelem = i  #Numero de elementos
    NPY = 1i#pontos internos na direção y
    NPX = NPY#pontos internos na direção x
    npg = 10    #apenas números pares

    ## Formatação dos dados ________________________________________________
    # println("1. Formatando os dados")
    # dad = format_dad(telastico1d(nelem, 3), NPX, NPY) # dados
    dad1 = format_dad(telasticosquare(nelem, 3), NPX, NPY) # dados
    dad2 = format_dad(telasticosquare(nelem, 3, planestress = true), NPX, NPY) # dados
    # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
    # u, ss = termoelasticidade(dad, npg, k=2e-6, θ=1)#temperatura constante
    EE1 = dad1.k.E
    v1 = dad1.k.nu
    k1 = dad1.k.k
    # fc(x, y) = [x, y]
    # df(x, y) = kchap * [0, 2c2 * y^1 + c1]
    # q2, dq2 = BEM.calc_forçacorpo(dad, fc, 10)
    # dMpe = BEM.Monta_dM_RIMd(dad, npg)
    # dq2t = dMpe * vcat(fc.(x, y)...)
    # [dq2 dq2t]

    # F, Fx, Fy = BEM.montaFs([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])
    # dF = [zeros(size(Fx)); zeros(size(Fy))]
    # for i = 1:ni(dad)+nc(dad)
    #     dF[2i-1, :] = Fx[i, :]
    #     dF[2i, :] = Fy[i, :]
    # end
    tdad = @timed u, ss1 = termoelasticidade(dad1, npg, θ = 1, metodo = metodo)#função
    u, ss2 = termoelasticidade(dad2, npg, θ = 1, metodo = metodo)#função

    # x = [dad.NOS[:, 1]; dad.pontos_internos[:, 1]];
    # y = [dad.NOS[:, 2]; dad.pontos_internos[:, 2]]
    # ssa = -k / 0.7 * f.(0, y) * dad.k.E
    # ua = uana.(y)
    # eu = nrmse(ua, u[:, 2])
    # es = nrmse(ssa, ss[:, 1])
    # [nc(dad) ni(dad) eu es tdad.time]
    # @infiltrate
    estrain = nrmse(-6.3e6, ss1[:, 1])
    estress = nrmse(-3.6e6, ss2[:, 1])
    [nc(dad1) ni(dad1) ss1[1] ss2[1] estrain estress tdad.time]
end
function testabar(i, metodo)
    nelem = i  #Numero de elementos
    NPY = 2i#pontos internos na direção y
    NPX = 2NPY#pontos internos na direção x
    npg = 10    #apenas números pares

    ## Formatação dos dados ________________________________________________
    # println("1. Formatando os dados")
    # dad = format_dad(telastico1d(nelem, 3), NPX, NPY) # dados
    dad = format_dad(telasticobar(nelem, 3, planestress = true), NPX, NPY) # dados
    # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
    # u, ss = termoelasticidade(dad, npg, k=2e-6, θ=1)#temperatura constante
    EE = 210e9
    # if caso == 1
    # carga(x, y) = [1, 0]
    # ux(x) = x * (8 - x) / 2 / EE
    # sx(x) = (4 - x)
    # elseif caso == 2
    carga(x, y) = [x, 0]
    ux(x) = x * (4^2 - x^2 / 3) / 2 / EE
    sx(x) = (4^2 - x^2) / 2
    # end
    # end
    tdad = @timed u, ss = termoelasticidade(dad, npg, θ = 0, carga = carga, metodo = metodo)#função

    x = [dad.NOS[:, 1]; dad.pontos_internos[:, 1]]

    ua = ux.(x)
    sa = sx.(x)
    eu = nrmse(ua, u[:, 1])
    emu = nme(ua, u[:, 1])
    es = nrmse(sa, ss[:, 1])
    ems = nme(sa, ss[:, 1])
    # @infiltrate
    [nc(dad) ni(dad) eu es emu ems tdad.time]
end
# resdibem = [testapatch(i, "DIBEM") for i in 1:3]
# lines(r.(x, y), ssa)
# lines!(r.(x, y), ss[:, 1])
# resdibem = [testagao(i, "DIBEM") for i in 1:6]
# resrim = [testagao(i, "RIM") for i in 1:4]
# resdibem = [testacil(i, "DIBEM") for i in 1:6]
# resrim = [testacil(i, "RIM") for i in 1:4]
# res = testabar(8, "DIBEM")
# resrim = [testabar(i, "DIBEM") for i in 1:6]
# resrim = reduce(vcat, resrim)
# # cliparray(resdibem)
# testacil(10, "DIBEM")
# testagao(5, "DIBEM")
# cliparray(resrim)

resdibem = [testacil(i, "DIBEM") for i = 1:6]
resdibem = reduce(vcat, resdibem)
resrim = [testacil(i, "RIM") for i = 1:3]
resrim = reduce(vcat, resrim)

# lines(y, -k / 0.7 * f.(0, y) * dad.k.E)
# lines!(y, ss[:, 1])
# lines(y, uana.(y) * 1e5)
# lines!(y, u[:, 2] * 1e5)
