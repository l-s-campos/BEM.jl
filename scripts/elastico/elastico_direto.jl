## Início da análise
# using Pkg
# Pkg.activate(pwd())
# Pkg.instantiate()
using DataFrames
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
# nelem = 10  #Numero de elementos
# order = 3
NPX = 10 #pontos internos na direção x
NPY = 10 #pontos internos na direção y
npg = 16    #apenas números pares
# ## Formatação dos dados ________________________________________________

prob = [cilindro, viga]
ana = [ana_cilindro, ana_viga]
res = DataFrame(
    prob = String[],
    n = Int[],
    order = Int[],
    t = Float64[],
    td = Float64[],
    speedup = Float64[],
    tsolve = Float64[],
    tdsolve = Float64[],
    speedupsolve = Float64[],
    e1 = Float64[],
    e2 = Float64[],
    e3 = Float64[],
    e1d = Float64[],
    e2d = Float64[],
    e3d = Float64[],
)
# res = zeros(0, 4)
# testeinterno = [zeros(10) (1 .- (0.8 .^ (1:5:50)))]
for p = 1:2, i = 1:7, j = 3:3
    # p = 2
    # i = 1
    # j = 2
    nelem = 2 * 2^i  #Numero de elementos
    order = j

    dad = format_dad(prob[p](nelem, order), NPX, NPY) # dados
    # include(datadir("dadpotencial.jl"))

    if prob[p] == cilindro
        corrige_CDC_cilindro(dad)
    elseif prob[p] == viga
        corrige_CDC_viga(dad)
    elseif prob[p] == placa_furo
        corrige_CDC_placa_furo(dad)
    end
    # dad.pontos_internos = testeinterno
    # function compara_potencial(dad, npg)
    reset_timer!()
    tHeG = @timed H, G = calc_HeG(dad, npg)  #importante
    tHeGd = @timed Hd, Gd = BEM.calc_HeGd(dad, 3)  #importante
    A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
    Ad, bd = aplicaCDC(Hd, Gd, dad) # Calcula a matriz A e o vetor b
    tsolve = @timed x = A \ b
    tsolved = @timed xd = Ad \ bd
    u, t = separa(dad, x) #importante
    ud, td = separa(dad, xd) #importante

    tens_cont, tens_nt = calc_tens_cont(dad, u, t)
    tens_contd, tens_ntd = calc_tens_cont(dad, ud, td)
    tens_int = calc_tens_int(dad, u, t)
    tens_intd = calc_tens_int(dad, ud, td)

    ts = [tHeG.time tHeGd.time tsolve.time tsolved.time]

    if prob[p] == cilindro
        ana_u, ana_tr, ana_tt = ana[p](dad)
        ur = -sqrt.(u[:, 1] .^ 2 + u[:, 2] .^ 2)
        urd = -sqrt.(ud[:, 1] .^ 2 + ud[:, 2] .^ 2)
        e1 = nme(ana_u, ur)
        e2 = nme(ana_tr[1:order*nelem], tens_nt[1:order*nelem, 2])
        e3 = nme(ana_tt[1:order*nelem], tens_nt[1:order*nelem, 1])

        e1d = nme(ana_u, urd)
        e2d = nme(ana_tr[1:order*nelem], tens_ntd[1:order*nelem, 2])
        e3d = nme(ana_tt[1:order*nelem], tens_ntd[1:order*nelem, 1])

    elseif prob[p] == viga
        ux, uy, σx, τxy, ux_i, uy_i, σx_i, τxy_i = ana[p](dad)
        e1 = nme(ux, u[:, 1])
        e2 = nme(uy, u[:, 2])
        e3 = nme(σx_i, tens_int[:, 1])
        e1d = nme(ux, ud[:, 1])
        e2d = nme(uy, ud[:, 2])
        e3d = nme(σx_i, tens_intd[:, 1])
    elseif prob[p] == placa_furo
        σx, σy, τ = ana[p](dad)
        e1 = nme(σx, tens_cont[:, 1])
        e2 = nme(σy, tens_cont[:, 2])
        e3 = nme(τ, tens_cont[:, 3])
        e1d = nme(σx, tens_contd[:, 1])
        e2d = nme(σy, tens_contd[:, 2])
        e3d = nme(τ, tens_contd[:, 3])
    end
    # lines(σx_analitico)
    # lines!(tens_cont[:, 1])
    # lines!(tens_contd[:, 1])
    # BEM.current_figure()
    # lines(τxy)
    # lines!(tens_cont[:, 3])
    # lines!(tens_contd[:, 3])
    # BEM.current_figure()
    # # lines(ux_analitico)
    # lines!(u[:, 1])
    # lines!(ud[:, 1])
    # BEM.current_figure()
    # lines(uy_analitico)
    # lines!(u[:, 2])
    # lines!(ud[:, 2])
    # BEM.current_figure()


    push!(
        res,
        [
            string(prob[p]),
            size(u, 1),
            order,
            ts[1],
            ts[2],
            ts[1] / ts[2],
            ts[3],
            ts[4],
            ts[3] / ts[4],
            e1,
            e2,
            e3,
            e1d,
            e2d,
            e3d,
        ],
    )
end
using CSV
for i in prob
    CSV.write(string(i, ".csv"), filter(row -> (row.prob == string(i)), res))
end