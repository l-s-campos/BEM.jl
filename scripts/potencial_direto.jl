## Início da análise
# using Pkg
# Pkg.activate(pwd())
# Pkg.instantiate()
using DataFrames, CSV
using DrWatson
# @quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 20  #Numero de elementos
order = 3
NPX = 10 #pontos internos na direção x
NPY = 10 #pontos internos na direção y
npg = 16    #apenas números pares
# ## Formatação dos dados ________________________________________________
# println("1. Formatando os dados");
# dad = format_dad(potencial1d(nelem, 2), NPX, NPY) # dados
# # dad0 = format_dad(potencial1d(nelem),NPX,NPY,0) # dados
# # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
# res = zeros(5, 4, 4, 10)
prob = [quarto_circ, placa_moulton, potencial1d, dad_laquini1, dad_laquini2, dad_laquini3]
ana = [ana_quarto_circ, T_ana_moulton, T_ana_1d, fa1, fa2, fa3]
res = DataFrame(
    prob = String[],
    n = Int[],
    order = Int[],
    t = Float64[],
    td = Float64[],
    speedup = Float64[],
    tsolve = Float64[],
    ti = Float64[],
    tid = Float64[],
    speedupi = Float64[],
    erms = Float64[],
    ermsd = Float64[],
    em1 = Float64[],
    em2 = Float64[],
)
# res = zeros(0, 4)
# testeinterno = [zeros(10) (1 .- (0.8 .^ (1:5:50)))]
for p = 1:1, i = 1:5, j = 2:5
    nelem = 4 * 2^i  #Numero de elementos
    order = j
    # p = 1

    dad = format_dad(prob[p](nelem, order), NPX, NPY) # dados
    # include(datadir("dadpotencial.jl"))

    if prob[p] == placa_moulton
        dad = corrigeCDC_moulton(dad)
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
    T, q = separa(dad, x) #importante
    Td, qd = separa(dad, xd) #importante
    # [tHeG.time tHeGd.time tsolve.time tsolved.time], T, Td, q, qd
    # end
    ti = @timed Ti = calc_Ti(dad, T, q, npg)
    tid = @timed Tid = BEM.calc_Tid(dad, Td, qd, 3)

    t = [tHeG.time tHeGd.time tsolve.time tsolved.time ti.time tid.time]
    # t, T, Td, q, qd = compara_potencial(dad, npg)
    Tana = ana[p](dad.NOS)
    e1 = norm(T - Tana) / norm(Tana)
    e2 = norm(Td - Tana) / norm(Tana)
    em1 = sum(abs, T - Tana) / sum(abs, Tana)
    em2 = sum(abs, Td - Tana) / sum(abs, Tana)

    Tiana = ana[p](dad.pontos_internos)
    ei1 = norm(Ti - Tiana) / norm(Tiana)
    ei2 = norm(Tid - Tiana) / norm(Tiana)
    eim1 = sum(abs, Ti - Tiana) / sum(abs, Tiana)
    eim2 = sum(abs, Tid - Tiana) / sum(abs, Tiana)
    # global res = abs.([res; [testeinterno[:, 2] .- 1 Tiana (Ti - Tiana) ./ Tiana (Tid - Tiana) ./ Tiana]])
    # @show cond(A) , cond(Ad)
    # @infiltrate
    # res[p, :] = [length(T) order t[1] t[2] t[3] t[1] / t[2] e1 e2 em1 em2]
    push!(
        res,
        [
            string(prob[p]),
            length(T),
            order,
            t[1],
            t[2],
            t[1] / t[2],
            t[3],
            t[5],
            t[6],
            t[5] / t[6],
            ei1,
            ei2,
            eim1,
            eim2,
        ],
    )
end
for i in prob, order = 2:5
    CSV.write(
        string(i, "_", order, ".csv"),
        filter(row -> (row.prob == string(i)) && (row.order == order), res),
    )
end
# p = lines(dad.NOS[:, 1], T)
# lines!(dad.NOS[:, 1], Td)
# lines!(dad.NOS[:, 1], Tana)
# p
