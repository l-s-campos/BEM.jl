## Início da análise
# using Pkg
# Pkg.activate(pwd())
# Pkg.instantiate()
using DataFrames, CSV
include(scriptsdir("includes.jl"))
@quickactivate "BEM"
nelem = 20  #Numero de elementos
order = 2
NPX = 10 #pontos internos na direção x
NPY = 10 #pontos internos na direção y
npg = 10    #apenas números pares
# ## Formatação dos dados ________________________________________________
# println("1. Formatando os dados");
# dad = format_dad(potencial1d(nelem, 2), NPX, NPY) # dados
# # dad0 = format_dad(potencial1d(nelem),NPX,NPY,0) # dados
# # dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
# res = zeros(5, 4, 4, 10)
prob = [placa_moulton, potencial1d, dad_laquini1, dad_laquini2, dad_laquini3]
ana = [T_ana_moulton, T_ana_1d, fa1, fa2, fa3]
res = DataFrame(
    prob = String[],
    n = Int[],
    order = Int[],
    t = Float64[],
    td = Float64[],
    tHd = Float64[],
    tsolve = Float64[],
    tsolved = Float64[],
    tsolveHd = Float64[],
    erms = Float64[],
    ermsd = Float64[],
    ermsHd = Float64[],
    em1 = Float64[],
    em2 = Float64[],
    em3 = Float64[],
)
for p = 1:1, i = 3:6, j = 2:2
    nelem = 4 * 2^i  #Numero de elementos
    order = j
    # p = 1

    dad = format_dad(prob[p](nelem, order), NPX, NPY) # dados
    # include(datadir("dadpotencial.jl"))

    if prob[p] == placa_moulton
        dad = corrigeCDC_moulton(dad)
    end

    # function compara_potencial(dad, npg)
    reset_timer!()
    tHeG = @timed H, G = calc_HeG(dad, npg,interno=true)  #importante
    tHeGd = @timed Hd, Gd = BEM.calc_HeGd(dad, 3,interno=true)  #importante
    tHeGHd = @timed HH, HG = BEM.calc_HeG_Hd(dad, atol = 1e-6)

    A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
    Ad, bd = aplicaCDC(Hd, Gd, dad) # Calcula a matriz A e o vetor b
    bH = aplicaCDC(HH, HG, dad)

    tsolve = @timed x = A \ b
    tsolved = @timed xd = Ad \ bd

    tsolve = @timed x, f = gmres(A, b, rtol = 1e-5, itmax = 100) #GMRES nas matrizes
    tsolved = @timed xd, f = gmres(Ad, bd, rtol = 1e-5, itmax = 100) #GMRES nas matrizes
    tsolveHd = @timed xH, f = gmres(HH, bH, rtol = 1e-5, itmax = 100) #GMRES nas matrizes



    T, q = separa(dad, x) #importante
    Td, qd = separa(dad, xd) #importante
    THd, qHd = separa(dad, xH) #importante
    # [tHeG.time tHeGd.time tsolve.time tsolved.time], T, Td, q, qd
    # end
    ti = @timed Ti = calc_Ti(dad, T, q, npg)
    tid = @timed Tid = BEM.calc_Tid(dad, Td, qd, 3)

    t =
        [tHeG.time  tHeGd.time  tHeGHd.time tsolve.time tsolved.time tsolveHd.time]
    # t, T, Td, q, qd = compara_potencial(dad, npg)
    Tana = ana[p](dad.NOS)
    e1 = norm(T - Tana) / norm(Tana)
    e2 = norm(Td - Tana) / norm(Tana)
    e3 = norm(THd[1:nc(dad)] - Tana) / norm(Tana)
    em1 = sum(abs, T - Tana) / sum(abs, Tana)
    em2 = sum(abs, Td - Tana) / sum(abs, Tana)
    em3 = sum(abs, THd[1:nc(dad)] - Tana) / sum(abs, Tana)

    Tiana = ana[p](dad.pontos_internos)
    ei1 = norm(Ti - Tiana) / norm(Tiana)
    ei2 = norm(Tid - Tiana) / norm(Tiana)
    # @infiltrate
    ei3 = norm(xH[nc(dad)+1:end] - Tiana) / norm(Tiana)
    eim1 = sum(abs, Ti - Tiana) / sum(abs, Tiana)
    eim2 = sum(abs, Tid - Tiana) / sum(abs, Tiana)
    eim3 = sum(abs, xH[nc(dad)+1:end] - Tiana) / sum(abs, Tiana)

    # res = DataFrame(prob=String[], n=Int[], order=Int[], t=Float64[], td=Float64[], tHd=Float64[], tsolve=Float64[], tsolved=Float64[], tsolveHd=Float64[], erms=Float64[], ermsd=Float64[], ermsHd=Float64[], em1=Float64[], em2=Float64[], em3=Float64[])

    # res[p, :] = [length(T) order t[1] t[2] t[3] t[1] / t[2] e1 e2 em1 em2]
    push!(
        res,
        [
            string(prob[p]),
            length(T),
            order,
            t[1],
            t[2],
            t[3],
            t[4],
            t[5],
            t[6],
            ei1,
            ei2,
            ei3,
            eim1,
            eim2,
            eim3,
        ],
    )
end
# for i in prob, order = 2:5
#     CSV.write(
#         string(i, "_", order, ".csv"),
#         filter(row -> (row.prob == string(i)) && (row.order == order), res),
#     )
# end
# p = lines(dad.NOS[:, 1], T)
# lines!(dad.NOS[:, 1], Td)
# lines!(dad.NOS[:, 1], Tana)
# p
