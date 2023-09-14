## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
# include(scriptsdir("placa.jl"))
nelem = [5]  #Numero de elementos
NPX = [20] #pontos internos na direção x
npg = [20]    #apenas números pares
## Formatação dos dados
# problema = ["sladek03_apoiado"]
# problema = ["redonda"]
problema = ["placa_furo_retangular"]
# problema = ["ex1"]
# carrega = ["comprex", "comprexy", "shear"]
# carrega = [-1, 0, 0]
# bc = ["FSCS"]
# bc = ["CCCC"]
bc = ["SSSS"]
# bc = ["SSSS", "SSSC", "CSSS", "SCSC", "CSCS", "FSSS", "SSSF", "FSCS", "SCSF", "FSFS", "SFSF", "CCCC"]
# anacomprex = Dict(bc .=> [4.0000, 4.8471, 5.7401, 6.7431, 7.6911, 1.4014, 2.3639, 1.6522, 2.3901, 0.9522, 2.0413, 10.0737])

# simula_placa_flambagem(dicts[1], [-1, 0, 0])
# simula_placa_flambagem(dicts[1])
# for i = 1:1
#     nelem = nelems[i]
#     NPX = NPXs[i]
params = @strdict nelem NPX npg problema bc
dicts = dict_list(params)
#     f, H, G, q, It, dNs, dMd, Md, dM, M = BEM.Matrizes_placa_flambagem(dicts[1])
#     for (i, d) in enumerate(dicts)
#         @show savename(d)
#         f = BEM.flambagem(f, d, H, G, It, dNs, dMd, Md, dM, M)
#         wsave(datadir("simulations\\flamba", savename(d, "jld2")), f)
#         # @show f
#     end
#     # MATLAB -> 44_49 -> 1718.375603segundos
# end

@unpack nelem, NPX, npg, problema, bc = dicts[1]
NPY = NPX
entrada = getfield(Main, Symbol(problema))(nelem, 3)
# entrada = getfield(Main, Symbol(problema))(nelem, 3, bc)

tdad = @timed dad = format_dad(entrada[1], NPX, NPY, canto=true) # dados
tdadpe = @timed dadpe = format_dad(entrada[2], NPX, NPY) # dados
tHeG = @timed H, G, q, It, dNs = calc_HeGeIt(dad, npg)  #importante
tHeGpl = @timed Hpe, Gpe = calc_HeG(dadpe, npg)
# dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
# println("2. Montando a matriz A e o vetor b")

# nosrestritos = [1 1]
A, b = aplicaCDC(Hpe, Gpe, dadpe) # Calcula a matriz A e o vetor b
# @infiltrate
x = A \ b
u, t = separa(dadpe, x)
tens_cont, tens_nt = calc_tens_cont(dadpe, u, t)
tens_int = calc_tens_int(dadpe, u, t, 30)


tdMd = @timed dMd = BEM.Monta_dM_RIMd(dad, npg)
# tdM = @timed dMd = BEM.Monta_dM_RIM(dad, npg)
# tMd = @timed Md = BEM.Monta_M_RIMd(dad, npg)
# tM = @timed Md = BEM.Monta_M_RIM(dad, npg)
# @infiltrate
# @show dad.pontos_internos
nosrestritos = [floor(Int, nelem / 2)+2 1
    floor(Int, nelem / 2)+2+nelem*3 2
    floor(Int, nelem / 2)+2+nelem*6 1
    floor(Int, nelem / 2)+2+nelem*9 2]

# ddMit = BEM.aplicaT(G, It, tens_nt, dNs, dMd[6], [tens_cont; tens_int])
dMit = BEM.aplicaT(dad, G, tens_nt, dNs, dMd[1], [tens_cont; tens_int])
# Mit = BEM.aplicaT(dad, Md, [tens_cont; tens_int])

ns = nc(dad) + ni(dad)
# ddMit = BEM.aplicaT(G, It, tens_nt, dNs, dMd[6], [-ones(ns) zeros(ns) zeros(ns)])
# dMit = BEM.aplicaT(dad, G, tens_nt, dNs, dMd[1], [-ones(ns) zeros(ns) zeros(ns)])
# Mit = BEM.aplicaT(dad, Md, [-ones(ns) zeros(ns) zeros(ns)])
# @infiltrate
# a, v = BEM.autovalor(i["H"], i["G"], i["Ib"], dad)
# a2, v2 = BEM.autovalor(H, G, ddMit, dad)
# a1, v1 = BEM.autovalor(H, G, dMit, dad)
a1, v1 = BEM.autovalor_num(H, G, dMit, dad)
# @time a1, v1 = BEM.autovalor(H, G, dMit, dad, num=true)
# a, v = BEM.autovalor(H, G, Mit, dad)
# a, v = BEM.autovalor(H, G, i["Ib"], dad)
# lambda2 = minimum(abs.(a2))
# lambda1 = minimum(abs.(a1))
# lambda = minimum(abs.(a))

# k2 = a2 / dad.k.D / pi^2
k1 = a1 / dad.k.D / pi^2
# k = a / dad.k.D / pi^2

# [k k1 k2]
# =#

