## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
# include(scriptsdir("placa.jl"))
nelem = [3]  #Numero de elementos
NPX = [5] #pontos internos na direção x
npg = [10]    #apenas números pares
# nt = [20, 40, 80, 160]    #apenas números pares
nt = [80]    #apenas números pares
## Formatação dos dados
problema = ["sladek03_apoiado"]
# problema = ["sladek03", "sladek03_apoiado"]
metodo = ["DRM"]
# metodo = ["Monta_M_RIMd", "Monta_M_RIM", "DRM"]
# metodo = ["Monta_M_RIMd"]
params = @strdict nelem NPX npg problema metodo nt

dicts = dict_list(params)
simula_placa_tempo(dicts[1])
for (i, d) in enumerate(dicts)
    @show savename(d)
    f = simula_placa_tempo(d)
    wsave(datadir("simulations", savename(d, "jld2")), f)
end
