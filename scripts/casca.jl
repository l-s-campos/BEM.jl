## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))

println("em desenvolvimento")
# include(scriptsdir("placacopy.jl"))
NPX = 3 #pontos internos na direção x
nelem = 3
npg = 16    #apenas números pares

entradaiso = bhattaiso(nelem, 3, "CCCC")
dad1 = format_dad(entradaiso[1], NPX, NPY, canto=true)
dadpe1 = format_dad(entradaiso[2], NPX, NPY, canto=true)
R11 = 100.0
R22 = 100.0
Casca = casca(dad1, dadpe1, R11, R22)
H, G = calc_HeG(Casca, npg)

M, P = Monta_M_RIMd(Casca, npg)


