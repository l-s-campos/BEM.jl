## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
# include(scriptsdir("placacopy.jl"))
NPX = 7#pontos internos na direção x
nelem = 7
npg = 16   #apenas números pares
## Formatação dos dados 
# problema = "sladek03"
# problema = "large1"
# metodo = "Monta_M_RIM"
# metodo = "DRM"
# metodo = ["Monta_M_RIMd"]
NPY = NPX
# entrada = bhatta2(nelem, 3, "CCCC")
# entrada1 = bhatta3(nelem, 3, "CCCC")
# entradaiso = bhattaiso(nelem, 3, "CCCC")
# entrada = reddy(nelem, 3, "SSSS")
# entrada = marc(nelem, 3, "CCCC")
# entrada = chia(nelem, 3, "SSSS")
# entrada = Putcha(nelem, 3, "SSSS")
entrada = large1(nelem, 3, "SSSS")

tdad = @timed dad = format_dad(entrada[1], NPX, NPY, canto=true) # dados
tdadpe = @timed dadpe = format_dad(entrada[2], NPX, NPY, canto=true) # dados

# tdad = @timed dad1 = format_dad(entrada1[1], NPX, NPY, canto=true) # dados
# tdadpe = @timed dadpe1 = format_dad(entrada1[2], NPX, NPY, canto=true) # dados

# dadiso = format_dad(entradaiso[1], NPX, NPY, canto=true) # dados
# dadpeiso = format_dad(entradaiso[2], NPX, NPY, canto=true) # dados
# dad = format_dad(vibra(nelem,3),NPX,NPY) # dados
# println("2. Montando a matriz A e o vetor b")
nt = 50
# ws = BEM.placa_grande(dad, dadpe, nt, npg, 0.5)
# ws1 = BEM.placa_grande(dad1, dadpe1, nt, npg, 0.5)



mats = BEM.matrizes_grande(dad, dadpe, npg)
ws = BEM.placa_grande(mats, nt, dad, dadpe, 0.5)

# mats1 = BEM.matrizes_grande(dad1, dadpe1, npg)
# ws1 = BEM.placa_grande(mats1, nt, dad, dadpe1, 0.5)

# matsiso = BEM.matrizes_grande(dadiso, dadpeiso, npg)
# wsiso = BEM.placa_grande(matsiso, nt, dadiso, dadpeiso, 0.5)

no_meio = ceil(Int, 2 * size(dad.NOS, 1) + size(dad.pontos_internos, 1) / 2)

h = 1
# a = 16

# ws_meio = -ws[no_meio, :]/h
ws_meio = -ws[no_meio, :]
ws_meio[abs.(ws_meio).>2] .= NaN;

# ws_meio1 = -ws1[no_meio, :] / h
# ws_meioiso = -wsiso[no_meio, :] / h

# E2 = dad.k.Material[1, 3]
# E1 = dad.k.Material[1, 2]
# nult = dad.k.Material[5]
# nutl = nult * E2 / E1
# D0 = E2 * h^3 / 12 / (1 - nult * nutl)


# ws_meio = -ws[no_meio, :] / dad.k.h
# Q = (1:nt) / nt * -dad.k.carga[3] / D0 / h
# Q = (1:nt) / nt * -dad.k.carga[3] * (1 - nult * nutl) / dad.k.D22 / h
# Q = (1:nt) / nt * -dad.k.carga[3]
# Q = (1:nt) / nt * -dad.k.carga[3] / E2 / h^4
# Q = (1:nt) / nt * -dad.k.carga[3] * (1 - nult * nutl) / E2 / h^4
# Q = (1:nt) / nt * -dad.k.carga[3] / E1 / h^4
Q = (1:nt) / nt * -dad.k.carga[3] / dadpe.k.E / dad.k.h^4

fig, ax, p = BEM.lines(Q, ws_meio)
# BEM.lines!(ds.qb, ds.wb)
# BEM.lines!(Q, ws_meio1)
# BEM.lines!(Q, ws_meioiso)
# BEM.scatter!(ds.x, ds.w12)
# BEM.scatter!(Q, ws_meio1)
BEM.DataInspector(fig)
ws_meio[end]

# qpala = [5, 10, 15, 20, 25]
# wpala = [0.5214
#     0.8521
#     1.080
#     1.255
#     1.399]
# qpala = [50, 100, 150, 200, 250]
# wpala = [0.1878
#     0.3576
#     0.5035
#     0.6380
#     0.7355]
# BEM.scatter!(qpala, wpala)
# fig