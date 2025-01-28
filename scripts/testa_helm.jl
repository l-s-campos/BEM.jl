using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))

nelem = 1000 #Numero de elementos
NPX = 2 #pontos internos na direção x
NPY = 2 #pontos internos na direção y
npg = 10    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
dad = format_dad(potencial1d(nelem), NPX, NPY) # dados
@time HH, HG = BEM.calc_HeG_Hd(dad, atol = 1e-6)
@time Hd, Gd = BEM.calc_HeGd(dad, 3)
# b = aplicaCDC(HH, HG, dad)
Ad, bd = aplicaCDC(Hd, Gd, dad)
xd = Ad \ bd

b = aplicaCDC(HH, HG, dad)

# tipoCDC = BEM.tipoCDC(dad)
# valoru, valorq = BEM.valoresCDC(dad)
# b = HG * valorq - HH * valoru
# @show norm(b[1:end-4] - bd)

@time x, f = gmres(HH, b, rtol = 1e-5, itmax = 100) #GMRES nas matrizes
@time xd, f = gmres(Ad, bd, rtol = 1e-5, itmax = 100) #GMRES nas matrizes
# @show norm(x[1:end-4] - xd)

HM = BEM.Monta_M_Hd(dad, 10)
M = BEM.Monta_M_RIMd(dad, npg)
size(Ad, 1)
