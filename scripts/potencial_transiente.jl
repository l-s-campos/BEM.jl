## Início da análise
# using Pkg
# Pkg.activate(pwd())
# Pkg.instantiate()
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 10  #Numero de elementos
order = 3
NPX = 10 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 12    #apenas números pares

Δt = 0.01 #step de tempo
tf = 5
t = 0:Δt:tf
pontos_t = length(t)

## Formatação dos dados ________________________________________________
dad = format_dad(difusao_placa(nelem, order), NPX, NPY) # dados
# dad0 = format_dad(potencial1d(nelem),NPX,NPY,0) # dados
# dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados

Ht, Gt = calc_HeG(dad, npg, interno = true)  #importante
M = BEM.Monta_M_RIMd(dad, npg)# calc_HeG_potencial linha 310
A, b = BEM.aplicaCDC(Ht - 11 * M / (6 * Δt), Gt, dad)
T = zeros(length(A[:, 1]), Int(pontos_t) + 4)

for i in range(4, length(T[1, :]) - 1)
    x = A \ (b + M * (-18 * T[:, i] + 9 * T[:, i-1] - 2 * T[:, i-2]) / (6 * Δt))

    T_aux, q = separa(dad, x) #format 479
    Ti = x[nc(dad)+1:end]

    T[:, i+1] = [T_aux; Ti]
end

it = 5
lines(T[1:nc(dad), 5])
lines!(ana_difusao_placa.(t[it-3], dad.NOS[:, 2]))
BEM.current_figure()


