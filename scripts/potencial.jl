## Início da análise
# using Pkg
# Pkg.activate(pwd())
# Pkg.instantiate()
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 20  #Numero de elementos
order = 2
NPX = 10 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 20    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
dad = format_dad(potencial1d(nelem, order), NPX, NPY) # dados
# dad0 = format_dad(potencial1d(nelem),NPX,NPY,0) # dados
# dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados

println("2. Montando a matriz A e o vetor b")
# H, G = BEM.calc_HeGd(dad, 3)
H, G = calc_HeG(dad, npg)  #importante
# H1,G1 = calc_HeG(dad,10*npg)  #importante


A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A \ b
println("4. Separando fluxo e temperatura")
T, q = separa(dad, x) #importante
# A1, b1 = calc_Aeb(dad, npg)  #importante
# x1 = A1 \ b1
println("5. Calculando nos pontos internos")
Ti = calc_Ti(dad, T, q, npg);

println("6. Gerando mapa de cor")
# @show erro = norm([T - dad.NOS[:, 1]; Ti - dad.pontos_internos[:, 1]])
# @show erro = norm(T - dad.NOS[:, 1])
# @show erro = norm( Ti - dad.pontos_internos[:, 1])

# Au, bu = BEM.calc_Aeb_emu(dad)
# Tu = Au \ bu

@show length(T), order
@show sum(abs, T - dad.NOS[:, 1]) / length(T)
# @show sum(abs, Tu - [dad.NOS[:, 1]; dad.pontos_internos[:, 1]]) / length(Tu)
# @show sum(abs, [T; Ti] - [dad.NOS[:, 1]; dad.pontos_internos[:, 1]]) / length(Tu)
@show sum(abs, Ti - dad.pontos_internos[:, 1]) / length(Ti)

# t, T, Td, q, qd=BEM.compara_potencial(potencial1d,5,2,10,1,1)
# geo=mostra_geometria(dad);
# mapa=mostra_resultado(dad,[T;Ti]);
# mean(abs, T - dad.NOS[:, 1])
#=
@time include(scriptsdir("potencial.jl"))
@time include(scriptsdir("potencial_direto.jl"))

=#