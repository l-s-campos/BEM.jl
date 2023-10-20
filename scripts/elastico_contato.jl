## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 1  #Numero de elementos
npi = 80
di = 1e-1

npg = 10    #apenas números pares

println("1. Formatando os dados");

dad = format_dad(Subregiao(nelem, 3), pontosi) # dados

println("2. Montando a matriz A e o vetor b")
H, G = calc_HeG(dad, npg, subcheia=false)  #importante


A, b = BEM.aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b

println("3. Resolvendo o sistema linear")
x = A \ b
println("4. Separando fluxo e temperatura")
u, t = separa(dad, x) #importante
tens_int = calc_tens_int(dad, u, t, fill(2, size(dad.pontos_internos)), 30)


