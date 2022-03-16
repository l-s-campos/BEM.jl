## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 20  #Numero de elementos
NPX = 2 #pontos internos na direção x
NPY = 2 #pontos internos na direção y
npg = 100    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
dad = format_dad_iga( placa482_iga(nelem),NPX,NPY) # dados
dad = aplicaCDCplaca_iga(dad)

# dad = format_dad_iga(potencial1diso(nelem),NPX,NPY) # dados
# dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados

println("2. Montando a matriz A e o vetor b")
H,G = calc_HeG(dad,npg)  #importante


A,b = aplicaCDC(H,G,dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A\b
println("4. Separando fluxo e temperatura")
T,q = separa(dad,x) #importante
println("6. Gerando mapa de cor")

Ti = calc_Ti(dad,T,q,npg);

Treal=placa482treal(dad)



# @show erro=norm([T-dad.NOS[:,1];Ti-dad.pontos_internos[:,1]])

# geo=mostra_geometria(dad)
# mapa=mostra_resultado(dad,[T;Ti])