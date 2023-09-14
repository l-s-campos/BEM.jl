## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 4  #Numero de elementos
NPX = 2 #pontos internos na direção x
NPY = 2 #pontos internos na direção y
npg = 10    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
dad = format_dad(elastico1d(nelem, 3), NPX, NPY) # dados

println("2. Montando a matriz A e o vetor b")
H, G = calc_HeG(dad, npg)  #importante


A, b = BEM.aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
# A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A \ b
println("4. Separando fluxo e temperatura")
u, t = separa(dad, x) #importante
tens_int = calc_tens_int(dad, u, t)
tens_cont, tens_nt = calc_tens_cont(dad, u, t)

S, D = calc_SeD(dad)
fatorcontorno = [fill(2, nc(dad)); ones(ni(dad))]
ss = reshape(D * t'[:] - S * u'[:], 3, :)' .* fatorcontorno



# println("6. Gerando mapa de cor")
# erro=norm([T-dad.NOS[:,1];Ti-dad.pontos_internos[:,1]])

# # geo=mostra_geometria(dad)
# # mapa=mostra_resultado(dad,[T;Ti])

