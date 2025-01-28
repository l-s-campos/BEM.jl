## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 20  #Numero de elementos
NPX = 1 #pontos internos na direção x
NPY = 1 #pontos internos na direção y
npg = 10   #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
entrada = elastico1d_inclusao(nelem, 3)
dad = format_dad(entrada[1], NPX, NPY) # dadosH
dad1 = format_dad(entrada[2], NPX, NPY) # dados
dads = [dad, dad1]
# dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados

println("2. Montando a matriz A e o vetor b")
H, G = calc_HeG(dads, npg)  #importante
H1, G1 = calc_HeG(dad, npg)  #importante
H2, G2 = calc_HeG(dad1, npg)  #importante


A, b = aplicaCDC(H, G, dads) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A \ b
println("4. Separando fluxo e temperatura")
u, t = separa(dads, x) #importante
A1, b1 = aplicaCDC(H1, G1, dad) # Calcula a matriz A e o vetor b
x1 = A1 \ b1
u1, t1 = separa(dad, x1) #importante
[u1 u[1:12nelem, :]]
# plot([u1 u[1:12nelem, :]])
scatter(u)
# # A1,b1 = calc_Aeb(dad,npg)  #importante
# # x1=A1\b1
# # println("5. Calculando nos pontos internos")
# # Ti = calc_Ti(dad,T,q,npg);

# # println("6. Gerando mapa de cor")
# # erro=norm([T-dad.NOS[:,1];Ti-dad.pontos_internos[:,1]])

# # # geo=mostra_geometria(dad)
# # # mapa=mostra_resultado(dad,[T;Ti])
