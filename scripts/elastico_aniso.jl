## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 10  #Numero de elementos
NPX = 2 #pontos internos na direção x
NPY = 2 #pontos internos na direção y
npg = 30    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
dad = format_dad(elastico_aniso_1d(nelem, 3), NPX, NPY) # dados
# dad = format_dad_iga(elastico_aniso_1d_iga(0,0),NPX,NPY) # dados


#println("2. Montando a matriz A e o vetor b")
H, G = calc_HeG(dad, npg)  #importante

A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A \ b
println("4. Separando deslocamento e forças de superficie")
u, t = separa(dad, x) #importante
sum(H, dims=2)
# lines(u[:,1])
# A1,b1 = calc_Aeb(dad,npg)  #importante
# x1=A1\b1
# println("5. Calculando nos pontos internos")
# Ti = calc_Ti(dad,T,q,npg);

# println("6. Gerando mapa de cor")
# erro=norm([T-dad.NOS[:,1];Ti-dad.pontos_internos[:,1]])

# # geo=mostra_geometria(dad)
# # mapa=mostra_resultado(dad,[T;Ti])