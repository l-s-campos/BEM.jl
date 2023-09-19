#tem de mudar algo no GMRES para aceitar o HMAT
## Início da análise
# include("potencial_interpolado.jl")
using TimerOutputs
to = TimerOutput()
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 100  #Numero de elementos
NPX = 1 #pontos internos na direção x
NPY = 1 #pontos internos na direção y
npg = 16    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
@timeit to "dados" dad = format_dad_iga(potencial1diso(nelem), NPX, NPY) # dados
println("2. Montando a matriz A e o vetor b")
# Ht,Gt = calc_HeG(dad,1:nc(dad)+ni(dad),1:4*nelem,npg)  #importante
A1, b1 = @btime MatrizH(dad, BEM.Akmeans, 9);
H, G = @btime calc_HeG(dad, npg)  #importante
show(to)
Aicheia = BEM.montaAcheia(A1.A, A1.block, A1.Tree1, A1.Tree2, A1.cols, dad)
println("3. Resolvendo o sistema linear")
# x2,f = gmres(A,b,5,tol=1e-5,maxIter=100,out=0) #GMRES nas matrizes hierarquicas
H, G = calc_HeG(dad, npg)  #importante=
A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
x = A \ b
println("4. Separando fluxo e temperatura")
T, q = separa(dad, x) #importante
xi, f = gmres(A1, b1, 5, tol=1e-5, maxIter=100, out=0) #GMRES nas matrizes 
T1, q1 = separa(dad, xi) #importante
T1i = xi[nc(dad)+1:end]
println("5. Calculando nos pontos internos")
Ti = calc_Ti(dad, T, q, npg);

println("6. Gerando mapa de cor")
erro = norm([T - dad.NOS[:, 1]; Ti - dad.pontos_internos[:, 1]])
erro = norm([T1 - dad.NOS[:, 1]; T1i - dad.pontos_internos[:, 1]])
println("6. Gerando mapa de cor")
erro = norm([T - dad.NOS[:, 1]; Ti - dad.pontos_internos[:, 1]])
erroi = norm([T1 - dad.NOS[:, 1]; T1i - dad.pontos_internos[:, 1]])

# geo=mostra_geometria(dad)
# mapa=mostra_resultado(dad,[T;Ti])

# @btime matvec(Aaca,baca,block,Tree1,Tree2,dad)
# @btime Ac*baca