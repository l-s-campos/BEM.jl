## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 10  #Numero de elementos
NPX = 2 #pontos internos na direção x
NPY = 2 #pontos internos na direção y
npg = 20    #apenas números pares
tipo = 3
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
dad = format_dad(viga(nelem, tipo), NPX, NPY) # dados

println("2. Montando a matriz A e o vetor b")
H, G = calc_HeG(dad, npg, interno = false)  #importante
# H, G = calc_HeG(dad, npg, interno = true)  #importante
Hd, Gd = BEM.calc_HeGd(dad, 3)

# us, ts = BEM.gera_ut(dad)
# LD = H * us - G * ts

# BEM.corrige_diagonais!(dad, H, G)
# G

# LDd = Hd * us - Gd * ts
# [H * us G * ts]

Ad, bd = aplicaCDC(Hd, Gd, dad) # Calcula a matriz A e o vetor b
xd = Ad \ bd

A, b = BEM.aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
# A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A \ b
println("4. Separando fluxo e temperatura")
u, t = separa(dad, x) #importante
ud, t = separa(dad, xd) #importante


tens_int = calc_tens_int(dad, u, t)
tens_cont, tens_nt = calc_tens_cont(dad, u, t)

Hx, Gx, Hy, Gy = BEM.calc_dHedG(dad, 8)
ux = Hx * u'[:] - Gx * t'[:]
uy = Hy * u'[:] - Gy * t'[:]


ut = [u; x[2nc(dad)+1:2:end] x[2nc(dad)+2:2:end]]
~, Fx, Fy = BEM.montaFs([dad.NOS; dad.pontos_internos], smooth = 1e-8)
ux2 = [Fx * ut[:, 1] Fx * ut[:, 2]]'[:]
uy2 = [Fy * ut[:, 1] Fy * ut[:, 2]]'[:]

[ux ux2 uy uy2]

Mx, My = BEM.Monta_deriv_M_RIMd(dad, 8)


# S, D = calc_SeD(dad)
# fatorcontorno = [fill(2, nc(dad)); ones(ni(dad))]
# ss = reshape(D * t'[:] - S * u'[:], 3, :)' .* fatorcontorno



# println("6. Gerando mapa de cor")
# erro=norm([T-dad.NOS[:,1];Ti-dad.pontos_internos[:,1]])

# # geo=mostra_geometria(dad)
# # mapa=mostra_resultado(dad,[T;Ti])
