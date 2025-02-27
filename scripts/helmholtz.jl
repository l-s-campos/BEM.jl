## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 50  #Numero de elementos
NPX = 40 #pontos internos na direção x
NPY = 40 #pontos internos na direção y
npg = 10    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
prob = [helm1d, helmdirechlet, helmcirculo]#helmcirculoinfinito
ana = [ANA_helm1d, ANA_helmdirechlet, ANA_helmcirculo]
res = zeros(Float64, 0, 3)
# for fr = 0.7:1:6
# dad = format_dad(helmcirculoinfinito(nelem, 3), NPX, NPY) # dados
dad = format_dad(helm1d(nelem, 3, 3), NPX, NPY) # dados

dadpot, ω, c = potencial_correlato(dad::helmholtz)


# dad = format_dad(potencial1d(nelem), NPX, NPY) # dados
# dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
# corrigeCDC_helmdirechlet(dad)
# println("2. Montando a matriz A e o vetor b")
H, G = calc_HeG(dad, npg)  #importante
H1, G1 = BEM.calc_HeG_hiper(dad, npg)  #importante

A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
A1, b1 = aplicaCDC(H1, G1, dad) # Calcula a matriz A e o vetor b
# println("3. Resolvendo o sistema linear")
x = A \ b
x1 = A1 \ b1
xbm = (A + im / ω * A1) \ (b + im / ω * b1)


# println("4. Separando fluxo e temperatura")
T, q = separa(dad, x) #importante
Th, qh = separa(dad, x1) #importante
Tbm, qbm = separa(dad, xbm) #importante

T1 = real(T)
# q1 = real(q)

# Tana = sin.(dad.k.FR * dad.NOS[:, 2]) / dad.k.FR / cos(dad.k.FR)
# norm(T1 - Tana) / norm(Tana)


ana = (1 / ω * BEM.besselh(0, 2, ω) / BEM.besselh(1, 2, ω))
ω, ana, abs.(real(([mean(T) - ana mean(T1) - ana mean(Tbm) - ana])))
# global res = [res; [ω ana mean(T) mean(T1)]]


# tsf1 = zeros(eltype(H), nc(dad));
# qsf1 = zeros(eltype(H), nc(dad));
# for i = 1:nc(dad)
#     qsf1[i], tsf1[i] = BEM.calsolfund(dad.NOS[i, :], dad.normal[i, :], dad)
#     # @show i, tsf1[i]
# end
# # @show
# norm(H * tsf1 - G * qsf1), norm(H1 * tsf1 - G1 * qsf1)
# H1[1, 1:4]'
# @show H[1, 1], G[1, 1]

# corrigediag!(H, G, dad)
# @show H[1, 1], G[1, 1]

# end
Hpot, Gpot = calc_HeG(dadpot, 10, interno = true)  #importante

M = BEM.Monta_M_RIMd(dadpot, npg, tiporadial = "IQ")

Apot, bpot = aplicaCDC(Hpot + ((ω / c)^2) * M, Gpot, dadpot) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
xpot = Apot \ -bpot
println("4. Separando fluxo e temperatura")
T2, q2 = separa(dadpot, xpot) #importante
# 
Ta = ANA_helm1d(dad)

e1 = norm(T1 - Ta) / norm(Ta)
e2 = norm(T2 - Ta) / norm(Ta)

lines(T1, label = "Helmholtz $e1")
lines!(T2, label = "Poisson $e2")
lines!(Ta, label = "analitico")
BEM.axislegend(position = :lb)
BEM.current_figure()