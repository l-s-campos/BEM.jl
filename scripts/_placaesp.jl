## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
# include(scriptsdir("placacopy.jl"))
NPX = 3#pontos internos na direção x
nelem = 3
npg = 16 #apenas números pares
## Formatação dos dados 
# problema = "sladek03"
# problema = "large1"
# metodo = "Monta_M_RIM"
# metodo = "DRM"
# metodo = ["Monta_M_RIMd"]
NPY = NPX

# entrada = plespex1(nelem, 3, "CCCC")
entrada, entradaterm = termbuckl(nelem, 3, "SSSS")

dad = format_dad(entrada, NPX, NPY, canto=true) # dados
H, G, q = calc_HeG(dad, npg)

# A, B, bc_val = BEM.aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
# println("3. Resolvendo o sistema linear")
# x = A \ (B * bc_val + q)
# println("4. Separando fluxo e temperatura")
# u, t = separa(dad, x) #importante

# u[end] * 100 * dad.k.D
np = nc(dad) + ni(dad)
# tens = [-ones(np) zeros(np) zeros(np)];
# tens = [-ones(np) -ones(np) zeros(np)];
dadterm = format_dad(entradaterm, NPX, NPY); # dados
α = 2e-6
u, tens = termoelasticidade(dadterm, npg, θ=1);

tensa = [-ones(np) -ones(np) zeros(np)] * dad.k.h * dadterm.k.E * α / (1 - dadterm.k.nu);
f1 = 1
f2(x, y) = x
f3(x, y) = (x + y) / 2

# M = Monta_M_RIM(dad);
Md = Monta_M_RIMd(dad);

M2d = BEM.aplicaT(dad, Md, tens)
# M2 = BEM.aplicaT(dad, M, tens);
# M2 = BEM.aplicaT(dad, Md, tens);
# ad, v = BEM.autovalor(H, G, M2d, dad)
a, v = BEM.autovalor(H, G, M2, dad;)
a[1]
dad.k.D * pi^2 * (1 - dad.k.nu) / (4dadterm.k.E * α * dad.k.h) * 2
pi^2 * dad.k.h^2 / (48(1 + dad.k.nu) * α) * (2)
2dad.k.D * pi^2 / dad.k.h / dadterm.k.E / α * (1 - dadterm.k.nu)
k = a[1] / dad.k.D / pi^2
# kd = ad / dad.k.D / pi^2
# k[1], kd[1]