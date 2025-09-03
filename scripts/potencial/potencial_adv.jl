## Início da análise
# using Pkg
# Pkg.activate(pwd())
# Pkg.instantiate()
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
i = 2

npg = 12    #apenas números pares
nelem = 2 * i  #Numero de elementos
order = 3
NPX = 3 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
## Formatação dos dados ________________________________________________
dad = format_dad(potencial1d(nelem, order), NPX, NPY) # dados

Ht, Gt = BEM.calc_HeG(dad, npg, interno = true)

v1 = zeros(size(Ht, 1)) .+ 1
v2 = zeros(size(Ht, 1))

# Mh = BEM.Monta_M_RIMd_hermite(dad, npg);
M = BEM.Monta_M_RIMd_hermite(dad, npg);
# 
# M = BEM.Monta_M_RIMd(dad, npg, aug = true);
# M = BEM.Monta_M_RIMd(dad, npg, aug = false);
~, Fx, Fy = BEM.montaFs([dad.NOS; dad.pontos_internos], aug = true)
Ml = v1 .* M * Fx + v2 .* M * Fy
A, b = BEM.aplicaCDC(Ht - Ml, Gt, dad)

# N1 = diagm(dad.normal[:, 1])
# N2 = diagm(dad.normal[:, 2])
# Mq =
#     v1 .* (M * [N1; zeros(ni(dad), nc(dad))] / (-dad.k)) +
#     v2 .* (M * [N2; zeros(ni(dad), nc(dad))] / (-dad.k))

# Mu =
#     v1 .* (M * [N2.^2 * Fx[1:nc(dad), :]; Fx[nc(dad)+1:end, :]]) +
#     v2 .* (M * [N1.^2 * Fy[1:nc(dad), :]; Fy[nc(dad)+1:end, :]])
# A1, b1 = BEM.aplicaCDC(Ht - Mu, Gt + Mq, dad)
# [[-N2 * Fx[1:nc(dad), :]; Fx[nc(dad)+1:end, :]] * Tana (
#     M * [N1; zeros(ni(dad), nc(dad))] / (-dad.k)
# ) * q Fx * Tana Txana]
# # x = A1 \ b1

x = A \ b
Tc, q = separa(dad, x) #format 479
Ti = x[nc(dad)+1:end]
u(x, vx) = (exp(vx * x) - 1) / (vx * exp(vx))
Tana = u.([dad.NOS[:, 1]; dad.pontos_internos[:, 1]], v1)
udx(x, vx) = exp(vx * x - vx)
Txana = udx.([dad.NOS[:, 1]; dad.pontos_internos[:, 1]], v1)

e = nme(Tana, [Tc; Ti])
println(nc(dad) + ni(dad), ",", nc(dad), ",", ni(dad), ",", e)
#  [[abs.(N2) * Fx[1:nc(dad), :]; Fx[nc(dad)+1:end, :]]*Tana+ [N1; zeros(ni(dad), nc(dad))] / (-dad.k)*q Fx*Tana Txana]

# lines(Tana)
# lines!( [Tc; Ti])
# BEM.current_figure()
# nodes = [dad.NOS; dad.pontos_internos]
# Fx * nodes[:, 1]