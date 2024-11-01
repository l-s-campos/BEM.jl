## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 5  #Numero de elementos
tipo = 3
NPX = 2 #pontos internos na direção x
NPY = 2 #pontos internos na direção y
npg = 16    #apenas números pares
## Formatação dos dados ________________________________________________
# dad = format_dad(sapata_conforme(nelem, tipo), NPX, NPY) # dados
dad = format_dad(sapata(nelem, 3), NPX, NPY) # dados

H, G = calc_HeG(dad, npg)  #importante
# muda_nt!(dad, H, G)

nosrestritos = ceil.(Int, [nc(dad) - nelem * tipo * 2 - 1 1])
# nosrestritos = ceil.(Int, [nc(dad) - nelem * tipo * 1.5 1])
A, b = BEM.aplicaCDC(H, G, dad, nosrestritos) # Calcula a matriz A e o vetor b
# A, b = BEM.aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
# x = A \ b
# u, t = separa(dad, x, nosrestritos) #importante


h = calc_gap(dad)

A2 = [A; zeros(2 * size(h[1], 1), size(A, 1) + 2 * size(h[1], 1))]
b2 = [b; zeros(2 * size(h[1], 1), 1)]

# x0 = 0 * b2
# contato = verifica_contato_sem_atrito(x0, h)
# # Aplica a condição de contato em cada nó
# BEM.aplica_contato_sem_atrito!(h, contato, A2, b2, dad)
# # y0 = A2 * x0 - b2
# # x = x0 - A2 \ y0
# x = A2 \ b2
# # u, t = separa(dad, x, [], h)
# u, t = separa(dad, x, nosrestritos, h)
# mostra_deformação(dad, u, escala=10)


x0 = 0 * b2
x = BEM.Contato_sem_atrito_NL_newton(dad, x0, A2, b2, h, maxiter=100)
u, t = separa(dad, x, nosrestritos, h)
# # # prob = NonlinearProblem(Contato_sem_atrito_NL, x0, (A2, b2, h))

# # # sol = solve(prob, NewtonRaphson(), abstol=1e-8)

mostra_deformação(dad, u, escala=1)