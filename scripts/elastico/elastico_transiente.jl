## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 20  #Numero de elementos
NPX = 15 #pontos internos na direção x
NPY = 15#pontos internos na direção y
npg = 20    #apenas números pares
tipo = 2
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
dad = format_dad(elastico1d(nelem, tipo), NPX, NPY) # dados
ρ = 1


ind = ceil(Int, nelem * 1.5 * tipo)#indice metade do segmento 2, onde tem a força aplicada

println("2. Montando a matriz A e o vetor b")
H, G = calc_HeG(dad, npg, interno = true)  #importante
M = BEM.Monta_M_RIMd(dad, npg)# calc_HeG_potencial linha 310
Δt = 2e-3 #step de tempo
tf = 0.1
t = 0:Δt:tf
pontos_t = length(t)

M .= M / (Δt^2) * ρ
#estático
# A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
# x = A \ (b)
# u, t = separa(dad, x) #importante    
# u[ind, 1], (1 - .25^2) / 1e5


A, b = BEM.aplicaCDC(H + 2 * M, G, dad) #Houbolt
u = zeros((nc(dad) + ni(dad)) * 2, pontos_t)
du = zeros((nc(dad) + ni(dad)) * 2, pontos_t)
FA = lu(A)
for i = 4:pontos_t
    x = FA \ (b + M * (5 * u[:, i-1] - 4 * (u[:, i-2]) + 1 * (u[:, i-3])))
    ui, ti = separa(dad, x) #importante    
    u[1:2:2*nc(dad), i] = ui[:, 1]
    u[2:2:2*nc(dad), i] = ui[:, 2]
    u[2*nc(dad)+1:end, i] = x[2*nc(dad)+1:end]
    du[:, i] = (2u[:, i] - 5 * u[:, i-1] + 4 * (u[:, i-2]) - 1 * (u[:, i-3])) / (Δt)^2
end

lines(t, u[2ind-1, :])
# lines(t, du[2ind-1, :])

