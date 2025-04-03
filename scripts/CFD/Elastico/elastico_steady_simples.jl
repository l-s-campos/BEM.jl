## Início da análise
using DrWatson, Plots, CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
includet(srcdir("CFD.jl"))

#Dados de entrada

nelem = 20 #Numero de elementos
order = 2
nelem_circle = 20

NPX = 15 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

L = 1
Re = 1 # Na verdade é o 1/ν ou 1/μ

caso = "Cavidade"
λ = 10^6

dad = format_dad(cavity(nelem, order, L, 1, λ), NPX, NPY)


begin
    H, G = calc_HeG(dad, npg, interno = true)

    M = BEM.Monta_M_RIMd(dad, npg)
end

# ==============Discretização============#

begin
    nx = length(unique(dad.pontos_internos[:, 1]))
    ny = length(unique(dad.pontos_internos[:, 2]))

    x_order = sort(unique(dad.pontos_internos[:, 1]))
    y_order = sort(unique(dad.pontos_internos[:, 2]))

    dx = (maximum(x_order) - minimum(x_order)) / nx
    dy = (maximum(y_order) - minimum(y_order)) / ny

end

# =========Inicia variáveis===============#

begin

    n = ni(dad) + nc(dad)

    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)

    x = zeros(2 * nc(dad) + 2 * ni(dad))
    p = zeros(ni(dad) + nc(dad))


end
dNx, dNy = BEM.montaFs([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])[2:3]

# ==============Solução===============#

Re = 200.0

x1 = passo_não_linear(
    x,
    dad,
    H,
    G,
    M,
    dNx,
    dNy,
    Re = Re,
    relaxation = 0.05,
    erro_vel_min = 1e-6,
    maxiter = 1000,
)
x2 = passo_não_linear2(
    x,
    dad,
    H,
    G,
    M,
    dNx,
    dNy,
    Re = Re,
    relaxation = 0.05,
    erro_vel_min = 1e-6,
    maxiter = 1000,
)

# using BifurcationKit

# par = [dad, H, G, M, dNx, dNy, 800.0];

# função_não_linear(x, par)
# prob2 = BifurcationProblem(função_não_linear, x, par, 7)

# sol = BifurcationKit.solve(prob2, Newton(), NewtonPar(tol = 1e-6, verbose = true))
# x = sol.u

# br2 = continuation(
#     prob2,
#     PALC(),
#     ContinuationPar(p_min = 10.0, p_max = 100.0, newton_options = NewtonPar(tol = 1e-6)),
# )
# x = br2.sol[2].x



begin
    u = zeros(nc(dad) + ni(dad), 2)
    u[contorno, :], t = separa(dad, x)
    u[interno, 1] = x[2*nc(dad)+1:2:end]
    u[interno, 2] = x[2*nc(dad)+2:2:end]
    uc = u[contorno, :]
    ui = u[interno, :]

    dudx = dNx * u[:, 1]
    dudy = dNy * u[:, 1]
    dvdx = dNx * u[:, 2]
    dvdy = dNy * u[:, 2]

    p = -λ * (dudx + dvdy)

    # #=========Heatmap e Quiver=========#

    x_array = dad.pontos_internos[:, 1]
    y_array = dad.pontos_internos[:, 2]

    utotal = sqrt.(ui[:, 1] .^ 2 + ui[:, 2] .^ 2)

    BEM.heatmap(x_array, y_array, utotal, interpolate = true)
    escala = 10^(-1)
    BEM.quiver!(x_array, y_array, escala * ui[:, 1], escala * ui[:, 2], color = :black)
    BEM.current_figure()
end

