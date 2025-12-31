using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
E = 200e9
v = 0.3
E_dash = 2E / (1 - v^2)

n = 1000
x0 = -10.0
xf = 10.0
y = collect(range(x0, stop = xf, length = n + 1))
x = (y[1:end-1] + y[2:end]) / 2
node = hcat(collect(1:n), collect(2:n+1)) # Element connectivity
nodes = [Point2D(node[i, 1], node[i, 2]) for i = 1:size(node, 1)]
al = abs.(y[node[:, 2]] .- y[node[:, 1]])# Length of the element

RE = 1000
h0 = x .^ 2 / (2RE) # original separation between the cylinder and the half-plane surface


dad = dadHS2D(y, x, nodes, h0, al, E_dash)

Kcheia = calc_K_2d(x0, xf, n, E = E_dash)
KFMM = montaFMM(dad)
KHmat = monta_hmat(dad)

ph = xf / RE * E / 4
p = ph * sqrt.(1 .- ((x / xf) .^ 2))
W = sum(p .* al)

u_Hmat = KHmat * p
u_FMM = KFMM * p
u_cheio = Kcheia * p
[u_cheio u_Hmat u_FMM]

fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel = "x",
    ylabel = "u(x)",
    title = "Displacement along the rough contact interface",
)
lines!(ax, x, u, label = "H-matrix")
lines!(ax, x, u_cheio, label = "Full matrix")
lines!(ax, x, x .^ 2 / (2RE), label = "teste")
axislegend(ax)
fig

@time p_con, g = contact_pressure_force(
    dad,
    K,
    W / 2;
    p_min = 0,
    H = 5e10,
    err_tol = 1e-8,
    it_max = 1000,
    h_ref = norm(h0),
)


lines(x, g, label = "gap")
lines(x, p_con, label = "pressure")
# lines(x, h0, label = "gap")
Wsaida = sum(p_con .* al)

a = sqrt(4Wsaida * (1 - v^2) * RE / (pi * E))