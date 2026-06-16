using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
E = 210e9
v = 0.3
E_dash = 2 / ((1 - v^2) / E + (1 - v^2) / E) / 2


RE = 12e-3 #raio cilindro

n = 200
x0 = -.4e-3
xf = -x0

y = collect(range(x0, stop = xf, length = n + 1))
x = (y[1:end-1] + y[2:end]) / 2
node = hcat(collect(1:n), collect(2:n+1)) # Element connectivity
nodes = [Point2D(node[i, 1], node[i, 2]) for i = 1:size(node, 1)]
al = abs.(y[node[:, 2]] .- y[node[:, 1]])# Length of the element

h0 = x .^ 2 / (2RE) # original separation between the cylinder and the half-plane surface


dad = dadHS2D(y, x, nodes, h0, al, E_dash)

Kcheia = calc_K_2d(x0, xf, n, E = E_dash)
# KFMM = montaFMM(dad)
# KHmat = monta_hmat(dad)
W = 185 * 100

a = sqrt(4W * (1 - v^2) * RE / (pi * E))
a * 1e3

p_a(x) = abs(x) < a ? 2W / pi / a * sqrt(1 - x^2 / a^2) : 0
pa = p_a.(x)
maximum(pa) / 1e6
lines(x * 1e3, pa / 1e6)
p_con, g = contact_pressure_force(
    dad,
    Kcheia,
    W;
    H = 5e12,
    err_tol = 1e-8,
    it_max = 1000,
    h_ref = norm(h0),
);

Wsaida = sum(p_con .* al)

lines(x * 1e3, dad.h0, label = "gap")
lines(x * 1e3, g, label = "gap")
lines(x * 1e3, p_con / 1e6, label = "normal pressure")
# lines(x * 1e3, (pa-p_con), label = "gap")


w1, w2 =
    desgaste_2D(dad, Kcheia, W, k_ar1 = 2.e-14, k_ar2 = 3.e-14, δ = 100e-6, nsteps = 18000)




begin



    fig = Figure(size = (600, 800))
    ax = Axis(
        fig[1, 1],
        xlabel = "Horizontal Position (mm)",
        ylabel = "Vertical Position (μm)",
        title = "Cylinder wear",
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    # Plot lines for each n
    lines!(ax, x * 1e3, dad.h0 * 1e6, label = "n=0", linewidth = 3, color = :black)
    lines!(
        ax,
        x * 1e3,
        (w1[:, 1000] + dad.h0) * 1e6,
        label = "n=1000",
        linewidth = 3,
        color = :red,
    )
    lines!(
        ax,
        x * 1e3,
        (w1[:, 5000] + dad.h0) * 1e6,
        label = "n=5000",
        linewidth = 3,
        color = :blue,
    )
    lines!(
        ax,
        x * 1e3,
        (w1[:, 10000] + dad.h0) * 1e6,
        label = "n=10000",
        linewidth = 3,
        color = :green,
    )
    lines!(
        ax,
        x * 1e3,
        (w1[:, 18000] + dad.h0) * 1e6,
        label = "n=18000",
        linewidth = 3,
        color = :magenta,
    )



    ax2 = Axis(
        fig[2, 1],
        xlabel = "Horizontal Position (mm)",
        ylabel = "Vertical Position (μm)",
        title = "Plate wear",
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    # Plot lines for each n
    lines!(ax2, x * 1e3, 0 * dad.h0 * 1e6, label = "n=0", linewidth = 3, color = :black)
    lines!(
        ax2,
        x * 1e3,
        (-w2[:, 1000]) * 1e6,
        label = "n=1000",
        linewidth = 3,
        color = :red,
    )
    lines!(
        ax2,
        x * 1e3,
        (-w2[:, 5000]) * 1e6,
        label = "n=5000",
        linewidth = 3,
        color = :blue,
    )
    lines!(
        ax2,
        x * 1e3,
        (-w2[:, 10000]) * 1e6,
        label = "n=10000",
        linewidth = 3,
        color = :green,
    )
    lines!(
        ax2,
        x * 1e3,
        (-w2[:, 18000]) * 1e6,
        label = "n=18000",
        linewidth = 3,
        color = :magenta,
    )
    Legend(fig[1, 2], ax)
    Legend(fig[2, 2], ax)
    display(fig)
end
