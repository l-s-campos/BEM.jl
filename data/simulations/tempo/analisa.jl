using DataFrames
using DrWatson
using CairoMakie, MakiePublication
CairoMakie.activate!()
engastado =
    collect_results(datadir("simulations/tempo"); rexclude = [r"problema=sladek03_apoiado"])
apoiado =
    collect_results(datadir("simulations/tempo"); rinclude = [r"problema=sladek03_apoiado"])

function plot_dis(data; figure_padding = (2, 6, 1, 6), legend_margin = ((5, 0, 0, 0)))


    fig = Figure(figure_padding = figure_padding)
    ax = Axis(fig[1, 1], xlabel = L"t/t_o", ylabel = L"w/w_{sta}")
    for d in eachrow(data)
        # d=engastado[1,:]
        if d.metodo == "Monta_M_RIM"
            met = "RIM"
        elseif d.metodo == "Monta_M_RIMd"
            met = "DIM"
        elseif d.metodo == "DRM"
            met = "DRM"
        end
        if maximum(abs.(d.dis)) > 10
            continue
        end
        lines!(ax, d.t, d.dis)
        i = rand(0:10)
        CairoMakie.scatter!(
            ax,
            d.t[1+i:10:end],
            d.dis[1+i:10:end],
            label = "$met - Number of boundary nodes  = $(d.nelem*12)",
        )
        # CairoMakie.scatter!(ax, d.t[1+i:10:end], d.dis[1+i:10:end], label="$met - ni = $(d.NPX^2) - nb = $(d.nelem*12)")

    end
    include(datadir("simulations\\sladek.jl"))

    if data[1, :problema] == "sladek03_apoiado"
        # @infiltrate
        lines!(
            mlpgapoiado[:, 1],
            mlpgapoiado[:, 2],
            label = "MLPG - Sladek et al. (2006)",
            linestyle = ".",
        )
        lines!(
            femapoiado[:, 1],
            femapoiado[:, 2],
            label = "FEM - Sladek et al. (2006)",
            linestyle = "-",
        )
    else
        lines!(
            mlpgengastado[:, 1],
            mlpgengastado[:, 2],
            label = "MLPG - Sladek et al. (2006)",
            linestyle = ".",
        )
        lines!(
            femengastado[:, 1],
            femengastado[:, 2],
            label = "FEM - Sladek et al. (2006)",
            linestyle = "-",
        )
    end
    Legend(fig[1, 2], ax, merge = true)
    return fig
end
# myplot(apoiado, 10)
function filtrov(npx, metodo, nelem, npg, nt, v)::Bool
    tnpx = npx in v[1]
    tmetodo = metodo in v[2]
    tnelem = nelem in v[3]
    tnpg = npg in v[4]
    tnt = nt in v[5]
    # @infiltrate
    tnpx && tmetodo && tnelem && tnpg && tnt
end
vals = [[3 11], ["Monta_M_RIM" "Monta_M_RIMd" "DRM"], [1 3 5], [10], [80]]
# vals = [[1], ["Monta_M_RIM"], [2], [10], [80]]
filtro(npx, metodo, nelem, npg, nt) = filtrov(npx, metodo, nelem, npg, nt, vals)
apoiadof = filter([:NPX, :metodo, :nelem, :npg, :nt] => filtro, apoiado)
myplot_web() = plot_dis(apoiadof)
save("apoiadonb.pdf", with_theme(myplot_web, theme_web(width = 800)))
# engastadof = filter([:NPX, :metodo, :nelem, :npg, :nt] => filtro, engastado)
# myplot_web() = plot_dis(engastadof)
# save("engastadonb.pdf", with_theme(myplot_web, theme_web(width=800)))

# function plot_tempo(data; figure_padding=(2, 6, 1, 6), legend_margin=((5, 0, 0, 0)))
#     fig = Figure(figure_padding=figure_padding)
#     a1 = subset(data, :metodo => ByRow(==("Monta_M_RIM")))
#     dofs1 = a1[:, :NPX] .^ 2 + a1[:, :nelem] * 12
#     tm1 = a1[:, :timerM]
#     a2 = subset(data, :metodo => ByRow(==("Monta_M_RIMd")))
#     dofs2 = a2[:, :NPX] .^ 2 + a2[:, :nelem] * 12
#     tm2 = a2[:, :timerM]
#     a3 = subset(data, :metodo => ByRow(==("DRM")))
#     dofs3 = a3[:, :NPX] .^ 2 + a3[:, :nelem] * 12
#     tm3 = a3[:, :timerM]

#     ax = Axis(fig[1, 1], yscale=log10, xscale=log10, xlabel="Number of points", ylabel="Time(s)")
#     scatter!(dofs1, tm1, label="RIM")
#     # ax = Axis(scene, xlabel="Number of points", ylabel="Time(s)")
#     scatter!(dofs2, tm2, label="DIM")
#     scatter!(dofs3, tm3, label="DRM")
#     Legend(fig[1, 2], ax, merge=true)
#     return fig
# end
# myplot_web() = plot_tempo(apoiadof)

# save("tempo.pdf", with_theme(myplot_web, theme_web(width=800)))
