function mostra_geometria(dad)
    f = Figure(size = (800, 800), Aspect = 1)
    f[1, 1] = Axis(f)
    pos = f[1, 1]
    s = 0.1 * maximum(dad.NOS)
    for el in dad.ELEM
        r = range(-1, stop = 1, length = 10)
        x = dad.NOS[el.indices, :]   # Coordenada (x,y) dos nós geométricos
        N_geo = BEM.calc_fforma.(r, Ref(el))
        ps = hcat([N_geo[i][1] for i = 1:10]...)' * x
        # lines(ps[:, 1], ps[:, 2])
        lines!(pos, ps[:, 1], ps[:, 2], color = :black)
        xi = (ps[2, 1] - ps[1, 1]) / 2
        yi = (ps[2, 2] - ps[1, 2]) / 2
        lines!(
            pos,
            [ps[1, 1] - yi, ps[1, 1] + yi],
            [ps[1, 2] - xi, ps[1, 2] + xi],
            color = :black,
            linestyle = :dashdot,
        )

        xe = (ps[end, 1] - ps[end-1, 1]) / 2
        ye = (ps[end, 2] - ps[end-1, 2]) / 2
        lines!(
            pos,
            [ps[end, 1] - ye, ps[end, 1] + ye],
            [ps[end, 2] - xe, ps[end, 2] + xe],
            color = :black,
            linestyle = :dashdot,
        )
        if typeof(dad) == elastico
            if el.tipoCDC[1] == 0
                scatter!(
                    pos,
                    x[:, 1],
                    x[:, 2],
                    color = :blue,
                    marker = :rtriangle,
                    markersize = 30,
                )
                # arrows!(pos, x[:, 1], x[:, 2], fill(s, size(x, 1)), fill(0, size(x, 1)), color=:blue)
            elseif el.tipoCDC[1] == 1 && el.valorCDC[1] != 0
                # scatter!(pos, x[:, 1], x[:, 2], color=:red)
                arrows!(
                    pos,
                    x[:, 1],
                    x[:, 2],
                    sign(el.valorCDC[1]) * fill(s, size(x, 1)),
                    fill(0, size(x, 1)),
                    color = :red,
                )
            end
            if el.tipoCDC[2] == 0
                scatter!(
                    pos,
                    x[:, 1],
                    x[:, 2],
                    color = :blue,
                    marker = :utriangle,
                    markersize = 30,
                )
                # arrows!(pos, x[:, 1], x[:, 2], fill(0, size(x, 1)), fill(s, size(x, 1)), color=:blue)
            elseif el.tipoCDC[2] == 1 && el.valorCDC[2] != 0
                # scatter!(pos, x[:, 1], x[:, 2], color=:red)
                arrows!(
                    pos,
                    x[:, 1],
                    x[:, 2],
                    fill(0, size(x, 1)),
                    sign(el.valorCDC[2]) * fill(s, size(x, 1)),
                    color = :red,
                )
            end
            if el.tipoCDC == 0
                scatter!(pos, x[:, 1], x[:, 2], color = :blue)
            elseif el.tipoCDC == 1
                scatter!(pos, x[:, 1], x[:, 2], color = :red)
            end
        end
        # @infiltrate
    end
    ps = dad.pontos_internos
    if isempty(ps)
        return f
    end
    scatter!(pos, ps[:, 1], ps[:, 2], color = :black)
    ps = dad.NOS
    scatter!(pos, ps[:, 1], ps[:, 2], color = :black)

    DataInspector(f)
    f
end

function mostra_resultado(dad, Ts; levels = 10)
    pts = [dad.NOS; dad.pontos_internos]

    tri = triangulate(pts')
    # @infiltrate
    # triin.segmentlist=reshape(reduce(vcat, [el.indices for el in dad.ELEM]')',2,:)
    # (triout, vorout) = triangulate("pQ", triin)
    # mesh(fig[1, 1][1, 1], triout.pointlist', triout.trianglelist', color=Ts, shading=false)
    # Colorbar(fig[1, 1][1, 2])
    # Colorbar(scene)
    f, ax, tr = tricontourf(tri, Ts, levels = levels)
    Colorbar(f[1, 2], tr)
    f
end
function mostra_deformação(dad, u; escala = 1)
    f = Figure(size = (800, 800), Aspect = 1)
    f[1, 1] = Axis(f)
    pos = f[1, 1]
    for el in dad.ELEM
        r = range(-1, stop = 1, length = 10)
        x = dad.NOS[el.indices, :]   # Coordenada (x,y) dos nós geométricos
        N_geo = BEM.calc_fforma.(r, Ref(el))
        ps = hcat([N_geo[i][1] for i = 1:10]...)' * x
        # @infiltrate
        psd = hcat([N_geo[i][1] for i = 1:10]...)' * (x + escala * u[el.indices, :])
        # lines(ps[:, 1], ps[:, 2])
        lines!(pos, ps[:, 1], ps[:, 2], color = :black)
        lines!(pos, psd[:, 1], psd[:, 2], color = :red)
    end
    DataInspector(f)
    f
end
