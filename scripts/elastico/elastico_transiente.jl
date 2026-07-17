## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
begin

    # for nelem in [10, 20, 30, 40, 50]  #Numero de elementos
    nelem = 50  #Numero de elementos
    NPX = 70 #pontos internos na direção x
    NPY = NPX#pontos internos na direção y
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
end

begin
    Δt = 1e-4
    # for Δt in [1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5] #step de tempo
    println("nc(dad) = $(nc(dad)), ni(dad) = $(ni(dad)), dt = $Δt")

    tf = 0.1
    t = 0:Δt:tf
    pontos_t = length(t)

    Mt = M / (Δt^2) * ρ
    #estático
    # A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
    # x = A \ (b)
    # u, t = separa(dad, x) #importante    
    # u[ind, 1], (1 - .25^2) / 1e5


    A, b = BEM.aplicaCDC(H + 2 * Mt, G, dad) #Houbolt
    u = zeros((nc(dad) + ni(dad)) * 2, pontos_t)
    du = zeros((nc(dad) + ni(dad)) * 2, pontos_t)
    FA = lu(A)
    BEM.@showprogress "passo de tempo" for i = 4:pontos_t
        x = FA \ (b + Mt * (5 * u[:, i-1] - 4 * (u[:, i-2]) + 1 * (u[:, i-3])))
        ui, ti = separa(dad, x) #importante    
        if abs(ui[ind]) > 1e-3
            u .= NaN
            println("erro")
            break
        end
        u[1:2:(2*nc(dad)), i] = ui[:, 1]
        u[2:2:(2*nc(dad)), i] = ui[:, 2]
        u[(2*nc(dad)+1):end, i] = x[(2*nc(dad)+1):end]
        du[:, i] = (2u[:, i] - 5 * u[:, i-1] + 4 * (u[:, i-2]) - 1 * (u[:, i-3])) / (Δt)^2
    end
    #     end
    # end

    function c1(E, v, ρ)
        sqrt(E * (1 - v) / (ρ * (1 + v) * (1 - 2v)))
    end

    function omega(n, E, v, ρ, L)
        (2n - 1) * pi * c1(E, v, ρ) / (2L)
    end

    function u_ana(x2, t, L, P0, E, v, ρ; N = 100)
        s = 0.0
        for n = 1:N
            ωn = omega(n, E, v, ρ, L)
            s +=
                ((-1)^(n - 1)) / (2n - 1)^2 *
                sin((2n - 1) * pi * x2 / (2L)) *
                (1 - cos(ωn * t))
        end
        return (8 * L * P0) / (pi^2 * E) * s
    end

    function sigma_ana(x2, t, L, P0, E, v, ρ; N = 100)
        s = 0.0
        for n = 1:N
            ωn = omega(n, E, v, ρ, L)
            s +=
                ((-1)^(n - 1)) / (2n - 1) *
                cos((2n - 1) * pi * x2 / (2L)) *
                (1 - cos(ωn * t))
        end
        return (4 * P0) / pi * s
    end

    L = 4.0
    P0 = 1.0
    E = 1e5
    v = 0.25
    ρ = 1.0
    x2 = L

    w = [omega(n, E, v, ρ, L) for n = 1:10]

    tvals = range(0.0, tf, length = 2000)
    uvals = [u_ana(x2, t, L, P0, E, v, ρ, N = 30) for t in tvals]
    sigmavals = [sigma_ana(0, t, L, P0, E, v, ρ, N = 10000) for t in tvals]

    # lines(tvals,uvals, color = :blue, linewidth = 2, label = "Displacement u(x₂,t)")
    # lines(tvals,sigmavals, color = :red, linewidth = 2, label = "Stress σ(x₂,t)")


    fig = Figure(size = (900, 500))
    ax = Axis(
        fig[1, 1],
        xlabel = "t (s)",
        ylabel = "u(x₂,t)",
        title = "Δt = $Δt s nb=$(nc(dad)) ni=$(ni(dad))",
    )

    l1 = lines!(ax, tvals, uvals, color = :blue, linewidth = 2, label = "Analytical")
    l2 = lines!(ax, t[1:(end-2)], u[2*ind-1, 3:end], label = "BEM Houbolt")

    axislegend(ax)

    # fig
    save("elastico_transiente_Δt$(Δt)_nb$(nc(dad))_ni$(ni(dad)).png", fig)

end
# lines(t, u[2ind-1, :])
# lines(t, du[2ind-1, :])

