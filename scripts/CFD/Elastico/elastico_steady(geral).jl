## Início da análise
using DrWatson, Plots, CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(srcdir("CFD.jl"))

#Dados de entrada

nelem = 15 #Numero de elementos
order = 2
nelem_circle = 20

NPX = 15 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

L = 1
Re = 100 # Na verdade é o 1/ν ou 1/μ

caso = "Cavidade"
λ = 10^5


## Formatação dos dados ________________________________________________

begin

    if caso == "Cavidade"
        dad = format_dad(cavity(nelem, order, L, Re, λ), NPX, NPY)

    elseif caso == "Von Karman"
        r = 0.1
        L = 30 * r
        h = 10 * r
        xc = 10 * r

        dad = format_dad(von_karman(nelem, order, L, h, r, xc, Re), NPX, NPY)

    elseif caso == "Expansion"
        L1 = 0.6
        h1 = 0.5
        L2 = 2.5
        h2 = 1
        dad = format_dad(expansion(nelem, order, L1, h1, L2, h2, Re, λ))
    elseif caso == "Contração"
        L1 = 1.5
        h1 = 1
        L2 = 1.5
        h2 = 0.5
        dad = format_dad(contraction(nelem, order, L1, h1, L2, h2, Re, λ))

    elseif caso == "Facing step"
        h = L / 3
        y_step = h / 2
        x_step = L / 5
        dad = format_dad(facing_step(nelem, order, L, h, y_step, x_step, Re, λ), NPX, NPY)

    elseif caso == "Canal entrada parabólica"

        h = L
        dpdx = -10^3

        dad = format_dad(channel(nelem, order, h, L, Re, λ), NPX, NPY)

        #Corrige CDC

        for i = 3*nelem+1:nelem*4
            n_indices = length(dad.ELEM[i].indices)
            for j = 1:n_indices
                indice = dad.ELEM[i].indices[j]
                dad.ELEM[i].valorCDC[1, j] =
                    dpdx * h^2 / (2 * 1) * (((dad.NOS[indice, 2] - h / 2) / h)^2 - 1 / 4)
            end
        end

    elseif caso == "Gradiente de pressão"

        dpdx = -10^3
        u_ana(y) = dpdx * L^2 / (2 * 1) * ((y / L)^2 - (y / L))
        dad = format_dad(gradiente_pressao(nelem, order, L, Re, λ), NPX, NPY)

    elseif caso == "Canal"

        h = 1

        dad = format_dad(channel(nelem, order, h, L, Re, λ), NPX, NPY)

    elseif caso == "Canal Simetria"
        h = L / 3
        dad = format_dad(channel_simetria(nelem, order, h, L, Re, λ), NPX, NPY)

    elseif caso == "Placa plana"

        h = 10 * L
        Ls = L / 4
        dad = format_dad(placa_plana_2(nelem, order, L, Ls, h, Re, λ), NPX, NPY)

    elseif caso == "Couette"
        u_ana(y) = y / L
        dad = format_dad(couette(nelem, order, L, Re, λ), NPX, NPY)

    elseif caso == "Taylor Couette"

        ra = 0.5
        ω_a = -1
        ω_b = 1
        rb = 2

        u_ana(r) =
            (ω_b * rb^2 - ω_a * ra^2) / (rb^2 - ra^2) * r +
            ((ω_a - ω_b) * ra^2 * rb^2 / (rb^2 - ra^2)) / r

        dad = format_dad(Taylor_Couette(nelem, order, ra, rb, Re, λ), NPX, NPY)

        # Corrige CDC ____________________________

        tangencial = zeros(length(dad.normal), 2)

        for i = 1:length(dad.normal[:, 1])
            tangencial[i, 1] = -dad.normal[i, 2]
            tangencial[i, 2] = dad.normal[i, 1]
        end

        # Cilindro interno ____________________________

        for i = nelem*3+1:nelem*4
            n_indices = length(dad.ELEM[i].indices)
            for j = 1:n_indices
                indice = dad.ELEM[i].indices[j]
                dad.ELEM[i].valorCDC[1, j] = (ω_a * ra) * tangencial[indice, 1]
                dad.ELEM[i].valorCDC[2, j] = (ω_a * ra) * tangencial[indice, 2]
            end
        end

        # Cilindro externo ____________________________

        for i = nelem+1:nelem*2
            n_indices = length(dad.ELEM[i].indices)
            for j = 1:n_indices
                indice = dad.ELEM[i].indices[j]
                dad.ELEM[i].valorCDC[1, j] = -(ω_b * rb) * tangencial[indice, 1]
                dad.ELEM[i].valorCDC[2, j] = -(ω_b * rb) * tangencial[indice, 2]
            end
        end

        #______________

        # Pontos internos ____________________________

        ntheta = 10
        nr = 10

        theta = range(0, pi / 2, length = ntheta)
        radius = range(ra, rb, length = nr)

        dad.pontos_internos = zeros(nr * ntheta, 2)

        for i = 1:nr
            for j = 1:ntheta
                dad.pontos_internos[(i-1)*nr+j, 1] = radius[i] * cos(theta[j])
                dad.pontos_internos[(i-1)*nr+j, 2] = radius[i] * sin(theta[j])
            end
        end

        #______________

    elseif caso == "Gravidade"

        g = 9.81
        h = 0.1
        L = 5 * h

        u_ana(y) = g / 1 * h^2 * ((y / h) - 1 / 2 * (y / h)^2)
        τ_ana(y) = 1 * g / 1 * h^2 * (1 / h - y / h^2)

        dad = format_dad(gravidade(nelem, order, L, h, g, Re, λ), NPX, NPY)

    elseif caso == "Step"

        h = L / 2
        L_step = L / 8
        x_step = L / 4

        dad = format_dad(step(nelem, order, L, h, L_step, x_step, Re, λ), NPX, NPY)
    end
end

# Refinamento (opcional)

#=
δx = 0.5
δy = 0.05
dad.pontos_internos = refina_global(dad,δx,δy)
=#

#=
#Von Karman
retangulo = [0.7, -0.25, 1.4, 0.25]
dad.pontos_internos = refina_local(dad,0.025,0.025,retangulo,true,1:4*order*nelem)
=#

#__________________________

# ==============Matrizes===============#

begin
    H, G = calc_HeG(dad, npg, interno = true)

    M = BEM.Monta_M_RIMd(dad, npg)

    A, b = BEM.aplicaCDC(H, G, dad)
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

    have_nolinear = true

    # IMPORTANTE
    constant_term = zeros(2 * (nc(dad) + ni(dad))) # Sem termo constante 
    #constant_term = -dpdx*[i % 2 == 1 ? 1 : 0 for i in 1:2*(nc(dad)+ni(dad))] #Caso Grad de pressão
    #constant_term = -g*[i % 2 == 1 ? 1 : 0 for i in 1:2*(nc(dad)+ni(dad))] #Caso Gravidade

    iter = 0
    erro_vel = 1
    #C = zeros(size(M))

    u = zeros(ni(dad) + nc(dad), 2)
    t = zeros(nc(dad), 2)
    u_before = deepcopy(u)
    nolinear = zeros(2 * (ni(dad) + nc(dad)))
    p = zeros(ni(dad) + nc(dad))

    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)

    relaxation = 0.1

    dist_cont_cont = zeros(nc(dad))

    for i = 1:nc(dad)
        if i == nc(dad)
            dist_cont_cont[i] = sqrt(
                (dad.NOS[i, 1] - dad.NOS[i-1, 1])^2 + (dad.NOS[i, 2] - dad.NOS[i-1, 2])^2,
            )
        else
            dist_cont_cont[i] = sqrt(
                (dad.NOS[i+1, 1] - dad.NOS[i, 1])^2 + (dad.NOS[i+1, 2] - dad.NOS[i, 2])^2,
            )
        end
    end

    index_dist_cont_int = zeros(nc(dad))

    for i = 1:nc(dad)
        dist =
            sqrt.(
                (dad.NOS[i, 1] .- dad.pontos_internos[:, 1]) .^ 2 +
                (dad.NOS[i, 2] .- dad.pontos_internos[:, 2]) .^ 2
            )
        index_dist_cont_int[i] = findfirst(x -> x == minimum(dist[dist.>1e-3]), dist)
    end
    index_dist_cont_int = Int.(index_dist_cont_int)


    dNx, dNy =
        BEM.montaFs([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])[2:3]
    Deriva_x = Dₓ(nx, ny, dx)
    Deriva_y = Dy(nx, ny, dy)

    # turbulencia

    #=
    lm = zeros(nc(dad)+ni(dad))
    wall_1 = 1:order*nelem*2
    wall_2 = order*nelem*3+1:order*nelem*4
    indices_wall = vcat(wall_1,wall_2)
    d = zeros(length(indices_wall),ni(dad))

    for i=1:ni(dad)
        d[:,i] = sqrt.((dad.NOS[indices_wall,1] .- dad.pontos_internos[i,1]).^2 .+ (dad.NOS[indices_wall,2] .- dad.pontos_internos[i,2]).^2)
    end

    dmax = maximum(d)
    lm[contorno] .= 0.09*dmax

    for i=1:ni(dad)
        lm[nc(dad)+i] = min(0.419*minimum(d[:,i]),0.09*dmax)
    end

    global correction = zeros(size(nolinear))
    =#
end

# Se tiver genérico

#=
begin 
    index_forward = Array{Union{Float64, Nothing}}(undef, ni(dad), 2)
    index_backward = Array{Union{Float64, Nothing}}(undef, ni(dad), 2)

    for i=1:ni(dad)
        index_forward[i,1] = indice_forward(dad,i,"x",x_order,y_order)
        index_backward[i,1] = indice_back(dad,i,"x",x_order,y_order)
        index_forward[i,2] = indice_forward(dad,i,"y",x_order,y_order)
        index_backward[i,2] = indice_back(dad,i,"y",x_order,y_order)
    end

    index_dist_int_int = zeros(ni(dad))
    for i=1:ni(dad)
        dist = abs.(dad.pontos_internos[i,2] .- dad.pontos_internos[:,2])
        index_dist_int_int[i] = findfirst(x -> x == minimum(dist[dist .> 1e-3]), dist)
    end

end
=#


#______________________________


# ==============Solução===============#

begin

    if have_nolinear
        if iter == 0
            global iter, u, t, erro, u_before, nx, ny, nolinear
        end

        while erro_vel > 10^-10


            # Derivadas

            #=
            dudx_int = Deriva_x*u[interno,1]
            dvdx_int = Deriva_x*u[interno,2]
            dudy_int = Deriva_y*u[interno,1]
            dvdy_int = Deriva_y*u[interno,2]
            =#

            #=
            dudx_cont,dudy_cont = deriva_contorno(dad,u[contorno,1],t[:,1],dist_cont_cont) 
            dvdx_cont,dvdy_cont = deriva_contorno(dad,u[contorno,2],t[:,2],dist_cont_cont)
            =#

            #=
            dudx_int,dudy_int=deriva_interno_generico(dad,u[interno,1],dx,dy,index_backward,index_forward,index_dist_int_int)
            dvdx_int,dvdy_int=deriva_interno_generico(dad,u[interno,2],dx,dy,index_backward,index_forward,index_dist_int_int)
            =#

            #=
            dudx_cont,dudy_cont = deriva_contorno_generico(dad,u[contorno,1],t[:,1],dist_cont_cont,4*nelem_circle*order)
            dvdx_cont,dvdy_cont = deriva_contorno_generico(dad,u[contorno,2],t[:,2],dist_cont_cont,4*nelem_circle*order)
            =#

            #=
            dudx_cont,dudy_cont = deriva_contorno_com_interno(dad,u[contorno,1],u[interno,1],dudx_int,dudy_int,index_dist_cont_int,dist_cont_cont)
            dudx_cont,dudy_cont = deriva_contorno_com_interno(dad,u[contorno,2],u[interno,2],dudx_int,dudy_int,index_dist_cont_int,dist_cont_cont)
            =#


            dudx = dNx * u[:, 1]
            dudy = dNy * u[:, 1]

            dvdx = dNx * u[:, 2]
            dvdy = dNy * u[:, 2]

            nolinear[1:2:2*n] = u[:, 1] .* dudx + u[:, 2] .* dudy
            nolinear[2:2:2*n] = u[:, 1] .* dvdx + u[:, 2] .* dvdy

            # Turbulence (RANS)
            #=
            G = 2*(dudx.^2 + dvdy.^2) + (dudy + dvdx).^2
            vt = lm.^2 .* sqrt.(G)

            correction[1:2:2*n] = (1 ./(1 .+ vt*Re))
            correction[2:2:2*n] = (1 ./(1 .+ vt*Re))
            =#

            #_____________________________________________

            # Solution

            x = A \ (b - M * nolinear)

            #x = A \ (b - M*(nolinear.*correction)) #RANS (muda as matrizes A, b e M tbm)

            #=
            u[1:nc(dad),:,i],t = separa(dad,x)

            u[nc(dad)+1:end,1,i] =  x[2*nc(dad)+1:2:end]
            u[nc(dad)+1:end,2,i] =  x[2*nc(dad)+2:2:end]
            =#

            u[1:nc(dad), :] =
                u[1:nc(dad), :] * (1 - relaxation) + relaxation * separa(dad, x)[1]
            t = t * (1 - relaxation) + relaxation * separa(dad, x)[2]

            u[interno, 1] =
                u[interno, 1] * (1 - relaxation) + relaxation * x[2*nc(dad)+1:2:end]
            u[interno, 2] =
                u[interno, 2] * (1 - relaxation) + relaxation * x[2*nc(dad)+2:2:end]

            erro_vel = nrmse(u_before, u)

            println("Iteração: $iter; erro_vel = $erro_vel")

            u_before = deepcopy(u)
            iter = iter + 1
        end

    elseif constant_term == zeros(2 * (nc(dad) + ni(dad)))

        x = A \ (b)

        u[1:nc(dad), :], t = separa(dad, x)


        u[interno, 1] = x[2*nc(dad)+1:2:end]
        u[interno, 2] = x[2*nc(dad)+2:2:end]

    else
        x = A \ (b + M * (constant_term))

        u[1:nc(dad), :], t = separa(dad, x)

        u[interno, 1] = x[2*nc(dad)+1:2:end]
        u[interno, 2] = x[2*nc(dad)+2:2:end]
    end
end

#________________________

# Computando variáveis

begin
    uc = u[contorno, :]
    ui = u[interno, :]

    #=
    dudx,dudy=deriva_interno_generico(dad,u[interno,1],dx,dy,index_backward,index_forward,index_dist_int_int)
    dvdx,dvdy=deriva_interno_generico(dad,u[interno,2],dx,dy,index_backward,index_forward,index_dist_int_int)

    dudx_cont,dudy_cont = deriva_contorno_generico(dad,u[contorno,1],t[:,1],dist_cont_cont,4*nelem_circle*order)
    dvdx_cont,dvdy_cont = deriva_contorno_generico(dad,u[contorno,2],t[:,2],dist_cont_cont,4*nelem_circle*order)
    =#


    dudx = dNx * u[:, 1]
    dudy = dNy * u[:, 1]
    dvdx = dNx * u[:, 2]
    dvdy = dNy * u[:, 2]


    #=
    dudx_int,dudy_int=deriva_interno_generico(dad,u[interno,1],dx,dy,index_backward,index_forward,index_dist_int_int)
    dvdx_int,dvdy_int=deriva_interno_generico(dad,u[interno,2],dx,dy,index_backward,index_forward,index_dist_int_int)

    dudx_cont,dudy_cont = deriva_contorno_com_interno(dad,u[contorno,1],u[interno,1],dudx_int,dudy_int,index_dist_cont_int,dist_cont_cont)
    dudx_cont,dudy_cont = deriva_contorno_com_interno(dad,u[contorno,2],u[interno,2],dudx_int,dudy_int,index_dist_cont_int,dist_cont_cont)
    =#

    p = -λ * (dudx + dvdy)
end

# #=========Plots variando nelem=========#

# erro_average = zeros(8)
# erro_L2 = zeros(8)
# erro_max = zeros(8)

# nelem=5:5:40

# p = Plots.plot(layout = (1, 1), xaxis=:log10, yaxis=:log10, xlabel="N° de elem.", ylabel="Erro")

# for order = 2:2
#     for i = 1:8
#         erro_average[i], erro_L2[i], erro_max[i] = solve_elastico(nelem[i], order,NPX,NPY)
#     end

#     if order == 2
#         Plots.plot!(p[1], nelem, erro_average, label="Erro médio", marker=3)
#         Plots.plot!(p[1], nelem, erro_L2, label="Erro L2", marker=3)
#         Plots.plot!(p[1], nelem, erro_max, label="Erro máximo", marker=3)
#         title!(p[1], "Ordem $order")  
#     elseif order == 3
#         Plots.plot!(p[2], nelem, erro_average, label="Erro médio", marker=3)
#         Plots.plot!(p[2], nelem, erro_L2, label="Erro L2", marker=3)
#         Plots.plot!(p[2], nelem, erro_max, label="Erro máximo", marker=3)
#         title!(p[2], "Ordem $order")
#     else
#         Plots.plot!(p[3], nelem, erro_average, label="Erro médio", marker=3)
#         Plots.plot!(p[3], nelem, erro_L2, label="Erro L2", marker=3)
#         Plots.plot!(p[3], nelem, erro_max, label="Erro máximo", marker=3)
#         title!(p[3], "Ordem $order")    
#     end
# end

# display(p)

# #=
# Plots.plot(nelem, centrada, label="Erro médio central", marker=3,xlabel="N° de elem.", ylabel="Erro")
# Plots.plot!(nelem, up, label="Erro médio upwind", marker=3)
# Plots.plot!(nelem, back, label="Erro médio backwind", marker=3)
# =#

# #=========Heatmap e Quiver=========#

x_array = dad.pontos_internos[:, 1]
y_array = dad.pontos_internos[:, 2]

utotal = sqrt.(ui[:, 1] .^ 2 + ui[:, 2] .^ 2)

BEM.heatmap(x_array, y_array, utotal)
escala = 10^(-1)
BEM.quiver!(x_array, y_array, escala * ui[:, 1], escala * ui[:, 2], color = :black)
BEM.current_figure()
# #=========Contour=========#

# Plots.contour(x_unique, y_unique, utotal, color=:viridis, levels=10)

# #=======Plots=======#

# #Couette

# u_ana_array = u_ana.(dad.pontos_internos[:,2])
# u_num = ui[:,1]

# τ_ana_array = ones(ni(dad))
# τ_num = Dₓ_1D(ni(dad),dad.pontos_internos[2,2] - dad.pontos_internos[1,2])*u_num

# Plots.plot(dad.pontos_internos[:,2],u_ana_array, label="Analítico", xaxis="Velocidade [m/s]",yaxis="Posição vertical [m]")
# Plots.scatter!(dad.pontos_internos[:,2],u_num,marker=(3,:x,:red),label="Numérico")

# Plots.plot(dad.pontos_internos[:,2],τ_ana_array,color=:blue,label="Analítico", xaxis="τ [Pa]",yaxis="y [m]")
# Plots.scatter!(dad.pontos_internos[:,2],τ_num,marker=(3,:x,:red), label="Numérico")

# #Cavidade

# Plots.scatter(ui[Int((nx)/2):nx:end,1],dad.pontos_internos[Int((nx)/2):nx:end,2],
# marker=(3,:x,:red), label="Numérico")

# Plots.plot(dad.pontos_internos[Int((nx+1)/2):nx:end,1],-ui[Int((nx+1)/2):Int((nx+1)/2)+nx,2])

# # Dados no arquivo Cavidade_validation

# Plots.plot(cavidade_uy_ref[2:16,2],cavidade_uy_ref[2:16,1],title="Re = 100"
# ,label="Analítico", xaxis="Velocidade [m/s]",yaxis="Posição vertical [m]")

# Plots.plot(cavidade_uy_ref_400[2:16,2],cavidade_uy_ref_400[2:16,1],title="Re = 400"
# ,label="Analítico", xaxis="Velocidade [m/s]",yaxis="Posição vertical [m]")

# Plots.plot(cavidade_uy_ref_1000[2:16,2],cavidade_uy_ref_1000[2:16,1],title="Re = 1000"
# ,label="Analítico", xaxis="Velocidade [m/s]",yaxis="Posição vertical [m]")


# #Canal

# intervalo = 3*order*nelem+1:4*order*nelem
# y_plot= dad.NOS[intervalo,2]

# u_ana_array = 3/2*(1 .- 4*(y_plot/(2*h)).^2)
# u_num = uc[intervalo,1]


# Plots.plot(y_plot,u_ana_array, label="Analítico", yaxis="Velocidade [m/s]",xaxis="Posição vertical [m]")
# Plots.scatter(reverse(y_plot),u_num,label="Numérico",marker=(3,:x,:red)) # Velocidade na saída

# position = 5
# pos_frac = round(L*position/nx,digits=2)
# Plots.plot(ui[position:nx:end,1],label="x = $pos_frac",xlabel="Nó interno",ylabel="u [m/s]") # Velocidade na entrada


# #=======Comparação com Analítico==========#

# # Von Karman #Rodar o arquivo do FVM

# intervalo = 1:Int(4*nelem_circle*order)
# angle = atan.(dad.normal[intervalo,2],(dad.normal[intervalo,1]))

# intervalo_1 = 1:Int(2*nelem_circle*order)
# intervalo_2 = Int(2*nelem_circle*order)+1:Int(4*nelem_circle*order)

# name_plot = "BEM - ne = $nelem_circle"


# Plots.scatter(angle_VFM*180/pi,p_VFM[:,4]/maximum(p_VFM[:,4]),
# label="Volumes finitos",yaxis="Pressão [Pa]",xaxis="Ângulo [°]",marker=:o)

# Plots.scatter!(angle*180/pi,([reverse(pressure[intervalo_1]);reverse(pressure[intervalo_2])])/maximum(pressure[intervalo]),
# marker=:x,label=name_plot)


# #Plots.scatter(angle*180/pi,τ,label="Tensão cisalhante",yaxis="Tensão [Pa]",xaxis="Ângulo [°]")

# F_num =0
# for i in intervalo
#     F_num += τ[i]*δ[i]
# end
# F_num

# # Poseullie (Grad de pressão)

# intervalo = order*nelem+1:2*order*nelem
# y_plot= dad.NOS[intervalo,2]

# u_ana_array = u_ana.(y_plot)
# u_num = uc[intervalo,1]

# τ_ana_array = τ_ana.(y_plot)
# τ_num = Dₓ_1D(ni(dad),dad.pontos_internos[2,2] - dad.pontos_internos[1,2])*u_num

# Plots.plot(y_plot,u_ana_array, label="Analítico", yaxis="Velocidade [m/s]",xaxis="Posição vertical [m]")
# Plots.scatter!(y_plot,u_num/Re,marker=(3,:x,:red),label="Numérico")

# Plots.plot(y_plot,τ_ana_array,color=:blue,label="Analítico", yaxis="τ [Pa]",xaxis="y [m]")
# Plots.scatter!(y_plot,τ_num,marker=(3,:x,:red), label="Numérico")

# # Taylor couette

# intervalo = nelem*order*0+1:nelem*order*1

# x_num = dad.NOS[intervalo,1]
# u_ana_array = u_ana.(x_num)
# u_num = -uc[intervalo,2]

# Plots.plot(x_num,u_ana_array,xlabel="Raio",
# ylabel="u",label="Solução analítica", title="Taylor Couette")
# Plots.scatter!(x_num,u_num,label="Solução numérica",marker=(:x, 3),color=:red)

# # Gravidade

# intervalo = order*nelem+1:2*order*nelem
# y_plot= dad.NOS[intervalo,2]
# u_ana_array = u_ana.(y_plot)
# u_num = -uc[intervalo,1]

# Plots.plot(y_plot,u_ana_array,xlabel="y [m]",
# ylabel="u [m/s]",label="Solução analítica",color=:blue, title="Gravidade") # Velocidade na saída
# Plots.scatter!(y_plot,u_num,marker=(3,:x,:red),label="Solução numérica")

# #=======Erros==========#

# erro_average = nrmse(u_ana_array,reverse(u_num))
# erro_L2 = nme(u_ana_array,reverse(u_num))
# erro_max = maximum(abs.((reverse(u_num) .- u_ana_array)/(maximum(u_ana_array)-minimum(u_ana_array))))

# erro_average = nrmse(τ_ana_array,τ_num)
# erro_L2= nme(τ_ana_array,τ_num)
# erro_max = maximum(abs.((τ_num .- τ_ana_array)/(maximum(τ_ana_array)-minimum(τ_ana_array))))






