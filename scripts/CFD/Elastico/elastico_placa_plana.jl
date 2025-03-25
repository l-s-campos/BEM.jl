## Início da análise
using DrWatson, Plots, CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))

#CDC de derivadas zeros ou velocidade u=1 e v=0 no topo deu na mesma

#Dados de entrada

nelem = 20 #Numero de elementos
order = 2

NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

L = 0.1
Re = 1

λ = 10^6

## Formatação dos dados ________________________________________________

h=40*L
Ls=L/4
dad = format_dad(placa_plana_2(nelem, order,L,Ls,h,Re,λ), NPX, NPY)

#=
dx = 0.07
dy=dx*10
dad.pontos_internos = refina_global(dad,dx,dy)
=#


δx = 0.03*0.1
δy = 0.01
xi=0.26*0.1
yi=0.1
xf=1.1*0.1
yf=1
retangulo = [xi, yi, xf, yf]
dad.pontos_internos = refina_local(dad,δx,δy,retangulo,false,nothing)


#__________________________

# ==============Matrizes===============#

begin
    H, G = calc_HeG(dad, npg)  

    M = BEM.Monta_M_RIMd(dad, npg)

    A, b = BEM.aplicaCDC(H, G, dad)
end

# ==============Discretização com refinamento============#
#=
begin     
    nx_refinado = Int(round((xf-xi)/δx))+1
    ny_refinado = Int(round((yf - yi)/δy))+1

    inicio_refinado = ni(dad) - (nx_refinado)*(ny_refinado) + 1

    nx = length(unique(dad.pontos_internos[1:inicio_refinado-1,1]))
    ny = length(unique(dad.pontos_internos[1:inicio_refinado-1,2]))

    x_order = sort(unique(dad.pontos_internos[:,1]))
    y_order = sort(unique(dad.pontos_internos[:,2]))

    dx = (maximum(unique(dad.pontos_internos[1:inicio_refinado-1,1]))-minimum(unique(dad.pontos_internos[1:inicio_refinado-1,1])))/nx
    dy = (maximum(unique(dad.pontos_internos[1:inicio_refinado-1,2]))-minimum(unique(dad.pontos_internos[1:inicio_refinado-1,2])))/ny
end
=#

# ==============Discretização sem refinamento============#
begin     
    nx = length(unique(dad.pontos_internos[:,1]))
    ny = length(unique(dad.pontos_internos[:,2]))

    x_order = sort(unique(dad.pontos_internos[:,1]))
    y_order = sort(unique(dad.pontos_internos[:,2]))

    dx = (maximum(x_order)-minimum(x_order))/nx
    dy = (maximum(y_order)-minimum(y_order))/ny
    
end


# =========Inicia variáveis===============#

begin  

    have_nolinear = true

    # IMPORTANTE
    #constant_term = zeros(2*(nc(dad)+ni(dad))) # Sem termo constante 
    #constant_term = -dpdx*[i % 2 == 1 ? 1 : 0 for i in 1:2*(nc(dad)+ni(dad))] #Caso Grad de pressão
    #constant_term = -g*[i % 2 == 1 ? 1 : 0 for i in 1:2*(nc(dad)+ni(dad))] #Caso Gravidade

    iter = 0
    erro = 1
    C = zeros(size(M))

    u = zeros(nc(dad) + ni(dad),2)
    t = zeros(nc(dad),2)
    u_before = deepcopy(u)
    nolinear = zeros(2*(ni(dad)+nc(dad)))
    p = zeros(ni(dad)+nc(dad))

    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)
    relaxation = 1

    dist_cont_cont = zeros(nc(dad))

    for i=1:nc(dad)
        if i == nc(dad)
            dist_cont_cont[i] = sqrt((dad.NOS[i,1]-dad.NOS[i-1,1])^2 + (dad.NOS[i,2]-dad.NOS[i-1,2])^2)
        else
            dist_cont_cont[i] = sqrt((dad.NOS[i+1,1]-dad.NOS[i,1])^2 + (dad.NOS[i+1,2]-dad.NOS[i,2])^2)
        end
    end

end

# Se tiver genérico

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

#______________________________


# ==============Solução===============#

begin

    if have_nolinear

        while erro > 10^-10 && iter < 5000

            if iter == 0
                global iter,u,t,erro,u_before,nx,ny,nolinear,C
            end
        
            # Derivadas
            
            #dudx, dudy, dvdx, dvdy = deriva_ui(dad,u[interno,:],dx,dy)

            #=
            dudx = Dₓ(nx, ny, dx)*u[interno,1]
            dvdx = Dₓ(nx, ny, dx)*u[interno,2]
            dudy = Dy(nx, ny, dy)*u[interno,1]
            dvdy = Dy(nx, ny, dy)*u[interno,2]
            
            
            dudx_cont,dudy_cont = deriva_contorno(dad,u[contorno,1],t[:,1],dist_cont_cont) #Talvez seja vezes Re/3 no t, ou seja, divido por E
            dvdx_cont,dvdy_cont = deriva_contorno(dad,u[contorno,2],t[:,2],dist_cont_cont)
            =#

            
            dudx,dudy=deriva_interno_generico(dad,u[interno,1],dx,dy,index_backward,index_forward,index_dist_int_int)
            dvdx,dvdy=deriva_interno_generico(dad,u[interno,2],dx,dy,index_backward,index_forward,index_dist_int_int)

            dudx_cont,dudy_cont = deriva_contorno_generico(dad,u[contorno,1],t[:,1],dist_cont_cont)
            dvdx_cont,dvdy_cont = deriva_contorno_generico(dad,u[contorno,2],t[:,2],dist_cont_cont)
            

            #=
            dudx_cont,dudy_cont = deriva_contorno_com_interno(dad,u[contorno,1],u[interno,1],dudx,dudy)
            dvdx_cont,dvdy_cont = deriva_contorno_com_interno(dad,u[contorno,2],u[interno,2],dvdx,dvdy)
            =#
            
            
            nolinear[2*nc(dad)+1:2:2*(nc(dad)+ni(dad))] = u[interno,1].*dudx .+ u[interno,2].*dudy
            nolinear[2*nc(dad)+2:2:2*(nc(dad)+ni(dad))] = u[interno,1].*dvdx .+ u[interno,2].*dvdy

            nolinear[1:2:2*nc(dad)] = u[contorno,1].*dudx_cont .+ u[contorno,2].*dudy_cont
            nolinear[2:2:2*nc(dad)] = u[contorno,1].*dvdx_cont .+ u[contorno,2].*dvdy_cont
            

            #_____________________________________________
            
            # Solution

            x = A \ (b - M*nolinear)

            #=
            u[1:nc(dad),:,i],t = separa(dad,x)
            
            u[nc(dad)+1:end,1,i] =  x[2*nc(dad)+1:2:end]
            u[nc(dad)+1:end,2,i] =  x[2*nc(dad)+2:2:end]
            =#
        
            u[1:nc(dad),:] = u[1:nc(dad),:]*(1-relaxation) + relaxation*separa(dad,x)[1]
            t =  t*(1-relaxation) + relaxation*separa(dad,x)[2]
            
            u[interno,1] =  u[interno,1]*(1-relaxation) + relaxation*x[2*nc(dad)+1:2:end]
            u[interno,2] =  u[interno,2]*(1-relaxation) + relaxation*x[2*nc(dad)+2:2:end]

           
            erro = nrmse(u_before,u)
           

            println("Iteração: $iter; Erro = $erro")
        
            u_before = deepcopy(u)
            iter = iter+1
        end

    elseif constant_term == zeros(2*(nc(dad)+ni(dad)))

        x = A \ (b)
    
        u[1:nc(dad),:],t = separa(dad,x)
        
        u[interno,1] =  x[2*nc(dad)+1:2:end]
        u[interno,2] =  x[2*nc(dad)+2:2:end]

    else
        x = A \ (b + M*(constant_term))
    
        u[1:nc(dad),:],t = separa(dad,x)
        
        u[interno,1] =  x[2*nc(dad)+1:2:end]
        u[interno,2] =  x[2*nc(dad)+2:2:end]
    end
end

#________________________


#=========Heatmap e Quiver=========#

uc = u[contorno,:]
ui = u[interno,:]

x_array = dad.pontos_internos[:,1]
y_array = dad.pontos_internos[:,2]

utotal = sqrt.(ui[:,1].^2 + ui[:,2].^2)

int_refinado = inicio_refinado:ni(dad)

BEM.heatmap(x_array[int_refinado], y_array[int_refinado], utotal[int_refinado])
fig, ax, hm = BEM.heatmap(x_array[int_refinado], y_array[int_refinado], utotal[int_refinado])
BEM.Colorbar(fig[:, end+1], hm,label="Velocidade [m/s]")
fig

escala = 10^(-1)
BEM.quiver(x_array,y_array,escala*ui[:,1],escala*ui[:,2],color=utotal)

#=========Contour=========#

Plots.contour(x_unique, y_unique, utotal, color=:viridis, levels=10)

#=======Plots=======#

#Placa plana

espessura_polinomio(x) = sqrt(30*x/Re)

espessura_Blausius(x) = 5/sqrt(1/(x*1/Re))
τ_Blausius(x) = 0.332*1*sqrt(1/(Re*x))

espessura_num = zeros(nx_refinado)

for i=0:nx_refinado-1
    position_iter = inicio_refinado + i
    #diference = abs.(ui[position_iter:nx_refinado:end,1] .- 0.99)
    #index_diference = findfirst(x -> x == minimum(diference),diference)
    #espessura_num[i+1] = index_diference*δy
    for j=0:ny_refinado-1
        if ui[position_iter+nx_refinado*j,1] >= 0.99
            espessura_num[i+1] = j*δy
            break
        end
    end
end

x_array_placa = range(xi-Ls,xf-Ls,length=nx_refinado)
espessura_ana_polinomio = espessura_polinomio.(x_array_placa)
espessura_ana_Blausius = espessura_Blausius.(x_array_placa)

τ_ana_Blausius = τ_Blausius.(x_array_placa)
intervalo_placa = Int(nelem/2*order+1):Int((nelem/2+nelem*2)*order)
τ_num = t[intervalo_placa]

Plots.plot(x_array_placa,espessura_ana_polinomio/10, label="Polinomial aproximado", xaxis="Posição [m]",yaxis="δ [m]")
Plots.plot!(x_array_placa,espessura_ana_Blausius/10, label="Blausius exata", xaxis="Posição [m]",yaxis="δ [m]")
Plots.scatter!(x_array_placa,espessura_num,marker=(:x,:red),label="Numérico")


#Tensão de cisalhamento
Plots.plot(x_array_placa,τ_ana_Blausius, label="Blausius exata", xaxis="Posição [m]",yaxis="τ [Pa]")
Plots.scatter!(dad.NOS[intervalo_placa,1],-τ_num,marker=(:x,:red),label="Numérico")


Plots.plot(ui[inicio_refinado+5:nx_refinado:end,1])