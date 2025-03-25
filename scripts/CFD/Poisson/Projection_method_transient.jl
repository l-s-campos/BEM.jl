using DrWatson,Plots,CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))

nelem = 20  #Numero de elementos
order = 2

NPX = 30 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

#=============== CASO ===============#

caso = "Von Karman"
Re = 10

#Tempo

begin
    dt = 0.01
    t_i = 0
    t_f = 1

    t = range(t_i, stop=t_f, step=dt)
    nt = length(t)
end

# Tempo artificial do loop interno

tau = dt*10^-4

#_______________________________

# ======== INÍCIO ===============#

begin
    
    if caso == "Cavidade"

        dad_u = format_dad(cavity_poisson(nelem, order,Re,"u"), NPX, NPY)
        dad_v = format_dad(cavity_poisson(nelem, order,Re,"v"), NPX, NPY)
        dad_p = format_dad(cavity_poisson(nelem, order,1,"p"), NPX, NPY)

        println("Caso Cavidade, Reynolds = $Re")

    elseif caso == "Von Karman"

        r=0.01
        L=30*r
        h = 10*r
        xc = 10*r

        println("Caso Von Karman, Reynolds = $Re")

        dad_u = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"u"), NPX, NPY)
        dad_v = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"v"), NPX, NPY)
        dad_p = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"p"), NPX, NPY)
        
    end

    #________________________

    Ht = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),3)
    Gt = zeros(ni(dad_u)+nc(dad_u),nc(dad_u),3)
    M = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),3)

    #=
    A_Houbolt = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),3)
    b_Houbolt = zeros(ni(dad_u)+nc(dad_u),3)

    
    A_Adams = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Adams = zeros(ni(dad_u)+nc(dad_u),2)
    =#

    A_Euler = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Euler = zeros(ni(dad_u)+nc(dad_u),2) 
    
    # Velocidade u
    
    Ht[:,:,1], Gt[:,:,1] = BEM.calc_HeGt(dad_u)
    M[:,:,1] = BEM.Monta_M_RIMd(dad_u, npg)

    #=
    A_Houbolt[:,:,1], b_Houbolt[:,1] = BEM.aplicaCDC(Ht[:,:,1] - 11*M[:,:,1]/(6*tau), Gt[:,:,1], dad_u)
    
    A_Adams[:,:,1], b_Adams[:,1] = BEM.aplicaCDC(Ht[:,:,1] - (24/(55*tau))*M[:,:,1], Gt[:,:,1], dad_u)
    =#

    A_Euler[:,:,1], b_Euler[:,1] = BEM.aplicaCDC(Ht[:,:,1] - M[:,:,1]/tau, Gt[:,:,1], dad_u)
    
    # Velocidade v
    
    Ht[:,:,2], Gt[:,:,2] = BEM.calc_HeGt(dad_v)
    M[:,:,2] = BEM.Monta_M_RIMd(dad_v, npg)

    #=
    A_Houbolt[:,:,2], b_Houbolt[:,2] = BEM.aplicaCDC(Ht[:,:,2] - 11*M[:,:,2]/(6*tau), Gt[:,:,2], dad_v)
    
    A_Adams[:,:,2], b_Adams[:,2] = BEM.aplicaCDC(Ht[:,:,2] - (24/(55*tau))*M[:,:,2], Gt[:,:,2], dad_v)
    =#
    A_Euler[:,:,2], b_Euler[:,2] = BEM.aplicaCDC(Ht[:,:,2] - M[:,:,2]/tau, Gt[:,:,2], dad_v)
    

    # Pressão

    Ht[:,:,3], Gt[:,:,3] = BEM.calc_HeGt(dad_p)
    M[:,:,3] = BEM.Monta_M_RIMd(dad_p, npg)
    A, b = BEM.aplicaCDC(Ht[:,:,3], Gt[:,:,3], dad_p)
end

#_______________________________

# Discretization

begin
    nx = length(unique(dad_u.pontos_internos[:,1]))
    ny = length(unique(dad_u.pontos_internos[:,1]))

    x_order = sort(unique(dad_u.pontos_internos[:,1]))
    y_order = sort(unique(dad_u.pontos_internos[:,2]))

    dx = (maximum(x_order)-minimum(x_order))/nx
    dy = (maximum(y_order)-minimum(y_order))/ny

    Cₒ =  dt*1/dx
    println("Número de Courant:$Cₒ")
end

#_______________

# Iniciando variáveis

begin

    u = zeros(ni(dad_u)+nc(dad_u),2,nt)
    p = zeros(ni(dad_u)+nc(dad_u),nt)

    dpdx = zeros(nc(dad_u)+ni(dad_u))
    dpdy = zeros(nc(dad_u)+ni(dad_u))
    p_flux = zeros(nc(dad_u))
    t = zeros(nc(dad_u),2)

    u_new = zeros(ni(dad_u)+nc(dad_u),2)
    u_prev = zeros(ni(dad_u)+nc(dad_u),2)

    p_new = zeros(ni(dad_u)+nc(dad_u))
    p_prev = zeros(ni(dad_u)+nc(dad_u))

    interno = nc(dad_u)+1:nc(dad_u)+ni(dad_u)
    contorno = 1:nc(dad_u)

    nolinear = zeros(2*(ni(dad_u)+nc(dad_u)))
    unsteady = zeros(2*(ni(dad_u)+nc(dad_u)))
end

#_______________

# Se tiver genérico

begin 
    index_forward = Array{Union{Float64, Nothing}}(undef, ni(dad_u), 2)
    index_backward = Array{Union{Float64, Nothing}}(undef, ni(dad_u), 2)

    for i=1:ni(dad_u)
        index_forward[i,1] = indice_forward(dad_u,i,"x",x_order,y_order)
        index_backward[i,1] = indice_back(dad_u,i,"x",x_order,y_order)
        index_forward[i,2] = indice_forward(dad_u,i,"y",x_order,y_order)
        index_backward[i,2] = indice_back(dad_u,i,"y",x_order,y_order)
    end

    index_dist_int_int = zeros(ni(dad_u))
    for i=1:ni(dad_u)
        dist = abs.(dad_u.pontos_internos[i,2] .- dad_u.pontos_internos[:,2])
        index_dist_int_int[i] = findfirst(x -> x == minimum(dist[dist .> 1e-3]), dist)
    end

    dist_cont_cont = zeros(nc(dad_u))

    for i=1:nc(dad_u)
        if i == nc(dad_u)
            dist_cont_cont[i] = sqrt((dad_u.NOS[i,1]-dad_u.NOS[i-1,1])^2 + (dad_u.NOS[i,2]-dad_u.NOS[i-1,2])^2)
        else
            dist_cont_cont[i] = sqrt((dad_u.NOS[i+1,1]-dad_u.NOS[i,1])^2 + (dad_u.NOS[i+1,2]-dad_u.NOS[i,2])^2)
        end
    end
end

#Loop solução

for i=4:nt

    iter = 0
    erro=1

    if i==4
        global nolinear,p,p_flux,u,t,dpdx,dpdy,
        u_new,u_prev,p_new,p_prev,unsteady
    end

    while erro > 10^-6

        iter+=1

        # Derivadas do campo de velocidade

        #=
        dudx = Dₓ(nx,ny,dx)*u_prev[interno,1]
        dudy = Dy(nx,ny,dx)*u_prev[interno,1]
        dvdx = Dₓ(nx,ny,dy)*u_prev[interno,2]
        dvdy = Dy(nx,ny,dy)*u_prev[interno,2]

        dudx_cont,dudy_cont = deriva_contorno(dad_u,u_prev[contorno,1],t[:,1])
        dvdx_cont,dvdy_cont = deriva_contorno(dad_v,u_prev[contorno,2],t[:,2])
        =#

        
        dudx,dudy=deriva_interno_generico(dad_u,u_prev[interno,1],dx,dy,index_backward,index_forward,index_dist_int_int)
        dvdx,dvdy=deriva_interno_generico(dad_v,u_prev[interno,2],dx,dy,index_backward,index_forward,index_dist_int_int)

        dudx_cont,dudy_cont = deriva_contorno_generico(dad_u,u_prev[contorno,1],t[:,1],dist_cont_cont)
        dvdx_cont,dvdy_cont = deriva_contorno_generico(dad_v,u_prev[contorno,2],t[:,2],dist_cont_cont)
        

        nolinear[1:2:2*nc(dad_u)] = u_prev[contorno,1].*dudx_cont .+ u_prev[contorno,2].*dudy_cont
        nolinear[2:2:2*nc(dad_u)] = u_prev[contorno,1].*dvdx_cont .+ u_prev[contorno,2].*dvdy_cont

        nolinear[2*nc(dad_u)+1:2:2*(nc(dad_u)+ni(dad_u))] = u_prev[interno,1].*dudx .+ u_prev[interno,2].*dudy
        nolinear[2*nc(dad_u)+2:2:2*(nc(dad_u)+ni(dad_u))] = u_prev[interno,1].*dvdx .+ u_prev[interno,2].*dvdy
        
        
        unsteady[1:2:end] = (u_prev[:,1] - u[:,1,i-1])/dt
        unsteady[2:2:end] = (u_prev[:,2] - u[:,2,i-1])/dt 
        

        #=
        unsteady[1:2:end] = (u[:,1,i-1] - u[:,1,i-2])/dt
        unsteady[2:2:end] = (u[:,1,i-1] - u[:,2,i-2])/dt
        =#
        #_________________

        # Obtendo a velocidade

        u_new[contorno,1], u_new[interno,1], t[:,1] = 
        poisson_Euler_unsteady(dad_u, nolinear[1:2:end] + unsteady[1:2:end],u_prev[:,1],A_Euler[:,:,1],b_Euler[:,1],M[:,:,1])
    
        u_new[contorno,2], u_new[interno,2], t[:,2] = 
        poisson_Euler_unsteady(dad_v, nolinear[2:2:end] + unsteady[2:2:end],u_prev[:,2],A_Euler[:,:,2],b_Euler[:,2],M[:,:,2]) 

        # Obtendo a pressão

        # Derivadas atualizadas __________
        
        #=
        dudx = Dₓ(nx,ny,dx)*u_new[interno,1]
        dvdy = Dy(nx,ny,dy)*u_new[interno,2]
    
        dudx_cont,_ = deriva_contorno(dad_u,u_new[contorno,1],t[:,1])
        _,dvdy_cont = deriva_contorno(dad_v,u_new[contorno,2],t[:,2])
        =#

        
        dudx,dudy=deriva_interno_generico(dad_u,u_new[interno,1],dx,dy,index_backward,index_forward,index_dist_int_int)
        dvdx,dvdy=deriva_interno_generico(dad_v,u_new[interno,2],dx,dy,index_backward,index_forward,index_dist_int_int)

        dudx_cont,dudy_cont = deriva_contorno_generico(dad_u,u_new[contorno,1],t[:,1],dist_cont_cont)
        dvdx_cont,dvdy_cont = deriva_contorno_generico(dad_v,u_new[contorno,2],t[:,2],dist_cont_cont)
        

        # Achando a pressão

        p_new[contorno], p_new[interno], p_flux = poisson_steady(dad_p,1/tau*([dudx_cont;dudx] + [dvdy_cont;dvdy]),A,b,M[:,:,3])
        
        #=
        dpdx_int = Dₓ(nx,ny,dx)*p_new[interno]
        dpdy_int = Dy(nx,ny,dy)*p_new[interno]
        
        dpdx_cont,dpdy_cont = deriva_contorno(dad_u,p_new[contorno],p_flux)
        =#

        
        dpdx_int,dpdy_int=deriva_interno_generico(dad_u,p_new[interno],dx,dy,index_backward,index_forward,index_dist_int_int)

        dpdx_cont,dpdy_cont = deriva_contorno_generico(dad_u,p_new[contorno],t[:,1],dist_cont_cont)
        
    
        dpdx = [dpdx_cont;dpdx_int]
        dpdy = [dpdy_cont;dpdy_int]
        
        
        #Correção da velocidade
    
        u_new[interno,1] = u_new[interno,1] - tau*dpdx[interno]
        u_new[interno,2] = u_new[interno,2] - tau*dpdy[interno]

        erro = nrmse(u_prev,u_new)
        
        if iter % 10 == 0
            println("Erro: $erro")
        end

        p_prev = deepcopy(p_new)
        u_prev = deepcopy(u_new)

    end

    u[:,:,i] = u_new
    p[:,i] = p_new

    u_mean = mean(u[:,1,i])

    println("Iteração: $(i-3); u_mean = $u_mean")
end

#_______________

# Mapa de cores

x_array = dad_u.pontos_internos[:,1]
y_array = dad_u.pontos_internos[:,2]


uc = u[contorno,1,end]
ui = u[interno,1,end]
vi = u[interno,2,end]

pc = p[contorno,end]
p_int = p[interno,end]

utotal = sqrt.(ui.^2 + vi.^2)

BEM.heatmap(
    x_array,
    y_array,
    utotal,
    colormap = :viridis 
)

scale = 10^-2
BEM.quiver(x_array,y_array,scale*ui,scale*vi,color=utotal)


iter = 200
tempo = iter*dt
ui = u[interno,1,iter]
vi = u[interno,2,iter]

Plots.plot!(x_order,-vi[Int((nx)/2+1):Int((nx)/2)+nx],label="Tempo = $tempo s")

Plots.plot(x_order,-vi[Int((nx)/2+1):Int((nx)/2)+nx]
,yaxis="v [m/s]",xaxis="Posição horizontal [m]",label="Tempo = $tempo s") #Primeiro plot



Plots.plot(ui[Int((nx)/2):nx:end],dad_u.pontos_internos[Int((nx)/2):nx:end,2]
,xaxis="u [m/s]",yaxis="Posição vertical [m]",label="Tempo = $tempo s") #Primeiro plot


Plots.plot!(ui[Int((nx)/2):nx:end],dad_u.pontos_internos[Int((nx)/2):nx:end,2],label="Tempo = $tempo s")



# Filme
output_folder = "frames/"
isdir(output_folder) || mkdir(output_folder)

for i=1:4:nt

    ui = u[interno,1,i]
    vi = u[interno,2,i]

    utotal = sqrt.(ui.^2 + vi.^2)

    fig, ax = Figure(resolution = (800, 800)), Axis(fig[1, 1])
    BEM.quiver!(ax, x_array, y_array, scale .* ui, scale .* vi, 
        color = utotal, colormap = :viridis)

    save(joinpath(output_folder, "frame_$i.png"), fig)
end

movie(output_folder, "output_movie.mp4", framerate = 30)
