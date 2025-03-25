using DrWatson,Plots,CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))

nelem = 20  #Numero de elementos
order = 2

NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares


#Tempo

begin
    dt = 0.001
    t_i = 0
    t_f = 1

    t = range(t_i, stop=t_f, step=dt)
    nt = length(t)
end

#_______________________________

# ======== INÍCIO ===============

Re = 10
L=5
h=1
caso = "Canal"

begin

    println("Caso $caso, Reynolds = $Re")

    if caso == "Von Karman"

        r = 0.5 #Não pode ser muito pequeno se não dá errado
        L = 20*r
        h = 10*r
        xc = 5*r

        dad_u = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"u"), NPX, NPY)
        dad_v = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"v"), NPX, NPY)
        dad_p = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,1,"p"), NPX, NPY)

    elseif caso == "Cavidade"

        dad_u = format_dad(cavity_poisson(nelem, order,Re,"u"), NPX, NPY)
        dad_v = format_dad(cavity_poisson(nelem, order,Re,"v"), NPX, NPY)
        dad_p = format_dad(cavity_poisson(nelem, order,1,"p"), NPX, NPY)

    elseif caso == "Canal"
        dad_u = format_dad(channel_poisson(nelem, order,L,h,Re,"u"), NPX, NPY)
        dad_v = format_dad(channel_poisson(nelem, order,L,h,Re,"v"), NPX, NPY)
        dad_p = format_dad(channel_poisson(nelem, order,L,h,1,"p"), NPX, NPY) 
    end

    dad = dad_u 


    # Refinamento local

    #=
    retangulo = [0.3, 1.2, 0.7, 1.8]
    intervalo_objeto = 1:2*order*nelem

    dad_u.pontos_internos = refina_local(dad_u,0.01,0.01,retangulo,true,intervalo_objeto)
    dad_v.pontos_internos = refina_local(dad_v,0.01,0.01,retangulo,true,intervalo_objeto)
    dad_p.pontos_internos = refina_local(dad_p,0.01,0.01,retangulo,true,intervalo_objeto)
    =#

    #________________________

    Ht = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),3)
    Gt = zeros(ni(dad_u)+nc(dad_u),nc(dad_u),3)
    M = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),3)

    A_Houbolt = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Houbolt = zeros(ni(dad_u)+nc(dad_u),2)

    A_Adams = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Adams = zeros(ni(dad_u)+nc(dad_u),2)

    A_Euler = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Euler = zeros(ni(dad_u)+nc(dad_u),2) 

    A = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u))
    b = zeros(ni(dad_u)+nc(dad_u))

    # Velocidade u

    Ht[:,:,1], Gt[:,:,1] = BEM.calc_HeGt(dad_u)
    M[:,:,1] = BEM.Monta_M_RIMd(dad_u, npg)
    A_Houbolt[:,:,1], b_Houbolt[:,1] = BEM.aplicaCDC(Ht[:,:,1] - 11*M[:,:,1]/(6*dt), Gt[:,:,1], dad_u)
    A_Adams[:,:,1], b_Adams[:,1] = BEM.aplicaCDC(Ht[:,:,1] - (24/(55*dt))*M[:,:,1], Gt[:,:,1], dad_u)
    A_Euler[:,:,1], b_Euler[:,1] = BEM.aplicaCDC(Ht[:,:,1] - M[:,:,1]/dt, Gt[:,:,1], dad_u)

    # Velocidade v

    Ht[:,:,2], Gt[:,:,2] = BEM.calc_HeGt(dad_v)
    M[:,:,2] = BEM.Monta_M_RIMd(dad_v, npg)
    A_Houbolt[:,:,2], b_Houbolt[:,2] = BEM.aplicaCDC(Ht[:,:,2] - 11*M[:,:,2]/(6*dt), Gt[:,:,2], dad_v)
    A_Adams[:,:,2], b_Adams[:,2] = BEM.aplicaCDC(Ht[:,:,2] - (24/(55*dt))*M[:,:,2], Gt[:,:,2], dad_v)
    A_Euler[:,:,2], b_Euler[:,2] = BEM.aplicaCDC(Ht[:,:,2] - M[:,:,2]/dt, Gt[:,:,2], dad_v)

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

    dx = (maximum(dad_u.pontos_internos[:,1])-minimum(dad_u.pontos_internos[:,1]))/nx
    dy = (maximum(dad_u.pontos_internos[:,2])-minimum(dad_u.pontos_internos[:,2]))/ny

    x_order = sort(unique(dad_u.pontos_internos[:,1]))
    y_order = sort(unique(dad_u.pontos_internos[:,2]))

    Cₒ =  dt*1/dx
    println("Número de Courant:$Cₒ")
end

#_______________

# Iniciando variáveis

begin
    nolinear = zeros(2*(ni(dad_u)+nc(dad_u)))
    u = zeros(ni(dad_u)+nc(dad_u),2,nt)
    u_dot = zeros(nc(dad_u) + ni(dad_u),2,nt)
    p = zeros(ni(dad_u)+nc(dad_u))
    t = zeros(nc(dad_u),2)

    interno = nc(dad_u)+1:nc(dad_u)+ni(dad_u)
    contorno = 1:nc(dad_u)
    generico = false
    temporal_method = "Euler"

    relaxation = 0.001

    dist_cont_cont = zeros(nc(dad))

    for i=1:nc(dad)
        if i == nc(dad)
            dist_cont_cont[i] = sqrt((dad.NOS[i,1]-dad.NOS[i-1,1])^2 + (dad.NOS[i,2]-dad.NOS[i-1,2])^2)
        else
            dist_cont_cont[i] = sqrt((dad.NOS[i+1,1]-dad.NOS[i,1])^2 + (dad.NOS[i+1,2]-dad.NOS[i,2])^2)
        end
    end

    dNx,dNy = BEM.montaFs([dad_u.NOS;dad_u.pontos_internos])[2:3]
    Deriva_x = Dₓ(nx, ny, dx) 
    Deriva_y = Dy(nx, ny, dy)
end

#Se tiver genérico
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

#_______________

#Loop solução

if generico == false

    for i=4:nt

        if i==4
            global nolinear,p,u,t,u_dot
        end
        

        # Derivadas do campo de velocidade

        #=
        dudx = Deriva_x*u[interno,1]
        dvdx = Deriva_x*u[interno,2]
        dudy = Deriva_y*u[interno,1]
        dvdy = Deriva_y*u[interno,2]

        dudx_cont,dudy_cont = deriva_contorno(dad_u,u[contorno,1,i-1],t[:,1],dist_cont_cont)
        dvdx_cont,dvdy_cont = deriva_contorno(dad_v,u[contorno,2,i-1],t[:,2],dist_cont_cont)
        =#

        #=
        dudx,dudy=deriva_interno_generico(dad,u[interno,1],dx,dy,index_backward,index_forward,index_dist_int_int)
        dvdx,dvdy=deriva_interno_generico(dad,u[interno,2],dx,dy,index_backward,index_forward,index_dist_int_int)

        dudx_cont,dudy_cont = deriva_contorno_generico(dad,u[contorno,1],t[:,1],dist_cont_cont)
        dvdx_cont,dvdy_cont = deriva_contorno_generico(dad,u[contorno,2],t[:,2],dist_cont_cont)
        =#

        #=
        nolinear[2*nc(dad_u)+1:2:2*(nc(dad_u)+ni(dad_u))] = u[interno,1,i-1].*dudx .+ u[interno,2,i-1].*dudy
        nolinear[2*nc(dad_u)+2:2:2*(nc(dad_u)+ni(dad_u))] = u[interno,1,i-1].*dvdx .+ u[interno,2,i-1].*dvdy
    
        nolinear[1:2:2*nc(dad_u)] = u[contorno,1,i-1].*dudx_cont .+ u[contorno,2,i-1].*dudy_cont
        nolinear[2:2:2*nc(dad_u)] = u[contorno,1,i-1].*dvdx_cont .+ u[contorno,2,i-1].*dvdy_cont
        =#


        dudx = dNx*u[:,1,i-1]
        dudy = dNy*u[:,1,i-1]

        dvdx = dNx*u[:,2,i-1]
        dvdy = dNy*u[:,2,i-1]

        nolinear[1:2:end] = u[:,1,i-1].*dudx[:] + u[:,2,i-1].*dudy[:]
        nolinear[2:2:end] = u[:,1,i-1].*dvdx[:] + u[:,2,i-1].*dvdy[:]
    
        #_________________
    
        # Obtendo a velocidade
    
        #Houbolt

        #=
        u[contorno,1,i], u[interno,1,i], t[:,1], u_dot[:,1,:] = 
        poisson_Houbolt(dad_u,(nolinear[1:2:end]),i,u[:,1,:],u_dot[:,1,:],A_Houbolt[:,:,1],b_Houbolt[:,1],M[:,:,1])

        u[contorno,2,i], u[interno,2,i], t[:,2], u_dot[:,2,:] = 
        poisson_Houbolt(dad_v,(nolinear[2:2:end]),i,u[:,2,:],u_dot[:,2,:],A_Houbolt[:,:,2],b_Houbolt[:,2],M[:,:,2])     
        =#

        #Adams

        #=
        u[contorno,1,i], u[interno,1,i], t[:,1], u_dot[:,1,:] = 
        poisson_Adams(dad_u,(nolinear[1:2:end]),i,u[:,1,:],u_dot[:,1,:],A_Houbolt[:,:,1],A_Adams[:,:,1],b_Houbolt[:,1],b_Adams[:,1],M[:,:,1])

        u[contorno,2,i], u[interno,2,i], t[:,2], u_dot[:,2,:] = 
        poisson_Adams(dad_v,(nolinear[2:2:end]),i,u[:,2,:],u_dot[:,2,:],A_Houbolt[:,:,2],A_Adams[:,:,2],b_Houbolt[:,2],b_Adams[:,2],M[:,:,2]) 
        =#

        # Euler
        
        u[contorno,1,i], u[interno,1,i], t[:,1] = 
        poisson_Euler(dad_u,nolinear[1:2:end],i,u[:,1,:],A_Euler[:,:,1],b_Euler[:,1],M[:,:,1])

        u[contorno,2,i], u[interno,2,i], t[:,2] = 
        poisson_Euler(dad_v,nolinear[2:2:end],i,u[:,2,:],A_Euler[:,:,2],b_Euler[:,2],M[:,:,2])
        
        #=
        u_cont_x,u_int_x,t_x= 
        poisson_Euler(dad_u,(nolinear[1:2:end]),i,u[:,1,:],A_Euler[:,:,1],b_Euler[:,1],M[:,:,1])

        u_cont_y,u_int_y,t_y= 
        poisson_Euler(dad_v,(nolinear[2:2:end]),i,u[:,2,:],A_Euler[:,:,2],b_Euler[:,2],M[:,:,2])

        u[contorno,1,i] = u_cont_x*relaxation + u[contorno,1,i-1]*(1-relaxation) 
        u[interno,1,i] = u_int_x*relaxation +  u[interno,2,i-1]*(1-relaxation)
        t[:,1] = t_x*relaxation + t[:,1]*(1-relaxation) 

        u[contorno,2,i] = u_cont_y*relaxation + u[contorno,2,i-1]*(1-relaxation) 
        u[interno,2,i] = u_int_y*relaxation +  u[interno,2,i-1]*(1-relaxation)
        t[:,2] = t_y*relaxation + t[:,2]*(1-relaxation)
        =#

        # Obtendo a pressão e sua derivada
    
        #Derivadas atualizadas __________
        
        #=
        dudx = Deriva_x*u[interno,1,i]
        dvdy = Deriva_y*u[interno,2,i]
    
        dudx_cont,_ = deriva_contorno(dad_u,u[contorno,1,i],t[:,1])
        _,dvdy_cont = deriva_contorno(dad_v,u[contorno,2,i],t[:,2])
        =#

        dudx = dNx*u[:,1,i-1]
        dvdy = dNy*u[:,2,i-1]
        
        # Achando a pressão
    
        #p[contorno], p[interno], p_flux = poisson_steady(dad_p,1/dt*([dudx_cont;dudx] + [dvdy_cont;dvdy]),A,b,M[:,:,3])
        
        p[contorno], p[interno], p_flux = poisson_steady(dad_p,1/dt*(dudx + dvdy),A,b,M[:,:,3])

        #=
        dpdx_int = Deriva_x*p[interno]
        dpdy_int = Deriva_y*p[interno]

        dpdx_cont,dpdy_cont = deriva_contorno(dad_u,p[contorno],p_flux)
    
    
        dpdx = [dpdx_cont;dpdx_int]
        dpdy = [dpdy_cont;dpdy_int]
        =#

        dpdx = dNx*p
        dpdy = dNy*p
        
        #Correção da velocidade
        

        u[interno,1,i] = u[interno,1,i] - dt*dpdx[interno]
        u[interno,2,i] = u[interno,2,i] - dt*dpdy[interno]

        erro = nrmse(u[:,:,i-1],u[:,:,i])
    
        u_mean = mean(u[:,1,i])
    
        println("Iteração: $(i-3); Erro: $erro; u_mean = $u_mean")
    
    end

else

    for i=4:nt

        if i==4
            global nolinear,p,u,t,u_dot,
            A,b,A_Houbolt,A_Adams,b_Houbolt,b_Adams,M,A_Euleer,b_Euler
        end
    
        # Derivadas do campo de velocidade
    
        dudx,dudy = deriva_interno_generico(dad_u,u[interno,1,i-1],dx,dy)
        dvdx,dvdy = deriva_interno_generico(dad_u,u[interno,2,i-1],dx,dy)

        dudx_cont,dudy_cont = deriva_contorno_generico(dad_u,u[contorno,1,i-1],t[:,1])
        dvdx_cont,dvdy_cont = deriva_contorno_generico(dad_u,u[contorno,2,i-1],t[:,2])
    
        nolinear[2*nc(dad_u)+1:2:2*(nc(dad_u)+ni(dad_u))] = u[interno,1,i-1].*dudx .+ u[interno,2,i-1].*dudy
        nolinear[2*nc(dad_u)+2:2:2*(nc(dad_u)+ni(dad_u))] = u[interno,1,i-1].*dvdx .+ u[interno,2,i-1].*dvdy
    
        nolinear[1:2:2*nc(dad_u)] = u[contorno,1,i-1].*dudx_cont .+ u[contorno,2,i-1].*dudy_cont
        nolinear[2:2:2*nc(dad_u)] = u[contorno,1,i-1].*dvdx_cont .+ u[contorno,2,i-1].*dvdy_cont
        
    
        #_________________
    
        # Obtendo a velocidade
    
        if temporal_method == "Houbolt"
        
            u[contorno,1,i], u[interno,1,i], t[:,1], u_dot[:,1,:] = 
            poisson_Houbolt(dad_u,(nolinear[1:2:end]),i,u[:,1,:],u_dot[:,1,:],A_Houbolt[:,:,1],b_Houbolt[:,1],M[:,:,1])
    
            u[contorno,2,i], u[interno,2,i], t[:,2], u_dot[:,2,:] = 
            poisson_Houbolt(dad_v,(nolinear[2:2:end]),i,u[:,2,:],u_dot[:,2,:],A_Houbolt[:,:,2],b_Houbolt[:,2],M[:,:,2]) 
    
        elseif temporal_method == "Adams"
        
            u[contorno,1,i], u[interno,1,i], t[:,1], u_dot[:,1,:] = 
            poisson_Adams(dad_u,(nolinear[1:2:end]),i,u[:,1,:],u_dot[:,1,:],A_Houbolt[:,:,1],A_Adams[:,:,1],b_Houbolt[:,1],b_Adams[:,1],M[:,:,1])
    
            u[contorno,2,i], u[interno,2,i], t[:,2], u_dot[:,2,:] = 
            poisson_Adams(dad_v,(nolinear[2:2:end]),i,u[:,2,:],u_dot[:,2,:],A_Houbolt[:,:,2],A_Adams[:,:,2],b_Houbolt[:,2],b_Adams[:,2],M[:,:,2]) 
        else
    
            u[contorno,1,i], u[interno,1,i], t[:,1] = 
            poisson_Euler(dad_u,(nolinear[1:2:end]),i,u[:,1,:],A_Euler[:,:,1],b_Euler[:,1],M[:,:,1])
        
            u[contorno,2,i], u[interno,2,i], t[:,2] = 
            poisson_Euler(dad_v,(nolinear[2:2:end]),i,u[:,2,:],A_Euler[:,:,2],b_Euler[:,2],M[:,:,2]) 
    
        end
    
        # Obtendo a pressão e sua derivada
    
        #Derivadas atualizadas __________

        dudx,dudy = deriva_interno_generico(dad_u,u[interno,1,i],dx,dy)
        dvdx,dvdy = deriva_interno_generico(dad_u,u[interno,2,i],dx,dy)
    

        dudx_cont,dudy_cont = deriva_contorno_generico(dad_u,u[contorno,1,i],t[:,1])
        dvdx_cont,dvdy_cont = deriva_contorno_generico(dad_u,u[contorno,2,i],t[:,2])
        
        # Achando a pressão
    
        p[contorno], p[interno], p_flux = poisson_steady(dad_p,1/dt*([dudx_cont;dudx] + [dvdy_cont;dvdy]),A,b,M[:,:,3])

        dpdx_int,dpdy_int = deriva_interno_generico(dad_u,p[interno],dx,dy)
        dpdx_cont,dpdy_cont = deriva_contorno_generico(dad_u,p[contorno],p_flux)
    
        dpdx = [dpdx_cont;dpdx_int]
        dpdy = [dpdy_cont;dpdy_int]
        
        
        #Correção da velocidade
    
        u[interno,1,i] = u[interno,1,i] - dt*dpdx[interno]
        u[interno,2,i] = u[interno,2,i] - dt*dpdy[interno]
    
        erro = nrmse(u[:,:,i-1],u[:,:,i])
    
        u_mean = mean(u[:,1,i])
    
        println("Iteração: $(i-3); Erro: $erro; u_mean = $u_mean")
    
    end

end


#_______________

# Mapa de cores

x_array = dad_u.pontos_internos[:,1]
y_array = dad_u.pontos_internos[:,2]

uc = u[contorno,1,end]

ui = u[interno,1,end]
vi = u[interno,2,end]

utotal = sqrt.(ui.^2 + vi.^2)

BEM.heatmap(
    x_array,
    y_array,
    utotal,
    colormap = :viridis 
)

scale = 10^-0
BEM.quiver(x_array,y_array,scale*ui,scale*vi,color=utotal)


Plots.plot(uc)

Plots.plot(ui[Int((nx)/2):nx:end],dad_u.pontos_internos[Int((nx)/2):nx:end,2])
Plots.plot(vi[Int((nx)/2):Int((nx)/2)+nx])

# Erros


cavidade_uy_ref = [1.0000 1.00000
0.9766 0.84123
0.9688 0.78871
0.9609 0.73722
0.9531 0.68717
0.8516 0.23151
0.7344 0.00332
0.6172 -0.13641
0.5000 -0.20581
0.4531 -0.21090
0.2813 -0.15662
0.1719 -0.10150
0.1016 -0.06434
0.0703 -0.04775
0.0625 -0.04192
0.0547 -0.03717
0.0000 0.00000]

u_ana_array = cavidade_uy_ref[2:16,2]

u_num = ui[Int((nx+1)/2):nx:end,1]

erro_average = nrmse(u_ana_array,u_num)
erro_L2 = nme(u_ana_array,u_num)
erro_max = maximum(abs.((u_num .- u_ana_array)/(maximum(u_ana_array)-minimum(u_ana_array))))



