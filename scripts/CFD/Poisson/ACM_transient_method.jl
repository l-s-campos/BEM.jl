using DrWatson,Plots,CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))

nelem = 20  #Numero de elementos
order = 2

NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

#=============== CASO ===============#

caso = "Cavidade"

#Tempo

begin
    dt = 0.05
    t_i = 0
    t_f = 5

    t = range(t_i, stop=t_f, step=dt)
    nt = length(t)
end

# Tempo artificial do loop interno

tau = dt/50

#_______________________________

# ======== INÍCIO ===============#

begin
    
    if caso == "Cavidade"

        Re = 100

        dad_u = format_dad(cavity_poisson(nelem, order,Re,"u"), NPX, NPY)
        dad_v = format_dad(cavity_poisson(nelem, order,Re,"v"), NPX, NPY)

        println("Caso Cavidade, Reynolds = $Re")

    elseif caso == "Von Karman"

        r=0.1 #Não pode ser muito pequeno se não dá errado
        L=20*r
        h = 10*r
        xc = 5*r

        Re = (2*r)*u₀
        println("Caso Von Karman, Reynolds = $Re")

        dad_u = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,u₀,"u"), NPX, NPY)
        dad_v = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,u₀,"v"), NPX, NPY)
        
    end

    # Refinamento local

    #=
    retangulo = [0.3, 1.2, 0.7, 1.8]
    intervalo_objeto = 1:2*order*nelem

    dad_u.pontos_internos = refina_local(dad_u,0.01,0.01,retangulo,true,intervalo_objeto)
    dad_v.pontos_internos = refina_local(dad_v,0.01,0.01,retangulo,true,intervalo_objeto)
    =#

    #________________________

    Ht = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    Gt = zeros(ni(dad_u)+nc(dad_u),nc(dad_u),2)
    M = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)

    #=
    A_Houbolt = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Houbolt = zeros(ni(dad_u)+nc(dad_u),2)

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

end

#_______________________________

# Discretization

begin
    nx = length(unique(dad_u.pontos_internos[:,1]))
    ny = length(unique(dad_u.pontos_internos[:,1]))

    dx = (maximum(dad_u.pontos_internos[:,1])-minimum(dad_u.pontos_internos[:,1]))/nx
    dy = (maximum(dad_u.pontos_internos[:,2])-minimum(dad_u.pontos_internos[:,2]))/ny

    Cₒ =  dt*1/dx
    println("Número de Courant:$Cₒ")
end

#_______________

# Iniciando variáveis

begin

    u = zeros(ni(dad_u)+nc(dad_u),2,nt)
    p = zeros(ni(dad_u)+nc(dad_u),nt)

    #u_dot = zeros(nc(dad_u) + ni(dad_u),2,4)

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

    a = 1 #Compressibilidade aritificial

    
    dist_index = zeros(nc(dad))
    index_dist = zeros(nc(dad))

    for k=1:nc(dad)
        dist = sqrt.((dad.NOS[k,1] .- dad.pontos_internos[:,1]).^2 .+ (dad.NOS[k,2] .- dad.pontos_internos[:,2]).^2)
        dist_index[k] = minimum(dist)
        index_dist[k] = Int(findfirst(x -> x == dist_index[k], dist))
    end
end

#_______________

#Loop solução

for i=4:nt

    iter = 0
    erro=1

    if i==4
        global nolinear,p,p_flux,u,t,dpdx,dpdy,
        A_Houbolt,A_Adams,b_Houbolt,b_Adams,M,A_Euleer,b_Euler,
        u_new,u_prev,p_new,p_prev,unsteady,dist_index,index_dist
    end

    while erro > 10^-5

        iter+=1

        # Derivadas do campo de velocidade
        
        dudx = Dₓ(nx,ny,dx)*u_prev[interno,1]
        dudy = Dy(nx,ny,dx)*u_prev[interno,1]
        dvdx = Dₓ(nx,ny,dy)*u_prev[interno,2]
        dvdy = Dy(nx,ny,dy)*u_prev[interno,2]

        dudx_cont,dudy_cont = deriva_contorno(dad_u,u_prev[contorno,1],t[:,1])
        dvdx_cont,dvdy_cont = deriva_contorno(dad_v,u_prev[contorno,2],t[:,2])

        nolinear[1:2:2*nc(dad_u)] = u_prev[contorno,1].*dudx_cont .+ u_prev[contorno,2].*dudy_cont
        nolinear[2:2:2*nc(dad_u)] = u_prev[contorno,1].*dvdx_cont .+ u_prev[contorno,2].*dvdy_cont

        nolinear[2*nc(dad_u)+1:2:2*(nc(dad_u)+ni(dad_u))] = u_prev[interno,1].*dudx .+ u_prev[interno,2].*dudy
        nolinear[2*nc(dad_u)+2:2:2*(nc(dad_u)+ni(dad_u))] = u_prev[interno,1].*dvdx .+ u_prev[interno,2].*dvdy
        
        unsteady[1:2:end] = (u_prev[:,1] - u[:,1,i-1])/dt
        unsteady[2:2:end] = (u_prev[:,2] - u[:,2,i-1])/dt 

        #_________________

        # Obtendo a velocidade

        u_new[contorno,1], u_new[interno,1], t[:,1] = 
        poisson_Euler_unsteady(dad_u, nolinear[1:2:end] + dpdx + unsteady[1:2:end],u_prev[:,1],A_Euler[:,:,1],b_Euler[:,1],M[:,:,1])
    
        u_new[contorno,2], u_new[interno,2], t[:,2] = 
        poisson_Euler_unsteady(dad_v, nolinear[2:2:end] + dpdy + unsteady[2:2:end],u_prev[:,2],A_Euler[:,:,2],b_Euler[:,2],M[:,:,2]) 

        # Obtendo a pressão e sua derivada

        # Derivadas atualizadas __________
         
        dudx = Dₓ(nx,ny,dx)*u_new[interno,1]
        dvdy = Dy(nx,ny,dy)*u_new[interno,2]
        
        # Achando a pressão

        p_new[interno] = p_prev[interno] - a*tau*(dudx + dvdy)

        # Aplicando CDC

        for k=1:nc(dad)

            if k in 1:2*order*nelem || k in 3*order*nelem+1:4*order*nelem # Paredes com condição de fluxo zero
    
                p[k] = p[nc(dad) + Int(index_dist[k])]
    
            else # Pressão especificada
                p[k] = 0
                p_flux[k] = (p[k] - p[nc(dad) + Int(index_dist[k])])/dist_index[k]
            end
        end   

        dpdx_int = Dₓ(nx,ny,dx)*p_new[interno]
        dpdy_int = Dy(nx,ny,dy)*p_new[interno]
        dpdx_cont,dpdy_cont = deriva_contorno(dad_u,p_new[contorno],p_flux)

        dpdx = [dpdx_cont;dpdx_int]
        dpdy = [dpdy_cont;dpdy_int]

        erro = nrmse(u_prev,u_new)
        
        if iter % 100 == 0
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
pi = p[interno,end]

utotal = sqrt.(ui.^2 + vi.^2)

BEM.heatmap(
    x_array,
    y_array,
    ui,
    colormap = :viridis 
)

scale = 10^-1
BEM.quiver(x_array,y_array,scale*ui,scale*vi,color=utotal)


Plots.plot(vi[Int(nx/2):Int(nx/2)+nx])
Plots.plot(ui[Int(nx/2):nx:end])
Plots.plot(p[contorno,end])
