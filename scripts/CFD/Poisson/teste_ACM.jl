using DrWatson,Plots,CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))

nelem = 20  #Numero de elementos
order = 2

NPX = 25 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares



#Tempo

begin
    dt = 1e-3
    t_i = 0
    t_f = 5e-0

    t = range(t_i, stop=t_f, step=dt)
    nt = length(t)
end

#_______________________________

# ======== INÍCIO ===============#

caso = "Cavidade"
Re = 1

begin

    println("Caso $caso, Reynolds = $Re") 
    
    if caso == "Cavidade"

        dad_u = format_dad(cavity_poisson(nelem, order,Re,"u"), NPX, NPY)
        dad_v = format_dad(cavity_poisson(nelem, order,Re,"v"), NPX, NPY)

    elseif caso == "Von Karman"

        r = 0.5 #Não pode ser muito pequeno se não dá errado
        L=20*r
        h = 10*r
        xc = 5*r

        dad_u = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"u"), NPX, NPY)
        dad_v = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"v"), NPX, NPY)
        
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

    A_Houbolt = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Houbolt = zeros(ni(dad_u)+nc(dad_u),2)

    A_Adams = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Adams = zeros(ni(dad_u)+nc(dad_u),2)

    A_Euler = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Euler = zeros(ni(dad_u)+nc(dad_u),2) 

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
    nolinear = zeros(2*(ni(dad_u)+nc(dad_u)))
    u = zeros(ni(dad_u)+nc(dad_u),2,nt)
    u_dot = zeros(nc(dad_u) + ni(dad_u),2,nt)
    dpdx = zeros(nc(dad_u)+ni(dad_u))
    dpdy = zeros(nc(dad_u)+ni(dad_u))

    p = zeros(ni(dad_u)+nc(dad_u),nt)
    p_flux = zeros(nc(dad_u))
    t = zeros(nc(dad_u),2)

    interno = nc(dad_u)+1:nc(dad_u)+ni(dad_u)
    contorno = 1:nc(dad_u)

    a = 100 #Compressibilidade aritificial
    generico = false


    minimum_distance = zeros(nc(dad_u))
    index_dist = zeros(nc(dad_u))

    for k=1:nc(dad_u)
        dist = sqrt.((dad_u.NOS[k,1] .- dad_u.pontos_internos[:,1]).^2 .+ (dad_u.NOS[k,2] .- dad_u.pontos_internos[:,2]).^2)
        index_dist = findfirst(x -> x == minimum(dist), dist)
        minimum_distance[k] = dist[index_dist]
    end

end

#_______________

#Loop solução

for i=4:nt

    erro = 1
    
    if i==4
        global nolinear,p,p_flux,u,t,u_dot,dpdx,dpdy,
        A_Houbolt,A_Adams,b_Houbolt,b_Adams,M,A_Euleer,b_Euler
    end
    

    # Derivadas do campo de velocidade


    # Com termo compressivo, como no vídeo
    #=
    du2dx = Dₓ(nx,ny,dx)*(u[interno,1,i-1].^2)
    du2dy = Dy(nx,ny,dx)*(u[interno,1,i-1].^2)

    duvdx = Dₓ(nx,ny,dx)*(u[interno,1,i-1] .*u[interno,2,i-1])
    duvdy = Dy(nx,ny,dx)*(u[interno,1,i-1] .*u[interno,2,i-1])

    dv2dx = Dₓ(nx,ny,dx)*(u[interno,2,i-1].^2)
    dv2dy = Dy(nx,ny,dx)*(u[interno,2,i-1].^2)

    term_1_x_cont,_ = deriva_contorno_com_interno(dad_u,u[contorno,1,i-1].^2, u[interno,1,i-1].^2, du2dx,du2dy)
    _,term_2_x_cont = deriva_contorno_com_interno(dad_u,u[contorno,1,i-1] .*u[contorno,2,i-1], 
    u[interno,1,i-1] .*u[interno,2,i-1],duvdx, duvdy)

    _,term_1_y_cont = deriva_contorno_com_interno(dad_u,u[contorno,2,i-1].^2, u[interno,2,i-1].^2,dv2dx, dv2dy)
    term_2_y_cont = term_2_x_cont
    
    nolinear[2*nc(dad_u)+1:2:2*(nc(dad_u)+ni(dad_u))] = du2dx + duvdy
    nolinear[2*nc(dad_u)+2:2:2*(nc(dad_u)+ni(dad_u))] = dv2dy + duvdx

    nolinear[1:2:2*nc(dad_u)] = term_1_x_cont + term_2_x_cont
    nolinear[2:2:2*nc(dad_u)] = term_1_y_cont + term_2_y_cont
    =#

    dudx = Dₓ(nx,ny,dx)*u[interno,1,i-1]
    dudy = Dy(nx,ny,dx)*u[interno,1,i-1]
    dvdx = Dₓ(nx,ny,dy)*u[interno,2,i-1]
    dvdy = Dy(nx,ny,dy)*u[interno,2,i-1]

    dudx_cont,dudy_cont = deriva_contorno(dad_u,u[contorno,1,i-1],t[:,1])
    dvdx_cont,dvdy_cont = deriva_contorno(dad_v,u[contorno,2,i-1],t[:,2])

    nolinear[2*nc(dad_u)+1:2:2*(nc(dad_u)+ni(dad_u))] = u[interno,1,i-1].*dudx .+ u[interno,2,i-1].*dudy
    nolinear[2*nc(dad_u)+2:2:2*(nc(dad_u)+ni(dad_u))] = u[interno,1,i-1].*dvdx .+ u[interno,2,i-1].*dvdy

    nolinear[1:2:2*nc(dad_u)] = u[contorno,1,i-1].*dudx_cont .+ u[contorno,2,i-1].*dudy_cont
    nolinear[2:2:2*nc(dad_u)] = u[contorno,1,i-1].*dvdx_cont .+ u[contorno,2,i-1].*dvdy_cont

    #_________________

    # Obtendo a velocidade
    
    u[contorno,1,i], u[interno,1,i], t[:,1] = 
    poisson_Euler(dad_u,nolinear[1:2:end] + dpdx,i,u[:,1,:],A_Euler[:,:,1],b_Euler[:,1],M[:,:,1])

    u[contorno,2,i], u[interno,2,i], t[:,2] = 
    poisson_Euler(dad_v,nolinear[2:2:end] + dpdy,i,u[:,2,:],A_Euler[:,:,2],b_Euler[:,2],M[:,:,2])


    # Obtendo a pressão e sua derivada

    #Derivadas atualizadas __________

    dudx = Dₓ(nx,ny,dx)*u[interno,1,i]
    dvdy = Dy(nx,ny,dy)*u[interno,2,i]

    # Achando a pressão
    
    p[interno,i] = p[interno,i-1] - a*dt*(dudx + dvdy)

    #Aplicando CDC
    
    for k=1:nc(dad_u)
        if k in 1:2*order*nelem || k in 3*order*nelem+1:4*order*nelem # Paredes com condição de fluxo zero

            p[k,i] = p[nc(dad_u) + index_dist,i]
            #p_flux[k] = 0
        else # Pressão especificada
            p[k,i] = 0
            p_flux[k] = (p[k,i] - p[nc(dad_u) + index_dist,i])/minimum_distance[k]
        end
    end  
    
    
    dpdx_int = Dₓ(nx,ny,dx)*p[interno,i]
    dpdy_int = Dy(nx,ny,dy)*p[interno,i]
    dpdx_cont,dpdy_cont = deriva_contorno(dad_u,p[contorno,i],p_flux)

    #dpdx_cont = zeros(nc(dad_u))
    #dpdy_cont = zeros(nc(dad_u))

    dpdx = [dpdx_cont;dpdx_int]
    dpdy = [dpdy_cont;dpdy_int]
    
    erro = nrmse(u[:,:,i-1],u[:,:,i])

    u_mean = mean(u[:,1,i])

    println("Iteração: $(i-3); Erro: $erro; u_mean = $u_mean")
end

#_______________

# Mapa de cores

x_array = dad_u.pontos_internos[:,1]
y_array = dad_u.pontos_internos[:,2]

uc = u[contorno,1,end]

ui = u[interno,1,end]
vi = u[interno,2,end]

p_cont = p[contorno,end]
p_int = p[interno,end]

utotal = sqrt.(ui.^2 + vi.^2)

BEM.heatmap(x_array, y_array, utotal)
fig, ax, hm = BEM.heatmap(x_array, y_array, utotal)
BEM.Colorbar(fig[:, end+1], hm,label="Velocidade [m/s]")
fig

scale = 10^-1
BEM.quiver(x_array,y_array,scale*ui,scale*vi,color=utotal)


Plots.plot(vi[Int(nx/2):Int(nx/2)+nx])
Plots.plot(ui[Int(nx/2):nx:end])
Plots.plot(p[contorno,end])