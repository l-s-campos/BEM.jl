## Início da análise
using DrWatson, Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))


#  Malha
nelem = 20 
order = 2

NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

# ======== INÍCIO ===============

L = 1
h = 1
Re = 10

begin
    
    #Cavidade

    dad_u = format_dad(cavity_poisson(nelem, order,Re,"u"), NPX, NPY)
    dad_v = format_dad(cavity_poisson(nelem, order,Re,"v"), NPX, NPY)
    dad_p = format_dad(cavity_poisson(nelem, order,1,"p"), NPX, NPY)
    dad = dad_u

    #Canal plano
    #=
    dad_u = format_dad(channel_poisson(nelem, order, L, h, u₀,"u"), NPX, NPY)
    dad_v = format_dad(channel_poisson(nelem, order, L, h, u₀,"v"), NPX, NPY)
    dad_p = format_dad(channel_poisson(nelem, order, L, h, u₀,"p"), NPX, NPY)
    =#

    # Placa plana
    #=
    dad_u = format_dad(channel_poisson(nelem, order, L, h, u₀,"u"), NPX, NPY)
    dad_v = format_dad(channel_poisson(nelem, order, L, h, u₀,"v"), NPX, NPY)
    dad_p = format_dad(channel_poisson(nelem, order, L, h, u₀,"p"), NPX, NPY)
    =#

    Ht = zeros(ni(dad)+nc(dad),ni(dad)+nc(dad),3)
    Gt = zeros(ni(dad)+nc(dad),nc(dad),3)
    M = zeros(ni(dad)+nc(dad),ni(dad)+nc(dad),3)

    A = zeros(ni(dad)+nc(dad),ni(dad)+nc(dad),3)
    b = zeros(ni(dad)+nc(dad),3)

    # Velocidade u

    Ht[:,:,1], Gt[:,:,1] = BEM.calc_HeGt(dad_u)
    M[:,:,1] = BEM.Monta_M_RIMd(dad_u, npg)
    A[:,:,1], b[:,1] = BEM.aplicaCDC(Ht[:,:,1], Gt[:,:,1], dad_u) 

    # Velocidade v

    Ht[:,:,2], Gt[:,:,2] = BEM.calc_HeGt(dad_v)
    M[:,:,2] = BEM.Monta_M_RIMd(dad_v, npg)
    A[:,:,2], b[:,2] = BEM.aplicaCDC(Ht[:,:,2], Gt[:,:,2], dad_v)

    # Pressão

    Ht[:,:,3], Gt[:,:,3] = BEM.calc_HeGt(dad_p)
    M[:,:,3] = BEM.Monta_M_RIMd(dad_p, npg)
    A[:,:,3], b[:,3] = BEM.aplicaCDC(Ht[:,:,3], Gt[:,:,3], dad_p)
end

# Discretização

begin
    nx = length(unique(dad.pontos_internos[:,1]))
    ny = length(unique(dad.pontos_internos[:,2]))

    dx = (maximum((dad.pontos_internos[:,1]))-minimum((dad.pontos_internos[:,1])))/nx
    dy = (maximum((dad.pontos_internos[:,2]))-minimum((dad.pontos_internos[:,2])))/ny
end

# Variáveis

begin
    iter=0
    erro = 1

    u = zeros(nc(dad) + ni(dad),2)
    t = zeros(nc(dad),2)
    u_before = deepcopy(u)

    p = zeros(ni(dad) + nc(dad))
    nolinear = zeros(2*(nc(dad)+ni(dad)))
    dpdx = zeros(nc(dad)+ni(dad))
    dpdy = zeros(nc(dad)+ni(dad))

    relaxation = vec([1, 1, 1]) #fatores de relaxação para u,v e p, respectivamente

    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)

    dist_index = zeros(nc(dad))
    index_dist = zeros(nc(dad))

    for k=1:nc(dad)
        dist = sqrt.((dad.NOS[k,1] .- dad.pontos_internos[:,1]).^2 .+ (dad.NOS[k,2] .- dad.pontos_internos[:,2]).^2)
        dist_index[k] = minimum(dist)
        index_dist[k] = findfirst(x -> x == dist_index[k], dist)
    end
end

while erro > 10^-8

    if iter == 0
        global iter,erro,u,t,u_before,u_int,uc,nx,ny,dx,dy,p,nolinear,relaxation,
        dpdx,dpdy,dist_index,index_dist
    end

    ## Derivadas das velocidades
    
    # Derivadas Internas

    dudx = Dₓ(nx,ny,dx)*u[interno,1]
    dudy = Dy(nx,ny,dx)*u[interno,1]
    dvdx = Dₓ(nx,ny,dy)*u[interno,2]
    dvdy = Dy(nx,ny,dy)*u[interno,2]


    nolinear[2*nc(dad_u)+1:2:2*(nc(dad_u)+ni(dad_u))] = u[interno,1].*dudx + u[interno,2].*dudy
    nolinear[2*nc(dad_u)+2:2:2*(nc(dad_u)+ni(dad_u))] = u[interno,1].*dvdx + u[interno,2].*dvdy
    
    # Derivadas no Contorno

    dudx_cont,dudy_cont = deriva_contorno(dad_u,u[contorno,1],t[:,1])
    dvdx_cont,dvdy_cont = deriva_contorno(dad_v,u[contorno,2],t[:,2])


    nolinear[1:2:2*nc(dad_u)] = u[contorno,1].*dudx_cont + u[contorno,2].*dudy_cont
    nolinear[2:2:2*nc(dad_u)] = u[contorno,1].*dvdx_cont + u[contorno,2].*dvdy_cont

    # Termo compressivel (No contorno deve ser zero pelas CDCs)

    #=
    dudx2_int = Dₓ(nx,ny,dx)*dudx
    dvdxdy_int = Dₓ(nx,ny,dx)*dvdy
    dvdy2_int = Dy(nx,ny,dy)*dvdy
    dudxdy_int = Dy(nx,ny,dy)*dudx
    =#

    ## Obtendo o campo de velocidades
    #=
    u[contorno,1], u[interno,1], t[:,1] = poisson_steady(dad_u,nolinear[1:2:end] + dpdx,A[:,:,1],b[:,1],M[:,:,1])
    u[contorno,2], u[interno,2], t[:,2] = poisson_steady(dad_v,nolinear[2:2:end] + dpdy,A[:,:,2],b[:,2],M[:,:,2])
    =#

    
    u_num_c, u_num_int, t[:,1] = poisson_steady(dad_u,nolinear[1:2:end]+dpdx,A[:,:,1],b[:,1],M[:,:,1])
    v_num_c, v_num_int, t[:,2] = poisson_steady(dad_v,nolinear[2:2:end]+dpdy,A[:,:,2],b[:,2],M[:,:,2])

    u[:,1] = [u_num_c;u_num_int]*relaxation[1] + (1-relaxation[1])*u[:,1]
    u[:,2] = [v_num_c;v_num_int]*relaxation[2] + (1-relaxation[2])*u[:,2]
    

    ## Poisson para a pressão

    #Encontro as novas derivadas
    
    dudx = Dₓ(nx,ny,dx)*u[interno,1]
    dudy = Dy(nx,ny,dx)*u[interno,1]
    dvdx = Dₓ(nx,ny,dy)*u[interno,2]
    dvdy = Dy(nx,ny,dy)*u[interno,2]

    dudx_cont,dudy_cont = deriva_contorno(dad_u,u[contorno,1],t[:,1])
    dvdx_cont,dvdy_cont = deriva_contorno(dad_v,u[contorno,2],t[:,2])

    f = [dudx_cont.^2 + 2*dudy_cont .* dvdx_cont + dvdy_cont.^2;
     dudx.^2 + 2*dudy.*dvdx + dvdy.^2]
    
    
    p_num_cont,p_num_int,p_flux = poisson_steady(dad_p,-f,A[:,:,3],b[:,3],M[:,:,3])
    p = [p_num_cont;p_num_int]*relaxation[3] + (1-relaxation[3])*p
    

    #p[contorno], p[interno], p_flux = poisson_steady(dad_p,-f,A[:,:,3],b[:,3],M[:,:,3])
    

    #Resolvendo poisson por dif finitas
    
    #=
    p[interno] = solve_poisson(nx, ny, dx, dy,-f,p[interno])

    #Aplicando CDC

    for k=1:nc(dad)

        if k in 1:2*order*nelem || k in 3*order*nelem+1:4*order*nelem # Paredes com condição de fluxo zero

            p[k] = p[nc(dad) + Int(index_dist[k])]

        else # Pressão especificada
            p[k] = 0
            p_flux[k] = (p[k] - p[nc(dad) + Int(index_dist[k])])/dist_index[k]
        end
    end 
    =#
    


    # Deriva a pressão
    
    dpdx = Dₓ(nx,ny,dx)*p[interno]
    dpdy = Dy(nx,ny,dy)*p[interno]

    dpdx_cont,dpdy_cont = deriva_contorno(dad_p,p[contorno],p_flux)

    dpdx = [dpdx_cont;dpdx]
    dpdy = [dpdy_cont;dpdy]

    # Corrige velocidades

    
    u_num_c, u_num_int, t[:,1] = poisson_steady(dad_u,nolinear[1:2:end]+dpdx,A[:,:,1],b[:,1],M[:,:,1])
    v_num_c, v_num_int, t[:,2] = poisson_steady(dad_v,nolinear[2:2:end]+dpdy,A[:,:,2],b[:,2],M[:,:,2])

    u[:,1] = [u_num_c;u_num_int]*relaxation[1] + (1-relaxation[1])*u[:,1]
    u[:,2] = [v_num_c;v_num_int]*relaxation[2] + (1-relaxation[2])*u[:,2]
    

    #=
    u[contorno,1], u[interno,1], t[:,1] = poisson_steady(dad_u,nolinear[1:2:end] + dpdx,A[:,:,1],b[:,1],M[:,:,1])
    u[contorno,2], u[interno,2], t[:,2] = poisson_steady(dad_v,nolinear[2:2:end] + dpdy,A[:,:,2],b[:,2],M[:,:,2])
    =#

    erro = nrmse(u_before,u)

    println("Erro = $erro")

    u_before = deepcopy(u)
    iter = iter+1
end

uc = u[contorno,:]
ui = u[interno,:]

x_array = dad.pontos_internos[:,1]
y_array = dad.pontos_internos[:,2]
y_unique = unique(y_array)

nx = length(unique(dad.pontos_internos[:,1]))
ny = length(unique(dad.pontos_internos[:,2]))

utotal = sqrt.(ui[:,1].^2 + ui[:,2].^2)

BEM.heatmap(x_array, y_array, utotal)
fig, ax, hm = BEM.heatmap(x_array, y_array, utotal)
BEM.Colorbar(fig[:, end+1], hm,label="Velocidade [m/s]")
fig

escala = 10^-1
BEM.quiver(x_array,y_array,escala*ui[:,1],escala*ui[:,2],color=utotal)