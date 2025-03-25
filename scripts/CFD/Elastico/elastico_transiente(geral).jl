## Início da análise
using DrWatson, Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))


#Dados de entrada

L = 1
Re = 1
λ = 10^5

nelem = 20 #Numero de elementos
nelem_circle=20
order = 2

NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

## Formatação dos dados ________________________________________________

#dad = format_dad(cavity(nelem, order,L, Re), NPX, NPY) # Cavidade
#dad = format_dad(channel_simetria(nelem, order,h,L,Re), NPX, NPY) # Canal com simetria
#dad = format_dad(channel(nelem, order,h,L,Re), NPX, NPY) # Canal
#dad = format_dad(placa_plana(nelem, order,L,Ls,h,Re), NPX, NPY) # Placa plana


r = 0.1
L = 30*r
h = 10*r
xc = 10*r

dad = format_dad(von_karman(nelem, order,L,h,r,xc,Re), NPX, NPY)


#______________________________

H, G = calc_HeG(dad, npg,interno=true)  #importante

M = BEM.Monta_M_RIMd(dad, npg)# calc_HeG_potencial linha 310

#Tempo

begin
    dt = 1e-3
    t_i = 0
    t_f = 1

    t = range(t_i, stop=t_f, step=dt)
    nt = length(t)
end

#A, b = BEM.aplicaCDC(H + M/dt, G, dad) #Euler
A, b = BEM.aplicaCDC(H + 11*M/(6*dt), G, dad) #Houbolt
#A, b = BEM.aplicaCDC(H, G, dad) 

#Discretization

begin
    nx = length(unique(dad.pontos_internos[:,1]))
    ny = length(unique(dad.pontos_internos[:,1]))

    x_order = sort(unique(dad.pontos_internos[:,1]))
    y_order = sort(unique(dad.pontos_internos[:,2]))

    dx = (maximum(dad.pontos_internos[:,1])-minimum(dad.pontos_internos[:,1]))/nx
    dy = (maximum(dad.pontos_internos[:,2])-minimum(dad.pontos_internos[:,2]))/ny

    Cₒ =  dt*1/dx
    println("Número de Courant:$Cₒ")
end

# Variaveis

begin
    n = ni(dad)+nc(dad)

    global u = zeros(nc(dad) + ni(dad),2,nt)
    global t = zeros(nc(dad),2)
    global u_n_corrigido = zeros((nc(dad)+ni(dad))*2)
    global u_dot = zeros((nc(dad) + ni(dad))*2,nt)
    global u_prev = deepcopy(u[:,:,1])

    global nolinear = zeros(2*(nc(dad)+ni(dad)))

    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)
    global iter=0
    global erro = 1

    relaxation = 0.1

    dist_cont_cont = zeros(nc(dad))

    for i=1:nc(dad)
        if i == nc(dad)
            dist_cont_cont[i] = sqrt((dad.NOS[i,1]-dad.NOS[i-1,1])^2 + (dad.NOS[i,2]-dad.NOS[i-1,2])^2)
        else
            dist_cont_cont[i] = sqrt((dad.NOS[i+1,1]-dad.NOS[i,1])^2 + (dad.NOS[i+1,2]-dad.NOS[i,2])^2)
        end
    end

    dNx,dNy = BEM.montaFs([dad.NOS;dad.pontos_internos])[2:3]

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


# Solution

for i=4:nt

    # u_dot 

    # Houbolt

    u_dot[:,i] = ((-18*u[:,:,i-1] + 9*(u[:,:,i-2]) - 2*(u[:,:,i-3]))/(6*dt))'[:]

    # Euler simples
    
    #u_dot[:,i] = (-u[:,:,i-1]/dt)'[:]
    

    while erro > 10^-8

        #Derivadas
        
        #=
        dudx_int = Dₓ(nx, ny, dx)*u[interno,1,i]
        dvdx_int = Dₓ(nx, ny, dx)*u[interno,2,i]

        dudy_int = Dy(nx, ny, dy)*u[interno,1,i]
        dvdy_int = Dy(nx, ny, dy)*u[interno,2,i]

        dudx_cont,dudy_cont = deriva_contorno(dad,u[contorno,1,i],t[:,1],dist_cont_cont)
        dvdx_cont,dvdy_cont = deriva_contorno(dad,u[contorno,2,i],t[:,2],dist_cont_cont)
        
        dudx = [dudx_cont;dudx_int]
        dudy = [dudy_cont;dudy_int]
        dvdx = [dvdx_cont;dvdx_int]
        dvdy = [dvdy_cont;dvdy_int]
        =#

        #=
        dudx,dudy=deriva_interno_generico(dad,u[interno,1,i],dx,dy,index_backward,index_forward,index_dist_int_int)
        dvdx,dvdy=deriva_interno_generico(dad,u[interno,2,i],dx,dy,index_backward,index_forward,index_dist_int_int)

        dudx_cont,dudy_cont = dNx*u[:,1,i],dNy*u[:,1,i]
        dvdx_cont,dvdy_cont = dNx*u[:,2,i],dNy*u[:,2,i]


        nolinear[2*nc(dad)+1:2:2*(nc(dad)+ni(dad))] = u[interno,1,i].*dudx + u[interno,2,i].*dudy
        nolinear[2*nc(dad)+2:2:2*(nc(dad)+ni(dad))] = u[interno,1,i].*dvdx + u[interno,2,i].*dvdy
        
        nolinear[1:2:2*nc(dad)] = u[contorno,1,i].*dudx_cont .+ u[contorno,2,i].*dudy_cont
        nolinear[2:2:2*nc(dad)] = u[contorno,1,i].*dvdx_cont .+ u[contorno,2,i].*dvdy_cont
        =#

        #=
        erro_derivative=1
        while erro_derivative > 10^-5

            dudx = ux[1:2:end]
            dvdx = ux[2:2:end]

            dudy = uy[1:2:end]
            dvdy = uy[2:2:end]
            
            nolinear[1:2:2*n] = u[:,1].*dudx + u[:,2].*dudy
            nolinear[2:2:2*n] = u[:,1].*dvdx + u[:,2].*dvdy

            ux = Hx * u[1:nc(dad),:]'[:] - Gx * t[1:nc(dad),:]'[:] + Mx * nolinear
            uy = Hy * u[1:nc(dad),:]'[:] - Gy * t[1:nc(dad),:]'[:] + My * nolinear

            erro_derivative = nrmse(ux_before,ux)
            println("Erro derivadas: $erro_derivative")

            ux_before = deepcopy(ux)
        end
        =#
        
        nolinear[1:2:2*n] = u[:,1,i].*(dNx*u[:,1,i]) + u[:,2,i].*(dNy*u[:,1,i])
        nolinear[2:2:2*n] = u[:,1,i].*(dNx*u[:,2,i]) + u[:,2,i].*(dNy*u[:,2,i])
        

        #Solution

        x = A \ (b - M*(nolinear + u_dot[:,i]))

        u[contorno,:,i] = u[contorno,:,i]*(1-relaxation) + separa(dad,x)[1]*relaxation
        t = t*(1-relaxation) + separa(dad,x)[2]*relaxation
        
        u[interno,1,i] =  x[2*nc(dad)+1:2:end]*relaxation + (1-relaxation)*u[interno,1,i]
        u[interno,2,i] =  x[2*nc(dad)+2:2:end]*relaxation + (1-relaxation)*u[interno,2,i]
        
        erro = nrmse(u_prev,u[:,:,i])

        u_prev = u[:,:,i]

        println("Iteração $iter; Erro = $erro")
        iter+=1

    end

    iter = 0
    erro = 1

    u_mean = mean(u[:,:,i])

    println("Iteração temporal $(i-3); u médio = $u_mean")
    
end

tempo = 50
seconds = round(tempo*dt,digits=5)
uc = u[contorno,:,tempo]
ui = u[interno,:,tempo]

x_array = dad.pontos_internos[:,1]
y_array = dad.pontos_internos[:,2]

utotal = sqrt.(ui[:,1].^2 + ui[:,2].^2)

BEM.heatmap(x_array, y_array, utotal)
fig, ax, hm = BEM.heatmap(x_array, y_array, utotal)
BEM.Colorbar(fig[:, end+1], hm)
fig

escala = 10^-1
BEM.quiver(x_array,y_array,escala*ui[:,1],escala*ui[:,2],color=utotal)

Plots.plot(ui[Int(nx/2):nx:end,1],y_order,xlabel="u [m/s]",ylabel="Posição vertical [m]", label="t=$seconds s")

Plots.plot!(x_order,ui[Int((nx)/2)+1:Int((nx)/2)+nx,2],
ylabel="v [m/s]",xlabel="Posição horizontal [m]",label="Tempo = $seconds s")