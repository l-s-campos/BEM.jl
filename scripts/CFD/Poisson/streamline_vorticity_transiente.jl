using DrWatson,Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))

#Dados de entrada

begin
    L=1
    #=
    h=L/2
    xc=L*0.3
    r = 0.1*L
    =#
end

nelem = 25  #Numero de elementos
order = 2


NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares


## Formatação dos dados ________________________________________________

dad = format_dad(cavity_ω(nelem, order), NPX, NPY) #Cavidade
#dad = format_dad(von_karman_ω(nelem, order,L,h,r,xc), NPX, NPY) #Von_Karman

# Discretization

begin     
    nx = length(unique(dad.pontos_internos[:,1]))
    ny = length(unique(dad.pontos_internos[:,2]))

    x_order = sort(unique(dad.pontos_internos[:,1]))
    y_order = sort(unique(dad.pontos_internos[:,2]))

    dx = (maximum(x_order)-minimum(x_order))/nx
    dy = (maximum(y_order)-minimum(y_order))/ny
    
end

#Tempo

#Para escolher dt, podemos considerar um tempo de convecção e de advecção
#dt advecção = dx/u e dt difusão = dx^2/ν, aqui, escolhi a advecção

begin
    u_wall = 500
    dt = 0.01 

    t_i = 0
    t_f = 10
    t = range(t_i, stop=t_f, step=dt)
    nt = length(t)

    Cₒ = dt*u_wall/dx # Número de Courant
    println("Número de Courant = $Cₒ")
end

# Iniciar variáveis

begin
    ω = zeros(ni(dad)+nc(dad),nt)
    ψ = zeros(ni(dad)+nc(dad))

    ω_nolinear = zeros(ni(dad)+nc(dad))
    ω_flux = zeros(nc(dad))

    #ψ_flux = [zeros(2*order*nelem);ones(order*nelem);zeros(order*nelem)]
    ψ_flux = zeros(nc(dad))

    erro = 1
    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)

    relaxation = 0.2
end


for i = 4:nt

    if i == 4
        global ψ,ω,erro,iter,ω_nolinear,ω_flux,ψ_before
    end

    # Calculando ψ e seu termo não linear 

    #_, ψ[1+nc(dad):ni(dad)+nc(dad)], _ =  poisson_ψ(-ω[:,i-1],ψ_flux)

    ψ[interno] = ψ_dif_finitas(ψ[interno], ω[interno,i-1], nx, ny, dx, dy) # Problema ao aplicar em geometria qualquer

    # Calculando ω e seu termo não linear (envolve derivadas de psi)

    
    u_int = Dy(nx, ny, dy)*ψ[interno]
    v_int = -Dₓ(nx, ny, dx)*ψ[interno]
    
    #u_int,v_int = deriva_interno_generico(dad,ψ[interno],dx,dy)

    v_cont = zeros(1:nc(dad))
    u_cont = ψ_flux

    
    dωdx_int = Dₓ(nx, ny, dx)*ω[interno,i-1]
    dωdy_int = Dy(nx, ny, dy)*ω[interno,i-1]
    dωdx_cont,dωdy_cont = deriva_contorno(dad,ω[contorno],ω_flux)
    

    #dωdx_int,dωdy_int = deriva_interno_generico(dad,ω[interno,i-1],dx,dy)
    

    ω_nolinear = [u_cont.*dωdx_cont + v_cont.*dωdy_cont; u_int.*dωdx_int + v_int.*dωdy_int]

    ω_contorno, ω_interno, ω_fluxo =  poisson_ω_transiente(dad,ω_nolinear,ψ[interno],ω,i,u_wall)
    ω[contorno,i] = (1-relaxation)*ω[contorno,i] + relaxation*ω_contorno
    ω[interno,i] = (1-relaxation)*ω[interno,i] + relaxation*ω_interno
    ω_flux = (1-relaxation)*ω_flux + relaxation*ω_fluxo

    #ω[1:nc(dad),i], ω[interno,i], ω_flux =  poisson_ω_transiente_von_karman(dad,ω_nolinear,ψ[interno],ω,i)
    
    if i>4
        erro = nrmse(ω[:,i-1],ω[:,i])
        println("Iteração: $(i-3); Erro: $erro")
    end

end

#Mapa de cores

x_array = dad.pontos_internos[:,1]
y_array = dad.pontos_internos[:,2]

ui = Dy(nx, ny, dy)*ψ[interno]
vi = -Dₓ(nx, ny, dx)*ψ[interno]

utotal = sqrt.(ui.^2 + vi.^2)

BEM.heatmap(
    x_array,
    y_array,
    utotal,
    colormap = :viridis 
)

scale = 10^(-3)
BEM.quiver(x_array,y_array,scale*ui,scale*vi,color=utotal)


#Validation cavity


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

Plots.plot!(cavidade_uy_ref[2:16,2],cavidade_uy_ref[2:16,1],title="Re = 100"
,label="Analítico", xaxis="Velocidade [m/s]",yaxis="Posição vertical [m]",color=:blue)

Plots.plot(dad.pontos_internos[Int((nx)/2):nx:end,1],-vi[Int((nx)/2):Int((nx)/2)+nx])

Plots.scatter(ui[Int((nx)/2):nx:end]/150,y_order,
marker=(3,:x,:red), label="Numérico")