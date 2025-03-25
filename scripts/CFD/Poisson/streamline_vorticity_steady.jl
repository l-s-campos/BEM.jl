using DrWatson,Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))

nelem = 20  #Numero de elementos
order = 2

NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

#Só para dar início

dad = format_dad(cavity_ω(nelem, order), NPX, NPY)
M = BEM.Monta_M_RIMd(dad, npg)

# Discretization

begin
    nx = length(unique(dad.pontos_internos[:,1]))
    ny = length(unique(dad.pontos_internos[:,1]))

    dx = (maximum(dad.pontos_internos[:,1])-minimum(dad.pontos_internos[:,1]))/nx
    dy = (maximum(dad.pontos_internos[:,2])-minimum(dad.pontos_internos[:,2]))/ny
end

#Iniciar variáveis

begin
    n = ni(dad)+nc(dad)
    
    global ω = zeros(ni(dad)+nc(dad))
    global ψ = zeros(ni(dad)+nc(dad))
    global ψ_before = deepcopy(ψ)

    global nolinear = zeros(ni(dad)+nc(dad))
    global ω_flux = zeros(nc(dad))

    # Contorno especificado em fluxo e valor 0
    ψ_flux = [zeros(2*order*nelem);ones(order*nelem);zeros(order*nelem)]

    global erro = 1
    global iter = 0
    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)

    relaxation = [1, 1]

    v_cont = zeros(contorno)
    u_cont = ψ_flux

    dNx,dNy = BEM.montaFs([dad.NOS;dad.pontos_internos])[2:3]
end

#vetores distancia

begin
    j_int = Vector{Int}(undef, nc(dad))
    j_normal = Vector{Int}(undef, nc(dad))

    for i =1:nelem*4
        indices = dad.ELEM[i].indices
        for (index, j) in enumerate(indices)

            dist = (dad.NOS[j,1] .- dad.pontos_internos[:,1]).^2 .+ (dad.NOS[j,2] .- dad.pontos_internos[:,2]).^2
            j_int[j] = Int(findfirst(x -> x == minimum(dist),dist))
        end
    end

    for i=1:nc(dad)
        normal = dad.normal[i,:]
        dist_normal = zeros(ni(dad))
        for k=1:ni(dad)
            dist_normal[k] = dot(dad.NOS[i,:],normal) - dot(dad.pontos_internos[k,:],normal)
        end
        j_normal[i] = Int(findfirst(x -> x == minimum(abs.(dist_normal)),dist_normal))
    end
end

while erro > 10^-10 && iter < 1000


    # Calculando ψ por diferenças finitas (apenas internas pois contorno é especificado)

    ψ[interno] = solve_poisson(nx, ny, dx, dy,-ω[interno],ψ[interno])*relaxation[1] + (1-relaxation[1])*ψ[interno]
    
    # Calculando ω e seu termo não linear (envolve derivadas de psi)

    u_int = Dy(nx, ny, dy)*ψ[interno]
    v_int = -Dₓ(nx, ny, dx)*ψ[interno]

    u_total = [u_cont;u_int]
    v_total = [v_cont;v_int]

    # Derivadas de ω

    #=
    dωdx_int = Dₓ(nx, ny, dx)*ω[interno]
    dωdy_int = Dy(nx, ny, dy)*ω[interno]

    dωdx_cont,dωdy_cont = deriva_contorno(dad,ω[contorno],ω_flux)
    =#

    dωdx = dNx*ω
    dωdy = dNy*ω

    nolinear = u_total.*dωdx + v_total.*dωdy

    #nolinear = [u_cont.*dωdx_cont + v_cont.*dωdy_cont; u_int.*dωdx_int + v_int.*dωdy_int]

    ω, ω_flux =  poisson_ω_steady(dad,nolinear,ψ[interno],1,j_int,j_normal)

    ω = ω*(1-relaxation[2]) + ω*relaxation[2]

    if iter>0
        erro = nrmse(ψ_before,ψ)
        println("Iteração: $iter; Erro: $erro")
    end

    iter = iter+1

    ψ_before = deepcopy(ψ)

end

#Mapa de cores

x_array = dad.pontos_internos[:,1]
y_array = dad.pontos_internos[:,2]

ui = Dy(nx, ny, dy)*ψ[interno]
vi = -Dₓ(nx, ny, dx)*ψ[interno]

utotal = sqrt.(ui.^2 + vi.^2)

fig, ax, hm = BEM.heatmap(x_array, y_array, utotal,colormap=:viridis )
BEM.Colorbar(fig[:, end+1], hm,label="Velocidade [m/s]")
fig

scale = 10^(-2)
BEM.quiver(x_array,y_array,scale*ui,scale*vi,color=utotal,colormap = :viridis)


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

Plots.plot(cavidade_uy_ref[2:16,2],cavidade_uy_ref[2:16,1],title="Re = 100"
,label="Analítico", xaxis="Velocidade [m/s]",yaxis="Posição vertical [m]",color=:blue)


Plots.scatter!(-ui[Int((nx)/2):nx:end,1]/6,dad.pontos_internos[Int((nx)/2):nx:end,2],
marker=(3,:x,:red), label="Numérico")

Plots.plot(dad.pontos_internos[Int((nx)/2):nx:end,1],-ui[Int((nx)/2):Int((nx)/2)+nx,2])
