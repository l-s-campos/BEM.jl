using DrWatson,Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include("my_functions.jl")

nelem = 20  #Numero de elementos
order = 2


NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

#Só para dar início

dad = format_dad(cavity_ψ(nelem, order), NPX, NPY)

# Discretization

nx = length(unique(dad.pontos_internos[:,1]))
ny = length(unique(dad.pontos_internos[:,1]))

dx = (maximum(dad.pontos_internos[:,1])-minimum(dad.pontos_internos[:,1]))/nx
dy = (maximum(dad.pontos_internos[:,2])-minimum(dad.pontos_internos[:,2]))/ny

#Iteração

ω = zeros(ni(dad)+nc(dad))
ψ = zeros(ni(dad)+nc(dad))
ψ_before = deepcopy(ψ)

ω_nolinear = zeros(ni(dad)+nc(dad))
ω_flux = zeros(nc(dad))

# Contorno especificado em fluxo e valor 0
ψ_flux = [zeros(2*order*nelem);ones(order*nelem);zeros(order*nelem)]

erro = 1
iter = 0

while erro > 10^-8 && iter < 20

    if iter == 0
        global ψ,ω,erro,iter,ω_nolinear,ω_flux,ψ_before
    end

    # Calculando ψ e seu termo não linear 


    _, ψ[1+nc(dad):ni(dad)+nc(dad)], _ =  poisson_ψ(-ω,ψ_flux)
    
    # Calculando ω e seu termo não linear (envolve derivadas de psi)

    u_int = Dy(nx, ny, dy)*ψ[1+nc(dad):ni(dad)+nc(dad)]
    v_int = -Dₓ(nx, ny, dx)*ψ[1+nc(dad):ni(dad)+nc(dad)]

    v_cont = zeros(1:nc(dad))
    u_cont = ψ_flux


    dωdx_int = Dₓ(nx, ny, dx)*ω[1+nc(dad):ni(dad)+nc(dad)]
    dωdy_int = Dy(nx, ny, dy)*ω[1+nc(dad):ni(dad)+nc(dad)]
    dωdx_cont,dωdy_cont = deriva_contorno(dad,ω[1:nc(dad)],ω_flux)

    ω_nolinear = [u_cont.*dωdx_cont + v_cont.*dωdy_cont; u_int.*dωdx_int + v_int.*dωdy_int]

    ω[1:nc(dad)], ω[1+nc(dad):ni(dad)+nc(dad)], ω_flux =  poisson_ω(ω_nolinear,ψ[1+nc(dad):ni(dad)+nc(dad)])
    
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

ui = Dy(nx, ny, dy)*ψ[1+nc(dad):ni(dad)+nc(dad)]
vi = -Dₓ(nx, ny, dx)*ψ[1+nc(dad):ni(dad)+nc(dad)]


u_total = sqrt.(ui.^2 + vi.^2)

BEM.heatmap(
    x_array,
    y_array,
    u_total,
    colormap = :viridis 
)

BEM.quiver(x_array,y_array,ui,vi)
