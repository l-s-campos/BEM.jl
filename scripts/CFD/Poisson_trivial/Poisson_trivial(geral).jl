## Início da análise
using DrWatson, Plots,SparseArrays
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))


#  Malha
nelem = 20 
order = 2

NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

# Couette

function couette_poisson(nelem,order,NPX,NPY)

    L = 1
    u₀ = 1

    dad = format_dad(Couette_poisson(nelem, order,L,u₀), NPX, NPY)

    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)

    Ht, Gt = BEM.calc_HeGt(dad)
    A, b = BEM.aplicaCDC(Ht, Gt, dad) 

    M = zeros(ni(dad)+nc(dad),ni(dad)+nc(dad))
    f = zeros(ni(dad)+nc(dad))

    u=zeros(ni(dad)+nc(dad))
    t = zeros(nc(dad))

    u[contorno], u[interno], t = poisson_steady(dad,f,A,b,M)

    u_ana_array = dad.pontos_internos[:,2]
    u_num = u[interno,1]

    erro_average = nrmse(u_ana_array,u_num)
    erro_L2 = nme(u_ana_array,u_num)
    erro_max = maximum(abs.((u_num .- u_ana_array)/(maximum(u_ana_array)-minimum(u_ana_array))))

    return erro_average,erro_L2,erro_max

end


# Poiseullie

function poiseullie_poisson(nelem,order,NPX,NPY)

    L = 1
    dpdx = -10^3
    μ = 1

    u_ana(y) = dpdx*L^2/(2*μ)*(((y-L/2)/L)^2-1/4) 

    dad = format_dad(Poiseullie_poisson(nelem, order,L), NPX, NPY)

    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)


    Ht, Gt = BEM.calc_HeGt(dad)
    M = BEM.Monta_M_RIMd(dad, npg)
    A, b = BEM.aplicaCDC(Ht, Gt, dad)

    f = ones(nc(dad)+ni(dad))*dpdx
    
    u=zeros(ni(dad)+nc(dad))
    t = zeros(nc(dad))

    u[contorno], u[interno], t = poisson_steady(dad,f,A,b,M)
    
    uc = u[contorno,:]

    intervalo = order*nelem+1:2*order*nelem
    y_plot= dad.NOS[intervalo,2]

    dpdx = -10^3
    u_ana_array = dpdx*L^2/(2*1)*((y_plot/L).^2 - (y_plot/L))
    u_num = uc[intervalo,1]

    τ_ana_array = 1/L .+ L*dpdx*(y_plot/L .- 1/2)
    τ_num = Dₓ_1D(length(u_num),dad.NOS[order*nelem+2,2] - dad.NOS[order*nelem+1,2])*u_num
    
    #=
    erro_average = nrmse(u_ana_array,u_num)
    erro_L2 = nme(u_ana_array,u_num)
    erro_max = maximum(abs.((u_num .- u_ana_array)/(maximum(u_ana_array)-minimum(u_ana_array))))
    =#
    
    erro_average = nrmse(τ_ana_array,τ_num)
    erro_L2 = nme(τ_ana_array,τ_num)
    erro_max = maximum(abs.((τ_num .- τ_ana_array)/(maximum(τ_ana_array)-minimum(τ_ana_array))))
    
    
    println(erro_average)

    return erro_average,erro_L2,erro_max

end

erro_average = zeros(8)
erro_L2 = zeros(8)
erro_max = zeros(8)

nelem=5:5:40

p = Plots.plot(layout = (1, 2), xaxis=:log10, yaxis=:log10, xlabel="N° de elem.", ylabel="Erro")

for order = 2:3
    for i = 1:8
        erro_average[i], erro_L2[i], erro_max[i] = couette_poisson(nelem[i], order,NPX,NPY)
    end

    if order == 2
        Plots.plot!(p[1], nelem, erro_average, label="Erro médio", marker=3)
        Plots.plot!(p[1], nelem, erro_L2, label="Erro L2", marker=3)
        Plots.plot!(p[1], nelem, erro_max, label="Erro máximo", marker=3)
        title!(p[1], "Ordem $order")  
    elseif order == 3
        Plots.plot!(p[2], nelem, erro_average, label="Erro médio", marker=3)
        Plots.plot!(p[2], nelem, erro_L2, label="Erro L2", marker=3)
        Plots.plot!(p[2], nelem, erro_max, label="Erro máximo", marker=3)
        title!(p[2], "Ordem $order")
    else
        Plots.plot!(p[3], nelem, erro_average, label="Erro médio", marker=3)
        Plots.plot!(p[3], nelem, erro_L2, label="Erro L2", marker=3)
        Plots.plot!(p[3], nelem, erro_max, label="Erro máximo", marker=3)
        title!(p[3], "Ordem $order")    
    end
end

display(p)

# Poiseullie + Couette 
begin

    L = 1
    dpdx = -10^3
    μ = 1
    u₀ = 80

    u_ana(y) = u₀*y/L + L^2/(2*μ)*dpdx*((y/L)^2 - (y/L))
    τ_ana(y) = μ*u₀/L + L*dpdx*(y/L - 1/2)

    dad = format_dad(Poiseullie_Couette_poisson(nelem, order,L,u₀), NPX, NPY)

    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)


    Ht, Gt = BEM.calc_HeGt(dad)
    M = BEM.Monta_M_RIMd(dad, npg)
    A, b = BEM.aplicaCDC(Ht, Gt, dad)

    f = ones(nc(dad)+ni(dad))*dpdx

    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)

    u=zeros(ni(dad)+nc(dad))
    t = zeros(nc(dad))

    u[contorno], u[interno], t = poisson_steady(dad,f,A,b,M)
end


#Gravidade

begin

    g = 9.81
    h = 0.1
    L = 1
    v = 1
    μ = 1

    u_ana(y) = g/v*h^2*((y/h) - 1/2*(y/h)^2)
    τ_ana(y) = μ*g/v*h^2*(1/h - y/h^2)

    dad = format_dad(gravidade_poisson(nelem, order,h,L), NPX, NPY)

    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)
    
    Ht, Gt = BEM.calc_HeGt(dad)
    M = BEM.Monta_M_RIMd(dad, npg)
    A, b = BEM.aplicaCDC(Ht, Gt, dad)

    f = ones(nc(dad)+ni(dad))*(-g)

    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)

    u=zeros(ni(dad)+nc(dad))
    t = zeros(nc(dad))

    u[contorno], u[interno], t = poisson_steady(dad,f,A,b,M)

end

# Taylor Couette (coordenadas cilíndricas)

begin 
    ra=1
    rb=10
    ω_a=0
    ω_b=10

    u_ana(r) = (ω_b * rb^2 - ω_a * ra^2) / (rb^2 - ra^2)*r + ((ω_a - ω_b) * ra^2 * rb^2 / (rb^2 - ra^2))/r
    
    dad = format_dad(Taylor_Couette(nelem, order,ra,rb,ω_a,ω_b), NPX, NPY)

    interno = nc(dad)+1:ni(dad)+nc(dad)
    contorno = 1:nc(dad)


    Ht, Gt = BEM.calc_HeGt(dad)
    A, b = BEM.aplicaCDC(Ht, Gt, dad)

    M = zeros(nc(dad)+ni(dad),nc(dad)+ni(dad))
    f = zeros(nc(dad)+ni(dad))
    
    u=zeros(ni(dad)+nc(dad))
    t = zeros(nc(dad))

    u[contorno], u[interno], t = poisson_steady(dad,f,A,b,M)

    radius = sqrt.(dad.NOS[:,1].^2 + dad.NOS[:,2].^2)
    theta = atan.(dad.NOS[:,2]./dad.NOS[:,1])
end

# Pós processamento

intervalo = 1:order*nelem

uc = u[contorno]
ui = u[interno]

x_array = dad.pontos_internos[:,1]
y_array = dad.pontos_internos[:,2]

x_unique = sort(unique(x_array))
y_unique = sort(unique(y_array))

dy_cont = (rb-ra)/(order*nelem) #Verificar se está pegando do intervaloc certo!

r_contorno = radius[intervalo]

nx = length(unique(dad.pontos_internos[:,1]))
ny = length(unique(dad.pontos_internos[:,2]))

# Mapa de calor
fig, ax, hm = BEM.heatmap(x_array, y_array, ui)
BEM.Colorbar(fig[:, end+1], hm)
fig

# Campo vetorial
escala = 10^-2
BEM.quiver(x_array,y_array,escala*ui,zeros(length(ui)),color=ui)

# Plots

u_ana_array = u_ana.(r_contorno)
τ_ana_array = τ_ana.(r_contorno)

u_num = uc[intervalo]
τ_num = Dₓ_1D(length(u_num),dy_cont)*u_num*μ

Plots.plot(r_contorno,u_ana_array,color=:blue,label="Analítico", xaxis="Velocidade [m/s]",yaxis="y [m]")
Plots.scatter!(r_contorno,u_num,marker=(3,:x,:red), label="Numérico")

Plots.plot(r_contorno,τ_ana_array,color=:blue,label="Analítico", xaxis="τ [Pa]",yaxis="y [m]")
Plots.scatter!(r_contorno,τ_num,marker=(3,:x,:red), label="Numérico")

#Erros 

erro_average = nrmse(u_ana_array,u_num)
erro_L2= nme(u_ana_array,u_num)
erro_max = maximum(abs.((u_num .- u_ana_array)/(maximum(u_ana_array)-minimum(u_ana_array))))

erro_average = nrmse(τ_ana_array,τ_num)
erro_L2= nme(τ_ana_array,τ_num)
erro_max = maximum(abs.((τ_num .- τ_ana_array)/(maximum(τ_ana_array)-minimum(τ_ana_array))))






