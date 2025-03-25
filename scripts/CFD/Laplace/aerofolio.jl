using DrWatson, Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))


nelem = 5  #Numero de elementos
order = 2
NPX = 10 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

## Formatação dos dados ________________________________________________

#Usando dominio retangular

dad = format_dad(NACA_2(nelem, order), NPX, NPY) # dados

#Usando dominio circular
#dad = format_dad(NACA_potencial(nelem, order), NPX, NPY) # dados


# Corrige CDC _______________________________ Se for Naca_potencial
#=
begin
    for i = nelem*(n) + (nelem*4) + 1 : nelem*(n) + (nelem*4)*2
        n_indices = length(dad.ELEM[i].indices)
        for j = 1:n_indices
            indice = dad.ELEM[i].indices[j]
            dad.ELEM[i].valorCDC[j] = -dad.NOS[indice,1]
        end
    end 

    for i = nelem*(n) + (nelem*4)*2 + 1 : nelem*(n) + (nelem*4)*3
        n_indices = length(dad.ELEM[i].indices)
        for j = 1:n_indices
            indice = dad.ELEM[i].indices[j]
            dad.ELEM[i].valorCDC[j] = -dad.NOS[indice,1]
        end
    end 
end
=#

# Corrige CDC _______________________________ Se for Naca_2

begin

    for i = nelem*n + nelem*4 + 1 : nelem*n + nelem*4*2
        n_indices = length(dad.ELEM[i].indices)
        for j = 1:n_indices
            indice = dad.ELEM[i].indices[j]
            dad.ELEM[i].valorCDC[j] = dad.NOS[indice,1]
        end
    end 

    for i = nelem*n + 3*nelem*4 +1: nelem*n + 4*nelem*4
        n_indices = length(dad.ELEM[i].indices)
        for j = 1:n_indices
            indice = dad.ELEM[i].indices[j]
            dad.ELEM[i].valorCDC[j] = dad.NOS[indice,1]
        end
    end 

end


dad.pontos_internos = refina_local(dad,0.02,0.02,[-0.2,-0.3,1.2,0.3], true, 1:(n)*order*nelem)

#________________________________

H, G = calc_HeG(dad, npg)  #importante

A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b

x = A \ b

T, q = separa(dad, x) #importante

Ti = calc_Ti(dad, T, q, npg);


# Mapa de calor

# Discretização

begin

    x_array = dad.pontos_internos[:,1]
    y_array = dad.pontos_internos[:,2]

    x_order = sort(unique(x_array))
    y_order = sort(unique(y_array))

    nx = length(unique(x_order))
    ny = length(unique(y_order))

    dx = (maximum((x_order))-minimum((x_order)))/nx
    dy = (maximum((y_order))-minimum((y_order)))/ny

end

ui, vi = deriva_interno_generico(dad,Ti,dx,dy)
#ui, vi = deriva_Ti_potencial(dad,Ti,dx,dy)

utotal = sqrt.(ui.^2 + vi.^2)

fig, ax, hm = BEM.heatmap(x_array, y_array, utotal)
BEM.Colorbar(fig[:, end+1], hm)
fig

escala = 4e-2
BEM.quiver(x_array,y_array,escala*ui,escala*vi,color=utotal)


# Derivada em linha

intervalo = 1:(n)*order*nelem

u_cont = deriva_T(dad,T[intervalo],1)
u_int = deriva_T(dad,Ti[:],2)

x_chord = range(0,1,length=Int((length(u_cont[:]))/2))

low = reverse(u_cont[1:Int((intervalo[end])/2)].^2)
up = - u_cont[Int((intervalo[end])/2)+1:end].^2


Plots.scatter(x_chord,up ,xlabel="Chord line [m]",ylabel="Coeficiente de pressão",
label="superior")
Plots.scatter!(x_chord,low,label="inferior")

pressure_up = sum(up)
pressure_low = sum(low)

