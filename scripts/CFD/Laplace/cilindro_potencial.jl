using DrWatson, Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))

function indice_forward(ponto,dad,arg)

    x_ordenado = sort(unique(dad.pontos_internos[:,1]))
    y_ordenado = sort(unique(dad.pontos_internos[:,2])) 

    if arg=='y'
        indice_ordenado = zeros(2)

        indice_ordenado[1] = Int(findfirst(x -> x == dad.pontos_internos[ponto,1], x_ordenado))
        indice_ordenado[2] = Int(findfirst(x -> x == dad.pontos_internos[ponto,2], y_ordenado))

        indice_procurado = findfirst(row -> row[1] == x_ordenado[Int(indice_ordenado[1])] && row[2] == y_ordenado[Int(indice_ordenado[2]+1)],
        eachrow(dad.pontos_internos))
    else
        indice_ordenado = zeros(2)

        indice_ordenado[1] = Int(findfirst(x -> x == dad.pontos_internos[ponto,1], x_ordenado))
        indice_ordenado[2] = Int(findfirst(x -> x == dad.pontos_internos[ponto,2], y_ordenado))

        indice_procurado = findfirst(row -> row[1] == x_ordenado[Int(indice_ordenado[1]+1)] && row[2] == y_ordenado[Int(indice_ordenado[2])],
        eachrow(dad.pontos_internos))
    end

    return indice_procurado

end
function indice_back(ponto,dad,arg)

    x_ordenado = sort(unique(dad.pontos_internos[:,1]))
    y_ordenado = sort(unique(dad.pontos_internos[:,2]))

    if arg == 'y'
        indice_ordenado = zeros(2)

        indice_ordenado[1] = Int(findfirst(x -> x == dad.pontos_internos[ponto,1], x_ordenado))
        indice_ordenado[2] = Int(findfirst(x -> x == dad.pontos_internos[ponto,2], y_ordenado))

        indice_procurado = findfirst(row -> row[1] == x_ordenado[Int(indice_ordenado[1])] && row[2] == y_ordenado[Int(indice_ordenado[2]-1)],
        eachrow(dad.pontos_internos))
    else
        indice_ordenado = zeros(2)

        indice_ordenado[1] = Int(findfirst(x -> x == dad.pontos_internos[ponto,1], x_ordenado))
        indice_ordenado[2] = Int(findfirst(x -> x == dad.pontos_internos[ponto,2], y_ordenado))

        indice_procurado = findfirst(row -> row[1] == x_ordenado[Int(indice_ordenado[1]-1)] && row[2] == y_ordenado[Int(indice_ordenado[2])],
        eachrow(dad.pontos_internos))
    end
    
    return indice_procurado

end

function deriva_Ti_potencial(dad,Ti,dx,dy)

    n = length(Ti)

    ui = zeros(n)
    vi = zeros(n)

    for i=1:n

        if dad.pontos_internos[i,2] == minimum(dad.pontos_internos[:,2])
            vi[i] = (Ti[indice_forward(i,dad,'y')] - Ti[i])/dy

        elseif dad.pontos_internos[i,2] == maximum(dad.pontos_internos[:,2])
            vi[i] = (Ti[i] - Ti[indice_back(i,dad,'y')])/dy

        else
            if indice_back(i,dad,'y') !== nothing && indice_forward(i,dad,'y') !== nothing
                vi[i] = (Ti[indice_forward(i,dad,'y')] - Ti[indice_back(i,dad,'y')])/(2*dy)
               
            elseif indice_forward(i,dad,'y') !== nothing
                vi[i] = (Ti[indice_forward(i,dad,'y')] - Ti[i])/dy

            else
                vi[i] = (Ti[i] - Ti[indice_back(i,dad,'y')])/dy
            end
        end



        if dad.pontos_internos[i,1] == minimum(dad.pontos_internos[:,1])
            ui[i] = (Ti[i+1] - Ti[i])/dx

        elseif dad.pontos_internos[i,1] == maximum(dad.pontos_internos[:,1])
            ui[i] = (Ti[i] - Ti[i-1])/dx

        else
            if indice_back(i,dad,'x') !== nothing && indice_forward(i,dad,'x') !== nothing
                ui[i] = (Ti[i+1] - Ti[i-1])/(2*dx)

            elseif indice_forward(i,dad,'x') !== nothing
                ui[i] = (Ti[i+1] - Ti[i])/dx

            else
                ui[i] = (Ti[i] - Ti[i-1])/dx
            end
        end
    end
    return ui,vi
end


nelem = 20  #Numero de elementos
order = 2

NPX = 60 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 20    #apenas números pares

## Formatação dos dados ________________________________________________

dad = format_dad(circular_cylinder(nelem, order), NPX, NPY) # dados

# Corrige CDC ____________________________

begin
    for i = nelem*5+1:nelem*6
        n_indices = length(dad.ELEM[i].indices)
        for j = 1:n_indices
            indice = dad.ELEM[i].indices[j]
            dad.ELEM[i].valorCDC[j] = dad.NOS[indice,1]
        end
    end 

    for i = nelem*7+1:nelem*8
        n_indices = length(dad.ELEM[i].indices)
        for j = 1:n_indices
            indice = dad.ELEM[i].indices[j]
            dad.ELEM[i].valorCDC[j] = dad.NOS[indice,1]
        end
    end 
end
#_______________________________

H, G = calc_HeG(dad, npg)  #importante

A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b

x = A \ b

T, q = separa(dad, x) #importante

Ti = calc_Ti(dad, T, q, npg);

# Pós processamento ________________________

begin
    x_array = dad.pontos_internos[:,1]
    y_array = dad.pontos_internos[:,2]

    x_ordenado = sort(unique(x_array))
    y_ordenado = sort(unique(y_array))

    nx = length(unique(x_ordenado))
    ny = length(unique(y_ordenado))

    dx = (maximum((dad.pontos_internos[:,1]))-minimum((dad.pontos_internos[:,1])))/nx
    dy = (maximum((dad.pontos_internos[:,2]))-minimum((dad.pontos_internos[:,2])))/ny
end

#=
vi, ui = deriva_Ti_potencial(dad,Ti,dx,dy) 
vi = -vi
=#

ui, vi = deriva_Ti_potencial(dad,Ti,dx,dy) 

utotal = sqrt.(ui.^2 + vi.^2)

fig, ax, hm = BEM.heatmap(x_array, y_array, utotal)
BEM.Colorbar(fig[:, end+1], hm)
fig

escala = 0.05
BEM.quiver(x_array,y_array,escala*ui,escala*vi,color=utotal)