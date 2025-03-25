#Fazendo por linhas de corrente, lembrar que a vx = dfi/dy
# e vy = -di/dx

using DrWatson, Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))

r=1
L=5
h=2

nelem = 10  #Numero de elementos
order = 2

NPX = 20 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 20    #apenas números pares

## Formatação dos dados ________________________________________________

dad = format_dad(potencial_one_quarter(nelem, order), NPX, NPY) # dados

# Corrige CDC ____________________________


for i = nelem*4+1:nelem*5
    n_indices = length(dad.ELEM[i].indices)
    for j = 1:n_indices
        indice = dad.ELEM[i].indices[j]
        dad.ELEM[i].valorCDC[j] = dad.NOS[indice,2]
    end
end 


#____________________________
# Pontos internos ____________________________

#=
ntheta = 10
nr = 10

theta = range(0,2*pi,length=ntheta)
radius = range(1.5,4.5,length=nr)

dad.pontos_internos = zeros(nr*ntheta,2)

for i = 1:nr
    for j = 1:ntheta
        dad.pontos_internos[(i-1)*nr+j,1] = radius[i]*cos(theta[j])
        dad.pontos_internos[(i-1)*nr+j,2] = radius[i]*sin(theta[j])
    end
end
=#


n_int=50

dad.pontos_internos = zeros(n_int,2)

dad.pontos_internos[:,2] = range(0.2,h-0.2,length=n_int)
dad.pontos_internos[:,1] .= 4*L/5


#=
dad.pontos_internos = zeros(length(dad.NOS[nelem*order+1:nelem*order*2]),2)
dad.pontos_internos = dad.NOS[nelem*order+1:nelem*order*2,:]*0.75
=#

#________________
dx = zeros(length(dad.pontos_internos[:,1]))

for i =1:length(dad.pontos_internos[:,1])
    if i == 1

        dx[i] = sqrt((dad.pontos_internos[i+1,2] - dad.pontos_internos[i,2])^2 +
    (dad.pontos_internos[i+1,1] - dad.pontos_internos[i,1])^2)

    else

        dx[i] = sqrt((dad.pontos_internos[i,2] - dad.pontos_internos[i-1,2])^2 +
    (dad.pontos_internos[i,1] - dad.pontos_internos[i-1,1])^2)
    
    end
end 

#____________________________

H, G = calc_HeG(dad, npg)  #importante

A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b

x = A \ b

T, q = separa(dad, x) #importante

Ti = calc_Ti(dad, T, q, npg);

function deriva_Ti(Ti,dx)
    n_Ti = length(Ti)
    dTi = zeros(n_Ti)

    for i = 1:n_Ti
        if i == 1
            dTi[i] = (Ti[i+1] - Ti[i])/(dx[i])
        else
            dTi[i] = (Ti[i] - Ti[i-1])/(dx[i])
        end
    end
    return dTi
end

xi = dad.pontos_internos[:,1]
yi = dad.pontos_internos[:,2]

u = deriva_Ti(Ti,dx)