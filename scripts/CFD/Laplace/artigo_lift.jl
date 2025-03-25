using DrWatson, Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))


nelem = 4  #Numero de elementos
order = 2
NPX = 5 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

## Formatação dos dados ________________________________________________

dad_v = format_dad(NACA_artigo_v(nelem, order), NPX, NPY)
dad_w = format_dad(NACA_artigo_w(nelem, order), NPX, NPY)

final_NACA = (n-1)*order*nelem

#Corrige_CDC para w

for i = nelem*(n) + 4*nelem + 1 : nelem*(n) + 4*nelem*(4)
    n_indices = length(dad_w.ELEM[i].indices)
    for j = 1:n_indices
        indice = dad_w.ELEM[i].indices[j]
        dad_w.ELEM[i].valorCDC[j] = 1*dad_w.normal[indice,1]
    end
end 


#________________________________

function calc_potencial(dad)

    H, G = calc_HeG(dad, npg)  #importante

    A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b

    if dad == dad_v
        line_CDC = zeros(nc(dad))
        position = 3
        line_CDC[1+position] = 1
        line_CDC[final_NACA-position] = -1

        A[1+position,:] = line_CDC'
        b[1+position] = 1
    end

    x = A \ b

    T, q = separa(dad, x) #importante

    Ti = calc_Ti(dad, T, q, npg)

    return T,Ti,q
end

function deriva_T(dad,T,arg=1)
    nT = length(T)
    dT = zeros(nT)

    if arg == 1
        dx = sqrt((dad.NOS[10,1] - dad.NOS[9,1])^2 + (dad.NOS[10,2] - dad.NOS[9,2])^2)
    else
        dx = sqrt((dad.pontos_internos[10,1] - dad.pontos_internos[9,1])^2 + (dad.pontos_internos[10,2] - dad.pontos_internos[9,2])^2)
    end

    for i = 1:nT
        if i == 1
            dT[i] = (T[i+1] - T[i])/dx
        elseif i == nT
            dT[i] = (T[i] - T[i-1])/dx
        else
            dT[i] = (T[i+1] - T[i-1])/(dx)
        end
    end
    return dT
end

T_v,Ti_v,q_v = calc_potencial(dad_v)

T_w,Ti_w, q_w = calc_potencial(dad_w)


intervalo = 1:Int(final_NACA)

u_v = deriva_T(dad_v,T_v[intervalo])
u_w = deriva_T(dad_w,T_w[intervalo])

x_chord = range(0,1,length=Int((length(intervalo))/2))

ponto = 20

lower = Int(length(x_chord)*2 - (ponto-1))
upper = Int(ponto)

# aerofolio no sentido anti horario

tau = ((u_w[upper]) - (u_w[lower]))/
((u_v[upper]) - (-u_v[lower]))

lift = 2*tau


fi = (-tau)*T_v + T_w
u_cont = (-tau)*u_v + u_w


Plots.scatter(x_chord,reverse(u_cont[100:length(x_chord)].^2) ,xlabel="Chord line [m]",ylabel="Coeficiente de pressão",
label="inferior")
Plots.scatter!(x_chord,-(u_cont[length(x_chord)+2:length(x_chord)*2].^2),label="superior")



#Plots.scatter(dad_v.NOS[1:lower,1],dad_v.NOS[1:lower,2])

#Plots.scatter(dad_v.NOS[1:upper,1],dad_v.NOS[1:upper,2])

#=
alphas = 0:2:18
lift_artigo = collect(alphas) .* 0.9 ./ 8
lift_array  = [0.0003161381112205804,
0.23493297665536664,
 0.47027996326633065,
 0.7052255916562015,
 0.9402599547409319,
 1.174053107106474,
 1.409670176607804,
 1.6445091644700975,
 1.8777266810284106,
 2.11280995562874]

erro_average = nrmse(lift_artigo,lift_array)
erro_L2= nme(lift_artigo,lift_array)
erro_max = maximum(abs.((lift_array .- lift_artigo)/(maximum(lift_artigo)-minimum(lift_artigo))))

Plots.plot(alphas,lift_artigo,label="Artigo", color=:blue,yaxis="CL", xaxis="Ângulo de ataque [°]")
Plots.scatter!(alphas,lift_array,label="BEM Master",marker=(3, :x, :red))
annotate!(15, 0.75, text("Erro médio: 2.59e-2", :black, 8))
annotate!(15, 0.6, text("Erro L2: 2.20e-2", :black, 8))
annotate!(15, 0.45, text("Erro Máximo: 4.34e-2", :black, 8))
=#

