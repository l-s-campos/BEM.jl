## Início da análise
using DrWatson, Plots
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
#include(scriptsdir("my_functions.jl"))


#Dados de entrada

nelem = 10 #Numero de elementos
order = 2

NPX = 10 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

L = 1
Re = 100 # Na verdade é o 1/ν ou 1/μ

caso = "Canal entrada parabólica"
λ = 10^5


## Formatação dos dados ________________________________________________

begin

    if caso == "Cavidade"
        dad = format_dad(cavity(nelem, order,L,Re,λ), NPX, NPY)

    elseif caso == "Von Karman"
        nelem_circle = 15
        r = 0.1
        L = 30*r
        h = 10*r
        xc = 10*r

        dad = format_dad(von_karman2(nelem, order,L,h,r,xc,Re), NPX, NPY)

    elseif caso == "Expansion"
        L1 = 0.6
        h1 = 0.5
        L2= 2.5
        h2 = 1
        dad = format_dad(expansion(nelem, order, L1, h1, L2, h2,Re,λ))
    elseif caso == "Contração"
        L1 = 1.5
        h1 = 1
        L2= 1.5
        h2 = 0.5
        dad = format_dad(contraction(nelem, order, L1, h1, L2, h2,Re,λ))

    elseif caso == "Facing step"
        h=L/3
        y_step=h/2
        x_step=L/5
        dad = format_dad(facing_step(nelem, order,L,h,y_step,x_step,Re,λ), NPX, NPY)

    elseif caso == "Canal entrada parabólica"
        h = 0.1*L
        dpdx = -10^3   

        dad = format_dad(channel(nelem, order,h,L,Re,λ), NPX, NPY)

        #Corrige CDC

        for i = 3*nelem+1:nelem*4
            n_indices = length(dad.ELEM[i].indices)
            for j = 1:n_indices
                indice = dad.ELEM[i].indices[j]
                dad.ELEM[i].valorCDC[1,j] = dpdx*h^2/(2*1)*(((dad.NOS[indice,2]-h/2)/h)^2-1/4)
            end
        end 

    elseif caso == "Gradiente de pressão"

        dpdx = -10^3
        u_ana(y) = dpdx*L^2/(2*1)*((y/L)^2 - (y/L))
        dad = format_dad(gradiente_pressao(nelem, order,L,Re,λ), NPX, NPY)

    elseif caso == "Canal"

        h = 1
        L = 2*h

        dad = format_dad(channel(nelem, order,h,L,Re,λ), NPX, NPY)

    elseif caso == "Canal Simetria"
        h=L/3
        dad = format_dad(channel_simetria(nelem, order,h,L,Re,λ), NPX, NPY)

    elseif caso == "Placa plana"

        h=10*L
        Ls=L/4
        dad = format_dad(placa_plana_2(nelem, order,L,Ls,h,Re,λ), NPX, NPY)

    elseif caso == "Couette"
        u_ana(y) = y/L
        dad = format_dad(couette(nelem, order,L,Re,λ), NPX, NPY)

    elseif caso == "Taylor Couette"

        ra = 1
        ω_a = 0
        ω_b = 2
        rb = 2

        u_ana(r) = (ω_b * rb^2 - ω_a * ra^2) / (rb^2 - ra^2)*r + ((ω_a - ω_b) * ra^2 * rb^2 / (rb^2 - ra^2))/r

        dad = format_dad(Taylor_Couette(nelem, order,ra,rb,Re,λ), NPX, NPY)

        # Corrige CDC ____________________________
        
        tangencial = zeros(length(dad.normal),2)

        for i =1:length(dad.normal[:,1])
            tangencial[i,1] = -dad.normal[i,2]
            tangencial[i,2] = dad.normal[i,1]
        end

        # Cilindro interno ____________________________

        for i = nelem*3+1:nelem*4
            n_indices = length(dad.ELEM[i].indices)
            for j = 1:n_indices
                indice = dad.ELEM[i].indices[j]
                dad.ELEM[i].valorCDC[1,j] = (ω_a*ra)*tangencial[indice,1]
                dad.ELEM[i].valorCDC[2,j] = (ω_a*ra)*tangencial[indice,2]
            end
        end 

        # Cilindro externo ____________________________

        for i = nelem+1:nelem*2
            n_indices = length(dad.ELEM[i].indices)
            for j = 1:n_indices
                indice = dad.ELEM[i].indices[j]
                dad.ELEM[i].valorCDC[1,j] = -(ω_b*rb)*tangencial[indice,1]
                dad.ELEM[i].valorCDC[2,j] = -(ω_b*rb)*tangencial[indice,2]
            end
        end 

        #______________

        # Pontos internos ____________________________

        ntheta = 10
        nr = 10

        offset_theta = 0.1
        offset_radius= 0.1

        theta = range(offset_theta,pi/2-offset_theta,length=ntheta)
        radius = range(ra+offset_radius,rb-offset_radius,length=nr)

        dad.pontos_internos = zeros(nr*ntheta,2)

        for i = 1:nr
            for j = 1:ntheta
                dad.pontos_internos[(i-1)*nr+j,1] = radius[i]*cos(theta[j])
                dad.pontos_internos[(i-1)*nr+j,2] = radius[i]*sin(theta[j])
            end
        end

        #______________

    elseif caso == "Gravidade"

        g = 9.81
        h = 0.1
        L=5*h

        u_ana(y) = g/1*h^2*((y/h) - 1/2*(y/h)^2)
        τ_ana(y) = 1*g/1*h^2*(1/h - y/h^2)

        dad = format_dad(gravidade(nelem, order,L,h,g,Re,λ), NPX, NPY)

    elseif caso == "Step"

        h=L/2
        L_step = L/8
        x_step = L/4

        dad = format_dad(step(nelem, order, L, h, L_step, x_step, Re, λ), NPX, NPY)
    end
end

#=
#Von Karman
retangulo = [0.7, -0.25, 1.4, 0.25]
dad.pontos_internos = refina_local(dad,0.025,0.025,retangulo,true,1:4*order*nelem)
=#

#__________________________

# ==============Matrizes===============#

begin
    H, G = calc_HeG(dad, npg, interno = true)  #importante
    Hx, Gx, Hy, Gy = BEM.calc_dHedG(dad, 8)

    M = BEM.Monta_M_RIMd(dad, npg) # Para obter uma nova em outro Re, lembrar que ela é diretamente proporcional a Re
    Mx, My = BEM.Monta_deriv_M_RIMd(dad, 8)                                                                                                                                                                                                                                                                   

    A, b = BEM.aplicaCDC(H, G, dad)
end

# ==============Discretização============#

begin     
    nx = length(unique(dad.pontos_internos[:,1]))
    ny = length(unique(dad.pontos_internos[:,2]))

    x_order = sort(unique(dad.pontos_internos[:,1]))
    y_order = sort(unique(dad.pontos_internos[:,2]))

    dx = (maximum(x_order)-minimum(x_order))/nx
    dy = (maximum(y_order)-minimum(y_order))/ny
    
end

# =========Inicia variáveis===============#

begin  

    have_nolinear = false

    # IMPORTANTE
    #constant_term = true # Sem termo constante 
    constant_term = -dpdx*[i % 2 == 1 ? 1 : 0 for i in 1:2*(nc(dad)+ni(dad))] #Caso Grad de pressão
    #constant_term = -g*[i % 2 == 1 ? 1 : 0 for i in 1:2*(nc(dad)+ni(dad))] #Caso Gravidade

    n = ni(dad)+nc(dad)

    global iter = 0
    global erro_vel = 1
    global erro_derivative = 1
    C = zeros(size(M))

    global u = zeros(n,2)
    global t = zeros(nc(dad),2)
    global u_before = deepcopy(u)
    global nolinear = zeros(2*(n))

    global dudx=zeros(n)
    global dudy=zeros(n)
    global dvdx=zeros(n)
    global dvdy=zeros(n)

    global ux = zeros(2*n)
    global uy = zeros(2*n)
    global ux_before = zeros(2*n)

    interno = nc(dad)+1:n
    contorno = 1:nc(dad)

    interno_x2 = nc(dad)*2+1:2*n
    contorno_x2 = 1:2*nc(dad)
    
    relaxation = 0.01
end


#______________________________


# ==============Solução===============#

dNx,dNy = BEM.montaFs([dad.NOS; dad.pontos_internos],[dad.NOS; dad.pontos_internos])[2:3]

#=
p = (A,b,M,dNx,dNy)

function nolinear_solve(x,p)
    # A,b,M,dNx,dNy = p
    
    u[1:nc(dad),:] = separa(dad,x)[1]
    
    u[interno,1] =  x[2*nc(dad)+1:2:end]
    u[interno,2] =  x[2*nc(dad)+2:2:end]

    nolinear[1:2:2*n] = u[:,1].*p[4]*u[:,1] + u[:,2].*p[5]*u[:,1]
    nolinear[2:2:2*n] = u[:,1].*p[4]*u[:,2] + u[:,2].*p[5]*u[:,2]

    R = p[1]*x-p[2]+p[3]*nolinear
    
end

prob = NonlinearProblem(nolinear_solve, x*0, p)
#sol = solve(prob, NewtonRaphson(autodiff =AutoFiniteDiff()),maxiters = 10, abstol = 1e-6)
sol = solve(prob, Broyden(),maxiters = 1000, abstol = 1e-8)
=#

#Solução com relaxação

begin

    if have_nolinear
        
        while erro_vel > 10^-10
        
            # Derivadas 
            
            
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
            
            
            dudx = dNx*u[:,1]
            dudy = dNy*u[:,1]
            dvdx = dNx*u[:,2]
            dvdy = dNy*u[:,2]
            
            
            nolinear[1:2:2*n] = u[:,1].*dudx + u[:,2].*dudy
            nolinear[2:2:2*n] = u[:,1].*dvdx + u[:,2].*dvdy
            
        
            #____________________________________
            
            # Solution

            x = A \ (b - M*nolinear)
        
            u[1:nc(dad),:] = u[1:nc(dad),:]*(1-relaxation) + relaxation*separa(dad,x)[1]
            t =  t*(1-relaxation) + relaxation*separa(dad,x)[2]
            
            u[interno,1] =  u[interno,1]*(1-relaxation) + relaxation*x[2*nc(dad)+1:2:end]
            u[interno,2] =  u[interno,2]*(1-relaxation) + relaxation*x[2*nc(dad)+2:2:end]

            erro_vel = nrmse(u_before,u)
            
            println("Iteração: $iter; erro_vel = $erro_vel")
        
            u_before = deepcopy(u)
            iter = iter+1
        end
        

    
    elseif constant_term == false

        x = A \ (b)
    
        u[1:nc(dad),:],t = separa(dad,x)
        
        u[interno,1] =  x[2*nc(dad)+1:2:end]
        u[interno,2] =  x[2*nc(dad)+2:2:end]

    else
        x = A \ (b + M*(constant_term))
    
        u[1:nc(dad),:],t = separa(dad,x)
        
        u[interno,1] =  x[2*nc(dad)+1:2:end]
        u[interno,2] =  x[2*nc(dad)+2:2:end]
    end
end



#________________________

# Computando variáveis

begin 
    #=
    u[1:nc(dad),:] = separa(dad,sol.u)[1]

    u[interno,1] =  sol.u[2*nc(dad)+1:2:end]
    u[interno,2] =  sol.u[2*nc(dad)+2:2:end]
    =#
    uc = u[contorno,:]
    ui = u[interno,:]

    # Pela eq integral de contorno

    #=
    ddx,ddy = BEM.calc_derivative_int_M(dad,u,t,Mx,My)

    dudx_solfund = ddx[1:2:2*(nc(dad)+ni(dad))]
    dudy_solfund = ddy[1:2:2*(nc(dad)+ni(dad))]
    dvdx_solfund = ddx[2:2:2*(nc(dad)+ni(dad))]
    dvdy_solfund = ddy[2:2:2*(nc(dad)+ni(dad))]
    =#

    pressure = -λ*(dudx + dvdy) 
end


#=========Heatmap e Quiver=========#

x_array = dad.pontos_internos[:,1]
y_array = dad.pontos_internos[:,2]

utotal = sqrt.(ui[:,1].^2 + ui[:,2].^2)

ux = Hx * u[1:nc(dad),:]'[:] - Gx * t[1:nc(dad),:]'[:] - Mx * constant_term
uy = Hy * u[1:nc(dad),:]'[:] - Gy * t[1:nc(dad),:]'[:] - My * constant_term




BEM.heatmap(x_array, y_array, utotal)
fig, ax, hm = BEM.heatmap(x_array, y_array, utotal)
BEM.Colorbar(fig[:, end+1], hm,label="Velocidade [m/s]")
fig

escala = 10^(-1)
BEM.quiver(x_array,y_array,escala*ui[:,1],escala*ui[:,2],color=utotal)

# Plots

Plots.scatter!(ui[Int((nx)/2):nx:end,1],dad.pontos_internos[Int((nx)/2):nx:end,2],
marker=(3,:x,:red), label="Numérico")

   dudx = dNx*u[:,1]
            dudy = dNy*u[:,1]
            dvdx = dNx*u[:,2]
            dvdy = dNy*u[:,2]