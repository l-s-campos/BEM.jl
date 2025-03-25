using DrWatson,Plots,CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(scriptsdir("my_functions.jl"))

nelem = 15  #Numero de elementos
order = 2

NPX = 25 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 10    #apenas números pares

#Tempo

begin
    dt = 1e-3
    t_i = 0
    t_f = 1e1

    t = range(t_i, stop=t_f, step=dt)
    nt = length(t)
end

#_______________________________

# ======== INÍCIO ===============#

caso = "Cavidade"
Re = 1

begin

    println("Caso $caso, Reynolds = $Re") 
    
    if caso == "Cavidade"

        dad_u = format_dad(cavity_poisson(nelem, order,Re,"u"), NPX, NPY)
        dad_v = format_dad(cavity_poisson(nelem, order,Re,"v"), NPX, NPY)

    elseif caso == "Von Karman"

        r = 0.5 #Não pode ser muito pequeno se não dá errado
        L=20*r
        h = 10*r
        xc = 5*r

        dad_u = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"u"), NPX, NPY)
        dad_v = format_dad(cilindro_chorin(nelem, order,r,L,h,xc,Re,"v"), NPX, NPY)
        
    end

    # Refinamento local

    #=
    retangulo = [0.3, 1.2, 0.7, 1.8]
    intervalo_objeto = 1:2*order*nelem

    dad_u.pontos_internos = refina_local(dad_u,0.01,0.01,retangulo,true,intervalo_objeto)
    dad_v.pontos_internos = refina_local(dad_v,0.01,0.01,retangulo,true,intervalo_objeto)
    =#

    #________________________

    dad = dad_u

    Ht = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    Gt = zeros(ni(dad_u)+nc(dad_u),nc(dad_u),2)
    M = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)

    A_Houbolt = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Houbolt = zeros(ni(dad_u)+nc(dad_u),2)

    A_Adams = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Adams = zeros(ni(dad_u)+nc(dad_u),2)

    A_Euler = zeros(ni(dad_u)+nc(dad_u),ni(dad_u)+nc(dad_u),2)
    b_Euler = zeros(ni(dad_u)+nc(dad_u),2) 

    # Velocidade u

    Ht[:,:,1], Gt[:,:,1] = BEM.calc_HeGt(dad_u)
    M[:,:,1] = BEM.Monta_M_RIMd(dad_u, npg)
    A_Houbolt[:,:,1], b_Houbolt[:,1] = BEM.aplicaCDC(Ht[:,:,1] - 11*M[:,:,1]/(6*dt), Gt[:,:,1], dad_u)
    A_Adams[:,:,1], b_Adams[:,1] = BEM.aplicaCDC(Ht[:,:,1] - (24/(55*dt))*M[:,:,1], Gt[:,:,1], dad_u)
    A_Euler[:,:,1], b_Euler[:,1] = BEM.aplicaCDC(Ht[:,:,1] - M[:,:,1]/dt, Gt[:,:,1], dad_u)

    # Velocidade v

    Ht[:,:,2], Gt[:,:,2] = BEM.calc_HeGt(dad_v)
    M[:,:,2] = BEM.Monta_M_RIMd(dad_v, npg)
    A_Houbolt[:,:,2], b_Houbolt[:,2] = BEM.aplicaCDC(Ht[:,:,2] - 11*M[:,:,2]/(6*dt), Gt[:,:,2], dad_v)
    A_Adams[:,:,2], b_Adams[:,2] = BEM.aplicaCDC(Ht[:,:,2] - (24/(55*dt))*M[:,:,2], Gt[:,:,2], dad_v)
    A_Euler[:,:,2], b_Euler[:,2] = BEM.aplicaCDC(Ht[:,:,2] - M[:,:,2]/dt, Gt[:,:,2], dad_v)

end

#_______________________________

# Discretization

begin
    nx = length(unique(dad_u.pontos_internos[:,1]))
    ny = length(unique(dad_u.pontos_internos[:,1]))

    dx = (maximum(dad_u.pontos_internos[:,1])-minimum(dad_u.pontos_internos[:,1]))/nx
    dy = (maximum(dad_u.pontos_internos[:,2])-minimum(dad_u.pontos_internos[:,2]))/ny

    x_order = sort(unique(dad_u.pontos_internos[:,1]))
    y_order = sort(unique(dad_u.pontos_internos[:,2]))

    Cₒ =  dt*1/dx
    println("Número de Courant:$Cₒ")
end

#_______________

# Iniciando variáveis

begin
    n = ni(dad)+nc(dad)

    global nolinear = zeros(2*(ni(dad)+nc(dad)))
    global u = zeros(ni(dad)+nc(dad),2,nt)
    global u_dot = zeros(nc(dad) + ni(dad),2,nt)
    global dpdx = zeros(nc(dad)+ni(dad))
    global dpdy = zeros(nc(dad)+ni(dad))

    global p = zeros(ni(dad)+nc(dad),nt)
    global p_flux = zeros(nc(dad))
    global t = zeros(nc(dad),2)

    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)

    a = 20 #Compressibilidade aritificial
    generico = false
    temporal_method = "Euler"

    dist_index = zeros(nc(dad))
    index_dist = zeros(nc(dad))

    for k=1:nc(dad)
        dist = sqrt.((dad.NOS[k,1] .- dad.pontos_internos[:,1]).^2 .+ (dad.NOS[k,2] .- dad.pontos_internos[:,2]).^2)
        dist_index[k] = minimum(dist)
        index_dist[k] = Int(findfirst(x -> x == dist_index[k], dist))
    end

    dNx,dNy = BEM.montaFs([dad.NOS;dad.pontos_internos])[2:3]
    Deriva_x = Dₓ(nx, ny, dx) 
    Deriva_y = Dy(nx, ny, dy)

    dist_cont_cont = zeros(nc(dad))

    for i=1:nc(dad)
        if i == nc(dad)
            dist_cont_cont[i] = sqrt((dad.NOS[i,1]-dad.NOS[i-1,1])^2 + (dad.NOS[i,2]-dad.NOS[i-1,2])^2)
        else
            dist_cont_cont[i] = sqrt((dad.NOS[i+1,1]-dad.NOS[i,1])^2 + (dad.NOS[i+1,2]-dad.NOS[i,2])^2)
        end
    end
end

#_______________

#Loop solução

for i=4:nt

    # Derivadas do campo de velocidade
    
    #=
    if generico

        dudx,dudy = deriva_interno_generico(dad,u[interno,1,i-1],dx,dy)
        dvdx,dvdy = deriva_interno_generico(dad_v,u[interno,2,i-1],dx,dy)

        dudx_cont,dudy_cont = deriva_contorno_generico(dad,u[contorno,1,i-1],t[:,1]) # Transição = 4*order*nelem
        dvdx_cont,dvdy_cont = deriva_contorno_generico(dad_v,u[contorno,2,i-1],t[:,2])
        
    else
        dudx = Deriva_x*u[interno,1,i-1]
        dudy = Deriva_y*u[interno,1,i-1]
        dvdx = Deriva_x*u[interno,2,i-1]
        dvdy = Deriva_y*u[interno,2,i-1]

        dudx_cont,dudy_cont = deriva_contorno(dad,u[contorno,1,i-1],t[:,1]) #Talvez multiplicar por Re/3 pra tirar o E
        dvdx_cont,dvdy_cont = deriva_contorno(dad_v,u[contorno,2,i-1],t[:,2])
    end
    =#

    #=
    dudx_int = Deriva_x*u[interno,1,i-1]
    dudy_int = Deriva_y*u[interno,1,i-1]
    dvdx_int = Deriva_x*u[interno,2,i-1]
    dvdy_int = Deriva_y*u[interno,2,i-1]
    
    
    dudx_cont,dudy_cont = deriva_contorno(dad,u[contorno,1,i-1],t[:,1]/dad.k,dist_cont_cont) 
    dvdx_cont,dvdy_cont = deriva_contorno(dad_v,u[contorno,2,i-1],t[:,2]/dad.k,dist_cont_cont)

    
    nolinear[2*nc(dad)+1:2:2*(nc(dad)+ni(dad))] = u[interno,1,i-1].*dudx_int + u[interno,2,i-1].*dudy_int
    nolinear[2*nc(dad)+2:2:2*(nc(dad)+ni(dad))] = u[interno,1,i-1].*dvdx_int + u[interno,2,i-1].*dvdy_int

    nolinear[1:2:2*nc(dad)] = u[contorno,1,i-1].*dudx_cont + u[contorno,2,i-1].*dudy_cont
    nolinear[2:2:2*nc(dad)] = u[contorno,1,i-1].*dvdx_cont + u[contorno,2,i-1].*dvdy_cont
    =#

    
    dudx = dNx*u[:,1,i-1]
    dudy = dNy*u[:,1,i-1]

    dvdx = dNx*u[:,2,i-1]
    dvdy = dNy*u[:,2,i-1]
    
    nolinear[1:2:2*n] = u[:,1,i-1].*dudx + u[:,2,i-1].*dudy
    nolinear[2:2:2*n] = u[:,1,i-1].*dvdx + u[:,2,i-1].*dvdy

    
    #_________________

    # Obtendo a velocidade

    #=
    if temporal_method == "Houbolt"
    
        u[contorno,1,i], u[interno,1,i], t[:,1], u_dot[:,1,:] = 
        poisson_Houbolt(dad,nolinear[1:2:end] + dpdx,i,u[:,1,:],u_dot[:,1,:],A_Houbolt[:,:,1],b_Houbolt[:,1],M[:,:,1])

        u[contorno,2,i], u[interno,2,i], t[:,2], u_dot[:,2,:] = 
        poisson_Houbolt(dad,nolinear[2:2:end] + dpdy,i,u[:,2,:],u_dot[:,2,:],A_Houbolt[:,:,2],b_Houbolt[:,2],M[:,:,2]) 

    elseif temporal_method == "Adams"
    
        u[contorno,1,i], u[interno,1,i], t[:,1], u_dot[:,1,:] = 
        poisson_Adams(dad,nolinear[1:2:end]+ dpdx,i,u[:,1,:],u_dot[:,1,:],A_Houbolt[:,:,1],A_Adams[:,:,1],b_Houbolt[:,1],b_Adams[:,1],M[:,:,1])

        u[contorno,2,i], u[interno,2,i], t[:,2], u_dot[:,2,:] = 
        poisson_Adams(dad,nolinear[2:2:end]+ dpdy,i,u[:,2,:],u_dot[:,2,:],A_Houbolt[:,:,2],A_Adams[:,:,2],b_Houbolt[:,2],b_Adams[:,2],M[:,:,2]) 
    else

        u[contorno,1,i], u[interno,1,i], t[:,1] = 
        poisson_Euler(dad,nolinear[1:2:end]+ dpdx,i,u[:,1,:],A_Euler[:,:,1],b_Euler[:,1],M[:,:,1])
    
        u[contorno,2,i], u[interno,2,i], t[:,2] = 
        poisson_Euler(dad,nolinear[2:2:end]+ dpdy,i,u[:,2,:],A_Euler[:,:,2],b_Euler[:,2],M[:,:,2]) 

    end
    =#
  
    #=
    u[contorno,1,i], u[interno,1,i], t[:,1], u_dot[:,1,:] = 
    poisson_Houbolt(dad,nolinear[1:2:end] + dpdx,i,u[:,1,:],u_dot[:,1,:],A_Houbolt[:,:,1],b_Houbolt[:,1],M[:,:,1])

    u[contorno,2,i], u[interno,2,i], t[:,2], u_dot[:,2,:] = 
    poisson_Houbolt(dad,nolinear[2:2:end] + dpdy,i,u[:,2,:],u_dot[:,2,:],A_Houbolt[:,:,2],b_Houbolt[:,2],M[:,:,2]) 
    =#

    
    u[contorno,1,i], u[interno,1,i], t[:,1] = 
    poisson_Euler(dad,nolinear[1:2:end]+ dpdx,i,u[:,1,:],A_Euler[:,:,1],b_Euler[:,1],M[:,:,1])

    u[contorno,2,i], u[interno,2,i], t[:,2] = 
    poisson_Euler(dad,nolinear[2:2:end]+ dpdy,i,u[:,2,:],A_Euler[:,:,2],b_Euler[:,2],M[:,:,2])
    

    # Obtendo a pressão e sua derivada

    #Derivadas atualizadas __________

    #=
    if generico
        dudx,dudy = deriva_interno_generico(dad,u[interno,1,i],dx,dy)
        dvdx,dvdy = deriva_interno_generico(dad,u[interno,2,i],dx,dy)
    else   
        dudx = Deriva_x*u[interno,1,i]
        dvdy =Deriva_y*u[interno,2,i]
    end
    =#

    #=
    dudx_int = Deriva_x*u[interno,1,i]
    dvdy_int = Deriva_y*u[interno,2,i]
    =#

    
    dudx = dNx*u[:,1,i]
    dvdy = dNy*u[:,2,i]
    

    # Achando a pressão
    #=
    if temporal_method == "Euler" || temporal_method == "Adams"
        p[interno,i] = p[interno,i-1] - a*dt*(dudx + dvdy)
    else
        p[interno,i] = (18*p[interno,i-1] - 9*p[interno,i-2] + 2*p[interno,i-3] - (6*dt)*a*(dudx + dvdy))/11
    end
    =#

    p[:,i] = p[:,i-1] - a*dt*(dudx[:] + dvdy[:])

    #Aplicando CDC
    #=
    for k=1:nc(dad)

        if k in 1:2*order*nelem || k in 3*order*nelem+1:4*order*nelem # Paredes com condição de fluxo zero

            p[k] = p[nc(dad) + Int(index_dist[k])]

        else # Pressão especificada
            p[k] = 0
            p_flux[k] = (p[k] - p[nc(dad) + Int(index_dist[k])])/dist_index[k]
        end
    end   
    =#

    dpdx = dNx*p[:,i]
    dpdy = dNy*p[:,i]
    
    #=
    if generico
        dpdx_int,dpdy_int = deriva_interno_generico(dad,p[interno,i],dx,dy)
        dpdx_cont,dpdy_cont = deriva_contorno_generico(dad,p[contorno,i],p_flux)
    else   
        dpdx_int = Deriva_x*p[interno,i]
        dpdy_int = Deriva_y*p[interno,i]
        dpdx_cont,dpdy_cont = deriva_contorno(dad,p[contorno,i],p_flux)
    end
    =#

    #=
    dpdx_int = Deriva_x*p[interno,i]
    dpdy_int = Deriva_y*p[interno,i]
    dpdx_cont,dpdy_cont = deriva_contorno(dad,p[contorno,i],p_flux)

    dpdx = [dpdx_cont;dpdx_int]
    dpdy = [dpdy_cont;dpdy_int]
    =#

    erro = nrmse(u[:,:,i-1],u[:,:,i])

    u_mean = mean(u[:,1,i])

    println("Iteração: $(i-3); Erro: $erro; u_mean = $u_mean")
end

#_______________

# Mapa de cores

x_array = dad.pontos_internos[:,1]
y_array = dad.pontos_internos[:,2]

uc = u[contorno,1,8400]

ui = u[interno,1,8400]
vi = u[interno,2,8400]

p_cont = p[contorno,8400]
p_int = p[interno,8400]

utotal = sqrt.(ui.^2 + vi.^2)

BEM.heatmap(x_array, y_array, utotal)
fig, ax, hm = BEM.heatmap(x_array, y_array, utotal)
BEM.Colorbar(fig[:, end+1], hm,label="Velocidade [m/s]")
fig

scale = 10^-1
BEM.quiver(x_array,y_array,scale*ui,scale*vi,color=utotal)


Plots.plot(vi[Int(nx/2):Int(nx/2)+nx])
Plots.plot(ui[Int((nx)/2):nx:end],y_order)
Plots.plot(p[contorno,end])

#Comparação com outros metodos

u_LBM =
[-2.09163243839448e-05 0
-8.70563713048693e-05 0.0101010101010101
-0.000151069067967581 0.0202020202020202
-0.000212901209139955 0.0303030303030303
-0.000272611262211442 0.0404040404040404
-0.000330292493051965 0.0505050505050505
-0.000386059604300867 0.0606060606060606
-0.000440039197720849 0.0707070707070707
-0.000492363156282164 0.0808080808080808
-0.000543164032681337 0.0909090909090909
-0.000592571910849309 0.101010101010101
-0.000640712279687244 0.111111111111111
-0.000687704637259365 0.121212121212121
-0.000733661592898752 0.131313131313131
-0.000778688319210126 0.141414141414141
-0.000822882238295596 0.151515151515152
-0.000866332863473924 0.161616161616162
-0.000909121741174019 0.171717171717172
-0.000951322449735451 0.181818181818182
-0.000993000631099295 0.191919191919192
-0.00103421403007322 0.202020202020202
-0.00107501253343074 0.212121212121212
-0.00111543819248142 0.222222222222222
-0.00115552522985085 0.232323232323232
-0.00119530001847059 0.242424242424242
-0.0012347810379434 0.252525252525253
-0.00127397879832676 0.262626262626263
-0.0013128957388228 0.272727272727273
-0.00135152609231176 0.282828282828283
-0.0013898557244557 0.292929292929293
-0.00142786193865783 0.303030303030303
-0.00146551325629195 0.313131313131313
-0.00150276916361842 0.323232323232323
-0.00153957983521319 0.333333333333333
-0.00157588582542403 0.343434343434343
-0.00161161773796722 0.353535353535354
-0.00164669586536273 0.363636363636364
-0.00168102980855106 0.373737373737374
-0.00171451806877005 0.383838383838384
-0.00174704762223073 0.393939393939394
-0.00177849347036085 0.404040404040404
-0.00180871817631464 0.414141414141414
-0.00183757138165271 0.424242424242424
-0.00186488931401935 0.434343434343434
-0.00189049428144945 0.444444444444444
-0.00191419416445766 0.454545454545455
-0.0019357819041908 0.464646464646465
-0.00195503500171823 0.474747474747475
-0.00197171505080189 0.484848484848485
-0.00198556696507483 0.494949494949495
-0.00199631459090902 0.505050505050505
-0.00200367125738855 0.515151515151515
-0.0020073371669071 0.525252525252525
-0.00200698668801539 0.535353535353535
-0.00200227928842191 0.545454545454545
-0.00199285662358643 0.555555555555556
-0.00197834261880802 0.565656565656566
-0.00195834338609782 0.575757575757576
-0.00193244753051629 0.585858585858586
-0.00190022629216129 0.595959595959596
-0.00186123409007973 0.606060606060606
-0.00181500892725123 0.616161616161616
-0.00176107321645498 0.626262626262626
-0.00169893449247948 0.636363636363636
-0.00162808656508101 0.646464646464646
-0.00154801058300778 0.656565656565657
-0.00145817655796997 0.666666666666667
-0.00135804482331621 0.676767676767677
-0.00124706796988437 0.686868686868687
-0.0011246927382699 0.696969696969697
-0.000990362402117625 0.707070707070707
-0.000843519126481349 0.717171717171717
-0.000683606825981004 0.727272727272727
-0.000510074012176578 0.737373737373737
-0.000322377142237194 0.747474747474747
-0.000119983964578532 0.757575757575758
9.76226430674795e-05 0.767676767676768
0.00033094083993171 0.777777777777778
0.00058044552000374 0.787878787878788
0.000846584028726109 0.797979797979798
0.00112977145808322 0.808080808080808
0.00143038599554646 0.818181818181818
0.00174876393515479 0.828282828282828
0.00208519481117641 0.838383838383838
0.00243991634050952 0.848484848484849
0.00281310961622913 0.858585858585859
0.00320489436950187 0.868686868686869
0.00361532472313704 0.878787878787879
0.00404438548507195 0.888888888888889
0.00449198939052372 0.898989898989899
0.00495797576539467 0.909090909090909
0.00544211101926778 0.919191919191919
0.00594409224786119 0.929292929292929
0.00646355435571347 0.939393939393939
0.00700008354588903 0.94949494949495
0.0075532372822733 0.95959595959596
0.00812257608750025 0.96969696969697
0.00870770055806972 0.97979797979798
0.0093082586045877 0.98989898989899
0.00992473339136288 1
]

Plots.scatter(u_ACM[Int(nx/2):nx:end],y_order, label="u_ACM", xaxis="u [m/s]",yaxis="Posição vertical [m]",marker=(3,:x,:red))
Plots.scatter!(u_elastico[Int(nx/2):nx:end],y_order, label="u_elástico",marker=(3,:x,:black))
Plots.scatter!(u_projecao_Euler[Int(nx/2):nx:end],y_order, label="u_projeção (Euler)",marker=(3,:x,:blue))
Plots.scatter!(u_projecao_Houbolt[Int(nx/2):nx:end],y_order, label="u_projeção (Houbolt)",marker=(3,:x,:green))
Plots.scatter!(-u_ωψ[Int(nx/2):nx:end]/maximum(abs.(u_ωψ[Int(nx/2):nx:end])),y_order, label="u_ψω",marker=(3,:x,:purple))
Plots.plot!(10^2*u_LBM[:,1],u_LBM[:,2], label="u_LBM",color=:orange)


Plots.scatter(ui[Int(nx/2):nx:end],y_order, label="u_projeção (Houbolt)",marker=(3,:x,:green))