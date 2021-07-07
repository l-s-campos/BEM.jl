function calc_HeG(dad::potencial,npg=8)
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS,1)
  H=zeros(n,n)
  G=zeros(n,n)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  contafonte=1
  for elem_i in dad.ELEM  #Laço dos pontos fontes
    for ind_elem = elem_i.indices
      pf = dad.NOS[ind_elem,:]   # Coordenada (x,y)  dos pontos fonte
      for elem_j in dad.ELEM  #Laço dos elementos
        x = dad.NOS[elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
        Δelem=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
        eet=(elem_j.ξs[end] -elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
        N_geo,~=calc_fforma(eet,elem_j)
        ps=N_geo'*x
        b=norm(ps'-pf)/norm(Δelem)
        eta,Jt=sinhtrans(qsi,eet,b)
        h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad)
        H[contafonte,elem_j.indices]=h
        G[contafonte,elem_j.indices]=g
      end
      contafonte+=1
    end
  end
  for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    H[i,i]+=-0.5
  end
  H,G
end

#__________________________________________________________________________________________________________
"Funcao para calcular a temperatura pelo metodo da integracao de contorno utilizada no estimador de erro recursivo"
function integraelem(pf,x,eta,w,elem,dad::potencial)
  h = zeros(Float64,size(elem))
  g = zeros(Float64,size(elem))
  
  for k = 1:size(w,1)
    N,dN = calc_fforma(eta[k],elem)
    pg = N'*x    # Ponto de gauss interpolador
    r = pg'-pf      # Distancia entre ponto de gauss e ponto fonte
    dxdqsi = dN'*x   # dx/dξ & dy/dξ
    dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
    sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
    sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
    Qast,Tast=calsolfund(r,[sy,-sx],dad)
    h+=N*Qast*dgamadqsi*w[k]
    g+=N*Tast*dgamadqsi*w[k]
    
  end
  h,g
end

function calc_Ti(dad::potencial,T,q,npg=8)
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.pontos_internos,1)
  Ti=zeros(n)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  
  for i= 1:n
    pf = dad.pontos_internos[i,:]   # Coordenada (x,y) dos pontos fonte
    for elem_j in dad.ELEM  #Laço dos elementos
      x = dad.NOS[elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
      Δelem=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
      eet=(elem_j.ξs[end] -elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
      N_geo,~=calc_fforma(eet,elem_j)
      ps=N_geo'*x
      b=norm(ps'-pf)/norm(Δelem)
      eta,Jt=sinhtrans(qsi,eet,b)
      h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad)
      Ti[i]+=h'*T[elem_j.indices]-g'*q[elem_j.indices]
    end
  end
  Ti
end
function calc_Ti(dad::potencial_iso,T,q,npg=8)
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.pontos_internos,1)
  Ti=zeros(n)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  
  for i= 1:n
    pf = dad.pontos_internos[i,:]   # Coordenada (x,y) dos pontos fonte
    for elem_j in dad.ELEM  #Laço dos elementos
      xf=elem_j.limites[:,2]
      x0=elem_j.limites[:,1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
      Δelem=xf-x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
      eet=2*dot(Δelem,pf.-x0)/norm(Δelem)^2-1
      pc=dad.pontos_controle[:,elem_j.indices]
      # @infiltrate
      cf=pc[1:2,:]./pc[4,:]'
      N,dN=calc_fforma(eet,elem_j,pc[4,:])        
      ps=cf*N
      b=norm(ps-pf)/norm(Δelem)

      eta,Jt=sinhtrans(qsi,eet,b)
      h,g= integrabezier(pf,cf,pc[4,:],eta,w.*Jt,elem_j,dad)
    
      Ti[i]+=h'*T[elem_j.indices]-g'*q[elem_j.indices]
    end
  end
  Ti
end

function calc_Aeb(dad::potencial,npg=8)
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS,1)
  ni = size(dad.pontos_internos,1)
  A=zeros(n+ni,n)
  B=zeros(n+ni)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  
  for i =1:n+ni  #Laço dos pontos fontes
    pf = [dad.NOS;dad.pontos_internos][i,:]   # Coordenada (x,y)  dos pontos fonte
    for elem_j in dad.ELEM  #Laço dos elementos
      x = dad.NOS[elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
      Δelem=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
      eet=(elem_j.ξs[end] -elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
      N_geo,~=calc_fforma(eet,elem_j)
      ps=N_geo'*x
      b=norm(ps'-pf)/norm(Δelem)
      eta,Jt=sinhtrans(qsi,eet,b)
      h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad)
      h[elem_j.indices.==i]=h[elem_j.indices.==i].-0.5
      if elem_j.tipoCDC==1
        A[i,elem_j.indices]=h
        B[i]+=dot(g,elem_j.valorCDC)
      else
        A[i,elem_j.indices]=-g
        B[i]+=-dot(h,elem_j.valorCDC)
      end
      
    end
  end
  [A [zeros(n,ni);-diagm(ones(ni))]],B
end


function calc_HeG(dad::potencial,b1,b2,npg=8)
  n1 = size(b1,1)
  n2 = size(b2,1)
  H=zeros(n1,0)
  G=zeros(n1,0)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM[b2]  #Laço dos elementos
    x = [dad.NOS;dad.pontos_internos][elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
    
    h1 = zeros(n1,size(elem_j));
    g1 = zeros(n1,size(elem_j));
    for i =1:n1
      ind= b1[i]  #Laço dos pontos fontes
      pf = [dad.NOS;dad.pontos_internos][ind,:]   # Coordenada (x,y)  dos pontos fonte
      Δelem=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
      eet=(elem_j.ξs[end] -elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
      N_geo,~=calc_fforma_gen(eet,elem_j.ξs)
      ps=N_geo'*x
      b=norm(ps'-pf)/norm(Δelem)
      eta,Jt=sinhtrans(qsi,eet,b)
      h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad)
      h[elem_j.indices.==ind]=h[elem_j.indices.==ind].-0.5
      
      h1[i,:] += h
      g1[i,:] += g
    end
    H= [H h1]
    G= [G g1]
  end
  H,G
end


function calc_HeG_interp(dad::potencial,b1,b2,npg=8,ninterp=3)
  collocCoord=[dad.NOS;dad.pontos_internos][b1,:]
  xmax=maximum(collocCoord,dims=1)
  xmin=minimum(collocCoord,dims=1)

  xs=criapontosinterp(ninterp)
  fontes,L,ninterp1,ninterp2=gera_interpolação(ninterp,collocCoord,xmax,xmin,xs)

  H=zeros(ninterp1*ninterp2,0)
  G=zeros(ninterp1*ninterp2,0)
  n1,n2=Nlinear(xs)
  xks=n1*xmin+n2*xmax
  
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM[b2]  #Laço dos elementos
    x = dad.NOS[elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
    
    h1 = zeros(ninterp1*ninterp2,size(elem_j));
    g1 = zeros(ninterp1*ninterp2,size(elem_j));
    ci=0
    for i2 =1:ninterp1
      for i1 =1:ninterp2
        ci+=1

        pf = [xks[i1,1],xks[i2,2]]   # Coordenada (x,y)  dos pontos fonte
        Δelem=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
        eet=(elem_j.ξs[end] -elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
        N_geo,~=calc_fforma_gen(eet,elem_j.ξs)
        ps=N_geo'*x
        b=norm(ps'-pf)/norm(Δelem)
        eta,Jt=sinhtrans(qsi,eet,b)
        h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad)

        h1[ci,:] += h
        g1[ci,:] += g
      end
    end
    H= [H h1]
    G= [G g1]
  end
  L,H,G
end

function gera_interpolação(ninterp,NOS,xmax,xmin,xs,ϵ=1e-6)
if(abs(xmax[1]-xmin[1])<ϵ)
  fontes=(2. .*(NOS[:,2] .-xmin[2])./(xmax[2]-xmin[2]).-1);
  L=lagrange(fontes,xs,ninterp);
  ninterp2=ninterp
  ninterp1=1
elseif(abs(xmax[2]-xmin[2])<ϵ)
  fontes=(2. .*(NOS[:,1] .-xmin[1])./(xmax[1]-xmin[1]).-1);
  L=lagrange(fontes,xs,ninterp);
  ninterp2=1
  ninterp1=ninterp
else
  fontes=[(2. .*(NOS[:,1] .-xmin[1])./(xmax[1]-xmin[1]).-1) (2. .*(NOS[:,2] .- xmin[2])./(xmax[2]-xmin[2]).-1)]
  L=lagrange(fontes,xs,ninterp,xs,ninterp)
  ninterp2=ninterp
  ninterp1=ninterp
end
fontes,L,ninterp1,ninterp2
end



function calc_HeG(dad::potencial_iso,npg=8)
  # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  nfonte = size(dad.NOS,1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS,1)
  H=zeros(n,n)
  G=zeros(n,n)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM  #Laço dos elementos
    xf=elem_j.limites[:,2]
    x0=elem_j.limites[:,1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
    Δelem=xf-x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
    pc=dad.pontos_controle[:,elem_j.indices]
    # @infiltrate
    cf=pc[1:2,:]./pc[4,:]'
    for contafonte =1 : nfonte
      pf = dad.NOS[contafonte,:]   # Coordenada (x,y)  dos pontos fonte
      eet=2*dot(Δelem,pf.-x0)/norm(Δelem)^2-1
        N,dN=calc_fforma(eet,elem_j,pc[4,:])        
        ps=cf*N
        b=norm(ps-pf)/norm(Δelem)
        eta,Jt=sinhtrans(qsi,eet,b)
        h,g= integrabezier(pf,cf,pc[4,:],eta,w.*Jt,elem_j,dad)
              
        H[contafonte,elem_j.indices]+=h
        G[contafonte,elem_j.indices]+=g
      
    end
  end
   H-dad.E/2,G
end
function integrabezier(pf,cf,we,eta,w,elem::bezier,prob::potencial_iso)
  h = zeros(Float64,size(elem))
  g = zeros(Float64,size(elem))
  for k = 1:size(w,1)
    N,dN = calc_fforma(eta[k],elem,we)
    pg = cf*N    # Ponto de gauss interpolador
    r = pg-pf      # Distancia entre ponto de gauss e ponto fonte
    dxdqsi = cf*dN   # dx/dξ & dy/dξ
    dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
    sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
    sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
    Qast,Tast=calsolfund(r,[sy,-sx],prob)
    # h+=N*dgamadqsi*w[k]
    # g+=N*dgamadqsi*w[k]
    h+=N*Qast*dgamadqsi*w[k]/2
    g+=N*Tast*dgamadqsi*w[k]/2
    
  end
  h,g
end


function calc_HeG_interp(dad::potencial_iso,b1,b2,npg=8,ninterp=3)
  collocCoord=[dad.NOS;dad.pontos_internos][b1,:]
  xmax=maximum(collocCoord,dims=1)
  xmin=minimum(collocCoord,dims=1)

  xs=criapontosinterp(ninterp)
  fontes,L,ninterp1,ninterp2=gera_interpolação(ninterp,collocCoord,xmax,xmin,xs)

  H=zeros(ninterp1*ninterp2,0)
  G=zeros(ninterp1*ninterp2,0)
  n1,n2=Nlinear(xs)
  xks=n1*xmin+n2*xmax
  
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM[b2]  #Laço dos elementos
    xf=elem_j.limites[:,2]
    x0=elem_j.limites[:,1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
    Δelem=xf-x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
    pc=dad.pontos_controle[:,elem_j.indices]
    # @infiltrate
    cf=pc[1:2,:]./pc[4,:]'
    
    h1 = zeros(ninterp1*ninterp2,size(elem_j));
    g1 = zeros(ninterp1*ninterp2,size(elem_j));
    ci=0
    for i2 =1:ninterp1
      for i1 =1:ninterp2
        ci+=1

        pf = [xks[i1,1],xks[i2,2]]   # Coordenada (x,y)  dos pontos fonte
        eet=2*dot(Δelem,pf.-x0)/norm(Δelem)^2-1
        N,dN=calc_fforma(eet,elem_j,pc[4,:])        
        ps=cf*N
        b=norm(ps-pf)/norm(Δelem)
        eta,Jt=sinhtrans(qsi,eet,b)
        h,g= integrabezier(pf,cf,pc[4,:],eta,w.*Jt,elem_j,dad) 
        h1[ci,:] += h
        g1[ci,:] += g
      end
    end
    H= [H h1]
    G= [G g1]
  end
  L,H,G
end

function calc_HeG(dad::potencial_iso,b1,b2,npg=8)
  n1 = size(b1,1)
  n2 = size(b2,1)
  H=zeros(n1,0)
  G=zeros(n1,0)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM[b2]  #Laço dos elementos
    xf=elem_j.limites[:,2]
    x0=elem_j.limites[:,1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
    Δelem=xf-x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
    pc=dad.pontos_controle[:,elem_j.indices]
    # @infiltrate
    cf=pc[1:2,:]./pc[4,:]'   
    h1 = zeros(n1,size(elem_j));
    g1 = zeros(n1,size(elem_j));
    for i =1:n1
      ind= b1[i]  #Laço dos pontos fontes
      pf = [dad.NOS;dad.pontos_internos][ind,:]   # Coordenada (x,y)  dos pontos fonte
      eet=2*dot(Δelem,pf.-x0)/norm(Δelem)^2-1
        N,dN=calc_fforma(eet,elem_j,pc[4,:])        
        ps=cf*N
        b=norm(ps-pf)/norm(Δelem)
        eta,Jt=sinhtrans(qsi,eet,b)
        h,g= integrabezier(pf,cf,pc[4,:],eta,w.*Jt,elem_j,dad)
      if ind in elem_j.sing
        h1[i,:] += h-0.5*dad.E[ind,elem_j.indices]
      else
        h1[i,:] += h
      end

      g1[i,:] += g
    end
    H= [H h1]
    G= [G g1]
  end
  H,G
end



function calsolfund(r,n,prob::potencial)
  R=norm(r)
  Qast=dot(r,n)/R^2/(2*π)       # Equação 4.36
  Tast=-log(R)/(2*π*prob.k)
  Qast,Tast
end

function calsolfund(r,n,prob::potencial_iso)
  R=norm(r)
  Qast=dot(r,n)/R^2/(2*π)       # Equação 4.36
  Tast=-log(R)/(2*π*prob.k)
  Qast,Tast
end