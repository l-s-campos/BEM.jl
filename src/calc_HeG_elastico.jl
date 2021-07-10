function calc_HeG(dad::Union{elastico, elastico_aniso},npg=8)
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS,1)
  H=zeros(2*n,2*n)
  G=zeros(2*n,2*n)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
      for i=1:n
            pf = dad.NOS[i,:]   # Coordenada (x,y)  dos pontos fonte
            for elem_j in dad.ELEM  #Laço dos elementos
              x = dad.NOS[elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
              Δelem=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
              eet=(elem_j.ξs[end] -elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
              N_geo=calc_fforma(eet,elem_j,false)
              ps=N_geo'*x
              b=norm(ps'-pf)/norm(Δelem)
              eta,Jt=sinhtrans(qsi,eet,b)
              # eta,Jt=telles(qsi,eet)
              h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad)
              cols=[2elem_j.indices.-1 2elem_j.indices]'[:]
              H[2i-1:2i,cols]=h
              G[2i-1:2i,cols]=g
          end
      end
  
  for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    H[2i-1:2i,2i-1:2i].=0
    H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
  end
H,G
end

#__________________________________________________________________________________________________________
"Funcao para calcular fazer integracao no contorno "
function integraelem(pf,x,eta,w,elem,dad::Union{elastico, elastico_aniso})
  h = zeros(Float64,2,2*size(elem))
  g = zeros(Float64,2,2*size(elem))
  Nm = zeros(Float64,2,2*size(elem))
  for k = 1:size(w,1)
    N,dN=calc_fforma(eta[k],elem)
      pg = N'*x    # Ponto de gauss interpolador
      dxdqsi = dN'*x   # dx/dξ & dy/dξ
      dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
      sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
      sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
      uast,tast=calsolfund(pg,pf,[sy,-sx],dad)
      Nm[1,1:2:end]=N
      Nm[2,2:2:end]=N
      h+=tast*Nm*dgamadqsi*w[k]
      g+=uast*Nm*dgamadqsi*w[k]
      # @infiltrate

end
h,g
end



function calc_Aeb(dad::Union{elastico, elastico_aniso},npg=8)
  nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS,1)
  A=zeros(2*n,2*n)
  B=zeros(2*n)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
      for i=1:n
        pf = [dad.NOS;dad.pontos_internos][i,:]   # Coordenada (x,y)  dos pontos fonte
        for elem_j in dad.ELEM  #Laço dos elementos
              x = dad.NOS[elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
              Δelem=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
              eet=(elem_j.ξs[end] -elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
              N=calc_fforma(eet,elem_j,false)
              ps=N_geo'*x
              b=norm(ps'-pf)/norm(Δelem)
              eta,Jt=sinhtrans(qsi,eet,b)
              # eta,Jt=telles(qsi,eet)
              h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad)
              # falta a integral singular
              for ii =1:2
              if elem_j.tipoCDC[ii]==1
                A[2i-1:2i,2elem_j.indices.-(2-ii)]=h[:,ii:2:end]
                B[2i-1:2i]=g[:,ii:2:end]*elem_j.valorCDC[ii,:]
              else
                  A[2i-1:2i,2elem_j.indices]=-g[:,ii:2:end]
                  B[2i-1:2i]+=-h[:,ii:2:end]*elem_j.valorCDC[ii,:]
              end
              end
          end
      end
A,B
end

function calc_HeG(dad::elastico_iso,npg=8)
  # nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS,1)
  H=zeros(2*n,2*n)
  G=zeros(2*n,2*n)
  qsi,w = gausslegendre(npg)    # Quadratura de gauss
  for elem_j in dad.ELEM  #Laço dos elementos
    xf=elem_j.limites[:,2]
    x0=elem_j.limites[:,1]     # Δx e Δy entre o primeiro e ultimo nó geometrico
    Δelem=xf-x0     # Δx e Δy entre o primeiro e ultimo nó geometrico
    pc=dad.pontos_controle[:,elem_j.indices]
    # @infiltrate
    cf=pc[1:2,:]./pc[4,:]'
    for i =1 : n
      pf = dad.NOS[i,:]   # Coordenada (x,y)  dos pontos fonte
      eet=2*dot(Δelem,pf.-x0)/norm(Δelem)^2-1
        N,dN=calc_fforma(eet,elem_j,pc[4,:])        
        ps=cf*N
        b=norm(ps-pf)/norm(Δelem)
        eta,Jt=sinhtrans(qsi,eet,b)
        h,g= integrabezier(pf,cf,pc[4,:],eta,w.*Jt,elem_j,dad)
              if i in elem_j.sing
                h=integrabeziersing(pf,cf,pc[4,:],eta,w.*Jt,elem_j,dad,eet)
              end
              cols=[2elem_j.indices.-1 2elem_j.indices]'[:]
              H[2i-1:2i,cols]=h
              G[2i-1:2i,cols]=g
          end
      end
  

  H[1:2:end,1:2:end]+=dad.E/2
  H[2:2:end,2:2:end]+=dad.E/2
H,G
end


function integrabeziersing(pf,cf,we,eta,w,elem::bezier,dad::elastico_iso,eet)
  h = zeros(Float64,2,2*size(elem))
  basisrc,dN=calc_fforma(eet,elem,we)       
  dxdqsi = cf*dN   # dx/dξ & dy/dξ
  dgamadqsif = norm(dxdqsi)/2  # dΓ/dξ = J(ξ) Jacobiano    
  matbasisrc=zeros(2,2*size(elem))
  matbasisrc[1,1:2:end]=basisrc
  matbasisrc[2,2:2:end]=basisrc

  E,ν=dad.Ev[1],dad.Ev[2]


  hterm=[0	-(1-2*ν)/(4*pi*(1-ν))
         (1-2*ν)/(4*pi*(1-ν)) 	0]
  htermMatrix=hterm*matbasisrc

  Nm = zeros(Float64,2,2*size(elem))


  for k = 1:size(w,1)
    N,dN = calc_fforma(eta[k],elem,we)
    pg = cf*N    # Ponto de gauss interpolador
    dxdqsi = cf*dN   # dx/dξ & dy/dξ
    dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
    sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
    sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
    uast,tast=calsolfund(pg,pf,[sy,-sx],dad)
    # h+=N*dgamadqsi*w[k]
    # g+=N*dgamadqsi*w[k]
    Nm[1,1:2:end]=N
    Nm[2,2:2:end]=N
    h+=(tast*Nm*dgamadqsi/2-htermMatrix/(eta[k]-eet))*w[k]
    # @infiltrate
  end
# @show h
  if abs(eet)==1
    beta_m=1/dgamadqsif
    h+=htermMatrix*log(abs(2/beta_m))*sign(-eet)
#        println("h = $(htermMatrix*log(abs(2/beta_m))*sign(-eet))")
else
    h+=htermMatrix*log(abs((1-eet)/(1+eet)));
end
  h
end



function integrabezier(pf,cf,we,eta,w,elem::bezier,dad::elastico_iso)
  h = zeros(Float64,2,2*size(elem))
  g = zeros(Float64,2,2*size(elem))
  Nm = zeros(Float64,2,2*size(elem))

  for k = 1:size(w,1)
    N,dN = calc_fforma(eta[k],elem,we)
    pg = cf*N    # Ponto de gauss interpolador
    r = pg-pf      # Distancia entre ponto de gauss e ponto fonte
    dxdqsi = cf*dN   # dx/dξ & dy/dξ
    dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
    sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
    sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
    uast,tast=calsolfund(r,[sy,-sx],dad)
    Nm[1,1:2:end]=N
    Nm[2,2:2:end]=N
    h+=tast*Nm*dgamadqsi*w[k]/2
    g+=uast*Nm*dgamadqsi*w[k]/2

  end
  h,g
end


function calsolfund(pg,pf,n,dad::Union{elastico, elastico_iso})
  # @infiltrate
  E,v=dad.Ev[1],dad.Ev[2]
  r = pg'-pf      # Distancia entre ponto de gauss e ponto fonte
  
    GE=E/(2*(1+v));
    # Distance of source and field points
    R = norm(r);
  
    # Components of the unity vector in the radial direction
    r1 = r/R
  
    # Plane elasticity fundamental solutions
    prod1 = 4*pi*(1-v);
    prod2 = (3-4*v)*log(1/R);
    prod3 = dot(r1,n)
  
    u11 = (prod2 + r1[1]^2)/(2*prod1*GE);
    u22 = (prod2 + r1[2]^2)/(2*prod1*GE);
    u12 = (r1[1]*r1[2])/(2*prod1*GE);
    u21=u12;
  
    t11 = -(prod3*((1-2*v)+2*r1[1]^2))/(prod1*R);
    t22 = -(prod3*((1-2*v)+2*r1[2]^2))/(prod1*R);
    t12 = -((prod3*2*r1[1]*r1[2])-(1-2*v)*(r1[1]*n[2]-r1[2]*n[1]))/(prod1*R)
    t21 = -((prod3*2*r1[1]*r1[2])-(1-2*v)*(r1[2]*n[1]-r1[1]*n[2]))/(prod1*R)
        uast = [u11 u12
           u21 u22];
  
         tast = [t11   t12
           t21   t22];
  
  uast,tast
  end


function calsolfund(pg,pf,n,dad::elastico_aniso)
  # @infiltrate
  mi=dad.k[1,:]
  A=dad.k[2:3,:]
  q=dad.k[4:5,:]
  g=dad.k[6:7,:]

  xcampo=pg[1]
  ycampo=pg[2]
  xf=pf[1]
  yf=pf[2]


    #Cálculo da distância do ponto fonte (xf,yf) ao ponto campo
    z1 = xcampo - xf+mi[1]*(ycampo - yf);
    z2 = xcampo - xf+mi[2]*(ycampo - yf);
    
    # Solução fundamental de deslocamento
    
    lns=[log(z1)     0
            0  log(z2)];
    
    
    uast = 2*real(A*lns*conj(q)');
    
    # Solução fundamental de forcas de superficies
    
    mi_n_z=[(mi[1]*n[1]-n[2])/z1         0
                   0          (mi[2]*n[1]-n[2])/z2];
    
    tast = 2*real(A*mi_n_z*conj(g)');
  
  uast,tast
  end
