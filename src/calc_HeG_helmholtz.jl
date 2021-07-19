function calsolfund(r,n,dad::helmholtz)
  R=norm(r)
  ZR = real(dad.k.FR*R/dad.k.CW);
  Z  = complex(0.,ZR);
  Tast = besselk(0,Z)/(2*π*dad.k.GE)
  # Qast = -dot(r,n)*Z/R/(2*π)*besselk(1,Z) # 
  Qast = -im*dot(r,n)/R*(dad.k.FR/(2*π*dad.k.CW))*besselk(1,Z)
  return Qast,Tast
  end

  function calc_HeG(dad::helmholtz,npg=8)
    nelem = size(dad.ELEM,1)    # Quantidade de elementos discretizados no contorno
    n = size(dad.NOS,1)
    H=zeros(Complex,n,n);
    G=zeros(Complex,n,n);
    qsi,w = gausslegendre(npg)    # Quadratura de gauss
    # contafonte=1
    for i in 1:n  #Laço dos pontos fontes
      # for ind_elem = elem_i.indices
        pf = dad.NOS[i,:]   # Coordenada (x,y)  dos pontos fonte
        for elem_j in dad.ELEM#[4:4]  #Laço dos elementos
          x = dad.NOS[elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
          Δelem=x[end,:]-x[1,:]     # Δx e Δy entre o primeiro e ultimo nó geometrico
          eet=(elem_j.ξs[end] -elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
          N_geo,~=calc_fforma(eet,elem_j)
          ps=N_geo'*x
          b=norm(ps'-pf)/norm(Δelem)
          eta,Jt=sinhtrans(qsi,eet,b)
          # eta,Jt=telles(qsi,eet)
          # eta,Jt=qsi,w*0 .+1
          # @show eet
          # @infiltrate
          h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad)
          # @show h
          H[i,elem_j.indices]=h
          G[i,elem_j.indices]=g
          
        end
        # contafonte+=1
      # end
    end
    for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
      H[i,i]+=0.5
    end
    H,G
  end
  