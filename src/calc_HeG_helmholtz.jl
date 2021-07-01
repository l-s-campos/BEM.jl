function calc_HeG(dad::escalar,npg=8)
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
        eta,wt=pontosintegra(x,elem_j,qsi,w)
        h,g= integraelem(pf,x,eta,wt,elem_j,dad.k,typeof(dad))
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
function integraelem(pf,x,eta,w,elem,kmat,prob::escalar)
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
    Qast,Tast=calsolfund(r,[sy,-sx],kmat,prob)
    h+=N*Qast*dgamadqsi*w[k]
    g+=N*Tast*dgamadqsi*w[k]    
  end
  h,g
end

function calsolfund(r,n,kmat,prob::helmholtz)
  FR,CW,GE=kmat
  R=norm(r)
  ZR = real(FR*R/CW);
  Z  = complex(0.,ZR);
  Tast = besselk(0,Z)/(2*π*GE)
  Qast = -im*dot(r,n)/R*(FR/(2*π*CW))*besselk(1,Z) # ERRAAAADO!!!!
  #  @show dot(r,n)  
  return Qast,Tast
  end
