function calc_HeG(dad::elastico,npg=8)
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
              h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad.Ev[1],dad.Ev[2])
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
function integraelem(pf,x,eta,w,elem,E,v)
  h = zeros(Float64,2,2*size(elem))
  g = zeros(Float64,2,2*size(elem))
  Nm = zeros(Float64,2,2*size(elem))
  for k = 1:size(w,1)
    N,dN=calc_fforma(eta[k],elem)
      pg = N'*x    # Ponto de gauss interpolador
      r = pg'-pf      # Distancia entre ponto de gauss e ponto fonte
      dxdqsi = dN'*x   # dx/dξ & dy/dξ
      dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
      sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
      sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
      uast,tast=calsolfund(r,[sy,-sx],E,v)
      Nm[1,1:2:end]=N
      Nm[2,2:2:end]=N
      h+=tast*Nm*dgamadqsi*w[k]
      g+=uast*Nm*dgamadqsi*w[k]
      # @infiltrate

end
h,g
end

function calsolfund(r,n,E,v)
# @infiltrate
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

function calc_Aeb(dad::elastico,npg=8)
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
              h,g= integraelem(pf,x,eta,w.*Jt,elem_j,dad.Ev[1],dad.Ev[2])
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



