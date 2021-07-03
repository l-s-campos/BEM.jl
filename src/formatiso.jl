
function calc_ncont(SEGMENTOS)
  # conta o número de contornos e o número de segmentos em cada contorno e 
  #     cada linha inicia com o primeiro segmento daquele contorno
  num_seg=size(SEGMENTOS,1);
  t=1;
  p2=0;
  p_ini = SEGMENTOS[1,2];
  icont=1; # índice de contornos
  contorno=zeros(Int32,1,2)
  contorno[1,1]=p_ini;
  while(t<num_seg)    # While over all lines
      while(p2!=p_ini)
          p1  = SEGMENTOS[t,2]
          p2  = SEGMENTOS[t,3]
          if p2 == p_ini
              if t < num_seg
                  p_ini = SEGMENTOS[t+1,2]
                  icont=icont+1;
                  contorno=[contorno; zeros(Int32,1,2)]
                  contorno[icont,1]=p_ini;
                  contorno[icont-1,2]=p_ini-contorno[icont-1,1]
              end;
          end;
          t=t+1
      end                             #end of while p2
  end
  if(icont>1)
      contorno[icont,2]=num_seg-sum(contorno[1:icont-1,2]);
  else
      contorno[1,2]=num_seg;
  end
  return contorno
  end

# function format_dad_iso(PONTOS,SEGMENTOS) 
  function  format_dad_iso(entrada,NPX=2,NPY=2)
    prob,PONTOS,SEGMENTOS,MALHA,CCSeg,k=entrada

#   println("PONTOS = $PONTOS")
#   println("SEGMENTOS = $SEGMENTOS")
  contorno=calc_ncont(SEGMENTOS);
#   println(" contorno = $contorno")
  ncont=size(contorno,1);   # num de contornos
#   @show ncont      # <------------                 emerinserted           emerinserted
#   @show contorno   # <------------                 emerinserted           emerinserted
  
#   println(" ncont = $ncont")
  icrv=0;
#   ncurves=Int64(0);
  ncurves=sum(contorno[:,2])  # num de segmentos ou de curvas
#   println("ncurves = $ncurves")
#   crv=Array{Curve,ncurves} #(ncurves)  # EmerEdited
  crv=Array{Curve}(undef,ncurves)    # EmerEdited
#   jj=0;
  for j=1:ncont  # looping pelos contornos um por um
      pini=contorno[j,1];  # Ponto onde começa o contorno j
      nseg=contorno[j,2];  # num de segmentos do contorno j
      for i=1:nseg   # looping por cada segmento do contorno j
          icrv=icrv+1;
          raio=SEGMENTOS[i+pini-1,4] #define valores para o raio
          np1=Int32(SEGMENTOS[i+pini-1,2]); # Define o primeir o ponto do segmento
          np2=Int32(SEGMENTOS[i+pini-1,3]); #Define o segundo/último ponto do segmento
          p1=[PONTOS[np1,2], PONTOS[np1,3],0]; #Define o primeiro ponto da curva
          p2=[PONTOS[np2,2], PONTOS[np2,3],0]; #Define o segundo/último ponto da curva/segmento
            
          if(raio==0)
              crv[icrv]=nrbline(p1,p2);
          else
              xc,yc=calcula_centro(p1[1],p1[2],p2[1],p2[2],raio);
  #            @show i,p1,p2,xc,yc,raio    #  <-----         emerinserted      emerinserted
              sang,eang = calcula_arco(p1[1],p1[2],p2[1],p2[2],xc,yc,raio);
              centro=[xc,yc,0];

 #             @show sang,eang    #  <-----         emerinserted      emerinserted
          #    println(raio)
          #    println(centro)
          #    println(sang)
          #    println(eang)
              crv[icrv]=nrbcirc(abs(raio),centro,sang,eang); 
          
          end
      end
  end
  ncrv=size(crv,1)

      for i=1:ncrv
        degree=crv[i].order-1
        if MALHA[i,3]>degree
          coefs,knots = bspdegelev(degree,crv[i].coefs,crv[i].knots,MALHA[i,3])
          crv[i] = nrbmak(coefs,knots)
      end
  end
  
  for i=1:ncrv
    
        h=MALHA[i,2]
        if h>0
        novosnos=range(0,stop=1,length=h+2)
      degree=crv[i].order-1
      coefs,knots = bspkntins(degree,crv[i].coefs,crv[i].knots,novosnos[2:end-1])
      crv[i] = nrbmak(coefs,knots)
      end
  end
  NOS = zeros(Float64,0,2)
  indfonte,indcoluna,indbezier,tipoCDC,valorCDC,E,NOS,collocPts,sing=indices(crv,CCSeg)
  pts_controle=hcat([crv[i].coefs for i =1:ncrv]...)
  cont_el=0
  num_elementos=size(indcoluna,1)
  ELEM = Vector{bezier}(undef,num_elementos)
# @infiltrate
for i = 1:ncrv # Corre as NURBS para variar o elemento
for j = 1:size(crv[i].conn, 1)  # Corre as curvas de Bézier
    cont_el+=1
    p=crv[i].order-1
    pts_controle_e = pts_controle[1:2,indcoluna[cont_el]]   # Coordenada (x,y) dos pontos de pontos_controle
    we = pts_controle[4,indcoluna[cont_el]]   # Coordenada (x,y) dos pontos de pontos_controle
    Ce=crv[i].C[:,:,j]
    Be=Bernsteins(p, -1)
    Wb=Ce'*we # Bézier weights (nodal values)
    # @infiltrate
    wb = dot(Be,Wb);   #  Bézier weight function
    # Shape function and derivatives
    R = diagm(we)*Ce*Be/wb    

    x0=(pts_controle_e./ [ we'; we'])*R # Calcula as coordenadas xy dos pontos de integração
    Be=Bernsteins(p, 1)
    # @infiltrate
    wb = dot(Be,Wb);   #  Bézier weight function
    # Shape function and derivatives
    R = diagm(we)*Ce*Be/wb    
    xf=(pts_controle_e./ [ we'; we'])*R # Calcula as coordenadas xy dos pontos de integração
    # @infiltrate
    ELEM[cont_el] = bezier(indcoluna[cont_el],Ce,p,[x0 xf],Wb,sing[cont_el])
end
end
# @infiltrate
prob(NOS,pts_controle,gera_p_in(NPX,NPY,PONTOS,SEGMENTOS),tipoCDC,valorCDC,ELEM,k,E)


end


function bezierExtraction(knot, p)
  uk = unique(knot)
  nb1 = length(uk) - 1
  conn = Vector(undef, nb1)
  conn[1] = 1:p + 1
  for i = 2:nb1
      conn[i] = conn[i - 1] .+ sum(knot .== uk[i])
  end
  m  = length(knot) - p - 1;
  a  = p + 1;
  b  = a + 1;
  nb1 = length(uk) - 1
  C = zeros(p + 1, p + 1, nb1)
  C[:,:,1] = Matrix(1.0I, p + 1, p + 1)
  nb = 1
  while b <= m
      C[:,:,nb + 1] = Matrix(1.0I, p + 1, p + 1)
      i = b;
      while b <= m && knot[b + 1] == knot[b]
          b = b + 1;
      end
  
      multiplicity = b - i + 1;
      if multiplicity < p
          numerator = knot[b] - knot[a];
          alphas = zeros(p)
          for j = p:-1:multiplicity + 1
              alphas[j - multiplicity] = numerator / (knot[a + j] - knot[a]);
          end
          r = p - multiplicity;
          for j = 1:r
              save = r - j + 1;
              s = multiplicity + j;
              for k = p + 1:-1:s + 1
                  alpha = alphas[k - s];
                  C[:,k,nb] = alpha * C[:,k,nb] + (1 - alpha) * C[:,k - 1,nb];
              end
              if b <= m
                  C[save:save + j,save,nb + 1] = C[p - j + 1:p + 1,p + 1,nb]
              end
          end
          nb = nb + 1;
          if b <= m
              a = b;
              b = b + 1;
          end
      elseif multiplicity == p
          if b <= m
                              nb = nb + 1; a = b; b = b + 1;
          end
      end
  end
  C, nb, conn, uk
end


function indices(crv,CCSeg)
    n = length(crv);    # Number of curves

    z=0;#ncollocpoints
    for k=1:n
        for i=1:crv[k].number
            z=z+1
        end
    end

    numcurva=zeros(Integer,z)
    collocPts=zeros(z)
    collocCoord=zeros(z,2)
    nnos=zeros(Integer,n)
    nbezier=zeros(Integer,n)

    if size(CCSeg,2)==3
    vCDC=zeros(0,1)
    tCDC=zeros(Int64,0)
elseif size(CCSeg,2)==5
    vCDC=zeros(0,2)
    tCDC=zeros(Int64,0,2)
    tCDC=zeros(Int64,0)
elseif size(CCSeg,2)==7
    vCDC=zeros(0,3)
    tCDC=zeros(Int64,0,3)
end

    for k=1:n
        p=crv[k].order-1;
        nnos[k]=crv[k].number;
        nbezier[k]=size(crv[k].conn,1)
        valorCDC=CCSeg[k,3:2:end];
        tipoCDC=round(Int,CCSeg[k,2:2:end]);
        for i=1:crv[k].number
            vCDC= [vCDC;valorCDC];
            tCDC = [tCDC;tipoCDC];
        end
    end

    nnos2=cumsum([0 nnos'],dims=2);
    nbezier2=cumsum([0 nbezier'],dims=2);


    indfonte = Array{Int}(undef,z, 2);
    indbezier = Array{Int}(undef,0, 2);
    indcoluna = Array{Any}(undef,0,1)


    for i = 1:n
        for j = 1:size(crv[i].conn, 1)
            indbezier=[indbezier;i j]
            indcoluna=[indcoluna
                       [(crv[i].conn[j] .+ nnos2[i])]]
        end
    end
    k = 1
    for k1 = 1:n
        for k2 = 1:crv[k1].number
            indfonte[k,:]=[k1 k2]
            k += 1
        end
    end




    E = spzeros(z,z);

    k = 1
    for k1 = 1:n
        uu = unique(crv[k1].knots)
        for k2 = 1:crv[k1].number
            ind = sum(crv[k1].fontes[k2].pts .> uu)
            E[k,crv[k1].conn[ind] .+ nnos2[k1]] += crv[k1].fontes[k2].basis
            k += 1
        end
    end
    cont=1
    sing=zeros(Int64,z)
ci=1
    for c in crv
        for f in c.fontes
            collocCoord[cont,:]=f.coords[1:2]
            collocPts[cont]=f.pts
            test=(collocPts[cont].>=crv[ci].range[:,1]) .& (collocPts[cont].<crv[ci].range[:,2])
            # @infiltrate
            sing[cont]=(1:size(crv[ci].conn,1))[test][1]+nbezier2[ci]
            cont+=1
        end
        ci+=1
    end
    sing2=[findall(sing.==i) for i =1:nbezier2[ci]]

    indfonte,indcoluna,indbezier,tCDC,vCDC,E,collocCoord,collocPts,sing2
end
