

function interpola(r)
    	if r == 0
        return 0
    end
    r^2 * log(r)
    # r

end
function int_interpolaρdρ(r)
    (4 * r^4 * log(r) - r^4) / 16
    # r^3/3
end


function Monta_M_RIMd(dad,npg)
    n_nos = size(dad.NOS, 1)
    nelem = size(dad.ELEM, 1)
    n_noi = size(dad.pontos_internos,1); #Number of internal nodes

    qsi, w = gausslegendre(npg)
    n_pontos = n_nos + n_noi;
    nodes = [dad.NOS;dad.pontos_internos]
    M = zeros(n_pontos);
    M1 = zeros(n_pontos);
    F,D=FeD(dad)
    M,M1=calcMs(dad,npg)
	# @show size(M)
	# @show length(M)
    A = ones(length(M)) * M' / F .* D
    for i = 1:n_pontos #Laço dos pontos radiais
        A[i,i] = 0
        A[i,i] = -sum(A[i,:])        
    end
    A+diagm(0 => M1)
    M,M1,F,D
end

function  calc_md(x,pf, k, qsi, w, elem)
    npg = length(w);
    m_el, m_el1 = 0, 0

    for i = 1:npg
       N = calc_fforma_gen(qsi[i],elem.ξs)
      pg = N'*x    # Ponto de gauss interpolador
      r = pg'-pf      # Distancia entre ponto de gauss e ponto fonte
      dN_geo = calc_dfforma_gen(qsi[i],elem.ξs) # calcula dN\dξ N1,N2 e N3
      dxdqsi = dN_geo'*x   # dx/dξ & dy/dξ
      dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
      sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
      sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ

     nx = sy; # Componente x do vetor normal unit�rio
     ny = -sx; # Componente y do vetor normal unit�rio
     # @infiltrate
     r = pg' - pf
     R = norm(r)
     m = int_interpolaρdρ(R)
     m1 = -(2 * R^2 * log(R) - R^2) / 4 / (2 * π * k)
    # calcula_Fd(pr, pf, pg, [nx,ny], k, qsi2, w2);
     m_el += dot([nx,ny], r) / norm(r)^2 * m * dgamadqsi * w[i]
     m_el1 += dot([nx,ny], r) / norm(r)^2 * m1 * dgamadqsi * w[i]
    end
    return m_el, m_el1
end
function Finterp(dad,Tree1,Tree2,block;ninterp=3,compressão=true,ϵ=1e-3)
    # arg = [NOS1,NOS_GEO1,tipoCDC,valorCDC,normal,ELEM1,k]
    #         1      2        3      4        5     6   7
    n = size(block,1)               # Quantidade de Submatrizes
    Faca = Array{Any}(undef,n,2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    Daca = Array{Any}(undef,n)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    nodes = [dad.NOS;dad.pontos_internos]
    M,M1=calcMs(dad,8)
    for i=1:n                       # Para cada Submatriz
        # @timeit to "Para cada Submatriz" begin
        b1 = Tree1[block[i,1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i,2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # Submatriz = Produto cartesiano I x J
        cols=vcat([dad.ELEM[ii].indices for ii in b2]...)

        if block[i,3]==0                # Se esses blocos não são admissiveis
            Faca[i,1],Daca[i,1] = FeD(dad,b1,cols)
            # Salva na linha i da 1º coluna da matriz Aaca a matriz H e salva G na matriz B
        else                              # Caso contrario (Se blocos são admissiveis)
            Faca[i,1],Faca[i,2],Daca[i] =FeDinterp(dad,b1,cols,8,ninterp)    
        end
        # end
        
    end
    return Faca,Daca
end

function FeD(dad,b1=0,b2=0)
    nodes = [dad.NOS;dad.pontos_internos]
    if b1==0
        n_pontos=size(nodes,1)
        b1=1:n_pontos
        b2=1:n_pontos    
    end
    n1,n2=size(b1,1),size(b2,1)
    F=zeros(n1,n2)
    D=zeros(n1,n2)
    
            for i =1:n1
            xi = nodes[b1[i],1];
            yi = nodes[b1[i],2];
            for j =1:n2
                if b1[i]==b2[j]
                    continue
                end
                xj = nodes[b2[j],1];
                yj = nodes[b2[j],2];
                r = sqrt((xi - xj)^2 + (yi - yj)^2);
                F[i,j] = interpola(r);      
                D[i,j] = -log(r) / (2 * π * dad.k)
            end
        end
        F,D
    end
    

    function calcMs(dad,npg)
        nodes = [dad.NOS;dad.pontos_internos]
        n_pontos=size(nodes,1)
        M = zeros(n_pontos);
        M1 = zeros(n_pontos);
        qsi, w = gausslegendre(npg)

        for i = 1:n_pontos #Laço dos pontos radiais
            pf = nodes[i,:]            
              for elem_j in dad.ELEM  #Laço dos elementos
                      x = dad.NOS[elem_j.indices,:]   # Coordenada (x,y) dos nós geométricos
                    m_el, m_el1 = calc_md(x,pf, dad.k, qsi, w, elem_j)
                    M[i] = M[i] + m_el
                    M1[i] = M1[i] + m_el1
                 end
        end
        M,M1
    end

    function FeDinterp(dad::potencial,b1,b2,npg=8,ninterp=3)
        nodes=[dad.NOS;dad.pontos_internos][b1,:]
        
        nodes2=[dad.NOS;dad.pontos_internos][b2,:]
        xmax=maximum(nodes,dims=1)
        xmin=minimum(nodes,dims=1)
      
        xs=criapontosinterp(ninterp)
        fontes,L,ninterp1,ninterp2=gera_interpolação(ninterp,nodes,xmax,xmin,xs)
      
        F=zeros(ninterp1*ninterp2,0)
        D=zeros(ninterp1*ninterp2,0)
        n1,n2=Nlinear(xs)
        xks=n1*xmin+n2*xmax
        nb2=size(b2,1)

        qsi,w = gausslegendre(npg)    # Quadratura de gauss
        for j =1:nb2

          f1 = zeros(ninterp1*ninterp2);
          d1 = zeros(ninterp1*ninterp2);
          ci=0
          for i2 =1:ninterp1
            for i1 =1:ninterp2
              ci+=1
              pf = [xks[i1,1],xks[i2,2]]   # Coordenada (x,y)  dos pontos fonte
              r = norm(pf-nodes2[j,:])
              f = interpola(r);      
              d = -log(r) / (2 * π * dad.k)
      
              f1[ci] += f
              d1[ci] += d
            end
          end
          F= [F f1]
          D= [D d1]
        end
        L,F,D
      end
