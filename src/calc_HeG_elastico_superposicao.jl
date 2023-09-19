function calc_HeG(dad::Vector{elastico}, npg=8)
  nc = size(dad, 1)    # Quantidade de elementos discretizados no contorno

  n1 = size(dad[1].NOS, 1)
  nt = 0
  for d in dad
    nt += size(d.NOS, 1)
  end
  n1 = size(dad[1].NOS, 1)

  H = zeros(2 * nt, 2 * nt)
  G = zeros(2 * nt, 2 * nt)

  qsi, w = gausslegendre(npg)    # Quadratura de gauss
  iold = 0
  for idad = 1:nc
    n = size(dad[idad].NOS, 1)
    for i = 1:n
      pf = dad[idad].NOS[i, :]   # Coordenada (x,y)  dos pontos fonte
      jold = 0
      for jdad = 1:nc

        for elem_j in dad[jdad].ELEM  #Laço dos elementos
          x = dad[jdad].NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
          Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
          eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
          N_geo = calc_fforma(eet, elem_j, false)
          ps = N_geo' * x
          b = norm(ps' - pf) / norm(Δelem)
          eta, Jt = sinhtrans(qsi, eet, b)
          # @infiltrate
          # eta,Jt=telles(qsi,eet)
          h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad[jdad])
          cols = [2elem_j.indices .- 1 2elem_j.indices]'[:]
          if jdad == 1
            mult = 1
          else
            mult = dad[jdad].k.E / dad[1].k.E
          end
          # @infiltrate jdad == 2
          H[(2i-1:2i).+iold, cols.+jold] = mult * h
          G[(2i-1:2i).+iold, cols.+jold] = g
        end
        jold += 2 * size(dad[jdad].NOS, 1)
      end
    end
    iold += 2n
  end

  for i = 1:n1                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    H[2i-1:2i, 2i-1:2i] .= 0
    H[2i-1:2i, 2i-1:2i] = -[sum(H[2i-1:2i, 1:2:2n1], dims=2) sum(H[2i-1:2i, 2:2:2n1], dims=2)]
  end
  for i = n1+1:nt                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    H[2i-1:2i, 2i-1:2i] .= 0
    H[2i-1:2i, 2i-1:2i] = -[sum(H[2i-1:2i, 1:2:end], dims=2) sum(H[2i-1:2i, 2:2:end], dims=2)]
  end
  # for i = 1:n1                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
  #   H[2i-1:2i, 2i-1:2i] .= 0
  #   H[2i-1:2i, 2i-1:2i] = -[sum(H[2i-1:2i, 1:2:2n1], dims=2) sum(H[2i-1:2i, 2:2:2n1], dims=2)]
  # end
  # for i = n1+1:nt                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
  #   H[2i-1:2i, 2i-1:2i] .= 0
  #   H[2i-1:2i, 2i-1:2i] = -[sum(H[2i-1:2i, 2n1+1:2:end], dims=2) sum(H[2i-1:2i, 2n1+2:2:end], dims=2)]
  #   H[2i-1, 2i-1] = 1
  #   H[2i, 2i] = 1
  # end
  H, G
end



function aplicaCDC(H, G, dad::Vector{elastico})
  nc = size(dad, 1)    # Quantidade de elementos discretizados no contorno
  n = 0
  for d in dad
    n += size(d.NOS, 1)
  end
  A = zeros(2 * n, 2 * n)
  b = zeros(2 * n)
  jold = 0
  for jdad = 1:nc
    for elem_i in dad[jdad].ELEM, i in 1:2  # Laço dos pontos fontes
      ind_elem = elem_i.indices
      # @infiltrate jdad == 2
      if elem_i.tipoCDC[i] == 0
        A[:, 2ind_elem.+(i-2).+jold] = -G[:, 2ind_elem.+(i-2).+jold]
        b += -H[:, 2ind_elem.+(i-2).+jold] * elem_i.valorCDC[i, :]
      elseif elem_i.tipoCDC[i] == 1
        A[:, 2ind_elem.+(i-2).+jold] = H[:, 2ind_elem.+(i-2).+jold]
        b += G[:, 2ind_elem.+(i-2).+jold] * elem_i.valorCDC[i, :]
      end
    end
    jold += 2 * size(dad[jdad].NOS, 1)
  end

  A, b
end


function separa(dad::Vector{elastico}, x)
  # Separa fluxo e temperatura

  # ncdc = número de linhas da matriz CDC
  # T = vetor que contêm as temperaturas nos nós
  # q = vetor que contêm o fluxo nos nós
  nc = size(dad, 1)    # Quantidade de elementos discretizados no contorno
  n = 0
  for d in dad
    n += size(d.NOS, 1)
  end
  T = zeros(2n)
  q = zeros(2n)
  jold = 0
  for jdad = 1:nc

    for elem_i in dad[jdad].ELEM, i in 1:2   # Laço dos pontos fontes
      ind_elem = elem_i.indices
      if elem_i.tipoCDC[i] == 0
        T[2ind_elem.+(i-2).+jold] = elem_i.valorCDC[i, :] # A temperatura é a condição de contorno
        q[2ind_elem.+(i-2).+jold] = x[2ind_elem.+(i-2).+jold] # O fluxo é o valor calculado
      elseif elem_i.tipoCDC[i] == 1
        T[2ind_elem.+(i-2).+jold] = x[2ind_elem.+(i-2).+jold] # 
        q[2ind_elem.+(i-2).+jold] = elem_i.valorCDC[i, :]
      end
    end
    jold += 2 * size(dad[jdad].NOS, 1)
  end

  [T[1:2:end] T[2:2:end]], [q[1:2:end] q[2:2:end]]
end