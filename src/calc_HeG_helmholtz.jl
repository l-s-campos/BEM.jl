function calsolfund(r, n, dad::helmholtz)
  R = norm(r)
  ZR = dad.k.FR * R / dad.k.CW

  # Tast = besselk(0, Z) / (2 * π * dad.k.GE)
  # Qast = -im * dot(r, n) / R * (dad.k.FR / (2 * π * dad.k.CW)) * besselk(1, Z)
  Tast = im / 4 * hankelh1(0, ZR)
  Qast = -dad.k.FR / dad.k.CW * im / 4 * hankelh1(1, ZR) * dot(r, n) / R

  return Qast, Tast
end
function calsolfund_hiper(r, n, nf, dad::helmholtz)
  R = norm(r)
  ZR = dad.k.FR * R / dad.k.CW

  Tast = dad.k.FR / dad.k.CW * im / 4 * hankelh1(1, ZR) * dot(r, nf) / R
  Qast = dad.k.FR / dad.k.CW * im / 4 / R * hankelh1(1, ZR) * dot(nf, n) - (dad.k.FR / dad.k.CW)^2 * im / 4 * hankelh1(2, ZR) * dot(r, n) / R * dot(r, nf) / R

  return Qast, Tast
end

function calc_HeG(dad::helmholtz, npg=8)
  nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS, 1)
  H = zeros(Complex, n, n)
  G = zeros(Complex, n, n)
  qsi, w = gausslegendre(npg)    # Quadratura de gauss
  # contafonte=1
  for i in 1:n  #Laço dos pontos fontes
    # for ind_elem = elem_i.indices
    pf = dad.NOS[i, :]   # Coordenada (x,y)  dos pontos fonte
    for elem_j in dad.ELEM#[4:4]  #Laço dos elementos
      x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
      Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
      eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
      N_geo, ~ = calc_fforma(eet, elem_j)
      ps = N_geo' * x
      b = norm(ps' - pf) / norm(Δelem)
      eta, Jt = sinhtrans(qsi, eet, b)
      # eta,Jt=telles(qsi,eet,b)
      # eta,Jt=qsi,w*0 .+1
      # @show i, eet
      # @infiltrate
      h, g = integraelem(pf, x, eta, w .* Jt, elem_j, dad)
      # @show h
      H[i, elem_j.indices] = h
      G[i, elem_j.indices] = g

    end
    # contafonte+=1
    # end
  end
  for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    H[i, i] += 0.5
  end
  H, G
end
function calc_HeG_hiper(dad::helmholtz, npg=8)
  nelem = size(dad.ELEM, 1)    # Quantidade de elementos discretizados no contorno
  n = size(dad.NOS, 1)
  H = zeros(Complex, n, n)
  G = zeros(Complex, n, n)
  qsi, w = gausslegendre(npg)    # Quadratura de gauss
  for i in 1:n  #Laço dos pontos fontes
    # for ind_elem = elem_i.indices
    pf = dad.NOS[i, :]   # Coordenada (x,y)  dos pontos fonte
    for elem_j in dad.ELEM#[4:4]  #Laço dos elementos
      x = dad.NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos

      nosing = elem_j.indices .== i
      if sum(nosing) == 1
        no_pf = findfirst(nosing)
        xi0 = elem_j.ξs[no_pf]
        # h, g = integraelemsing(pf, nf, x, qsi2, w2, elem_j, dad, xi0)
        # hn, gn = integraelemsing_num(pf, nf, x, elem_j, dad, pre, xi0, 30)
        h, g = integraelem_hiper_sing(pf, dad.normal[i, :], x, xi0, elem_j, dad)
      else
        Δelem = x[end, :] - x[1, :]     # Δx e Δy entre o primeiro e ultimo nó geometrico
        eet = (elem_j.ξs[end] - elem_j.ξs[1]) * dot(Δelem, pf .- x[1, :]) / norm(Δelem)^2 + elem_j.ξs[1]
        N_geo, ~ = calc_fforma(eet, elem_j)
        ps = N_geo' * x
        b = norm(ps' - pf) / norm(Δelem)
        eta, Jt = sinhtrans(qsi, eet, b)
        # eta,Jt=telles(qsi,eet,b)
        # eta,Jt=qsi,w*0 .+1
        # @show i, eet
        # @infiltrate
        h, g = integraelem_hiper(pf, dad.normal[i, :], x, eta, w .* Jt, elem_j, dad)
      end
      # @show h
      H[i, elem_j.indices] = h
      G[i, elem_j.indices] = g

    end
    # contafonte+=1
    # end
  end
  for i = 1:n                              #i=1:size(dad.NOS,1) #Laço dos pontos fontes
    G[i, i] += -0.5
  end
  # corrigediags!(H, G, dad, false)
  H, G
end


function corrigediag!(H, G, dad, interno=true)
  xmin = minimum(dad.NOS[:, 1])
  xmax = maximum(dad.NOS[:, 1])
  ymin = minimum(dad.NOS[:, 2])
  ymax = maximum(dad.NOS[:, 2])

  l = sqrt((xmax - xmin)^2 + (ymax - ymin)^2)# Largura do ret�ngulo que cont�m a geometria

  p1 = [xmin, ymin] .- 10l
  p2 = [xmax, ymax] .+ 10l

  tsf1 = zeros(eltype(H), nc(dad))
  qsf1 = zeros(eltype(H), nc(dad))
  tsf2 = zeros(eltype(H), nc(dad))
  qsf2 = zeros(eltype(H), nc(dad))
  for i = 1:nc(dad)
    qsf1[i], tsf1[i] = BEM.calsolfund(dad.NOS[i, :] - p1, dad.normal[i, :], dad)
    qsf2[i], tsf2[i] = BEM.calsolfund(dad.NOS[i, :] - p2, dad.normal[i, :], dad)
    H[i, i] = 0
    G[i, i] = 0
  end
  v1 = H * tsf1 - G * qsf1
  v2 = H * tsf2 - G * qsf2

  # @infiltrate
  for i = 1:nc(dad)
    if interno
      hgd = [tsf1[i] -qsf1[i]
        tsf2[i] -qsf2[i]] \ [-v1[i], -v2[i]]
    else
      hgd = [tsf1[i] -qsf1[i]
        tsf2[i] -qsf2[i]] \ ([-v1[i], -v2[i]] - [tsf1[i], tsf2[i]])
    end
    H[i, i] = hgd[1]
    G[i, i] = hgd[2]
  end
end

function corrigediags!(H, G, dad, interno=true)
  xmin = minimum(dad.NOS[:, 1])
  xmax = maximum(dad.NOS[:, 1])
  ymin = minimum(dad.NOS[:, 2])
  ymax = maximum(dad.NOS[:, 2])

  l = sqrt((xmax - xmin)^2 + (ymax - ymin)^2)# Largura do ret�ngulo que cont�m a geometria

  d = 3

  nf = length(dad.ELEM[1].indices) * 3
  np = length(dad.ELEM[1].indices)

  ps = [((xmax + xmin) / 2 .+ d * l * sin.(range(0, 2pi, length=nf + 1)[1:end-1])) ((ymax + ymin) / 2 .+ d * l * cos.(range(0, 2pi, length=nf + 1)[1:end-1]))]

  tsf = zeros(eltype(H), nc(dad), nf)
  qsf = zeros(eltype(H), nc(dad), nf)

  for e in dad.ELEM, ii = 1:np, ifonte = 1:nf
    qsf[e.indices[ii], ifonte], tsf[e.indices[ii], ifonte] = BEM.calsolfund(dad.NOS[e.indices[ii], :] - ps[ifonte, :], dad.normal[e.indices[ii], :], dad)
    H[e.indices[ii], e.indices] .= 0
    G[e.indices[ii], e.indices] .= 0
  end
  v = H * tsf - G * qsf
  # @infiltrate
  for e in dad.ELEM, ii = 1:np
    # e = dad.ELEM[1]
    #  ii = 1
    if interno
      hgd = conj.([tsf[e.indices, :]; -qsf[e.indices, :]]') \ -v[e.indices[ii], :]
    else
      hgd = conj.([tsf[e.indices, :]; -qsf[e.indices, :]]') \ (-v[e.indices[ii], :] - tsf[e.indices[ii], :])
    end
    H[e.indices[ii], e.indices] = hgd[1:np]
    G[e.indices[ii], e.indices] = hgd[np+1:end]
  end
end