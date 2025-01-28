function telles(gamm, eet, r = 0)

    eest = eet^2 - 1
    term1 = eet * eest + abs(eest)
    if term1 < 0
        term1 = (-term1)^(1 / 3)
        term1 = -term1
    else
        term1 = term1^(1 / 3)
    end

    term2 = eet * eest - abs(eest)
    if term2 < 0
        term2 = (-term2)^(1 / 3)
        term2 = -term2
    else
        term2 = term2^(1 / 3)
    end
    GAMM = term1 + term2 + eet

    if r < 0.05
        r = 0
    elseif r < 1.3
        r = 0.85 + 0.24 * log(r)
    elseif r < 3.618
        r = 0.893 + 0.0832 * log(r)
    else
        r = 1
    end

    Q = 1 + 3 * GAMM^2
    A = (1 - r) / Q
    B = -3 * GAMM * (1 - r) / Q
    C = (r + 3 * GAMM^2) / Q
    D = -B

    eta = A * gamm .^ 3 + B * gamm .^ 2 .+ C * gamm .+ D
    Jt = 3 * A * gamm .^ 2 .+ 2 * B * gamm .+ C
    return eta, Jt
end


function sinhtrans(u, a, b)
    #https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.1208
    # https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.2244 iterado
    # https://www.sciencedirect.com/science/article/pii/S0955799715001125 3d
    if b < 1e-12
        b = 1e-12
        # return telles(u, a, b)
        # return Monegato(u, a)
    end
    μ = 1 / 2 * (asinh((1 + a) / b) + asinh((1 - a) / b))
    η = 1 / 2 * (asinh((1 + a) / b) - asinh((1 - a) / b))

    # x=(asinh.((u.-a)./b).+η)./μ
    x = a .+ b * sinh.(μ * u .- η)
    J = b * μ * cosh.(μ * u .- η)
    x, J
end

function Ma(u, a, b)
    # https://link.springer.com/content/pdf/10.1007%2Fs00466-002-0340-0.pdf
    # recomenda dividir a integral em 2 na singularidade
    J = sqrt.(b .^ 2 .+ (u .- a) .^ 2)
    # x=log.(J.+(u.-a))
    # xi=log.(sqrt.(b.^2 .+(-1 .-a).^2).+(-1. -a))
    # xf=log.(sqrt.(b.^2 .+(1 .-a).^2).+(1. - a))
    x = 1 / 2 * (exp.(u) .- b^2 * exp.(-u)) .+ a
    xi = 1 / 2 * (exp.(-1) .- b^2 * exp.(1)) .+ a
    xf = 1 / 2 * (exp.(1) .- b^2 * exp.(-1)) .+ a
    x = 2 * (x .- (xf + xi) / 2) / (xf - xi)
    x, J
end


function Wenjing(r, a, m = 10)
    # https://link.springer.com/article/10.1007/s00466-008-0262-6
    r = (r .+ 1) / 2
    ro = 1 / (m - 1) * ((1 .- r .^ m) ./ (1 .- r) .- 1)
    dro = ((m * (ro .- 1) - ro) .* ro .^ m + ro) ./ ((m - 1) * (ro .- 1) .^ 2 .* ro)
    # dro = ((m * (r .- 1) - r) .* r .^ m + r) ./ ((m - 1) * (r .- 1) .^ 2 .* r)
    if a == -1
        return 2 * ro .- 1, dro
    elseif a == 1
        return 1 .- 2ro, dro
    else
        return [(a + 1) * ro .- 1; (1 - a) * ro .+ 0],
        [(a + 1) / 2 * dro; (1 - a) / 2 * dro]
    end

end
function Xie(u, a, b)
    # https://www.sciencedirect.com/science/article/pii/S0955799711000208



end
#
#=
f(x)=(1 .-x.^2)./sqrt.((x.-a).^2 .+b^2)
f1(x)=.5*log.((x.-a).^2 .+b^2)
f2(x)=.5 ./((x.-a).^2 .+b^2)
a=1
b=0.001

ana1=(-0.5 *(a - 1)* (log((a - 1)^2 + b^2) - 2) - b *atan((a - 1)/b))-(-0.5* (a + 1)* (log((a + 1)^2 + b^2) - 2) - b*atan((a + 1)/b))
ana2=-(0.5*atan((a - 1)/b))/b+(0.5*atan((a + 1)/b))/b
x, w = BEM.gausslegendre(20);
x1, j1 = BEM.sinhtrans(x, a, b);
x2, j2 = BEM.telles(x, a);
x3, j3 = BEM.Wenjing(x, a, 100);
x4, j4 = BEM.Monegato(x, a, 5.0);
[sum(f1(x) .* w)  sum(f1(x1) .* w .* j1) sum(f1(x2) .* w .* j2) sum(f1(x3) .* w .* j3) sum(f1(x4) .* w .* j4)] .- ana1
[sum(f2(x) .* w)  sum(f2(x1) .* w .* j1) sum(f2(x2) .* w .* j2) sum(f2(x3) .* w .* j3) sum(f2(x4) .* w .* j4)] .- ana2

=#
function pontosintegra(NOS, elem_j, ind_elem, qsi, w)

    nosing = elem_j.indices .== ind_elem

    if sum(nosing) == 1 #integração singular

        eet = (elem_j.ξs[nosing])[1]
        if eet^2 != 1
            x1 = qsi * (eet + 1) / 2 .+ (eet - 1) / 2
            x2 = qsi * (1 - eet) / 2 .+ (eet + 1) / 2
            qsis = [x1; x2]
            ws = [w * (1 + eet) / 2; w * (1 - eet) / 2]
            eta, Jt = Monegato(qsis, eet, 5.0)
            inds = eta .!= eet
            return eta[inds], ws[inds] .* Jt[inds]
        end
        # eta,Jt=sinhtrans(qsi,eet[1],0)
        eta, Jt = Monegato(qsi, eet, 5.0)
        inds = eta .!= eet
        return eta[inds], w[inds] .* Jt[inds]
    else
        xel = NOS[elem_j.indices, :]   # Coordenada (x,y) dos nós geométricos
        pf = NOS[ind_elem, :]

        if distancia(xel[1], xel, pf, elem_j) > 5 * elem_j.comprimento
            # @show distancia(xel[1],xel,pf,elem_j),elem_j.comprimento
            return qsi, w
        end
        res = optimize(
            x -> distancia(x[1], xel, pf, elem_j),
            [0.0],
            Newton(),
            Optim.Options(g_tol = 1e-3, iterations = 10),
        )
        # @infiltrate
        # eta,Jt=sinhtrans(qsi,Optim.minimizer(res)[1],Optim.minimum(res))
        eta, Jt =
            sinhtrans(qsi, Optim.minimizer(res)[1], Optim.minimum(res) / elem_j.comprimento)
        # @show Optim.minimizer(res)[1],Optim.minimum(res)/elem_j.comprimento
        return eta, w .* Jt
    end

end

function distancia(qsi, x, pf, elem)
    N_geo, ~ = calc_fforma(qsi, elem)
    ps = N_geo' * x
    norm(ps' - pf)
end


function Monegato(t, s0, q = 5.0)
    s0 = min(max(-1, s0), 1)

    δ = 2^(-q) * ((1 + s0)^(1 / q) + (1 - s0)^(1 / q))^q
    t0 = ((1 + s0)^(1 / q) - (1 - s0)^(1 / q)) / ((1 + s0)^(1 / q) + (1 - s0)^(1 / q))
    s = s0 .+ δ * (t .- t0) .^ q
    ds = q * δ * (t .- t0) .^ (q - 1)
    s, ds
end
