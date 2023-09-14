function mudavar2(x, p=5)
    t = x^(1 / p) / (x^(1 / p) + (1 - x)^(1 / p))
end
function mudavar1(t, p=5)
    x = t^(p) / (t^(p) + (1 - t)^(p))
end
function f_prime(t, p=5)
    numerator = p * (-(t - 1) * t)^(p - 1)
    denominator = (t^p + (1 - t)^p)^2
    return numerator / denominator
end
function corrige_x(x, b)
    for i in 1:size(x, 1)
        if x[i] > b
            x[i] = x[i] - b
        end
    end
    return x
end

"f->  funçao a ser integrada
m->ordem da singularidade
t->posição da singularidade 
n->Quantidade de pontos de integração
https://link.springer.com/article/10.1007/s10092-021-00446-1
t = 0.3
f1(x) = (1 + x - x^2) / (x - t)^1
f2(x) = (1 + x - x^2) / (x - t)^2
f3(x) = (1 + x - x^2) / (x - t)^3
eta, w = BEM.novelquad(2, (t - 0.5) * 2, 32)
F1(x)=-(t^2 - t - 1) log(x - t) - 1/2 x (2 t + x - 2)
dot(f1.(eta / 2 .+ 0.5), w) / 2-1.22523041106851637258923008288999
dot(f2.(eta / 2 .+ 0.5), w) / 2+6.42298561774988045927786175929650
dot(f3.(eta / 2 .+ 0.5), w) / 2-2.73546857952209343844408750481721"
function novelquad(m, t, n, p=5)
    a = 0
    b = 1
    h = (b - a) / n
    tau = mudavar2(t / 2 + 0.5, p)
    if m == 1
        xs = collect(tau .+ h * (1:n) .- h / 2)
        ws = ones(n) * h
    elseif m == 2 || m == 3
        xs = [tau .+ h * (1:n) .- h / 2; tau .+ h / 2 * (1:2n) .- h / 4]
        ws = [ones(n) * 2h; -ones(2n) * h / 2]
    elseif m == 4
        xs = [tau .+ h * (1:n) .- h / 2; tau .+ h / 2 * (1:2n) .- h / 4; tau .+ h / 4 * (1:4n) .- h / 8]
        ws = [ones(n) * 16h / 7; -ones(2n) * 5h / 7; -ones(4n) * h / 28]
    else
        return 0
    end
    return mudavar1.(corrige_x(xs, b)) * 2 .- 1, ws .* f_prime.(xs) * 2
end


function integraelemsing_num(pf, nf, x, elem, dad::Union{placa_fina,placa_fina_isotropica}, pre, xi0, npg=30)
    h = zeros(Float64, 2, 2 * size(elem))
    g = zeros(Float64, 2, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    eta, w = novelquad(2, xi0, npg)
    # @infiltrate

    for k = 1:size(w, 1)

        # @infiltrate
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        wast, vast, dast = calsolfund2(pg', pf, [sy, -sx], nf, dad, pre)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        h += vast * Nm * dgamadqsi * w[k]
        g += wast * Nm * dgamadqsi * w[k]

        # @infiltrate
    end
    h, g
end


function integradelemsing_num(pf, x, elem, dad::Union{elastico}, xi0, npg=30)
    d = zeros(Float64, 3, 2 * size(elem))
    s = zeros(Float64, 3, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    eta, w = novelquad(3, xi0, npg)
    # @infiltrate

    for k = 1:size(w, 1)

        # @infiltrate
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        D, S = caldsolfund(pg', pf, [sy, -sx], dad)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        D_2 = [D[1, 1, 1] D[2, 1, 1]
            D[1, 2, 2] D[2, 2, 2]
            D[1, 1, 2] D[2, 1, 2]]

        S_2 = [S[1, 1, 1] S[2, 1, 1]
            S[1, 2, 2] S[2, 2, 2]
            S[1, 1, 2] S[2, 1, 2]]
        # @infiltrate
        d += D_2 * Nm * dgamadqsi * w[k]
        s += S_2 * Nm * dgamadqsi * w[k]

        # @infiltrate
    end
    d, s
end

function integraelemsing_num(pf, x, elem, dad::Union{elastico}, xi0, npg=30)
    h = zeros(Float64, 2, 2 * size(elem))
    g = zeros(Float64, 2, 2 * size(elem))
    Nm = zeros(Float64, 2, 2 * size(elem))
    eta, w = novelquad(1, xi0, npg)
    # @infiltrate

    for k = 1:size(w, 1)

        # @infiltrate
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        uast, tast = calsolfund(pg', pf, [sy, -sx], dad, elem.regiao)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        h += tast * Nm * dgamadqsi * w[k]
        g += uast * Nm * dgamadqsi * w[k]
        # @infiltrate
    end
    h, g
end
