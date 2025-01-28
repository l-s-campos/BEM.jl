function integraelemsing(pf, nf, x, eta, w, elem, dad::placa_fina, eet)
    npg = size(w, 1)
    hs = zeros(Float64, 2 * 2 * size(elem), npg)
    gs = zeros(Float64, 2 * 2 * size(elem), npg)
    Nm = zeros(Float64, 2, 2 * size(elem))
    for k = 1:npg
        N, dN = calc_fforma(eta[k], elem)
        pg = N' * x # Ponto de gauss interpolador
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
        sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
        # @infiltrate
        uast, tast = calsolfund(pg', pf, [sy, -sx], nf, dad)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N
        hs[:, k] = (tast*Nm*dgamadqsi)[:]
        gs[:, k] = (uast*Nm*dgamadqsi)[:]
    end
    # @infiltrate
    d = (eta .- eet)
    hsd = [hs' gs'] .* d .^ 2
    as = hsd
    # @infiltrate
    Ne, dNe = calc_fforma_gen(eet, eta)
    past = Ne' * as
    dpast = dNe' * as
    # @infiltrate
    intgauss1 = ((hsd .- past .- dpast .* d) ./ (d .^ 2))' * w
    if abs(eet) == 1
        N, dN = calc_fforma(eet, elem)
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsif = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        N1, dN1 = calc_fforma(eet + 1e-10, elem)
        d2N = (dN1 - dN) / 1e-10
        d2xdqsi = d2N' * x # dx/dξ & dy/dξ
        gb2 = d2xdqsi * dxdqsi' / dgamadqsif^2#gamma/beta^2
        beta_m = 1 / dgamadqsif
        inte =
            intgauss1 +
            (dpast * log(2 / beta_m) * sign(-eet) - past * (gb2 * sign(-eet) + 1 / 2))'
    else
        inte =
            intgauss1 +
            (
                dpast * log(1 - eet) - dpast * log(1 + eet) +
                past * (-1 / (1 - eet) - 1 / (1 + eet))
            )'
    end
    # inte = intgauss1 + (dpast * log(1 - eet) - dpast * log(1 + eet) + past * (-1 / (1 - eet) - 1 / (1 + eet)))'
    intre = reshape(inte, 2, :)
    intre[:, 1:2*size(elem)], intre[:, 2*size(elem)+1:end]
end

function integraelemsing(pf, nf, x, eta, w, elem, dad::elastico, eet)
    npg = size(w, 1)
    hs = zeros(Float64, 2 * 2 * size(elem), npg)
    gs = zeros(Float64, 2 * 2 * size(elem), npg)
    Nm = zeros(Float64, 2, 2 * size(elem))
    ForwardDiff.derivative(g, pi / 4)integrando(eet)
    # @infiltrate
    d = (eta .- eet)
    hsd = [hs' gs'] .* d .^ 2
    as = hsd
    # @infiltrate
    Ne, dNe = calc_fforma_gen(eet, eta)
    past = Ne' * as
    dpast = dNe' * as
    # @infiltrate
    intgauss1 = ((hsd .- past .- dpast .* d) ./ (d .^ 2))' * w
    if abs(eet) == 1
        N, dN = calc_fforma(eet, elem)
        dxdqsi = dN' * x # dx/dξ & dy/dξ
        dgamadqsif = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
        N1, dN1 = calc_fforma(eet + 1e-10, elem)
        d2N = (dN1 - dN) / 1e-10
        d2xdqsi = d2N' * x # dx/dξ & dy/dξ
        gb2 = d2xdqsi * dxdqsi' / dgamadqsif^2#gamma/beta^2
        beta_m = 1 / dgamadqsif
        inte =
            intgauss1 +
            (dpast * log(2 / beta_m) * sign(-eet) - past * (gb2 * sign(-eet) + 1 / 2))'
    else
        inte =
            intgauss1 +
            (
                dpast * log(1 - eet) - dpast * log(1 + eet) +
                past * (-1 / (1 - eet) - 1 / (1 + eet))
            )'
    end
    # inte = intgauss1 + (dpast * log(1 - eet) - dpast * log(1 + eet) + past * (-1 / (1 - eet) - 1 / (1 + eet)))'
    intre = reshape(inte, 2, :)
    intre[:, 1:2*size(elem)], intre[:, 2*size(elem)+1:end]
end

function integrando(eet)
    N, dN = calc_fforma(eet, elem)
    pg = N' * x # Ponto de gauss interpolador
    dxdqsi = dN' * x # dx/dξ & dy/dξ
    dgamadqsi = norm(dxdqsi) # dΓ/dξ = J(ξ) Jacobiano
    sx = dxdqsi[1] / dgamadqsi # vetor tangente dx/dΓ
    sy = dxdqsi[2] / dgamadqsi # vetor tangente dy/dΓ
    # @infiltrate
    uast, tast = calsolfund(pg', pf, [sy, -sx], nf, dad)
    Nm[1, 1:2:end] = N
    Nm[2, 2:2:end] = N
    h = (tast*Nm*dgamadqsi)[:]
    g = (uast*Nm*dgamadqsi)[:]
    h, g
end
