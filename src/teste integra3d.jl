"""
    sinhtrans(u, a, b)

Perform a sinh transformation on the input parameter `u`.

This function implements the sinh transformation as described in:
- https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.1208
- https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.2244 (iterated version)
- https://www.sciencedirect.com/science/article/pii/S0955799715001125 (3D version)

# Arguments
- `u`: coordinate to be transformed
- `a`: position of the singularity
- `b`: distance of the singularity

# Returns
- `x`: Transformed coordinate
- `J`: Jacobian of the transformation

# Note
If `b` is less than 1e-12, it is set to 1e-12 to avoid numerical instability.
"""
function sinhtrans(u, a, b)
    # Ensure minimum value for b to avoid numerical instability
    b = max(b, 1e-12)

    # Calculate transformation parameters
    μ = 0.5 * (asinh((1 + a) / b) + asinh((1 - a) / b))
    η = 0.5 * (asinh((1 + a) / b) - asinh((1 - a) / b))

    # Perform transformation
    x = a + b * sinh(μ * u - η)
    J = b * μ * cosh(μ * u - η)
    return x, J
end

function calc_points_quadsing(par1::Float64, par2::Float64, xx1, ww1, xx2, ww2, b = 1e-3)
    X = [
        -1 -1
        1 -1
        1 1
        -1 1
        -1 -1
    ]

    npg1 = length(xx1)
    npg2 = length(xx2)

    pesos = zeros(npg1 * npg2 * 4, 3)
    cont = 1

    for kk = 1:4
        r1 = [X[kk, 1] - par1, X[kk, 2] - par2]
        r2 = [X[kk+1, 1] - par1, X[kk+1, 2] - par2]
        nr1 = norm(r1)
        nr2 = norm(r2)
        ang = acos((nr1^2 + nr2^2 - 4) / (2 * nr1 * nr2))
        ang1 = asin(nr2 * sin(ang) / 2)

        for i = 1:npg1
            for j = 1:npg2
                t = (xx2[j] + 1) / 2 * ang  # theta
                d = sin(t) * nr1 / sin(π - t - ang1)
                px = X[kk, 1] - (X[kk, 1] - X[kk+1, 1]) * d / 2

                py = X[kk, 2] - (X[kk, 2] - X[kk+1, 2]) * d / 2
                ri = [px - par1, py - par2]
                nri = norm(ri)

                nx, Jt = sinhtrans(xx1[i], -1, b / nri)

                pesos[cont, 1] = ri[1] * (nx + 1) / 2 + par1
                pesos[cont, 2] = ri[2] * (nx + 1) / 2 + par2

                # pesos[cont, 1] = ri[1] * (xx1[i] + 1) / 2 + par1
                # pesos[cont, 2] = ri[2] * (xx1[i] + 1) / 2 + par2

                # Je = ang / 2 * nri / 2 * nri * (xx1[i] + 1) / 2 * Jt
                Je = ang / 2 * nri / 2 * nri * (nx + 1) / 2 * Jt
                pesos[cont, 3] = ww1[i] * ww2[j] * Je

                cont += 1
            end
        end
    end
    pesos
end

using FastGaussQuadrature, LinearAlgebra, Plots
x1, w1 = gausslegendre(20);    # Quadratura de gauss
x2, w2 = gausslegendre(10);    # Quadratura de gauss
pesos = calc_points_quadsing(0.0, 0.0, x1, w1, x2, w2, sqrt(1e-3));
f(x, y, b) = x / sqrt((x - 5)^2 + (y - 5)^2)^b
f.(5pesos[:, 1] .+ 5, 5pesos[:, 2] .+ 5, 0)' * pesos[:, 3] * 25 - 500
f.(5pesos[:, 1] .+ 5, 5pesos[:, 2] .+ 5, 1)' * pesos[:, 3] * 25 - 176.2747

#maxima - quad_qags(integrate(1/(x^2+y^2+1/1000), x, -1, 1),y,-1,1);
f(x, y, b) = 1 / ((y^2 + x^2)^b + 1e-3)
f.(pesos[:, 1], pesos[:, 2], 1)' * pesos[:, 3] - 22.39523274905529
#0.0020785240507308345 20/10 - sem a transformada
#8.0866868756857e-11 20/10 - com a transformada
