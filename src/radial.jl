function radial(xj, xi)
    r = norm(xj - xi) + 1e-9

    dx = xj[1] - xi[1]
    dy = xj[2] - xi[2]

    f = r^3
    dfdx = 3r * dx
    dfdy = 3r * dy

    # f = r^2 * log(r)
    # dfdx = dx * (log(r) + 1)
    # dfdy = dy * (log(r) + 1)


    # f = r
    # dfdx = dx
    # dfdy = dy
    if r == 0
        dfdxx = 0
        dfdyy = 0
        dfdxy = 0
    else
        dfdxx = 3 * (2dx^2 + dy^2) / r
        dfdyy = 3 * (dx^2 + 2dy^2) / r
        dfdxy = 3 * (dx * dy) / r
    end
    f, dfdx, dfdy, dfdxx, dfdyy, dfdxy
end

function montaF(pa, pint)
    n_pontos2 = size(pa, 1)
    n_pontos1 = size(pint, 1)
    F = zeros(n_pontos1, n_pontos2)
    Fx = zeros(n_pontos1, n_pontos2)
    Fy = zeros(n_pontos1, n_pontos2)
    Fxx = zeros(n_pontos1, n_pontos2)
    Fyy = zeros(n_pontos1, n_pontos2)
    Fxy = zeros(n_pontos1, n_pontos2)
    for i = 1:n_pontos1
        xi = pint[i, :]
        for j = 1:n_pontos2

            xj = pa[j, :]
            # @infiltrate
            aux = radial(xj, xi)
            F[i, j] = aux[1]
            Fx[i, j] = aux[2]
            Fy[i, j] = aux[3]
            Fxx[i, j] = aux[4]
            Fyy[i, j] = aux[5]
            Fxy[i, j] = aux[6]
        end
    end
    F, Fx, Fy, Fxx, Fyy, Fxy
end

function montaP(pa, pint)
    #como farei apenas para funções contantes ou lineares
    n_pontos2 = 3
    # n_pontos2 = 6
    n_pontos1 = size(pa, 1)
    P = zeros(n_pontos1, n_pontos2)
    dPx = zeros(n_pontos1, n_pontos2)
    dPy = zeros(n_pontos1, n_pontos2)
    dPxx = zeros(n_pontos1, n_pontos2)
    dPyy = zeros(n_pontos1, n_pontos2)
    dPxy = zeros(n_pontos1, n_pontos2)
    for i = 1:n_pontos1
        P[i, 1] = 1
        P[i, 2] = pa[i, 1]
        P[i, 3] = pa[i, 2]

        # P[i, 4] = pa[i, 1] * pa[i, 2]
        # P[i, 5] = pa[i, 1]^2
        # P[i, 6] = pa[i, 2]^2

        dPx[i, 1] = 0
        dPx[i, 2] = 1
        dPx[i, 3] = 0
        # dPx[i, 4] = pa[i, 2]
        # dPy[i, 5] = 2pa[i, 1]
        # dPy[i, 6] = 0

        dPy[i, 1] = 0
        dPy[i, 2] = 0
        dPy[i, 3] = 1
        # dPy[i, 4] = pa[i, 1]
        # dPy[i, 5] = 0
        # dPy[i, 6] = 2pa[i, 2]

        # dPxy[i, 4] = 1
    end
    P, dPx, dPy, dPxx, dPyy, dPxy
end
using LinearAlgebra


function montaFs(nos, nosi = nos; smooth = 0.0)
    # F, dFdx, dFdy, dFdxx, dFdyy, dFdxy = montaF(nos, nosi)
    # P, dPx, dPy, dPxx, dPyy, dPxy = montaP(nos, nosi)
    F1, dFdx, dFdy, dFdxx, dFdyy, dFdxy = montaF(nosi, nosi)
    P1, dPx, dPy, dPxx, dPyy, dPxy = montaP(nosi, nosi)
    F, dFdx, dFdy, dFdxx, dFdyy, dFdxy = montaF(nos, nosi)
    P, dPx, dPy, dPxx, dPyy, dPxy = montaP(nos, nosi)
    In = Matrix{}(I, size(F1))
    F1 = F1 + I * smooth
    W = F1 \ P1 / (P1' / F1 * P1)
    # @infiltrate
    aux = (In - W * P1') / F1
    N = aux * F + W * P'
    dNx = aux * dFdx + W * dPx'
    dNy = aux * dFdy + W * dPy'
    # dNxx = aux * dFdxx + W * dPxx'
    # dNyy = aux * dFdyy + W * dPyy'
    # dNxy = aux * dFdxy + W * dPxy'

    N', dNx', dNy'
    # dNx', dNy', dNxx', dNyy', dNxy'
end

# pint = range(0, 1, length=10)
# pint = [(ones(10)*pint')[:] (ones(10) * pint')'[:]]

# F, dFdx, dFdy, dFdxx, dFdyy, dFdxy = montaF(pint, pint)
# P, dPx, dPy, dPxx, dPyy, dPxy = montaP(pint, pint)

# f(x) = x[1]^2 * x[2] + x[2]^3 + x[2]^2 / 2
# dfdx(x) = 2x[1] * x[2]
# dfdxx(x) = 2x[2]
# dfdy(x) = x[1]^2 + 3x[2]^2 + x[2] / 2

# In = Matrix{}(I, size(F))
# W = F \ P / (P' / F * P)
# # dNx = ((In - W * P') / F * dFdx + W * dPx')
# dNx = ((In - W * P') / F * dFdx + W * dPx')
# dNy = ((In - W * P') / F * dFdy + W * dPy')


# z = [f(pint[i, :]) for i = 1:100]
# zx = [dfdx(pint[i, :]) for i = 1:100]
# zxx = [dfdxx(pint[i, :]) for i = 1:100]
# zy = [dfdy(pint[i, :]) for i = 1:100]
# zax = dNx' * z
# zay = dNy' * z
# zax2 = dFdx' * (F \ z)
# zay2 = dFdy' * (F \ z)
# [zx zax zax2] .- zx
# @show sum(([zx zax zax2] .- zx) .^ 2, dims=1)
# @show sum(([zy zay zay2] .- zy) .^ 2, dims=1)
# zxxa1 = dNxx * z
# zxxa2 = dNx * (dNx * z)
# dNx, dNy = montaFs(pint)
