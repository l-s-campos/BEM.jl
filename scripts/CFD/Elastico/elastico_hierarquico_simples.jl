## Início da análise
using DrWatson,
    Plots,
    CairoMakie,
    HaltonSequences,
    DelimitedFiles,
    LinearMaps,
    ForwardDiff,
    Krylov,
    LinearOperators,
    PreallocationTools
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
include(srcdir("CFD.jl"))

#Dados de entrada

nelem = 20 #Numero de elementos
order = 2

NPX = 30 #pontos internos na direção x
NPY = NPX #pontos internos na direção y
npg = 14    #apenas números pares

L = 1
Re = 400 # Na verdade é o 1/ν ou 1/μ

λ = 10^8

## Formatação dos dados ________________________________________________

dad = format_dad(cavity(nelem, order, L, 1, λ), NPX, NPY)

# ==============Matrizes===============#

begin
    H, G = BEM.calc_HeG(dad, npg, interno = true)

    #Hd, Gd = calc_HeG(dad, npg, interno = true)

    #M = BEM.Monta_M_Hd(dad, npg)
    M = BEM.Monta_M_RIMd(dad, npg)
end

A, b = BEM.aplicaCDC(H, G * Re, dad)
#A, b = BEM.aplicaCDC(Hd, Gd*Re, dad)

F1, dFdx, dFdy, dFdxx, dFdyy, dFdxy =
    BEM.montaF([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])

dNx = -dFdx / F1
dNy = -dFdy / F1
# ==============Discretização============#

begin
    ncont = nc(dad)
    nint = ni(dad)

    n = nint + ncont

    interno = nc(dad)+1:n
    contorno = 1:nc(dad)

    interno2 = 2nc(dad)+1:2n
    contorno2 = 1:2nc(dad)
end

# ==============Matrizes hierarquicas===============#

#parametros
nmax = 50
eta = 3
rank = 20
rtol = 1e-6
atol = 1e-8

#define variaveis
elements = [
    [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
    [Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for i = 1:ni(dad)]
]
# elements = [elements elements]'[:]
struct elastMatrix <: AbstractMatrix{Float64}
    A::Matrix{Float64}
end
function Base.getindex(K::elastMatrix, i::Int, j::Int)
    return K.A[i, j]
    # return K.A[2i-1:2i, 2j-1:2j]
end
Base.size(K::elastMatrix) = size(K.A)

Xclt =
    Yclt = BEM.ClusterTree(
        repeat(elements, inner = 2),
        BEM.PrincipalComponentSplitter(nmax = nmax),
    )
Xclt1 = Yclt1 = BEM.ClusterTree(elements, BEM.PrincipalComponentSplitter(nmax = nmax))
# adm = WeakAdmissibilityStd()
adm = BEM.StrongAdmissibilityStd(eta = eta)
comp = BEM.PartialACA(atol = atol, rank = rank, rtol = rtol)


#monta matrizes

HA = BEM.assemble_hmatrix(
    elastMatrix(A),
    Xclt,
    Yclt;
    adm,
    comp,
    threads = false,
    distributed = false,
)

HM = BEM.assemble_hmatrix(
    elastMatrix(M),
    Xclt,
    Yclt;
    adm,
    comp,
    threads = false,
    distributed = false,
)

HdFdx = BEM.assemble_hmatrix(
    elastMatrix(dFdx),
    Xclt1,
    Yclt1;
    adm,
    comp,
    threads = false,
    distributed = false,
)

HdFdy = BEM.assemble_hmatrix(
    elastMatrix(dFdy),
    Xclt1,
    Yclt1;
    adm,
    comp,
    threads = false,
    distributed = false,
)

HF = BEM.assemble_hmatrix(
    elastMatrix(F1),
    Xclt1,
    Yclt1;
    adm,
    comp,
    threads = false,
    distributed = false,
)

HdNx = BEM.assemble_hmatrix(
    elastMatrix(dNx),
    Xclt1,
    Yclt1;
    adm,
    comp,
    threads = false,
    distributed = false,
)

HdNy = BEM.assemble_hmatrix(
    elastMatrix(dNy),
    Xclt1,
    Yclt1;
    adm,
    comp,
    threads = false,
    distributed = false,
)

BEM.compress!(HA, TSVD(atol = atol, rank = rank, rtol = rtol))
BEM.compress!(HM, TSVD(atol = atol, rank = rank, rtol = rtol))
BEM.compress!(HdFdx, TSVD(atol = atol, rank = rank, rtol = rtol))
BEM.compress!(HdFdy, TSVD(atol = atol, rank = rank, rtol = rtol))
BEM.compress!(HF, TSVD(atol = atol, rank = rank, rtol = rtol))
BEM.compress!(HdNx, TSVD(atol = atol, rank = rank, rtol = rtol))
BEM.compress!(HdNy, TSVD(atol = atol, rank = rank, rtol = rtol))


compress_HA = BEM.compression_ratio(HA)
compress_HM = BEM.compression_ratio(HM)
compress_HdFdx = BEM.compression_ratio(HdFdx)
compress_HdFdy = BEM.compression_ratio(HdFdy)
compress_HF = BEM.compression_ratio(HF)
compress_HdNx = BEM.compression_ratio(HdNx)
compress_HdNy = BEM.compression_ratio(HdNy)


#______________________________


# =========Inicia variáveis===============#

begin
    iter = 0
    erro_vel = 1

    u = zeros(ni(dad) + nc(dad), 2)
    t = zeros(nc(dad), 2)
    u_before = deepcopy(u)
    nolinear = zeros(2 * (ni(dad) + nc(dad)))

    interno = nc(dad)+1:nc(dad)+ni(dad)
    contorno = 1:nc(dad)

    relaxation = 1

    dnlv = zeros(2 * n)

    dnl = zeros(2n, 2n)

    contrib = zeros(ni(dad))
    dudx = zeros(n)
    dvdx = zeros(n)
    dudy = zeros(n)
    dvdy = zeros(n)
    #x = zeros(2 * (ni(dad) + nc(dad)))
    x = A \ b
end

# ==============Solução===============#


# Uma vez no início:
dNx_diag, dNy_diag, dNx_block, dNy_block =
    diag(dNx), diag(dNy), dNx[interno, interno], dNy[interno, interno]





n = length(x) ÷ 2

p = [
    DiffCache(zeros(n)),
    DiffCache(zeros(n)),
    DiffCache(zeros(n)),
    DiffCache(zeros(n)),
    DiffCache(zeros(n)),
    DiffCache(zeros(n)),
    DiffCache(zeros(2n)),
]

u_fix = similar(x, ncont);
u_fix .= u[1:ncont, 1];
v_fix = similar(x, ncont);
v_fix .= u[1:ncont, 2];

function F(x, p)
    T = eltype(x)

    # p[5] u = Vector{T}(undef, n)# temp u
    # p[6] v = Vector{T}(undef, n)# temp v
    u = get_tmp(p[5], x)
    v = get_tmp(p[6], x)
    @views u[1:ncont] .= convert.(T, u_fix)
    @views v[1:ncont] .= convert.(T, v_fix)
    @views u[ncont+1:end] .= x[2nc(dad)+1:2:end]
    @views v[ncont+1:end] .= x[2nc(dad)+2:2:end]

    dNx_u = get_tmp(p[1], x)
    dNy_u = get_tmp(p[2], x)
    dNx_v = get_tmp(p[3], x)
    dNy_v = get_tmp(p[4], x)

    dNx_u = dNx * u
    dNy_u = dNy * u
    dNx_v = dNx * v
    dNy_v = dNy * v

    # p[7] result = similar(x, 2n)
    result = get_tmp(p[7], x)
    @views result[1:2:end] .= u .* dNx_u .+ v .* dNy_u
    @views result[2:2:end] .= u .* dNx_v .+ v .* dNy_v

    return A * x - (b - Re * M * result)
end


# function F(x, dNx_u_pre, dNy_u_pre, dNx_v_pre, dNy_v_pre)
dNx_u = dNx * u[:, 1]
dNy_u = dNy * u[:, 1]
dNx_v = dNx * u[:, 2]
dNy_v = dNy * u[:, 2]

# dNx_u_pre = DiffCache(dNx_u)
# dNx_u_pre = zeros(BEM.Dual{Float64}, length(dNx_u))

# struct DualVector{T}
#     primal::T   # primeiro vetor
#     dual::T     # segundo vetor
# end

# dNx_u_pre = DualVector(dNx_u, dNx_u)
# dNy_u_pre = DualVector(dNy_u, dNy_u)
# dNx_v_pre = DualVector(dNx_v, dNx_v)
# dNx_v_pre = DualVector(dNy_v, dNy_v)

dNx_u_p = DiffCache(dNx_u)
dNy_u_p = DiffCache(dNy_u)
dNx_v_p = DiffCache(dNx_v)
dNy_v_p = DiffCache(dNy_v)
u_aux_p = DiffCache(dNx_u)
v_aux_p = DiffCache(dNx_u)
x_aux_p = DiffCache(x)

function F(x, dNx_u, dNy_u, dNx_v, dNy_v, u_aux, v_aux, x_aux, u_fix, v_fix, b)
    T = eltype(x)
    n = length(x) ÷ 2

    # u_aux = Vector{n,T}
    # v_aux = Vector{n,T}

    u_aux = get_tmp(u_aux, x)
    v_aux = get_tmp(v_aux, x)
    dNx_u = get_tmp(dNx_u, x)
    dNy_u = get_tmp(dNy_u, x)
    dNx_v = get_tmp(dNx_v, x)
    dNy_v = get_tmp(dNy_v, x)
    result = get_tmp(x_aux, x)
    # @show typeof(dNx_u), typeof(u_aux)


    @views u_aux[1:ncont] .= convert.(T, u_fix)
    @views v_aux[1:ncont] .= convert.(T, v_fix)
    @views u_aux[ncont+1:end] .= x[2ncont+1:2:end]
    @views v_aux[ncont+1:end] .= x[2ncont+2:2:end]

    # # dNx_u = dNx * u_aux
    # # dNy_u = dNy * u_aux
    # # dNx_v = dNx * v_aux
    # # dNy_v = dNy * v_aux

    mul!(dNx_u, dNx, u_aux)
    mul!(dNy_u, dNy, u_aux)
    mul!(dNx_v, dNx, v_aux)
    mul!(dNy_v, dNy, v_aux)

    # result = similar(x, 2n)
    @views result[1:2:2n] .= u_aux .* dNx_u .+ v_aux .* dNy_u
    @views result[2:2:2n] .= u_aux .* dNx_v .+ v_aux .* dNy_v

    # #@views result .= [nl_i nl_p]'[:]
    mul!(result, HM, result)
    mul!(result, HA, x, 1.0, Re)#mul!(C, A, B, α, β)
    # #Combined inplace matrix-matrix or matrix-vector multiply-add A B α + C β. 
    return result - b

    # return HA * x - (b - Re * HM * result)
    # result
end

# J(y, v) = ForwardDiff.derivative!(y, t -> F(x + t * v, p), 0)
J(y, v) = ForwardDiff.derivative!(
    y,
    t -> F(
        x + t * v,
        dNx_u_p,
        dNy_u_p,
        dNx_v_p,
        dNy_v_p,
        u_aux_p,
        v_aux_p,
        x_aux_p,
        u_fix,
        v_fix,
        b,
    ),
    0,
)
opH = LinearOperators.LinearOperator(Float64, 2n, 2n, true, true, (y, v) -> J(y, v))
@time J(x, ones(2n));
# @btime jacob_map!(x);

function jacob_map!(x)
    return HA * x + Re * HM * (dnl * x)
end

function jacob_map2!(x)
    contribu =
        u[:, 1] .* (HdNx * [zeros(nc(dad)); x[2nc(dad)+1:2:end]]) +
        u[:, 2] .* (HdNy * [zeros(nc(dad)); x[2nc(dad)+1:2:end]])
    contribv =
        u[:, 1] .* (HdNx * [zeros(nc(dad)); x[2nc(dad)+2:2:end]]) +
        u[:, 2] .* (HdNy * [zeros(nc(dad)); x[2nc(dad)+2:2:end]])

    dnl_i1 = contribu + Diagonal(dudx[:]) * [zeros(nc(dad)); x[2nc(dad)+1:2:end]]
    dnl_i2 = Diagonal(dudy[:]) * [zeros(nc(dad)); x[2nc(dad)+2:2:end]]
    dnl_p1 = Diagonal(dvdx[:]) * [zeros(nc(dad)); x[2nc(dad)+1:2:end]]
    dnl_p2 = contribv + Diagonal(dvdy[:]) * [zeros(nc(dad)); x[2nc(dad)+2:2:end]]

    dnlv[1:2:end] = dnl_i1 + dnl_i2 # u-u
    dnlv[2:2:end] = dnl_p1 + dnl_p2 # v-u
    #     dnl[2nc(dad)+1:2:end, 2nc(dad)+1:2:end] = dnl_i1  # u-u
    #     dnl[2nc(dad)+1:2:end, 2nc(dad)+2:2:end] = dnl_i2  # u-v  
    #     dnl[2nc(dad)+2:2:end, 2nc(dad)+1:2:end] = dnl_p1  # v-u
    #     dnl[2nc(dad)+2:2:end, 2nc(dad)+2:2:end] = dnl_p2  # v-v
    return HA * x + Re * HM * dnlv

end

tempo = @elapsed begin

    if iter == 0
        global iter, u, t, erro, u_before, nolinear, x, dnl
    end

    while erro_vel > 10^-6


        dudx = HdNx * u[:, 1]
        dudy = HdNy * u[:, 1]

        dvdx = HdNx * u[:, 2]
        dvdy = HdNy * u[:, 2]

        nolinear[1:2:2*n] = u[:, 1] .* dudx + u[:, 2] .* dudy
        nolinear[2:2:2*n] = u[:, 1] .* dvdx + u[:, 2] .* dvdy


        # Solution

        #calc dnl


        # tmonta = @elapsed begin
        #     contrib = u[nc(dad)+1:end, 1] .* dNx_block + u[nc(dad)+1:end, 2] .* dNy_block


        #     dnl_i1 = contrib + Diagonal(dudx[nc(dad)+1:end])
        #     dnl_i2 = Diagonal(dudy[nc(dad)+1:end])
        #     dnl_p1 = Diagonal(dvdx[nc(dad)+1:end])
        #     dnl_p2 = contrib + Diagonal(dvdy[nc(dad)+1:end])

        #     dnl[2nc(dad)+1:2:end, 2nc(dad)+1:2:end] = dnl_i1  # u-u
        #     dnl[2nc(dad)+1:2:end, 2nc(dad)+2:2:end] = dnl_i2  # u-v  
        #     dnl[2nc(dad)+2:2:end, 2nc(dad)+1:2:end] = dnl_p1  # v-u
        #     dnl[2nc(dad)+2:2:end, 2nc(dad)+2:2:end] = dnl_p2  # v-v
        # end

        #----krylov


        # opH = LinearOperators.LinearOperator(Float64, 2n, 2n, true, true, (y, v) -> J(y, v))


        # tresolve = @elapsed result = gmres(
        #     opH,
        #     F(
        #         x,
        #         dNx_u_p,
        #         dNy_u_p,
        #         dNx_v_p,
        #         dNy_v_p,
        #         u_aux_p,
        #         v_aux_p,
        #         x_aux_p,
        #         u_fix,
        #         v_fix,
        #         b,
        #     ),
        #     rtol = rtol,
        # )

        Residuo = HA * x - (b - Re * HM * nolinear)
        jacobiano = LinearMap(jacob_map2!, 2n, 2n)

        # d0=luA\(Residuo- Re * HM * dnlv)
        tresolve = @elapsed result = gmres(jacobiano, Residuo, d0, rtol = rtol)

        #tresolve = @elapsed result = gmres(opH, F(x, p), rtol = rtol)
        x = x - result[1]

        #println("tempo de resolver: $tresolve e tempo de montar: $tmonta")


        u[1:nc(dad), :] =
            u[1:nc(dad), :] * (1 - relaxation) + relaxation * separa(dad, x)[1]
        t = t * (1 - relaxation) + relaxation * separa(dad, x)[2]

        u[interno, 1] =
            u[interno, 1] * (1 - relaxation) + relaxation * x[2*nc(dad)+1:2:end]
        u[interno, 2] =
            u[interno, 2] * (1 - relaxation) + relaxation * x[2*nc(dad)+2:2:end]


        erro_vel = nrmse(u_before, u)


        println("Iteração: $iter; erro_vel = $erro_vel")

        u_before = deepcopy(u)
        iter = iter + 1
    end

end

#________________________

# Computando variáveis

begin
    uc = u[contorno, :]
    ui = u[interno, :]

    dudx = dNx * u[:, 1]
    dudy = dNy * u[:, 1]
    dvdx = dNx * u[:, 2]
    dvdy = dNy * u[:, 2]

    p = -λ * (dudx + dvdy)
end

