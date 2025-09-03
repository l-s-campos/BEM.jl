## Início da análise
using DrWatson, Plots, CairoMakie
@quickactivate "BEM"
include(scriptsdir("includes.jl"))

#Dados de entrada

nelem = 50 #Numero de elementos
order = 2
nelem_circle = 20

# NPX = 80
for NPX in [10, 20, 30, 40, 50, 60, 70, 80] #pontos internos na direção x
    NPY = NPX #pontos internos na direção y
    npg = 16    #apenas números pares

    L = 1
    Re = 1 # Na verdade é o 1/ν ou 1/μ

    caso = "Cavidade"
    λ = 10^6

    dad = format_dad(cavity(nelem, order, L, 1, λ), NPX, NPY)

    begin
        H, G = calc_HeG(dad, npg, interno = true)

        M = BEM.Monta_M_RIMd(dad, npg)
    end

    # ==============Discretização============#

    begin
        nx = length(unique(dad.pontos_internos[:, 1]))
        ny = length(unique(dad.pontos_internos[:, 2]))

        x_order = sort(unique(dad.pontos_internos[:, 1]))
        y_order = sort(unique(dad.pontos_internos[:, 2]))

        dx = (maximum(x_order) - minimum(x_order)) / nx
        dy = (maximum(y_order) - minimum(y_order)) / ny

    end

    # =========Inicia variáveis===============#

    begin

        n = ni(dad) + nc(dad)

        interno = nc(dad)+1:nc(dad)+ni(dad)
        contorno = 1:nc(dad)

        x = zeros(2 * nc(dad) + 2 * ni(dad))
        p = zeros(ni(dad) + nc(dad))


    end
    dNx, dNy =
        BEM.montaFs([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])[2:3]

    # ==============Solução===============#
    # =========HMAT===============#


    # begin
    # Re = 100
    for Re in [50, 100, 200, 300, 400]

        if Re < 200
            relax = 0.5
        elseif Re < 300
            relax = 0.1
        elseif Re < 600
            relax = 0.01
        else
            relax = 0.005
        end

        begin
            x = zeros(2 * nc(dad) + 2 * ni(dad))
            A, b = BEM.aplicaCDC(H, Re * G, dad)
            elements = [
                [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
                [
                    Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for
                    i = 1:ni(dad)
                ]
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
                Yclt = ClusterTree(
                    repeat(elements, inner = 2),
                    BEM.PrincipalComponentSplitter(nmax = 30),
                )
            Xclt1 = Yclt1 = ClusterTree(elements, BEM.PrincipalComponentSplitter(nmax = 30))
            # adm = WeakAdmissibilityStd()
            adm = StrongAdmissibilityStd(eta = 20.0)
            comp = PartialACA(atol = 1e-6, rank = 20, rtol = 1e-6)
            # comp = TSVD(atol = 1e-6, rank = 20, rtol = 1e-6)
            KA = elastMatrix(A)
            KM = elastMatrix(M)
            HA = assemble_hmatrix(
                KA,
                Xclt,
                Yclt;
                adm,
                comp,
                threads = false,
                distributed = false,
            )
            HM = assemble_hmatrix(
                KM,
                Xclt,
                Yclt;
                adm,
                comp,
                threads = false,
                distributed = false,
            )

            HdNx = assemble_hmatrix(
                elastMatrix(dNx),
                Xclt1,
                Yclt1;
                adm,
                comp,
                threads = false,
                distributed = false,
            )
            HdNy = assemble_hmatrix(
                elastMatrix(dNy),
                Xclt1,
                Yclt1;
                adm,
                comp,
                threads = false,
                distributed = false,
            )
            BEM.compress!(HA, TSVD(atol = 1e-6, rank = 20, rtol = 1e-6))
            BEM.compress!(HM, TSVD(atol = 1e-6, rank = 20, rtol = 1e-6))
            # tf1 = @timed FA = lu(A)#, atol = 1e-6, rank = 10, rtol = 1e-6)
            # tf2 = @timed FHA = lu(HA)#, atol = 1e-6, rank = 10, rtol = 1e-6)
            @show 2 * nc(dad) + 2 * ni(dad) compression_ratio(HM)
        end


        # t1 = @timed @profile passo_não_linear(
        #     x,
        #     dad,
        #     FA,
        #     b,
        #     M,
        #     dNx,
        #     dNy,
        #     Re = Re,
        #     relaxation = relax,
        #     erro_vel_min = 1e-6,
        #     maxiter = 1000,
        # )


        # t2 = @timed passo_não_linear_Hmat(
        #     x,
        #     dad,
        #     FHA,
        #     b,
        #     HM,
        #     HdNx,
        #     HdNy,
        #     Re = Re,
        #     relaxation = relax,
        #     erro_vel_min = 1e-6,
        #     maxiter = 1000,
        # )
        # println()
        # println(
        #     Re,
        #     ",",
        #     2 * nc(dad) + 2 * ni(dad),
        #     ",",
        #     t1.time,
        #     ",",
        #     t2.time,
        #     ",",
        #     t1.time / t2.time,
        #     ",",
        #     tf1.time,
        #     ",",
        #     tf2.time,
        #     ",",
        #     tf1.time / tf2.time,
        # )
    end
end
#     end
# end
# x = zeros(2 * nc(dad) + 2 * ni(dad));
# x2 = passo_não_linear2(
#     x,
#     dad,
#     H,
#     G,
#     M,
#     dNx,
#     dNy,
#     Re = Re,
#     relaxation = 0.1,
#     erro_vel_min = 1e-6,
#     maxiter = 1000,
# )
# x = zeros(2 * nc(dad) + 2 * ni(dad));
# x3 = passo_não_linear3(
#     x,
#     dad,
#     H,
#     G,
#     M,
#     dNx,
#     dNy,
#     Re = Re,
#     relaxation = 0.5,
#     erro_vel_min = 1e-6,
#     maxiter = 1000,
# )

# using NonlinearSolve, PreallocationTools
# A, b = BEM.aplicaCDC(H, Re * G, dad)
# nolinear = DiffCache(zeros(size(x)))
# u = DiffCache(zeros(nc(dad) + ni(dad), 2))
# t = DiffCache(zeros(nc(dad), 2))
# dudx = DiffCache(zeros(nc(dad) + ni(dad)))
# dudy = DiffCache(zeros(nc(dad) + ni(dad)))
# dvdx = DiffCache(zeros(nc(dad) + ni(dad)))
# dvdy = DiffCache(zeros(nc(dad) + ni(dad)))
# aux = DiffCache(zeros((nc(dad) + ni(dad)) * 2))
# Residuo = similar(x)

# Af = lu!(A)
# ps = (dad, A, b, M, dNx, dNy, Re, u, t, nolinear, dudx, dudy, dvdx, dvdy)

# prob = NonlinearProblem(ResiduoCFD!, x, ps);
# @time sol =
#     solve(prob, abstol = 1e-6, reltol = 1e-6, NewtonRaphson(; autodiff = AutoForwardDiff()));
# @time sol = solve(prob, LevenbergMarquardt(; autodiff = AutoForwardDiff()));
# @time sol = solve(prob, DFSane());
# @time sol = solve(prob, Broyden());
# @time sol = solve(prob, RobustMultiNewton());


# x = sol.u


# using BifurcationKit

# par = [dad, H, G, M, dNx, dNy, 800.0];

# função_não_linear(x, par)
# prob2 = BifurcationProblem(função_não_linear, x, par, 7)

# sol = BifurcationKit.solve(prob2, Newton(), NewtonPar(tol = 1e-6, verbose = true))
# x = sol.u

# br2 = continuation(
#     prob2,
#     PALC(),
#     ContinuationPar(p_min = 10.0, p_max = 100.0, newton_options = NewtonPar(tol = 1e-6)),
# )
# x = br2.sol[2].x



# begin
#     u = zeros(nc(dad) + ni(dad), 2)
#     u[contorno, :], t = separa(dad, x)
#     u[interno, 1] = x[2*nc(dad)+1:2:end]
#     u[interno, 2] = x[2*nc(dad)+2:2:end]
#     uc = u[contorno, :]
#     ui = u[interno, :]

#     dudx = dNx * u[:, 1]
#     dudy = dNy * u[:, 1]
#     dvdx = dNx * u[:, 2]
#     dvdy = dNy * u[:, 2]

#     p = -λ * (dudx + dvdy)

#     # #=========Heatmap e Quiver=========#

#     x_array = dad.pontos_internos[:, 1]
#     y_array = dad.pontos_internos[:, 2]

#     utotal = sqrt.(ui[:, 1] .^ 2 + ui[:, 2] .^ 2)

#     BEM.heatmap(x_array, y_array, utotal, interpolate = true)
#     escala = 10^(-1)
#     BEM.quiver!(x_array, y_array, escala * ui[:, 1], escala * ui[:, 2], color = :black)
#     BEM.current_figure()
# end

