function calc_K_2d(x0, xf, n; E = 200e9)
    y = range(x0, stop = xf, length = n + 1)
    x = (y[1:end-1] + y[2:end]) / 2
    node = hcat(collect(1:n), collect(2:n+1)) # Element connectivity
    al = abs.(y[node[:, 2]] .- y[node[:, 1]])# Length of the element

    A = zeros(n, n)  # Matrix A
    for i = 1:n
        for j = 1:n # Loop on field points (Column)
            # Compute parameters used in the formulas for the two integrals
            x1 = y[node[j, 1]] .- x[i]  # Distance from the source point to the beginning of the element
            x2 = y[node[j, 2]] .- x[i]
            r1 = abs.(x1) # Distance from the source point to the beginning of the element
            r2 = abs.(x2)# Distance from the source point to the end of the element

            A[i, j] = (al[j] .+ x1 .* log.(r1) .- x2 .* log.(r2)) * -4 / (pi * E)

        end
    end
    return A
end
ext_sqrt(p, q) = p + sqrt(p .^ 2 + q .^ 2)
function calc_K_3d(x0, xf, n1, y0 = x0, yf = xf, n2 = n1; E = 200e9)
    y1 = range(x0, stop = xf, length = n1 + 1)
    y2 = range(y0, stop = yf, length = n2 + 1)
    x1 = (y1[1:end-1] + y1[2:end]) / 2
    x2 = (y2[1:end-1] + y2[2:end]) / 2
    y = hcat(repeat(y1, 1, n2 + 1)[:], repeat(y2', n1 + 1, 1)[:])
    x = hcat(repeat(x1, 1, n2)[:], repeat(x2', n1, 1)[:])

    ind = 1:n1+1
    inds = (repeat(ind', n2 + 1, 1) .+ collect((0:n2) * (n1 + 1)))'
    node = zeros(Int, n1 * n2, 4)
    for j = 1:n2
        for i = 1:n1
            node[i+(j-1)*n1, :] = [inds[i, j] inds[i+1, j] inds[i, j+1] inds[i+1, j+1]]
        end
    end
    # Element connectivity
    al = [abs.(y[node[:, 2], 1] .- y[node[:, 1], 1]) abs.(
        y[node[:, 3], 2] .- y[node[:, 1], 2]
    )]# Length of the element
    n = n1 * n2
    A = zeros(n, n)  # Matrix A
    for i = 1:n
        for j = 1:n # Loop on field points (Column)
            # Compute parameters used in the formulas for the two integrals
            x1 = y[node[j, 1], 1] .- x[i, 1]  # Distance from the source point to the beginning of the element
            x2 = y[node[j, 1], 2] .- x[i, 2]  # Distance from the source point to the beginning of the element


            term_1 =
                (x1 + 0.5 * al[j, 1]) .* log(
                    ext_sqrt(x2 + 0.5 * al[j, 2], x1 + 0.5 * al[j, 1]) ./
                    ext_sqrt(x2 - 0.5 * al[j, 2], x1 + 0.5 * al[j, 1]),
                )
            term_2 =
                (x2 + 0.5 * al[j, 2]) .* log(
                    ext_sqrt(x1 + 0.5 * al[j, 1], x2 + 0.5 * al[j, 2]) ./
                    ext_sqrt(x1 - 0.5 * al[j, 1], x2 + 0.5 * al[j, 2]),
                )
            term_3 =
                (x1 - 0.5 * al[j, 1]) .* log(
                    ext_sqrt(x2 - 0.5 * al[j, 2], x1 - 0.5 * al[j, 1]) ./
                    ext_sqrt(x2 + 0.5 * al[j, 2], x1 - 0.5 * al[j, 1]),
                )
            term_4 =
                (x2 - 0.5 * al[j, 2]) .* log(
                    ext_sqrt(x1 - 0.5 * al[j, 1], x2 - 0.5 * al[j, 2]) ./
                    ext_sqrt(x1 + 0.5 * al[j, 1], x2 - 0.5 * al[j, 2]),
                )


            A[i, j] = 2 / (pi * E) * (term_1 + term_2 + term_3 + term_4)

        end
    end
    return A
end


@kwdef mutable struct dadHS2D <: AbstractMatrix{Float64}
    y::Vector{Float64}
    x::Vector{Float64}
    node::Vector{SVector{2,Int64}}
    h0::Vector{Float64}
    al::Vector{Float64}
    E::Float64 = 200e9
end
Base.size(K::dadHS2D) = length(K.x), length(K.x)

function Base.getindex(K::dadHS2D, i::Int, j::Int)
    x = K.y[K.node[j]] .- K.x[i]  # Distance from the source point to the  element extremes
    r = abs.(x) # Distance from the source point to the beginning of the element
    (K.al[j] .+ x[1] .* log.(r[1]) .- x[2] .* log.(r[2])) * -4 / (pi * K.E)
end


function monta_hmat(dad::dadHS2D)
    xe = [Point2D(e, e) for e in dad.x]


    splitter = BEM.PrincipalComponentSplitter(; nmax = 10)
    Xclt = ClusterTree(xe, splitter)

    adm = StrongAdmissibilityStd(eta = 3)
    comp = PartialACA(; atol = 1e-6)

    assemble_hmatrix(dad, Xclt, Xclt; adm, comp)
end
@kwdef mutable struct dadHS3D <: AbstractMatrix{Float64}
    y::Vector{SVector{2,Int64}}
    x::Vector{SVector{2,Int64}}
    node::Vector{SVector{4,Int64}}
    h0::Vector{Float64}
    al::Vector{SVector{2,Int64}}
    E::Float64 = 200e9
end
Base.size(K::dadHS3D) = size(K.x, 1), size(K.x, 1)

function Base.getindex(K::dadHS3D, i::Int, j::Int)
    x = K.y[K.node[j]] .- K.x[i]  # Distance from the source point to the the  element extremes


    term_1 =
        (x[1] + 0.5 * K.al[j][1]) .* log(
            ext_sqrt(x2 + 0.5 * K.al[j][2], x1 + 0.5 * K.al[j][1]) ./
            ext_sqrt(x2 - 0.5 * K.al[j][2], x1 + 0.5 * K.al[j][1]),
        )
    term_2 =
        (x2 + 0.5 * K.al[j][2]) .* log(
            ext_sqrt(x1 + 0.5 * K.al[j][1], x2 + 0.5 * K.al[j][2]) ./
            ext_sqrt(x1 - 0.5 * K.al[j][1], x2 + 0.5 * K.al[j][2]),
        )
    term_3 =
        (x1 - 0.5 * K.al[j][1]) .* log(
            ext_sqrt(x2 - 0.5 * K.al[j][2], x1 - 0.5 * K.al[j][1]) ./
            ext_sqrt(x2 + 0.5 * K.al[j][2], x1 - 0.5 * K.al[j][1]),
        )
    term_4 =
        (x2 - 0.5 * K.al[j][2]) .* log(
            ext_sqrt(x1 - 0.5 * K.al[j][1], x2 - 0.5 * K.al[j][2]) ./
            ext_sqrt(x1 + 0.5 * K.al[j][1], x2 - 0.5 * K.al[j][2]),
        )


    2 / (pi * K.E) * (term_1 + term_2 + term_3 + term_4)
end


"""
    contact_pressure(p_min, H, z, W_aim, Nx1, Nx2, dx1, dx2, fft2_Kernel_circ, err_tol, it_max, h_ref)

Calculates the contact pressure occurring when a rigid smooth surface is loaded against an elastic half-space, following a method similar to Polonsky and Keer (1999).

# Arguments
- `dad`: estrutura de dados do problema (dadHS2D ou dadHS3D).
- `K`: [m/Pa] influence matrix (can be a dense or H-matrix).
- `W_aim`: [N] desired normal load.
- `p_min`: [Pa] minimum pressure 
- `H`: [Pa] maximum pressure (hard-wall constraint).
- `err_tol`: [-] relative error tolerance.
- `it_max`: [-] maximum number of iterations.
- `h_ref`: [m] reference length of relative error.

# Returns
- `p_con`: [Pa] contact pressure field (Nx1 x Nx2 matrix).
- `g`: [m] residual of the gap height distribution (Nx1 x Nx2 matrix).
- `err`: [-] relative error history (vector).
"""
function contact_pressure_force(
    dad,
    K,
    W_aim;
    p_min = 0,
    H = 5e8,
    err_tol = 1e-9,
    it_max = 100,
    h_ref = 1e-6,
)
    n = size(K, 1)
    areatotal = sum(dad.al)
    # Initialize pressure field (use Float64 for calculations)
    p_con = fill(W_aim / areatotal, n) .* dad.al # [Pa] initial pressure field

    # Compute initial elastic deformation
    u = K * p_con # [m] elastic deformation 
    g = -u + dad.h0 # [m] residual of the gap height distribution

    # --- Initial Classification ---

    # Find linear indices of points that are in the elastic domain
    # Use findall for logical indexing
    A_el = findall((p_con .> p_min) .& (p_con .< H))

    # Set residual reference plane to zero in the elastic domain:
    # Use mean(g[A_el]) for the mean of the elements at indices A_el
    if !isempty(A_el)
        g_mean = mean(g[A_el])
    else
        # Fallback if no points are initially elastic (should be rare)
        g_mean = 0.0
    end
    g = g .- g_mean

    # Find indices for non-contact and plastic domains
    # Non-contact points
    A_nc_cr = findall((p_con .<= p_min) .& (g .<= 0)) # Correctly non-contact (p <= p_min AND g <= 0)
    A_nc_wr = findall((p_con .<= p_min) .& (g .> 0))  # Wrongly non-contact (p <= p_min BUT g > 0)

    # Plastic points
    A_pl_cr = findall((p_con .>= H) .& (g .>= 0))     # Correctly plastic (p >= H AND g >= 0)
    A_pl_wr = findall((p_con .>= H) .& (g .< 0))      # Wrongly plastic (p >= H BUT g < 0)

    # Within these points, the pressure distribution needs to be adjusted:
    # A_free is the set of indices where the pressure is allowed to change.
    A_free = vcat(A_el, A_nc_wr, A_pl_wr)

    G = sum(g[A_free] .^ 2) # [m^2] norm of the residual
    G_old = 1.0 # [m^2] previous norm of the residual
    delta = 0.0 # [-] flag whether to use conjugate gradient or steepest descend
    err = zeros(it_max) # [-] relative error history
    i_it = 0 # [-] iteration counter        
    t = zeros(n) # [Pa] search direction

    while i_it == 0 || (err[i_it] > err_tol && i_it < it_max)
        i_it += 1

        # --- Find search direction (Conjugate Gradient or Steepest Descent) ---

        # The term delta*(G/G_old) is the Polak-Ribière parameter beta,
        # which is zero for Steepest Descent (delta=0).
        t[A_free] = g[A_free] + delta * (G / G_old) * t[A_free]
        t[A_nc_cr] .= 0.0
        t[A_pl_cr] .= 0.0

        # --- Determine step length tau ---

        # Compute elastic deformation due to search direction t
        r = -(K * t) # r is the change in gap height due to the change in pressure t

        # Adjust reference plane:
        if !isempty(A_el)
            r_mean = mean(r[A_el])
        else
            r_mean = 0.0
        end
        r = r .- r_mean

        # Compute tau (Optimal step size)
        # Numerator: sum(g .* t) in A_free
        # Denominator: sum(r .* t) in A_free
        tau = sum(g[A_free] .* t[A_free]) / sum(r[A_free] .* t[A_free])

        # --- Update pressure and confine it ---

        p_con[A_free] = p_con[A_free] .- tau * t[A_free]

        # Confinement (Projection onto the convex set p_min <= p <= H)
        p_con[p_con.<p_min] .= p_min
        p_con[p_con.>H] .= H

        # --- Rescale pressure to meet the desired normal load W_aim ---

        # Recalculate A_unscalable based on the *updated* p_con, 
        # as confinement might have moved points out of A_free.
        A_fixed = (p_con .<= p_min) .| (p_con .>= H)

        # The set to be rescaled should be the points *not* fixed to p_min or H
        W_unscalable = sum((p_con.*dad.al)[A_fixed]) # Total load from fixed points
        W_scalable = sum((p_con.*dad.al)[.!A_fixed]) # Total load from points that can be scaled

        # Rescale only the 'free' points (A_free_updated)
        p_con[.!A_fixed] = (W_aim - W_unscalable) / W_scalable * p_con[.!A_fixed]

        # Find indices of points that are in the elastic domain due to the pressure
        # distribution (after update and rescaling):
        A_el = findall((p_con .> p_min) .& (p_con .< H))

        # Compute elastic deformation to find new residual of the gap height distribution
        u = K * p_con # [m] elastic deformation
        g = -u + dad.h0
        # g = -u + z

        # Set residual reference plane to zero in the elastic domain:
        if !isempty(A_el)
            g_mean = mean(g[A_el])
        else
            g_mean = 0.0
        end
        g = g .- g_mean

        # Find new indices for non-contact and plastic domains
        A_nc_cr = findall((p_con .<= p_min) .& (g .>= 0)) # Correctly non-contact
        A_nc_wr = findall((p_con .<= p_min) .& (g .< 0))  # Wrongly non-contact
        A_pl_cr = findall((p_con .>= H) .& (g .>= 0))     # Correctly plastic
        A_pl_wr = findall((p_con .>= H) .& (g .< 0))      # Wrongly plastic

        # Within these points, the pressure distribution needs to be adjusted:
        A_free = vcat(A_el, A_nc_wr, A_pl_wr)

        # --- Determine whether to use conjugate gradient or steepest descend ---

        # Use Conjugate Gradient (delta=1) only if all non-contact and plastic points are 'correct'
        if isempty(A_nc_wr) && isempty(A_pl_wr)
            delta = 1.0
        else
            delta = 0.0
        end

        # Save G for the next iteration:
        G_old = G

        # Compute norm of the residual:
        G = sum(g[A_free] .^ 2)

        # Compute relative error
        err[i_it] = sqrt(G) / h_ref
        # @infiltrate i_it == 100
    end
    err = err[1:i_it] # Trim error history to actual number of iterations
    # @show err

    return p_con, g
end
"""
    contact_pressure_disp(dad, K, g_aim; ε = 1e-8, maxiter = 100)
Calculates the contact pressure occurring when a rigid smooth surface is loaded against an elastic half-space, following Normal adhesive contact on rough surfaces:
efficient algorithm for FFT-based BEM resolution
    https://hal.science/hal-01755724v1/document

    # Arguments
- `dad`: estrutura de dados do problema (dadHS2D ou dadHS3D).
- `K`: [m/Pa] influence matrix (can be a dense or H-matrix).
- `g_aim`: [m] desired gap height distribution.
- `p_min`: [Pa] minimum pressure 
- `H`: [Pa] maximum pressure (hard-wall constraint).
- `err_tol`: [-] relative error tolerance.
- `it_max`: [-] maximum number of iterations.
- `h_ref`: [m] reference length of relative error.

# Returns
- `p_con`: [Pa] contact pressure field (Nx1 x Nx2 matrix).
- `g`: [m] residual of the gap height distribution (Nx1 x Nx2 matrix).
- `err`: [-] relative error history (vector).
    """
function contact_pressure_disp(dad, K, g_aim; ε = 1e-8, maxiter = 100)
    n = size(K, 1)

    t = 0.0
    R_old = 1.0
    g = fill(g_aim, n) # uniform init
    Klu = lu(K)
    for iter = 1:maxiter
        g_old = copy(g)

        # Elastic + total pressure
        q = Klu \ (-g + dad.h0)
        # Residual
        ∂Ep∂u = q # placeholder, depends on formulation
        ∂Ep∂ub = mean(∂Ep∂u)

        # Recenter q
        ql = q .- ∂Ep∂ub
        R = norm(ql)
        t = ql .+ (R / R_old) .* t
        R_old = R

        # Update with conjugate gradient
        r = Klu \ (t)
        r̄ = mean(r)
        rl = r .- r̄

        τ = sum(ql .* t) / sum(rl .* t)
        g .-= τ .* t

        # Truncate negatives
        g = max.(g, 0.0)

        # Internal loop for active set
        In_a = findall(g .== 0.0 && ql .< 0.0)
        while !isempty(In_a)
            if iter == 1
                g[In_a] .-= τ .* q[In_a]   # first entry
            else
                g[In_a] .-= α .* q[In_a]   # subsequent entries
            end
            # Compute current mean gap
            G = mean(g)
            g = g_aim / G * g

            g_old = copy(g)
            q = Klu \ (-g + dad.h0)
            In_c = findall(g .> 0.0)
            # Recenter q
            q̄ = mean(q[In_c])
            ql = q .- q̄

            # Search direction
            t = zeros(length(g))
            t[In_c] .= q[In_c]
            r = Klu \ (t)
            r̄ = mean(r)
            rl = r .- r̄

            α = sum(ql .* t) / sum(rl .* t)
            g .-= α .* t
            g = max.(g, 0.0)
            In_a = findall(g .== 0.0 && ql .< 0.0)

        end
        G = mean(g)
        g = g_aim / G * g
        # Convergence check
        if dot(g, q) < ε
            break
        end
    end
    q = Klu \ (-g + dad.h0)
    return q, g
end

function montaFMM(dad::dadHS2D; rtol = 1e-5)
    n = length(dad.x)
    sources = Matrix{Float64}(undef, 2, n)

    for j = 1:n
        sources[:, j] = @SVector [dad.x[j], 0.0]  # 2D points in the xy-plane
    end

    charges = Vector{Float64}(undef, n)
    MFMM = LinearMaps.LinearMap{Float64}(n, n) do y, x
        # multiply by weights and constant
        @. charges = 4 / (pi * dad.E) * dad.al * x
        out = FMM2D.rfmm2d(; sources = sources, charges = charges, eps = rtol, pg = 1)
        return copyto!(y, out.pot)
    end


    return MFMM
    # return MFMM + LinearMap(Mcorrige(dad))
end

function Mcorrige(dad::dadHS2D)
    n = length(dad.x)
    M = spzeros(n, n)
    prox = Ponto_prox(dad, 3)
    for i = 1:n
        for j in prox[i] # Loop on field points (Column)
            if i == j
                x = dad.y[dad.node[j]] .- dad.x[i]  # Distance from the source point to the  element extremes
                r = abs.(x) # Distance from the source point to the beginning of the element
                M[i, j] =
                    (dad.al[j] .+ x[1] .* log.(r[1]) .- x[2] .* log.(r[2])) * -4 /
                    (pi * dad.E)
            else
                x = dad.y[dad.node[j]] .- dad.x[i]  # Distance from the source point to the  element extremes
                r = abs.(x) # Distance from the source point to the beginning of the element
                M[i, j] =
                    (dad.al[j] .+ x[1] .* log.(r[1]) .- x[2] .* log.(r[2])) * -4 /
                    (pi * dad.E) -
                    4 / (pi * dad.E) * dad.al[j] * log(abs(dad.x[i] - dad.x[j]))


            end
        end
    end
    return M
end


function montaFMM(dad::dadHS3D; rtol = 1e-5)
    n = length(dad.x)
    sources = Matrix{Float64}(undef, 3, n)

    for j = 1:n
        sources[:, j] = @SVector [dad.x[j], 0.0]  # 3D points in the xy-plane
    end

    charges = Vector{Float64}(undef, n)
    MFMM = LinearMaps.LinearMap{Float64}(n, n) do y, x
        # multiply by weights and constant
        @. charges = -4 / (pi * dad.E) * dad.al * x
        out = FMM2D.rfmm3d(; sources = sources, eps = rtol, pg = 1)
        return copyto!(y, out.pot)
    end


    return MFMM + LinearMap(Mcorrige)
end

function Mcorrige(dad::dadHS3D)
    n = length(dad.x)
    M = spzeros(n, n)
    prox = Ponto_prox(dad, 2)
    for i = 1:n
        for j in prox[i] # Loop on field points (Column)
            if i == j
                x = dad.y[dad.node[j]] .- dad.x[i]  # Distance from the source point to the the  element extremes


                term_1 =
                    (x[1] + 0.5 * dad.al[j][1]) .* log(
                        ext_sqrt(x2 + 0.5 * dad.al[j][2], x1 + 0.5 * dad.al[j][1]) ./
                        ext_sqrt(x2 - 0.5 * dad.al[j][2], x1 + 0.5 * dad.al[j][1]),
                    )
                term_2 =
                    (x2 + 0.5 * dad.al[j][2]) .* log(
                        ext_sqrt(x1 + 0.5 * dad.al[j][1], x2 + 0.5 * dad.al[j][2]) ./
                        ext_sqrt(x1 - 0.5 * dad.al[j][1], x2 + 0.5 * dad.al[j][2]),
                    )
                term_3 =
                    (x1 - 0.5 * dad.al[j][1]) .* log(
                        ext_sqrt(x2 - 0.5 * dad.al[j][2], x1 - 0.5 * dad.al[j][1]) ./
                        ext_sqrt(x2 + 0.5 * dad.al[j][2], x1 - 0.5 * dad.al[j][1]),
                    )
                term_4 =
                    (x2 - 0.5 * dad.al[j][2]) .* log(
                        ext_sqrt(x1 - 0.5 * dad.al[j][1], x2 - 0.5 * dad.al[j][2]) ./
                        ext_sqrt(x1 + 0.5 * dad.al[j][1], x2 - 0.5 * dad.al[j][2]),
                    )


                2 / (pi * dad.E) * (term_1 + term_2 + term_3 + term_4)

                M[i, j] = 2 / (pi * dad.E) * (term_1 + term_2 + term_3 + term_4)
            else
                x = dad.y[dad.node[j]] .- dad.x[i]  # Distance from the source point to the the  element extremes


                term_1 =
                    (x[1] + 0.5 * dad.al[j][1]) .* log(
                        ext_sqrt(x2 + 0.5 * dad.al[j][2], x1 + 0.5 * dad.al[j][1]) ./
                        ext_sqrt(x2 - 0.5 * dad.al[j][2], x1 + 0.5 * dad.al[j][1]),
                    )
                term_2 =
                    (x2 + 0.5 * dad.al[j][2]) .* log(
                        ext_sqrt(x1 + 0.5 * dad.al[j][1], x2 + 0.5 * dad.al[j][2]) ./
                        ext_sqrt(x1 - 0.5 * dad.al[j][1], x2 + 0.5 * dad.al[j][2]),
                    )
                term_3 =
                    (x1 - 0.5 * dad.al[j][1]) .* log(
                        ext_sqrt(x2 - 0.5 * dad.al[j][2], x1 - 0.5 * dad.al[j][1]) ./
                        ext_sqrt(x2 + 0.5 * dad.al[j][2], x1 - 0.5 * dad.al[j][1]),
                    )
                term_4 =
                    (x2 - 0.5 * dad.al[j][2]) .* log(
                        ext_sqrt(x1 - 0.5 * dad.al[j][1], x2 - 0.5 * dad.al[j][2]) ./
                        ext_sqrt(x1 + 0.5 * dad.al[j][1], x2 - 0.5 * dad.al[j][2]),
                    )


                2 / (pi * dad.E) * (term_1 + term_2 + term_3 + term_4)

                M[i, j] =
                    2 / (pi * dad.E) * (term_1 + term_2 + term_3 + term_4) -
                    2 / (pi * dad.E) * dad.al[j][1] * dad.al[j][2] /
                    (sqrt((dad.x[i] - dad.x[j])^2))
            end
        end
    end
    return M
end



function Ponto_prox(dad, dist = 2)
    balltree = BallTree(dad.x')
    prox = Vector{Vector{Int64}}(undef, length(dad.x))
    for i = 1:length(dad.x)
        # Find all points within the tolerance distance from point i    
        prox[i] = inrange(balltree, [dad.x[i]], dad.al[i] * dist)
    end
    prox
end

export calc_K_2d, dadHS2D, dadHS3D, monta_hmat, montaFMM
export calc_K_3d, contact_pressure_force, contact_pressure_disp