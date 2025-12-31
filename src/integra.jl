function mudavar2(x, p = 5)
    t = x^(1 / p) / (x^(1 / p) + (1 - x)^(1 / p))
end
function mudavar1(t, p = 5)
    x = t^(p) / (t^(p) + (1 - t)^(p))
end
function f_prime(t, p = 5)
    numerator = p * (-(t - 1) * t)^(p - 1)
    denominator = (t^p + (1 - t)^p)^2
    return numerator / denominator
end
function corrige_x(x, b)
    for i = 1:size(x, 1)
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
function novelquad(m, t, n, p = 5)
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
        xs = [
            tau .+ h * (1:n) .- h / 2
            tau .+ h / 2 * (1:2n) .- h / 4
            tau .+ h / 4 * (1:4n) .- h / 8
        ]
        ws = [ones(n) * 16h / 7; -ones(2n) * 5h / 7; -ones(4n) * h / 28]
    else
        return 0
    end
    return mudavar1.(corrige_x(xs, b)) * 2 .- 1, ws .* f_prime.(xs) * 2
end


function integraelemsing_num(
    pf,
    nf,
    x,
    elem,
    dad::Union{placa_fina,placa_fina_isotropica},
    pre,
    xi0,
    npg = 30,
)
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

function integra_deriv_elemsing_num(pf, x, elem, dad::Union{elastico}, xi0, npg = 30)
    hx = zeros(Float64, 2, 2 * size(elem))
    gx = zeros(Float64, 2, 2 * size(elem))

    hy = zeros(Float64, 2, 2 * size(elem))
    gy = zeros(Float64, 2, 2 * size(elem))
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
        Ux, Tx, Uy, Ty = calc_deriv_solfund(pg', pf, [sy, -sx], dad, elem.regiao)
        Nm[1, 1:2:end] = N
        Nm[2, 2:2:end] = N

        hx += Tx * Nm * dgamadqsi * w[k]
        gx += Ux * Nm * dgamadqsi * w[k]
        hy += Ty * Nm * dgamadqsi * w[k]
        gy += Uy * Nm * dgamadqsi * w[k]

    end
    hx, gx, hy, gy
end


function integradelemsing_num(pf, x, elem, dad::Union{elastico}, xi0, npg = 30)
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
        D_2 = [
            D[1, 1, 1] D[2, 1, 1]
            D[1, 2, 2] D[2, 2, 2]
            D[1, 1, 2] D[2, 1, 2]
        ]

        S_2 = [
            S[1, 1, 1] S[2, 1, 1]
            S[1, 2, 2] S[2, 2, 2]
            S[1, 1, 2] S[2, 1, 2]
        ]
        # @infiltrate
        d += D_2 * Nm * dgamadqsi * w[k]
        s += S_2 * Nm * dgamadqsi * w[k]

        # @infiltrate
    end
    d, s
end

function integraelemsing_num(pf, x, elem, dad::Union{elastico}, xi0, npg = 30)
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



# -------------------------
# Geometry and shape helpers
# -------------------------
# Basic Lagrange shape functions for arbitrary order nodes (xi_nodes in [0,1])
# Evaluate N_i(ξ) and N_i'(ξ) for i=1..m
export lagrange_basis, lagrange_basis_and_derivative

function lagrange_basis(xi_nodes::AbstractVector{<:Real}, ξ::Real)
    m = length(xi_nodes)
    N = zeros(Float64, m)
    for i in 1:m
        prod = 1.0
        for j in 1:m
            if i != j
                prod *= (ξ - xi_nodes[j])/(xi_nodes[i] - xi_nodes[j])
            end
        end
        N[i] = prod
    end
    return N
end

function lagrange_basis_and_derivative(xi_nodes::AbstractVector{<:Real}, ξ::Real)
    m = length(xi_nodes)
    N = zeros(Float64, m)
    dN = zeros(Float64, m)
    for i in 1:m
        # basis
        Ni = 1.0
        for j in 1:m
            if i != j
                Ni *= (ξ - xi_nodes[j])/(xi_nodes[i] - xi_nodes[j])
            end
        end
        N[i] = Ni
        # derivative via sum over factors
        s = 0.0
        for k in 1:m
            if k == i
                continue
            end
            term = 1.0/(xi_nodes[i] - xi_nodes[k])
            for j in 1:m
                if j != i && j != k
                    term *= (ξ - xi_nodes[j])/(xi_nodes[i] - xi_nodes[j])
                end
            end
            s += term
        end
        dN[i] = s
    end
    return N, dN
end

# Evaluate geometry x(ξ), y(ξ) and derivatives given control node coordinates
# xp, yp are vectors of control coordinates in global order corresponding to xi_nodes
export eval_geometry

function eval_geometry(xp::AbstractVector{<:Real}, yp::AbstractVector{<:Real},
                       xi_nodes::AbstractVector{<:Real}, ξ::Real)
    N, dN = lagrange_basis_and_derivative(xi_nodes, ξ)
    x = dot(N, xp)
    y = dot(N, yp)
    xξ = dot(dN, xp)
    yξ = dot(dN, yp)
    # Jacobian magnitude |J| = sqrt(xξ^2 + yξ^2)
    J = hypot(xξ, yξ)
    return x, y, xξ, yξ, J
end

# -------------------------
# r^2(ξ) helper and derivative (complex allowed)
# -------------------------
export r2_and_derivative

# Evaluate r^2(ξ) = (x(ξ)-xs)^2 + (y(ξ)-ys)^2 and derivative wrt ξ.
# xs, ys can be real or complex (ComplexF64).
function r2_and_derivative(xp, yp, xi_nodes, ξ, xs, ys)
    # xi is possibly complex; we evaluate shape functions with complex arithmetic if necessary
    if isa(ξ, Complex)
        # convert nodes to Complex for basis evaluation
        xi_c = Complex.(xi_nodes)
        m = length(xi_nodes)
        N = zeros(ComplexF64, m)
        dN = zeros(ComplexF64, m)
        for i in 1:m
            Ni = complex(1.0)
            for j in 1:m
                if i != j
                    Ni *= (ξ - xi_c[j])/(xi_c[i] - xi_c[j])
                end
            end
            N[i] = Ni
            s = complex(0.0)
            for k in 1:m
                if k == i; continue; end
                term = 1/(xi_c[i] - xi_c[k])
                for j in 1:m
                    if j != i && j != k
                        term *= (ξ - xi_c[j])/(xi_c[i] - xi_c[j])
                    end
                end
                s += term
            end
            dN[i] = s
        end
        x = sum(N .* Complex.(xp))
        y = sum(N .* Complex.(yp))
        xξ = sum(dN .* Complex.(xp))
        yξ = sum(dN .* Complex.(yp))
        r2 = (x - xs)^2 + (y - ys)^2
        dr2 = 2*(x - xs)*xξ + 2*(y - ys)*yξ
        return r2, dr2
    else
        x, y, xξ, yξ, _ = eval_geometry(xp, yp, xi_nodes, ξ)
        r2 = (x - xs)^2 + (y - ys)^2
        dr2 = 2*(x - xs)*xξ + 2*(y - ys)*yξ
        return r2, dr2
    end
end

# -------------------------
# Appendix A: Newton–Raphson to find complex root ξs = a ± b i s.t. r^2(ξs) = 0
# (see Appendix A and A.2). We do Newton on complex ξ.
# -------------------------
export find_complex_root

function find_complex_root(xp::Vector{<:Real}, yp::Vector{<:Real}, xi_nodes::Vector{<:Real},
                           xs::Real, ys::Real; ξ0::Complex = 0.5 + 0.1im,
                           maxiter::Int=30, tol::Real=1e-12)
    ξ = ξ0
    for k in 1:maxiter
        r2, dr2 = r2_and_derivative(xp, yp, xi_nodes, ξ, complex(xs), complex(ys))
        if abs(r2) < tol
            return ξ, true, k
        end
        # Newton step: ξ_new = ξ - r2 / r2'
        if dr2 == 0
            return ξ, false, k
        end
        ξ = ξ - r2 / dr2
    end
    return ξ, false, maxiter
end

# -------------------------
# Appendix B: log correction for complex quasi-singularity
#    C_Gc = ∫_0^1 ln(w(ξ)) N(ξ) dξ - GL∫ ln(w) N(ξ) dξ  (B.12)
# We compute the exact integral with BigFloat quadrature (many points)
# and subtract the ordinary double-precision Gauss-Legendre result to form the correction.
# -------------------------
export compute_CGc_numeric

function compute_CGc_numeric(shape_node_xi::Vector{<:Real}, Nindex::Int,
                             xp::Vector{<:Real}, yp::Vector{<:Real},
                             xi_s::Complex; ngl::Int=16, high_precision_ng::Int=200)
    # shape_node_xi: natural nodes of element (size = o_e+1)
    # Nindex: which shape function index (1-based) — we integrate N_index(ξ) * ln(w)
    # xi_s is complex root a + i b
    # w(ξ) = (ξ - a)^2 + b^2  (Eq. A.1 + B.11-B.13)
    a = real(xi_s)
    b = imag(xi_s)
    # prepare GL on [0,1]
    ξ_gl, ω_gl = gausslegendre(ngl)
    # prepare high-precision composite quadrature nodes (uniform)
    # We use BigFloat arithmetic for the "exact" integral
    setprecision(256)
    # Evaluate N(ξ) with Float64 for GL and BigFloat for precise integral
    # GL integral (double precision)
    function N_at(ξ)
        Ns = lagrange_basis(shape_node_xi, ξ)
        return Ns[Nindex]
    end
    function lnw(ξ)
        return log((ξ - a)^2 + b^2)
    end
    GL_val = 0.0
    for i in eachindex(ξ_gl)
        GL_val += ω_gl[i] * N_at(ξ_gl[i]) * lnw(ξ_gl[i])
    end

    ξ_hp, ω_hp = gausslegendre(high_precision_ng)
    # Evaluate N(ξ) in BigFloat by evaluating Lagrange polynomials in BigFloat
    function N_at_big(ξb::BigFloat)
        m = length(shape_node_xi)
        xi_b = BigFloat.(shape_node_xi)
        Ns = zeros(BigFloat, m)
        for i in 1:m
            prod = BigFloat(1)
            for j in 1:m
                if i != j
                    prod *= (ξb - xi_b[j])/(xi_b[i] - xi_b[j])
                end
            end
            Ns[i] = prod
        end
        return Ns[Nindex]
    end
    a_b = BigFloat(a); b_b = BigFloat(b)
    GL_big = BigFloat(0)
    for i in eachindex(ξ_hp)
        Nval = N_at_big(ξ_hp[i])
        lnw_big = log((ξ_hp[i] - a_b)^2 + b_b^2)
        GL_big += ω_hp[i] * Nval * lnw_big
    end
    # Correction C_Gc = precise integral - (double-precision GL)
    C_Gc = Float64(GL_big) - GL_val
    return C_Gc
end

# -------------------------
# Appendix C: principal-value / finite-part integrals for real singularities.
# We'll compute PV ∫_0^1 N(ξ)/(ξ - a) dξ numerically using symmetric exclusion
# around a and extrapolating eps -> 0 (adaptive).
# This numerically delivers the Cauchy principal value (Appendix C).
# -------------------------
export pv_integral_over_xi

function pv_integral_over_xi(shape_node_xi::Vector{<:Real}, Nindex::Int, a::Real;
                             tol::Real=1e-12, maxiter::Int=12, ng_per_interval::Int=32)
    # If a is exactly at an endpoint (0 or 1), the PV may be treated via limits — we follow the same exclusion approach.
    # We'll integrate [0, a-eps] and [a+eps, 1] with GL on each subinterval and extrapolate.
    function integrate_excluding_eps(eps)
        val = 0.0
        # left interval
        if a - eps > 0
            ξl, ωl = gausslegendre(ng_per_interval)
            # map to [0, a-eps]
            for (ξi, wi) in zip(ξl, ωl)
                ξm = ( (a - eps)*ξi )
                Nval = lagrange_basis(shape_node_xi, ξm)[Nindex]
                val += wi * Nval / (ξm - a) * (a - eps) # Jacobian scaling
            end
        end
        # right interval
        if a + eps < 1
            ξr, ωr = gausslegendre(ng_per_interval)
            for (ξi, wi) in zip(ξr, ωr)
                # map ξi from [0,1] to [a+eps, 1] -> ξm = a+eps + (1 - (a+eps)) * ξi
                L = 1 - (a + eps)
                ξm = a + eps + L*ξi
                Nval = lagrange_basis(shape_node_xi, ξm)[Nindex]
                val += wi * Nval / (ξm - a) * L
            end
        end
        return val
    end
    # Extrapolate by halving eps until convergence
    eps = 1e-2
    prev = integrate_excluding_eps(eps)
    for k in 1:maxiter
        eps /= 2
        cur = integrate_excluding_eps(eps)
        if abs(cur - prev) < tol
            return cur, true, k
        end
        prev = cur
    end
    return prev, false, maxiter
end

# -------------------------
# Appendix D: Complex quasi-singularity corrections for 1/r^2, 1/r^4, 1/r^6
# Implementation note: we compute the necessary R arrays numerically following the expressions
# in D.1–D.2 by evaluating limits with high precision when needed (use r2_and_derivative with complex ξ).
# -------------------------
export complex_quasi_r_corrections

function complex_quasi_r_corrections(xp, yp, xi_nodes, xi_s::Complex)
    # returns dictionaries of R_G and R_H arrays approximated numerically
    a = real(xi_s); b = imag(xi_s)
    # Evaluate r2 and derivatives at xi_s via analytic continuation (we compute limits)
    # We compute evaluation by small complex offset to estimate needed derivatives numerically
    h = 1e-8 + 1e-8im
    r20, _ = r2_and_derivative(xp, yp, xi_nodes, xi_s, 0.0 + 0.0im, 0.0 + 0.0im)
    # but we need actual arrays from the paper: instead we compute R_G(1..2) and R_H(1..4) numerically
    # We'll compute sample values following the formulas references in the paper (D.4...D.15).
    # For brevity here we produce numerical approximations per each required R by perturbation.
    R = Dict{String,Any}()
    # A simple approach (practical) is to compute directional real/imag parts of
    #  (some combination) — user can refine with closed-form formulas from Appendices D if desired.
    # We'll return placeholders computed from finite differences:
    function numeric_R_expr(exprfun)
        # exprfun is a function ξ -> complex vector or matrix
        δ = 1e-7 + 1e-7im
        val = exprfun(xi_s)
        # no extrapolation; return value
        return val
    end
    # Example: R_G(1) approximated as 2*real(g/(r2')|ξs) per Eq. (24) surrogate
    # User note: these are numerical approximations and should be cross-checked if you require exact
    # closed-form polynomial coefficients from the paper.
    R["RG1"] = numeric_R_expr(ξ -> begin
        # compute g/(r2)' where g = [x^2 x*y; x*y y^2] * N_l (we pick N_l sum weight 1)
        x, y, xξ, yξ, J = eval_geometry(xp, yp, xi_nodes, real(ξ)) # careful: approximate via real ξ
        # fallback numeric complex eval:
        r2, dr2 = r2_and_derivative(xp, yp, xi_nodes, ξ, complex(0.0), complex(0.0))
        # build g using N_l = 1
        g = [x^2 x*y; x*y y^2]
        return 2*real.(g ./ dr2)
    end)
    # For now return this basic set; a full implementation would map all R_G(i) and R_H(i) per D.1, D.2.
    return R
end


# -------------------------
# Simple demonstration / example function
# -------------------------
export demo_example

function demo_example()
    println("CBEM Appendices demo")
    # Example cubic element nodes in xi = [0, 1/3, 2/3, 1] with example control points
    xi_nodes = [0.0, 1/3, 2/3, 1.0]
    xp = [0.0, 0.2, 0.7, 1.0]
    yp = [0.0, 0.5, 0.6, 0.0]
    # Example source point near element
    xs, ys = 0.25, 0.1
    # Initial complex guess for root
    ξ0 = 0.2 + 0.15im
    ξs, conv, iters = find_complex_root(xp, yp, xi_nodes, xs, ys; ξ0=ξ0)
    @printf("Found complex root ξs = %.12f %+.12fi  (conv=%s, iter=%d)\n", real(ξs), imag(ξs), string(conv), iters)

    # Compute numeric log correction for shape function 1
    C_Gc = compute_CGc_numeric(xi_nodes, 1, xp, yp, ξs; ngl=16, high_precision_ng=200)
    @printf("Numeric C_Gc (approx) = %.6e\n", C_Gc)

    # PV integral example for 1/(ξ - a)
    a = 0.35
    pvval, ok, its = pv_integral_over_xi(xi_nodes, 1, a; tol=1e-10)
    @printf("PV ∫ N1/(ξ - %.3f) dξ ≈ %.6e (converged=%s in %d iters)\n", a, pvval, string(ok), its)

    return Dict(
        :ξs => ξs,
        :C_Gc => C_Gc,
        :pvval => pvval,
    )
end
