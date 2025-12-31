# [Examples](@id examples)

## Thermal

This example demonstrates the application of the Boundary Element Method (BEM) to solve thermal conduction problems. It shows how to model steady-state heat transfer in a domain by discretizing only the boundary, reducing computational complexity compared to domain discretization methods.

### Initialization

Set up the environment and define problem parameters.

```julia
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 20  # Number of elements
order = 2   # Element order
NPX = 10    # Internal points in x direction
NPY = NPX   # Internal points in y direction
npg = 20    # Number of Gauss points (even numbers only)
```

### Data Formatting

Format the input data using the predefined geometry and boundary conditions.

```julia
println("1. Formatting the data")
dad = format_dad(potencial1d(nelem, order), NPX, NPY)  # Data structure
```

### Matrix Assembly

Compute the influence matrices H and G, then apply boundary conditions to form the system matrix A and vector b.

```julia
println("2. Assembling matrix A and vector b")
H, G = calc_HeG(dad, npg)  # Compute H and G matrices
A, b = aplicaCDC(H, G, dad)  # Apply boundary conditions
```

### Solving the System

Solve the linear system and separate the results into temperature and flux.

```julia
println("3. Solving the linear system")
x = A \ b
println("4. Separating temperature and flux")
T, q = separa(dad, x)
```

## Linear Elasticity

This example illustrates the use of BEM for solving linear elasticity problems. It demonstrates how to compute displacements and stresses in elastic solids under various loading conditions, leveraging the boundary-only discretization to efficiently handle 2D or 3D elastic problems.

### Initialization

Set up the environment for the elasticity analysis.

```julia
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 10  # Number of elements
NPX = 2     # Internal points in x direction
NPY = 2     # Internal points in y direction
npg = 20    # Number of Gauss points
tipo = 3    # Element type
```

### Data Formatting

Prepare the data for a beam problem.

```julia
println("1. Formatting the data")
dad = format_dad(viga(nelem, tipo), NPX, NPY)  # Beam data
```

### Matrix Assembly

Assemble the matrices for the elasticity problem.

```julia
println("2. Assembling matrix A and vector b")
H, G = calc_HeG(dad, npg, interno = false)  # Compute matrices (boundary only)
A, b = BEM.aplicaCDC(H, G, dad)  # Apply boundary conditions
```

### Solving the System

Solve for displacements and tractions.

```julia
println("3. Solving the linear system")
x = A \ b
println("4. Separating displacements and tractions")
u, t = separa(dad, x)
```
## Thermal Elasticity

This example demonstrates coupled thermal-elasticity analysis, where temperature fields induce thermal stresses in elastic materials. It combines heat conduction and structural mechanics to solve problems involving thermal expansion and stress generation.

### Benchmark Test Function

The benchmark function `testagao` evaluates the accuracy and performance of the thermal-elasticity solver by comparing numerical results with analytical solutions for varying mesh densities.

```julia
function testagao(i, metodo)
    nelem = i     # Number of elements
    NPY = 2i      # Internal points in y-direction
    NPX = 2NPY    # Internal points in x-direction
    npg = 10      # Number of Gauss points

    ## Data Formatting
    dad = format_dad(telasticogao(nelem, 3), NPX, NPY)

    # Define temperature distribution
    c0 = 0
    c1 = -60
    c2 = 40
    f(x, y) = c2 * y^2 + c1 * y + c0

    # Extract material properties
    EE = dad.k.E
    v = dad.k.nu
    k = dad.k.k

    # Analytical displacement solution
    uana(y) = (c1 / 2 * (y^2 - 0.5^2) + c2 / 3 * (y^3 + 0.5^3) + c0 * (y + 0.5)) * 1.3 / 0.7 * k

    # Solve thermal-elasticity problem
    tdad = @timed u, ss = termoelasticidade(dad, npg, θ = f, metodo = metodo)

    # Post-process results
    y = [dad.NOS[:, 2]; dad.pontos_internos[:, 2]]
    ssa = -k / 0.7 * f.(0, y) * dad.k.E
    ua = uana.(y)
    eu = nrmse(ua, u[:, 2])
    es = nrmse(ssa, ss[:, 1])

    return [nc(dad) ni(dad) eu es tdad.time]
end
```

### Thermal-Elasticity Solver Function

The core function `termoelasticidade` implements the coupled thermal-elastic BEM formulation. It handles temperature-induced deformations and computes displacements and stresses.

```julia
function termoelasticidade(dad, npg; θ = 1, carga = 0, metodo = "dibem")
    # Compute influence matrices
    H, G = calc_HeG(dad, npg, interno = true)
    A, b = BEM.aplicaCDC(H, G, dad)

    # Material properties
    EE = dad.k.E
    v = dad.k.nu
    k = dad.k.k
    kchap = EE * k / (1 - 2 * v)

    # Process temperature input
    np = nc(dad) + ni(dad)
    if θ isa Number
        theta = fill(θ, np)
    elseif θ isa Vector
        theta = θ
    elseif θ isa Function
        theta = [
            θ.(dad.NOS[:, 1], dad.NOS[:, 2])
            θ.(dad.pontos_internos[:, 1], dad.pontos_internos[:, 2])
        ]
    end

    # Compute thermal loads
    F, Fx, Fy = BEM.montaFs([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])
    dF = [zeros(size(Fx)); zeros(size(Fy))]
    for i = 1:np
        dF[2i-1, :] = Fx[i, :]
        dF[2i, :] = Fy[i, :]
    end

    # Thermal expansion terms
    if metodo == "DIBEM"
        Mpe = BEM.Monta_M_RIMd(dad, npg)
        dMpe = BEM.Monta_dM_RIMd(dad, npg)
    else
        Mpe = BEM.Monta_M_RIM(dad, npg)
        dMpe = BEM.Monta_dM_RIM(dad, npg)
    end
    q2 = Mpe * kchap * dF * theta
    dq2 = dMpe * kchap * dF * theta

    # Body forces (if any)
    if carga != 0
        carga_nodal = zeros(2np)
        for i = 1:nc(dad)
            carga_nodal[2i-1:2i] = carga(dad.NOS[i, 1], dad.NOS[i, 2])
        end
        for i = nc(dad)+1:np
            carga_nodal[2i-1:2i] = carga(dad.pontos_internos[i-nc(dad), 1], dad.pontos_internos[i-nc(dad), 2])
        end
        qc = Mpe * carga_nodal
        dqc = dMpe * carga_nodal
    else
        qc = 0
        dqc = 0
    end

    # Assemble thermal load vector
    normal_fonte = dad.normal
    q1 = G * (normal_fonte .* theta[1:nc(dad)])'[:] * kchap

    # Solve system
    x = A \ (b + q1 .- q2 .+ qc)

    # Separate results
    u, t, uint = separa(dad, x)

    # Compute stresses
    S, D = calc_SeD(dad)
    dq1 = D * (normal_fonte .* theta[1:nc(dad)])'[:] * kchap
    dq3 = [kchap * theta kchap * theta 0 * theta]
    fatorcontorno = [fill(2, nc(dad)); ones(ni(dad))]
    sigma = reshape(D * t'[:] - S * u'[:] + dq1 .- dq2 .+ dqc, 3, :)' .* fatorcontorno - dq3

    return [u; uint], sigma
end
```

This implementation supports different solution methods (RIM or DIBEM) and can handle various temperature distributions and mechanical loads.
