function calc_K(x0,xf,n;E=200e9)
y= range(x0, stop=xf, length=n+1)
x=(y[1:end-1]+y[2:end])/2
node=hcat(collect(1:n),collect(2:n+1)) # Element connectivity
al = abs.(y[ node[ :,2]] .- y[ node[ :,1]])# Length of the element

    A = zeros(n, n)  # Matrix A
for i = 1:n
    for j = 1:n # Loop on field points (Column)
        # Compute parameters used in the formulas for the two integrals
        x1 = y[node[j,1]] .- x[i]  # Distance from the source point to the beginning of the element
        x2 = y[node[j,2]] .- x[i]
        r1 = abs.(x1) # Distance from the source point to the beginning of the element
        r2 = abs.(x2)# Distance from the source point to the end of the element

        A[i, j]  = ( al[j] .+ x1 .* log.(r1) .- x2 .* log.(r2)) * -4/(pi*E)

    end
  end
    return A
end

@kwdef mutable struct kernelHS2D <: AbstractMatrix{Float64}
y::Vector{Float64}
x::Vector{Float64}
node=::Matrix{Int64}
al ::Vector{Float64}
E::Float64=200e9
end

function Base.getindex(K::kernelHS2D, i::Int, j::Int)
        x1 = K.y[K.node[j,1]] .- K.x[i]  # Distance from the source point to the beginning of the element
        x2 = K.y[node[j,2]] .- K.x[i]
        r1 = abs.(x1) # Distance from the source point to the beginning of the element
        r2 = abs.(x2)# Distance from the source point to the end of the element

        ( K.al[j] .+ x1 .* log.(r1) .- x2 .* log.(r2)) * -4/(pi*K.E)
end

