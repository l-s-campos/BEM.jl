using HaltonSequences
using CairoMakie
using StaticArrays

# https://arxiv.org/pdf/2206.01885
# https://arxiv.org/pdf/2102.05215
# https://arxiv.org/pdf/2212.12674
# https://arxiv.org/pdf/1705.05443

function anchor_net_with_selection(X::Vector{SVector{2,Float64}}, m::Int)
    n = length(X)
    d = 2
    # Step 1: Construct the low discrepancy set T in bounding box B0
    mins = reduce((a, b) -> min.(a, b), X)
    maxs = reduce((a, b) -> max.(a, b), X)
    B0 = (mins, maxs)
    s = m # number of anchor base points

    # Generate s low discrepancy points in B0. Example: Halton sequence approximation
    T = low_discrepancy_points(s, d, B0)

    # Step 2: Initialize groups Gi
    G = [Int[] for _ = 1:s]

    # Step 3-6: Assign each xj in X to closest anchor tk by infinity norm
    for j = 1:n
        xj = X[j]
        dists = zeros(Float64, s)
        for i = 1:s
            dists[i] = maximum(abs.(xj .- T[i]))
        end
        i_closest = argmin(dists)
        push!(G[i_closest], j)
    end

    # Step 7-10: Find nonempty Gi's, their boxes, and measures
    nonempty_indices = filter(i -> !isempty(G[i]), 1:s)
    Q = length(nonempty_indices)
    B_boxes = Vector{Tuple{AbstractVector{Float64},AbstractVector{Float64}}}(undef, Q)
    λ = zeros(Float64, Q)
    for q = 1:Q
        idx = nonempty_indices[q]
        points_q = X[G[idx]]
        mins_q = reduce((a, b) -> min.(a, b), points_q)
        maxs_q = reduce((a, b) -> max.(a, b), points_q)
        B_boxes[q] = (mins_q, maxs_q)
        λ[q] = prod(maxs_q .- mins_q)
    end

    # Step 11: Choose anchor points in each box in proportion to Lebesgue measures
    sum_λ = sum(λ)
    AX = Vector{SVector{2,Float64}}()
    for q = 1:Q
        n_points_q = round(Int, m * λ[q] / (sum_λ + eps()))
        AX_parts = low_discrepancy_points(n_points_q, d, B_boxes[q])
        append!(AX, AX_parts)
    end
    # Step 12: Full anchor net AX
    s_ax = length(AX)

    # Preallocate S for indices of closest points in X to anchor points
    S = zeros(Int, s_ax)
    for i = 1:s_ax
        y = AX[i]
        dists_x = zeros(Float64, n)
        for j = 1:n
            dists_x[j] = maximum(abs.(X[j] .- y))
        end
        S[i] = argmin(dists_x)
    end

    # Step 11: Choose anchor points in each box in proportion to Lebesgue measures
    sum_λ = sum(λ)
    AX = Vector{Vector{Float64}}(undef, 0)
    for q = 1:Q
        n_points_q = round(Int, m * λ[q] / (sum_λ + eps()))
        AX_parts = low_discrepancy_points(n_points_q, d, B_boxes[q])
        for p in AX_parts
            push!(AX, p)
        end
    end
    # Step 12: Full anchor net AX
    s_ax = length(AX)

    # Preallocate S for indices of closest points in X to anchor points
    S = zeros(Int, s_ax)
    for i = 1:s_ax
        y = AX[i]
        dists_x = zeros(Float64, n)
        for j = 1:n
            dists_x[j] = maximum(abs.(X[j] .- y))
        end
        S[i] = argmin(dists_x)
    end

    return S, AX
end

# Supporting function for low discrepancy points (Halton sequence)
function low_discrepancy_points(
    n::Int,
    d::Int,
    box::Tuple{AbstractVector{Float64},AbstractVector{Float64}},
)
    mins, maxs = box
    halton_seq = HaltonPoint(d, length = n)
    scaled_points = Vector{SVector{d,Float64}}(undef, n)
    for i = 1:n
        scaled_points[i] = halton_seq[i] .* (maxs - mins) .+ mins
    end
    return scaled_points
end

# Example usage
function example_smash()
    # Generate random 2D points
    n = 200
    X = [SVector{2,Float64}(rand() * 10, rand() * 10) for _ = 1:n]  # points in [0,10] x [0,10]
    m = 20  # number of anchor points

    # Run the anchor net selection
    S, AX = anchor_net_with_selection(X, m)

    # Plot
    fig = Figure()
    ax = Axis(fig[1, 1], title = "Anchor Net Selection Example", xlabel = "X", ylabel = "Y")

    # Plot all points
    scatter!(
        ax,
        [p[1] for p in X],
        [p[2] for p in X],
        color = :blue,
        label = "Data Points",
        markersize = 4,
    )

    # Plot anchor points
    ax_points = [p[1] for p in AX]
    ay_points = [p[2] for p in AX]
    scatter!(
        ax,
        ax_points,
        ay_points,
        color = :red,
        label = "Anchor Points",
        markersize = 8,
        marker = :star5,
    )

    # Plot selected points (closest to anchors)
    scatter!(
        ax,
        [X[i][1] for i in S],
        [X[i][2] for i in S],
        color = :green,
        label = "Closest Data Points to Anchors",
        markersize = 6,
        marker = :circle,
    )

    axislegend(ax)
    save("example_smash_plot.png", fig)
    println("Plot saved to example_smash_plot.png")
end

