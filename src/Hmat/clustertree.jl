"""
    mutable struct ClusterTree{N,T}

Tree structure used to cluster poitns of type `SVector{N,T}` into [`HyperRectangle`](@ref)s.

# Fields:
- `_elements::Vector{SVector{N,T}}` : vector containing the sorted elements.
- `container::HyperRectangle{N,T}` : container for the elements in the current node.
- `index_range::UnitRange{Int}` : indices of elements contained in the current node.
- `loc2glob::Vector{Int}` : permutation from the local indexing system to the
  original (global) indexing system used as input in the construction of the
  tree.
- `glob2loc::Vector{Int}` : inverse of `loc2glob` permutation.
- `children::Vector{ClusterTree{N,T}}`
- `parent::ClusterTree{N,T}`
"""
mutable struct ClusterTree{N,T}
    _elements::Vector{SVector{N,T}}
    container::HyperRectangle{N,T}
    index_range::UnitRange{Int}
    loc2glob::Vector{Int}
    glob2loc::Vector{Int}
    children::Vector{ClusterTree{N,T}}
    parentnode::ClusterTree{N,T}
    depth::Int
    Xstar::Vector{Int}
    function ClusterTree(
        els::Vector{SVector{N,T}},
        container,
        loc_idxs,
        loc2glob,
        glob2loc,
        children,
        parentnode,
        depth::Int = 0,
        Xstar = Vector{Int}(undef, 0),
    ) where {N,T}
        clt = new{N,T}(els, container, loc_idxs, loc2glob, glob2loc)
        clt.children = isnothing(children) ? Vector{typeof(clt)}() : children
        clt.parentnode = isnothing(parentnode) ? clt : parentnode
        clt.Xstar = Xstar
        clt.depth = depth
        return clt
    end
end

# convenience functions and getters
"""
    root_elements(clt::ClusterTree)

The elements contained in the root of the tree to which `clt` belongs.
"""
root_elements(clt::ClusterTree) = clt._elements

"""
    index_range(clt::ClusterTree)

Indices of elements in `root_elements(clt)` which lie inside `clt`.
"""
index_range(clt::ClusterTree) = clt.index_range

children(clt::ClusterTree) = clt.children

"""
    parentnode(clt::ClusterTree)

The node's parent. If `t` is a root, then `parent(t)==t`.
"""
parentnode(clt::ClusterTree) = clt.parentnode

"""
    container(clt::ClusterTree)

Return the object enclosing all the elements of the `clt`.
"""
container(clt::ClusterTree) = clt.container

"""
    elements(clt::ClusterTree)

Iterable list of the elements inside `clt`.
"""
elements(clt::ClusterTree) = view(root_elements(clt), index_range(clt))

"""
    loc2glob(clt::ClusterTree)

The permutation from the (local) indexing system of the elements of the `clt` to
the (global) indexes used upon the construction of the tree.
"""
loc2glob(clt::ClusterTree) = clt.loc2glob

"""
    glob2loc(clt::ClusterTree)

The inverse of [`loc2glob`](@ref).
"""
glob2loc(clt::ClusterTree) = clt.glob2loc

isleaf(clt::ClusterTree) = isempty(clt.children)
isroot(clt::ClusterTree) = parentnode(clt) == clt

depth(clt::ClusterTree) = clt.depth

diameter(node::ClusterTree) = diameter(container(node))
radius(node::ClusterTree) = radius(container(node))

"""
    distance(X::ClusterTree, Y::ClusterTree)

Distance between the containers of `X` and `Y`.
"""
function distance(node1::ClusterTree, node2::ClusterTree)
    return distance(container(node1), container(node2))
end

Base.length(node::ClusterTree) = length(index_range(node))

"""
    ClusterTree(elements,splitter;[copy_elements=true, threads=false])

Construct a `ClusterTree` from the  given `elements` using the splitting
strategy encoded in `splitter`. If `copy_elements` is set to false, the
`elements` argument are directly stored in the `ClusterTree` and are permuted
during the tree construction.
"""
function ClusterTree(
    elements,
    splitter = GeometricSplitter();
    copy_elements = true,
    threads = false,
)
    copy_elements && (elements = deepcopy(elements))
    bbox = bounding_box(elements)
    n = length(elements)
    irange = 1:n
    loc2glob = collect(irange)
    glob2loc = collect(irange) # used as buffer during the construction
    children = nothing
    parent = nothing
    Xstar = Vector{Int}(undef, 0)
    #build the root, then recurse
    root =
        ClusterTree(elements, bbox, irange, loc2glob, glob2loc, children, parent, 0, Xstar)
    _build_cluster_tree!(root, splitter, threads)
    # inverse the loc2glob permutation
    glob2loc .= invperm(loc2glob)
    # finally, permute the elements so as to use the local indexing
    copy!(elements, elements[loc2glob]) # faster than permute!
    # compute selected points after tree is built
    # compute_selected!(root)
    return root
end

function _build_cluster_tree!(current_node, splitter, threads, depth = 0)
    if should_split(current_node, depth, splitter)
        split!(current_node, splitter)
        if threads
            Threads.@threads for child in children(current_node)
                _build_cluster_tree!(child, splitter, threads, depth + 1)
            end
        else
            for child in children(current_node)
                _build_cluster_tree!(child, splitter, threads, depth + 1)
            end
        end
    end
    return current_node
end

function select_leaves!(tree; m = 10)

    dv = [i.depth for i in nodes(tree)]
    I = sortperm(-dv)
    nodes_ = nodes(tree)

    for i in I
        if isleaf(nodes_[i])
            filter_selected_leaf!(nodes_[i], m)
        else
            filter_selected_node!(nodes_[i], m)
        end
    end


end

function filter_selected_leaf!(current_node, m)
    els = index_range(current_node)
    points = current_node.loc2glob[els]
    S, AX = anchor_net_with_selection(root_elements(current_node)[points], m)
    current_node.Xstar = S
end
function filter_selected_node!(current_node, m)
    select = Int[]
    for node in current_node.children
        # @show node.depth, index_range(node), [node.selected_indices]
        append!(select, index_range(node)[node.Xstar])
    end

    # els = index_range(current_node)
    points = current_node.loc2glob[select]
    S, AX = anchor_net_with_selection(root_elements(current_node)[points], m)
    current_node.Xstar = select[S] .- index_range(current_node)[1] .+ 1
    # if isroot(current_node)
    #     return
    # end
    # if isempty(parentnode(current_node).selected_indices)
    #     current_node.parentnode.selected_indices = deepcopy(current_node.selected_indices)
    # else
    #     current_node.parentnode.selected_indices =
    #         [current_node.parentnode.selected_indices; current_node.selected_indices]
    #     # filter_selected_node!(parentnode(current_node), m)
    # end
end

function compute_selected!(current_node)


    # if isleaf(current_node)
    #     # els = elements(current_node)
    #     els = index_range(current_node)
    #     n = length(els)
    #     if n > 0
    #         # X = hcat(els...)'
    #         m = min(10, n)  # number of anchor points
    #         points = current_node._elements[els]
    #         S, AX = anchor_net_with_selection(points, m)
    #         # S are local indices in X (1 to n), map to global indices
    #         current_node.selected_indices = current_node.loc2glob[els[S]]
    #     end
    # else
    #     # First, compute for children
    #     for child in children(current_node)
    #         compute_selected!(child)
    #     end
    #     # Then, run on union of selected from children
    #     all_selected = Int[]
    #     for child in children(current_node)
    #         if !isnothing(child.selected_indices)
    #             append!(all_selected, child.selected_indices)
    #         end
    #     end
    #     # unique_selected = unique(all_selected)
    #     if !isempty(all_selected)
    #         points = current_node._elements[all_selected]
    #         n = length(points)
    #         m = min(10, n)
    #         S, AX = anchor_net_with_selection(points, m)
    #         current_node.selected_indices = all_selected[S]
    #     end
    # end
end

function Base.show(io::IO, tree::ClusterTree{N,T}) where {N,T}
    return print(io, "ClusterTree with $(length(tree.index_range)) points.")
end

function Base.summary(clt::ClusterTree)
    @printf "Cluster tree with %i elements" length(clt)
    nodes_ = nodes(clt)
    @printf "\n\t number of nodes: %i" length(nodes_)
    leaves_ = leaves(clt)
    @printf "\n\t number of leaves: %i" length(leaves_)
    points_per_leaf = map(length, leaves_)
    @printf "\n\t min number of elements per leaf: %i" minimum(points_per_leaf)
    @printf "\n\t max number of elements per leaf: %i" maximum(points_per_leaf)
    depth_per_leaf = map(depth, leaves_)
    @printf "\n\t min depth of leaves: %i" minimum(depth_per_leaf)
    @printf "\n\t max depth of leaves: %i" maximum(depth_per_leaf)
end
function ClusterTree(
    dad::potencial,
    splitter = CardinalitySplitter();
    threads = false,
    dist = 1e3,
)
    elements = [
        [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
        [Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for i = 1:ni(dad)]
    ]
    tipoCDC = BEM.tipoCDC(dad)
    elementsCDC = [
        [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
        [
            Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) .+ dist^2 for
            i = 1:ni(dad)
        ]
    ]

    # for i = 1:length(tipoCDC)
    #     if tipoCDC[i] == 1
    #         elementsCDC[i] = elementsCDC[i] .+ dist
    #     end
    # end
    for i = 1:length(elementsCDC)
        elementsCDC[i] += dist * [1, 1] .* .!tipoCDC[i, :]
    end


    yclt = ClusterTree(elementsCDC, splitter)
    corrige_yclt!(yclt, elements, dist)
    # @infiltrate
    # return ClusterTree(elements, splitter), ClusterTree(elements, splitter)
    return ClusterTree(elements, splitter), yclt
end
function corrige_yclt!(current_node, elements, dist)
    current_node._elements = elements
    irange = index_range(current_node)
    l2g = loc2glob(current_node)
    l2g[irange]
    current_node.container = bounding_box(elements[l2g[irange]])
    for child in children(current_node)
        # @infiltrate
        corrige_yclt!(child, elements, dist)
    end
end
function ClusterTree(
    dad::elastico,
    splitter = CardinalitySplitter();
    threads = false,
    dist = 1e6,
)
    elements = repeat(
        [
            [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
            [
                Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for
                i = 1:ni(dad)
            ]
        ],
        inner = 2,
    )
    tipoCDC = BEM.tipoCDC(dad)
    elementsCDC = repeat(
        [
            [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
            [
                Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for
                i = 1:ni(dad)
            ]
        ],
        inner = 2,
    )

    # for i = 1:length(tipoCDC)
    #     if tipoCDC[i] == 1
    #         elementsCDC[i] = elementsCDC[i] .+ dist
    #     end
    # end
    # for i = 1:nc(dad)
    #     elementsCDC[2i-1] .+= dist *  tipoCDC[2i-1]
    #     elementsCDC[2i] .+= dist *  tipoCDC[2i]
    # end
    elementsCDC[1:2:end] .+= dist * tipoCDC[1:2:end]
    elementsCDC[2:2:end] .-= dist * tipoCDC[2:2:end]
    yclt = ClusterTree(elementsCDC, splitter)
    corrige_yclt!(yclt, elements, dist)
    # @infiltrate
    # return ClusterTree(elements, splitter), ClusterTree(elements, splitter)
    return ClusterTree(elements, splitter), yclt
end