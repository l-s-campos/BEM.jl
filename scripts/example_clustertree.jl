using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
# Generate random 2D points
n = 1000
points = [SVector{2,Float64}(rand(), rand()) for _ = 1:n]

# Build ClusterTree
tree = ClusterTree(points, CardinalitySplitter(50))
BEM.select_leaves!(tree, m = 10)

println("Fields: ", fieldnames(typeof(tree)))

println("Root selected: ", tree.selected_indices)

# Function to print selected for nodes
function print_selected(node, level = 0)
    indent = "  "^level
    if isempty(BEM.children(node))
        println(
            "$(indent)Leaf with $(length(node)) points, selected: $(length(node.selected_indices))",
        )
    else
        println(
            "$(indent)Node with $(length(node)) points, selected: $(length(node.selected_indices))",
        )
        for child in BEM.children(node)
            print_selected(child, level + 1)
        end
    end
end
print_selected(tree)

println("Tree structure with selected points:")

