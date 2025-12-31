# [Data](@id data)

This section explains how to enter data for BEM simulations in the BEM.jl package. Data entry involves defining the geometry, boundary conditions, material properties, and discretization parameters. The following example function `potencial1d` demonstrates the data structure for a 1D thermal potential problem.

## Geometry Definition

The geometry is defined using points and segments. Points are specified in a matrix where each row contains the point number and its coordinates.

```julia
function potencial1d(ne = 15, tipo = 2)
    # Define the points of the geometry
    # PONTOS = [point_number x y]
    PONTOS = [
        1 0 0  # Point 1 at (0, 0)
        2 1 0  # Point 2 at (1, 0)
        3 1 1  # Point 3 at (1, 1)
        4 0 1  # Point 4 at (0, 1)
    ]
```

## Segments

Segments connect the points and define the boundaries. Each segment is defined by its number, start and end points, and curvature (radius). A radius of 0 indicates a straight line.

```julia
    # Define segments: [segment_number start_point end_point radius]
    # Radius: 0 = straight line, >0 = left curve, <0 = right curve
    SEGMENTOS = [
        1 1 2 0  # Straight segment from point 1 to 2
        2 2 3 0  # Straight segment from point 2 to 3
        3 3 4 0  # Straight segment from point 3 to 4
        4 4 1 0  # Straight segment from point 4 to 1
    ]
```

## Mesh Discretization

The mesh specifies how many elements to place on each segment and their type.

```julia
    # Define mesh: [segment_number number_of_elements element_type]
    MALHA = [
        1 ne tipo  # Segment 1: ne elements of type tipo
        2 ne tipo  # Segment 2: ne elements of type tipo
        3 ne tipo  # Segment 3: ne elements of type tipo
        4 ne tipo  # Segment 4: ne elements of type tipo
    ]
```

## Boundary Conditions

Boundary conditions are applied to each segment, specifying whether temperature or heat flux is prescribed.

```julia
    # Define boundary conditions: [segment_number condition_type condition_value]
    # condition_type: 0 = prescribed temperature, 1 = prescribed flux
    CCSeg = [
        1 1 0    # Segment 1: flux = 0
        2 1 -1   # Segment 2: flux = -1
        3 1 0    # Segment 3: flux = 0
        4 0 0    # Segment 4: temperature = 0
    ]
```

## Material Properties

Material properties such as thermal conductivity are defined as scalars or functions.

```julia
    # Thermal conductivity of the material
    k = 1
```

The function returns the data structures needed for the BEM solver.

```julia
    # Return the data for the potential problem
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end
```
