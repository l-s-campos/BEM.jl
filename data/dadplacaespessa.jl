function plespex1(ne=15, tipo=2, bc="SSSS")
    a = 1
    theta = 0

    POINTS = [1 0 0
        2 a 0
        3 a*(1+sind(theta)) a*cosd(theta)
        4 a*sind(theta) a*cosd(theta)]

    # Lines (how the points are joined)
    # 
    SEGMENTS = [1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0]


    # Discretization (elements per line)

    MESH = [1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0 1 0]
    engastada = [0 0 0 0 0 0]
    apoiada = [1 0 1 0 0 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    E = 200e9
    nu = 0.30
    h = 0.01
    # λ = sqrt(10) / h
    λ = 5 / 6

    D = E * h^3 / (12 * (1 - nu^2))


    A = 0
    B = 0
    C = 1
    placa_espessa_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, λ=λ, nu=nu, carga=[A, B, C])
end

function termbuckl(ne=15, tipo=2, bc="SSSS")
    a = 1
    theta = 0

    POINTS = [1 0 0
        2 a 0
        3 a*(1+sind(theta)) a*cosd(theta)
        4 a*sind(theta) a*cosd(theta)]

    # Lines (how the points are joined)
    # 
    SEGMENTS = [1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0]


    # Discretization (elements per line)

    MESH = [1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0 1 0]
    engastada = [0 0 0 0 0 0]
    apoiada = [1 0 1 0 0 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    CCSeg = [1 1 0 0 0
        2 0 0 1 0
        3 1 0 0 0
        4 0 0 1 0]

    E = 1
    nu = 0.30
    h = 0.01
    # λ = sqrt(10) / h
    λ = 5 / 6

    D = E * h^3 / (12 * (1 - nu^2))


    A = 0
    B = 0
    C = 0
    [[placa_espessa_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, λ=λ, nu=nu, carga=[A, B, C])], [elastico, POINTS, SEGMENTS, MESH, CCSeg, (E=E, nu=nu)]]
end