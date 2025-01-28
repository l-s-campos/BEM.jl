function sladek03(ne = 15, tipo = 2, n_dt = 100)

    POINTS = [
        1 0 0
        2 0.254 0
        3 0.254 0.254
        4 0 0.254
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0

    # #
    # BC_Segments = [1  2
    #                2  2
    #                3  2
    #                4  2];   #simplesmente apoiado
    BC_Segments = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
        'c' 0 0 0 0
    ]
    #engastado
    # Boundary condition of corners (thin plate)
    # cornerbc = 'f': Corners free
    # cornerbc = 'c': Corners clamped

    # type_cornerbc = 'c';
    aprfun = 3 # Type of approximation function
    dA = 32 / 2

    # Domain load
    A = 0
    B = 0
    C = 2.07e6

    # n_dt = 300;
    # dt =10e-3/n_dt;
    # thickness   = .0127;
    # rho = 7.166e3;
    #
    E2 = 0.6895e10
    E1 = 2 * E2
    nu12 = 0.3
    G12 = E2 / (2 * (1 + nu12))
    thickness = 0.0127
    h = thickness
    nu21 = nu12 * E2 / E1
    rho = 7.166e3


    D22 = E2 * h^3 / (12 * (1 - nu12 * nu21))
    a = 0.254
    #-------------------------------------------------------------------------#
    # Termo normalizador
    #-------------------------------------------------------------------------#
    to = a^2 * sqrt(rho * h / D22) / 4
    # to = 1.35e-2        #-------------------------------------------------------------------------#
    # N�mero(n_dt) de intervalos e Passos(dt) de Tempo
    #-------------------------------------------------------------------------#
    # n_dt = 100

    # dt = 0.9*to/n_dt; #simplesmente apoiado
    # dt =0.5*to/n_dt; #engastada
    dt = to / 2n_dt #

    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    Material = [1 E1 E2 G12 nu12 thickness / 2 0]
    dadmat = Compute_Material_Placa(Material)

    R = 1
    wn = [
        25.670
        45.090
        58.648
        71.211
        82.994
        100.929
    ] / a^2 * sqrt(dadmat.D11 / rho)#pg 565 berthelot

    placa_fina,
    POINTS,
    SEGMENTS,
    MESH,
    BC_Segments,
    (;
        dadmat...,
        carga = [A, B, C],
        rho = rho,
        thickness = thickness,
        dt = dt,
        n_dt = n_dt,
        to = to,
        w_st = 6.74e-3,
    )
end
function sladek03_apoiado(ne = 15, tipo = 2, n_dt = 100)

    POINTS = [
        1 0 0
        2 0.254 0
        3 0.254 0.254
        4 0 0.254
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0


    BC_Segments = [
        1 0 0 1 0
        2 0 0 1 0
        3 0 0 1 0
        4 0 0 1 0
        'c' 0 0 0 0
    ]
    Nx = -1
    Ny = 0
    Nxy = 0
    BC_Segments_pl = [
        1 1 Nxy 1 -Ny 2
        2 1 Nx 1 -Nxy 2
        3 1 -Nxy 1 Ny 2
        4 1 -Nx 1 Nxy 2
    ]
    #engastado
    # Boundary condition of corners (thin plate)
    # cornerbc = 'f': Corners free
    # cornerbc = 'c': Corners clamped

    # type_cornerbc = 'c';
    aprfun = 3 # Type of approximation function
    dA = 32 / 2

    # Domain load
    A = 0
    B = 0
    C = 2.07e6

    # n_dt = 300;
    # dt =10e-3/n_dt;
    # thickness   = .0127;
    # rho = 7.166e3;
    #
    E2 = 0.6895e10
    E1 = 2 * E2
    nu12 = 0.3
    G12 = E2 / (2 * (1 + nu12))
    thickness = 0.0127
    h = thickness
    nu21 = nu12 * E2 / E1
    rho = 7.166e3


    D22 = E2 * h^3 / (12 * (1 - nu12 * nu21))
    a = 0.254
    #-------------------------------------------------------------------------#
    # Termo normalizador
    #-------------------------------------------------------------------------#
    to = a^2 * sqrt(rho * h / D22) / 4
    # to = 1.35e-2
    #-------------------------------------------------------------------------#
    # N�mero(n_dt) de intervalos e Passos(dt) de Tempo
    #-------------------------------------------------------------------------#
    # n_dt = 100
    # @infiltrate
    dt = 0.9 * to / n_dt #simplesmente apoiado
    # dt =0.5*to/n_dt; #engastada
    # dt = to / 2n_dt #

    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    Material = [1 E1 E2 G12 nu12 thickness / 2 0]
    dadmat = Compute_Material_Placa(Material)
    # R = 1
    # ωmn = zeros(0)
    # for m in 1:3, n in 1:3
    #         ωmn = [ωmn; pi^2 / a^2 * sqrt(1 / rho) * sqrt(dadmat.D11 * m^4 + (dadmat.D12 + 2dadmat.D66) * m^2 * n^2 + dadmat.D22 * n^4)]
    # end
    # @show ωmn
    placa_fina,
    POINTS,
    SEGMENTS,
    MESH,
    BC_Segments,
    (;
        dadmat...,
        carga = [A, B, C],
        rho = rho,
        thickness = thickness,
        dt = dt,
        n_dt = n_dt,
        to = to,
        w_st = 23.42e-3,
    )
    # [[placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga=[A, B, C], rho=rho, thickness=thickness, dt=dt, n_dt=n_dt, to=to, w_st=23.42e-3)], [elastico_aniso, POINTS, SEGMENTS, MESH, BC_Segments_pl, Compute_Material(Material)]]
end
function large1(ne = 15, tipo = 2, bc = "SSSS")
    a = 1
    theta = 0

    POINTS = [
        1 0 0
        2 a 0
        3 a*(1+sind(theta)) a*cosd(theta)
        4 a*sind(theta) a*cosd(theta)
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0 2
        2 0 0 0 0 2
        3 0 0 0 0 2
        4 0 0 0 0 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    E = 1e7
    nu = 0.3
    h = 1

    D = E * h^3 / (12 * (1 - nu^2))

    A = 0
    B = 0
    C = -40 * h^4 * E
    [
        [
            placa_fina_isotropica,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (h = h, D = D, nu = nu, carga = [A, B, C]),
        ],
        [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E = E, nu = nu)],
    ]

    # Material = [1 1 1.0001 0.5 0.0 0.5 0]

    # dadmat = Compute_Material_Placa(Material)
    # A = 0
    # B = 0
    # C = 0
    # [[placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga=[A, B, C])], [elastico_aniso, POINTS, SEGMENTS, MESH, BC_Segments_pl, Compute_Material(Material)]]
end


function large1_aniso(ne = 15, tipo = 2, bc = "SSSS")
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a 1
        4 0 1
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0 2
        2 0 0 0 0 2
        3 0 0 0 0 2
        4 0 0 0 0 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    E1 = 4.5e6
    h = 1

    A = 0
    B = 0
    C = -100 * h^4 * E2 / a^4


    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    Material = [1 E1 E1 / 3 E1 / 6 0.25 h / 2 0]#pala1
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    # Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end

function Putcha(ne = 15, tipo = 2, bc = "CCCC")
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a a
        4 0 a
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    h = 1

    E1 = 60e6
    E2 = 1.5e6
    # % Ey=1.001e7;
    G12 = .9e6
    # % Gxy=3.7994e+006;
    nu12 = 0.25



    A = 0
    B = 0
    C = -250 * h^4 * E2

    D = E1 * h^3 / (12 * (1 - nu12^2))

    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    Material = [
        1 E1 E2 G12 nu12 h/8 0
        2 E1 E2 G12 nu12 h/8 45
        3 E1 E2 G12 nu12 h/8 -45
        4 E1 E2 G12 nu12 h/8 90
    ]
    dadmat = Compute_Material_Placa(Material)

    # [[placa_fina_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, nu=nu12, carga=[A, B, C])], [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E=E1, nu=nu12)]]

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end


function bhatta2(ne = 15, tipo = 2, bc = "CCCC")
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a 2a
        4 0 2a
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    h = 1

    E1 = 1.2e5
    E2 = 0.6e5
    G12 = 0.07e5
    nu12 = 0.071

    # E2 = 1.199999e5
    # G12 = E1 / 2 / (1 + nu12)

    # E1 = 1.e5
    # E2 = 0.042e5
    # G12 = 0.075e5
    # nu12 = 0.238

    # E1 = 1.e5#alwar
    # E2 = E1 / 4
    # nu12 = 0.3
    # nu21 = nu12 * E2 / E1
    # G12 = 0.35 * E2 / (1 - nu12 * nu21)

    # E1 = 60.3#marczak II
    # E2 = 20.1
    # nu12 = 0.25
    # nu21 = nu12 * E2 / E1
    # G12 = 12.06

    A = 0
    B = 0
    C = -230 * h^4 * E1 / a^4

    D = E1 * h^3 / (12 * (1 - nu12^2))

    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    # [[placa_fina_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, nu=nu12, carga=[A, B, C])], [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E=E1, nu=nu12)]]

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end
function chia(ne = 15, tipo = 2, bc = "CCCC")
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a a
        4 0 a
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    # h = 1
    # E1 = 4.5e+10
    # E2 = 1.5e+10
    # G12 = 7.5e+09
    # nu12 = 0.25

    h = 1
    E1 = 60.3e+09
    E2 = 20.1e+09
    G12 = 12.06e+09
    nu12 = 0.25

    # E1 = 1.0e7
    # E2 = 3 * E1
    # G12 = 0.6 * E1
    # nu12 = 0.25


    A = 0
    B = 0
    C = -60 * h^4 * E2 / a^4

    D = E1 * h^3 / (12 * (1 - nu12^2))

    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    # [[placa_fina_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, nu=nu12, carga=[A, B, C])], [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E=E1, nu=nu12)]]

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end
function ritz(ne = 15, tipo = 2, bc = "CCCC")
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a a
        4 0 a
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    # h = 1
    # E1 = 4.5e+10
    # E2 = 1.5e+10
    # G12 = 7.5e+09
    # nu12 = 0.25

    h = 1
    E1 = 3e6
    E2 = 1.28e6
    G12 = 0.37e6
    nu12 = 0.32

    # E1 = 1.0e7
    # E2 = 3 * E1
    # G12 = 0.6 * E1
    # nu12 = 0.25


    A = 0
    B = 0
    C = -100 * h^4 * E2 / a^4

    D = E1 * h^3 / (12 * (1 - nu12^2))

    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    # [[placa_fina_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, nu=nu12, carga=[A, B, C])], [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E=E1, nu=nu12)]]

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end
function marc(ne = 15, tipo = 2, bc = "CCCC")
    a = 1
    theta = 0

    POINTS = [
        1 0 0
        2 a 0
        3 a*(1+sind(theta)) a*cosd(theta)
        4 a*sind(theta) a*cosd(theta)
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#


    E = 1
    E1 = E

    E2 = E / 4
    nu12 = 0.30
    nu21 = nu12 * E2 / E1
    G12 = E2 / (1 - nu12 * nu21) * 0.35


    # E2 = 1.28
    # nu12 = 0.32
    # G12 = 0.37

    h = 1
    a = 1
    D22 = 0.0213
    Qbarra = 3200

    A = 0
    B = 0
    # C = -Qbarra * E2 * h / (a^4)
    C = -Qbarra * D22 * h / (a^4)


    D = E1 * h^3 / (12 * (1 - nu12^2))
    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    # [[placa_fina_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, nu=nu12, carga=[A, B, C])], [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E=E1, nu=nu12)]]

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end
function bhatta3(ne = 15, tipo = 2, bc = "CCCC")
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a a
        4 0 a
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    h = 1

    # E1 = 1.2e5
    # E2 = 0.06e5
    # G12 = 0.07e5
    # nu12 = 0.071

    # E2 = 1.199999e5
    # G12 = E1 / 2 / (1 + nu12)

    E1 = 1.e5
    E2 = 0.042e5
    G12 = 0.075e5
    nu12 = 0.238

    # E1 = 1.e5#alwar
    # E2 = E1 / 4
    # nu12 = 0.3
    # nu21 = nu12 * E2 / E1
    # G12 = 0.35 * E2 / (1 - nu12 * nu21)

    # E1 = 60.3#marczak II
    # E2 = 20.1
    # nu12 = 0.25
    # nu21 = nu12 * E2 / E1
    # G12 = 12.06

    A = 0
    B = 0
    C = -230 * h^4 * E1 / a^4

    D = E1 * h^3 / (12 * (1 - nu12^2))

    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    # [[placa_fina_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, nu=nu12, carga=[A, B, C])], [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E=E1, nu=nu12)]]

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end

function bhattaiso(ne = 15, tipo = 2, bc = "CCCC")
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a a
        4 0 a
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 1 0 0 0
        2 0 0 1 0
        3 1 0 0 0
        4 0 0 1 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    h = 1

    E1 = 1.2e5
    # E2 = 0.6e5
    E2 = 1.199999e5
    # G12 = 0.07e5
    nu12 = 0.3
    G12 = E1 / 2 / (1 + nu12)

    # E1 = 1.e5
    # E2 = 0.042e5
    # G12 = 0.075e5
    # nu12 = 0.238

    #E1 = 1.e5#alwar
    #E2 = E1 / 4
    #nu12 = 0.3
    #nu21 = nu12 * E2 / E1
    #G12 = 0.35 * E2 / (1 - nu12 * nu21)

    # E1 = 60.3#marczak II
    # E2 = 20.1
    # nu12 = 0.25
    # nu21 = nu12 * E2 / E1
    # G12 = 12.06

    A = 0
    B = 0
    C = -200 * h^4 * E1 / a^4

    D = E1 * h^3 / (12 * (1 - nu12^2))

    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    # Material = [1 E1 E2 G12 nu12 h / 2 0]
    # dadmat = Compute_Material_Placa(Material)

    [
        [
            placa_fina_isotropica,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (h = h, D = D, nu = nu12, carga = [A, B, C]),
        ],
        [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E = E1, nu = nu12)],
    ]

    # [[placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga=[A, B, C])], [elastico_aniso, POINTS, SEGMENTS, MESH, BC_Segments_pl, Compute_Material(Material)]]
end
function pala2(ne = 15, tipo = 2, bc = "SSSS")
    E1 = 60e6
    E2 = 1.5e6
    G12 = .9e6
    nu12 = 0.25
    a = 1
    b = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a b
        4 0 b
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0 2
        2 0 0 0 0 2
        3 0 0 0 0 2
        4 0 0 0 0 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    h = 1
    Az = 0
    Bz = 0
    Cz = -250 * h^4 * E2 / a^4

    Material = [
        1 E1 E2 G12 nu12 h/8 90
        2 E1 E2 G12 nu12 h/8 -45
        3 E1 E2 G12 nu12 h/8 45
        4 E1 E2 G12 nu12 h/8 0
    ]
    # Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    [
        [
            placa_fina,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (; dadmat..., carga = [Az, Bz, Cz]),
        ],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end


function ex1(ne = 15, tipo = 2, bc = "SSSS")
    a = 1
    theta = 0

    POINTS = [
        1 0 0
        2 a 0
        3 a*(1+sind(theta)) a*cosd(theta)
        4 a*sind(theta) a*cosd(theta)
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    Nx = -1
    Ny = 0
    Nxy = -0
    BC_Segments_pl = [
        1 1 Nxy 0 0 2
        2 1 Nx 1 -Nxy 2
        3 1 -Nxy 1 Ny 2
        4 0 0 1 Nxy 2
    ]
    # BC_Segments_pl = [1 1 Nxy 1 Ny 2
    # 2 1 Nx 1 -Nxy 2
    # 3 1 -Nxy 1 Ny 2
    # 4 0 0 1 Nxy 2]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    A = 0
    B = 0
    C = 0
    E = 1
    nu = 0.3
    h = 0.05
    D = E * h^3 / (12 * (1 - nu^2))
    [
        [
            placa_fina_isotropica,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (h = h, D = D, nu = nu, carga = [A, B, C]),
        ],
        [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E = E, nu = nu)],
    ]

    # Material = [1 1 1.0001 0.5 0.0 0.5 0]

    # dadmat = Compute_Material_Placa(Material)
    # A = 0
    # B = 0
    # C = 0
    # [[placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga=[A, B, C])], [elastico_aniso, POINTS, SEGMENTS, MESH, BC_Segments_pl, Compute_Material(Material)]]
end
function vibra(ne = 15, tipo = 2)

    POINTS = [
        1 0 0
        2 0.45 0
        3 0.45 0.35
        4 0 0.35
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0

    # #
    # BC_Segments = [1  2
    #                2  2
    #                3  2
    #                4  2];   #simplesmente apoiado
    BC_Segments = [
        1 0 0 1 0
        2 0 0 1 0
        3 0 0 1 0
        4 0 0 1 0
        'c' 0 0 0 0
    ]
    #engastado
    # Boundary condition of corners (thin plate)
    # cornerbc = 'f': Corners free
    # cornerbc = 'c': Corners clamped

    # type_cornerbc = 'c';
    aprfun = 3 # Type of approximation function
    dA = 32 / 2

    # Domain load
    A = 0
    B = 0
    C = 2.07e6

    # n_dt = 300;
    # dt =10e-3/n_dt;
    # thickness   = .0127;
    # rho = 7.166e3;
    #
    E2 = 12e10
    E1 = 1e10
    nu12 = 0.3
    G12 = 4.8e9
    thickness = 0.002
    h = thickness
    nu21 = nu12 * E2 / E1
    rho = 1510


    # D22 = E2 * h^3 / (12 * (1 - nu12 * nu21))
    a = 0.45
    #-------------------------------------------------------------------------#
    # Termo normalizador
    #-------------------------------------------------------------------------#
    to = a^2 * sqrt(rho * h / D22) / 4
    # to = 1.35e-2
    #-------------------------------------------------------------------------#
    # N�mero(n_dt) de intervalos e Passos(dt) de Tempo
    #-------------------------------------------------------------------------#
    n_dt = 100

    dt = 0.9 * to / n_dt #simplesmente apoiado
    # dt =0.5*to/n_dt; #engastada
    # dt = to / 2n_dt #

    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    Material = [1 E1 E2 G12 nu12 thickness / 2 0]
    dadmat = Compute_Material_Placa(Material)
    R = 450 / 350
    ωmn = zeros(0)
    for m = 1:3, n = 1:3
        ωmn = [
            ωmn
            pi^2 / a^2 *
            sqrt(1 / rho) *
            sqrt(
                dadmat.D11 * m^4 +
                (dadmat.D12 + 2dadmat.D66) * m^2 * n^2 * R^2 +
                dadmat.D22 * n^4 * R^4,
            )
        ]
    end
    @show ωmn
    placa_fina,
    POINTS,
    SEGMENTS,
    MESH,
    BC_Segments,
    (;
        dadmat...,
        carga = [A, B, C],
        rho = rho,
        thickness = thickness,
        dt = dt,
        n_dt = n_dt,
        to = to,
        w_st = 23.42e-3,
    )
end



function placa_com_furo(ne = 15, tipo = 2, bc = "SSSS")

    POINTS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
        5 0.5-b/2 0.5-c/2
        6 0.5+b/2 0.5-c/2
        7 0.5+b/2 0.5+c/2
        8 0.5-b/2 0.5+c/2
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2
        2 2 3
        3 3 4
        4 4 1
        5 5 6
        6 6 7
        7 7 8
        8 8 5
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        5 tipobc[bc[1]]
        6 tipobc[bc[2]]
        7 tipobc[bc[3]]
        8 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    Nx = -1
    Ny = 0
    Nxy = 0
    BC_Segments_pl = [
        1 1 Nxy 1 -Ny 2
        2 1 Nx 1 -Nxy 2
        3 1 -Nxy 1 Ny 2
        4 1 -Nx 1 Nxy 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    Material = [1 1 1.000 0.5 0.0 0.5 0]

    dadmat = Compute_Material_Placa(Material)
    A = 0
    B = 0
    C = 0
    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end

function placa_furo(ne = 15, tipo = 2, bc = "SSSS")
    d = 0.1
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a 1
        4 0 1
        5 0.5-d/2 0.5
        6 0.5+d/2 0.5
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
        5 5 6 -d/2
        6 6 5 -d/2
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        5 livre
        6 livre
        'c' 0 0 0 0
    ]
    Nx = -1
    Ny = 0
    Nxy = 0
    BC_Segments_pl = [
        1 1 Nxy 0 0 2
        2 1 Nx 1 -Nxy 2
        3 1 -Nxy 1 Ny 2
        4 0 0 1 Nxy 2
        5 1 0 1 0 2
        6 1 0 1 0 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    A = 0
    B = 0
    C = 0
    E = 1
    nu = 0.25
    h = 0.05
    D = E * h^3 / (12 * (1 - nu^2))
    [
        [
            placa_fina_isotropica,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (h = h, D = D, nu = nu, carga = [A, B, C]),
        ],
        [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E = E, nu = nu)],
    ]
end

function placa_furo_retangular(ne = 15, tipo = 2, bc = "SSSS")
    b = 0.375
    # b = d
    d = 0.250
    a = 1
    cy = 0.25
    cx = 0.0
    POINTS = [
        1 0 0
        2 a 0
        3 a 1
        4 0 1
        5 cx+0.5-b/2 cy+0.5-d/2
        6 cx+0.5-b/2 cy+0.5+d/2
        7 cx+0.5+b/2 cy+0.5+d/2
        8 cx+0.5+b/2 cy+0.5-d/2
    ]
    # POINTS = [1 -2.5 -2.5
    #         2 2.5 -2.5
    #         3 2.5 2.5
    #         4 -2.5 2.5
    #         5 -0.5 -0.5
    #         6 -0.5 0.5
    #         7 0.5 0.5
    #         8 0.5 -0.5]
    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 5 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        5 livre
        6 livre
        7 livre
        8 livre
        'f' 0 0 0 0
    ]
    Nx = -1
    Ny = 0
    Nxy = -0
    BC_Segments_pl = [
        1 1 Nxy 1 -Ny 2
        2 1 Nx 1 -Nxy 2
        3 1 -Nxy 1 Ny 2
        4 1 -Nx 1 Nxy 2
        5 1 0 1 0 2
        6 1 0 1 0 2
        7 1 0 1 0 2
        8 1 0 1 0 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    A = 0
    B = 0
    C = 0
    E = 1
    nu = 0.25
    h = 0.01
    D = E * h^3 / (12 * (1 - nu^2))
    [
        [
            placa_fina_isotropica,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (h = h, D = D, nu = nu, carga = [A, B, C]),
        ],
        [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E = E, nu = nu)],
    ]
end
function placacomfuro_4(ne = 15, tipo = 2, bc = "SSSS")
    d = 0.1
    a = 1
    # POINTS = [1 0 0
    #         2 a 0
    #         3 a 1
    #         4 0 1
    #         5 0.5-d/2 0.5
    #         6 0.5+d/2 0.5]
    POINTS = [
        1 d/2 0
        2 a 0
        3 a 1
        4 0 1
        5 0 d/2
    ]
    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 5 0
        5 5 1 -d
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        5 livre
        'c' 0 0 0 0
    ]
    Nx = -1
    Ny = 0
    Nxy = 0
    BC_Segments_pl = [
        1 1 Nxy 0 0 2
        2 1 Nx 1 -Nxy 2
        3 1 -Nxy 1 Ny 2
        4 0 0 1 Nxy 2
        4 1 0 1 0 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    A = 0
    B = 0
    C = 0
    E = 1
    nu = 0.25
    h = 0.05
    D = E * h^3 / (12 * (1 - nu^2))
    [
        [
            placa_fina_isotropica,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (h = h, D = D, nu = nu, carga = [A, B, C]),
        ],
        [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E = E, nu = nu)],
    ]
end
function redonda(ne = 15, tipo = 2, bc = "SS")
    d = 2
    POINTS = [
        1 -1 0
        2 0 -1
        3 1 0
        4 0 1
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 3 d/2
        2 3 1 d/2
    ]

    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        'c' 0 0 0 0
    ]
    Nx = -1
    Ny = 0
    Nxy = 0
    BC_Segments_pl = [
        1 3 Nx 0 0 2
        2 1 Nx 1 -Nxy 2
        3 1 Nx 1 -Nxy 2
        4 1 Nx 1 -Nxy 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    A = 0
    B = 0
    C = 0
    E = 1
    nu = 0.25
    h = 0.05
    D = E * h^3 / (12 * (1 - nu^2))
    [
        [
            placa_fina_isotropica,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (h = h, D = D, nu = nu, carga = [A, B, C]),
        ],
        [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E = E, nu = nu)],
    ]
end

function placa_furo_retangular4(ne = 15, tipo = 2)
    d = 0.1
    a = 1
    POINTS = [
        1 d/2 0
        2 a/2 0
        3 a/2 0.5
        4 0 0.5
        5 0 d/2
        6 d/2 d/2
    ]
    # POINTS = [1 -2.5 -2.5
    #         2 2.5 -2.5
    #         3 2.5 2.5
    #         4 -2.5 2.5
    #         5 -0.5 -0.5
    #         6 -0.5 0.5
    #         7 0.5 0.5
    #         8 0.5 -0.5]
    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 5 0
        5 5 6 0
        6 6 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 1 0 0 0
        2 apoiada
        3 apoiada
        4 1 0 0 0
        5 livre
        6 livre
        'f' 0 0 0 0
    ]
    Nx = -1
    Ny = 0
    Nxy = 0
    BC_Segments_pl = [
        1 1 Nxy 0 0 2
        2 1 Nx 1 -Nxy 2
        3 1 -Nxy 1 Ny 2
        4 0 0 1 Nxy 2
        5 1 0 1 0 2
        6 1 0 1 0 2
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#

    A = 0
    B = 0
    C = 0
    E = 1
    nu = 0.0
    h = 1
    D = E * h^3 / (12 * (1 - nu^2))
    [
        [
            placa_fina_isotropica,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments,
            (h = h, D = D, nu = nu, carga = [A, B, C]),
        ],
        [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E = E, nu = nu)],
    ]
end

function reddy(ne = 15, tipo = 2, bc = "CCCC")
    a = 12
    POINTS = [
        1 0 0
        2 a 0
        3 a a
        4 0 a
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    # BC_Segments_pl = [1 1 0 0 0
    #         2 0 0 1 0
    #         3 1 0 0 0
    #         4 0 0 1 0]
    # BC_Segments_pl = [1 0 0 1 0
    #         2 1 0 0 0
    #         3 0 0 1 0
    #         4 1 0 0 0]
    BC_Segments_pl = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#


    E1 = 3.e6
    E2 = 1.28e6
    nu12 = 0.32
    h = 0.138
    G12 = 0.37e6

    C = -2.0
    A = 0
    B = 0

    Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end
function bhatta3(ne = 15, tipo = 2, bc = "CCCC")
    a = 1
    POINTS = [
        1 0 0
        2 a 0
        3 a a
        4 0 a
    ]

    # Lines (how the points are joined)
    #
    SEGMENTS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]


    # Discretization (elements per line)

    MESH = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    #  Boundary conditions of lines
    # BC = [type value type value type value type value]
    # Se a CDC imposta for borda livre -> 1 0 1 0
    # Se a CDC imposta for borda engastada -> 0 0 0 0
    # Se a CDC imposta for borda apoiada -> 0 0 1 0
    livre = [1 0 1 0]
    engastada = [0 0 0 0]
    apoiada = [0 0 1 0]
    tipobc = Dict('F' => livre, 'C' => engastada, 'S' => apoiada)
    BC_Segments = [
        1 tipobc[bc[1]]
        2 tipobc[bc[2]]
        3 tipobc[bc[3]]
        4 tipobc[bc[4]]
        'c' 0 0 0 0
    ]
    BC_Segments_pl = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
    ]
    #-------------------------------------------------------------------------#
    # Matriz de propriedades do material
    #-------------------------------------------------------------------------#
    h = 1

    # E1 = 1.2e5
    # E2 = 0.06e5
    # G12 = 0.07e5
    # nu12 = 0.071

    # E2 = 1.199999e5
    # G12 = E1 / 2 / (1 + nu12)

    E1 = 1.e5
    E2 = 0.042e5
    G12 = 0.075e5
    nu12 = 0.238

    # E1 = 1.e5#alwar
    # E2 = E1 / 4
    # nu12 = 0.3
    # nu21 = nu12 * E2 / E1
    # G12 = 0.35 * E2 / (1 - nu12 * nu21)

    # E1 = 60.3#marczak II
    # E2 = 20.1
    # nu12 = 0.25
    # nu21 = nu12 * E2 / E1
    # G12 = 12.06

    A = 0
    B = 0
    C = -230 * h^4 * E1 / a^4

    D = E1 * h^3 / (12 * (1 - nu12^2))

    # Material = [1 E1 1.0001 * E1 E1 / 2 / (1 + 0.3) 0.3 h / 2 0]
    # Material = [1 E1 4 * E1 0.3581 * E1 0.3 h / 2 0]

    # Material = [1 E1 E1 / 3 E1 / 5 0.25 h / 2 0]

    Material = [1 E1 E2 G12 nu12 h / 2 0]
    dadmat = Compute_Material_Placa(Material)

    # [[placa_fina_isotropica, POINTS, SEGMENTS, MESH, BC_Segments, (h=h, D=D, nu=nu12, carga=[A, B, C])], [elastico, POINTS, SEGMENTS, MESH, BC_Segments_pl, (E=E1, nu=nu12)]]

    [
        [placa_fina, POINTS, SEGMENTS, MESH, BC_Segments, (; dadmat..., carga = [A, B, C])],
        [
            elastico_aniso,
            POINTS,
            SEGMENTS,
            MESH,
            BC_Segments_pl,
            Compute_Material(Material),
        ],
    ]
end
