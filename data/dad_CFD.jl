function Taylor_Couette(ne = 15, tipo = 2, ra = 1, rb = 2, Re = 50, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 ra 0
        2 rb 0
        3 0 rb
        4 0 ra
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 rb
        3 3 4 0
        4 4 1 -ra
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        1 0 0 1 0
        2 0 0 0 0
        3 1 0 0 0
        4 0 0 0 0
    ]

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function couette(ne = 15, tipo = 2, L = 1, Re = 10, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 L 0
        3 L L
        4 0 L
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        1 0 0 0 0
        2 1 0 0 0
        3 0 1 0 0
        4 1 0 0 0
    ]
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function cavity(ne = 15, tipo = 2, L = 1, Re = 1, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = μ * 2 * (1 + v)

    PONTOS = [
        1 0 0
        2 L 0
        3 L L
        4 0 L
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 1 0 0
        4 0 0 0 0
    ]

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function step(
    ne = 15,
    tipo = 2,
    L = 1,
    h = 0.4,
    L_step = 0.2,
    x_step = 0.2,
    Re = 1,
    λ = 10^6,
)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 x_step 0
        3 x_step L_step
        4 x_step+L_step L_step
        5 x_step+L_step 0
        6 L 0
        7 L h
        8 0 h
    ]

    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 5 0
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 1 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]
    CCSeg = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
        5 0 0 0 0
        6 1 0 0 0
        7 0 0 0 0
        8 0 1 0 0
    ]

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function expansion(ne = 15, tipo = 2, L1 = 1, h1 = 0.4, L2 = 1, h2 = 0.8, Re = 1, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 L1 0
        3 L1 -(h2 - h1)/2
        4 L1+L2 -(h2 - h1)/2
        5 L1+L2 h1+(h2-h1)/2
        6 L1 h1+(h2-h1)/2
        7 L1 h1
        8 0 h1
    ]

    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 5 0
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 1 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]
    CCSeg = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 1 0 0 0
        5 0 0 0 0
        6 0 0 0 0
        7 0 0 0 0
        8 0 1 0 0
    ]

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function contraction(
    ne = 15,
    tipo = 2,
    L1 = 1,
    h1 = 0.8,
    L2 = 2,
    h2 = 0.4,
    Re = 1,
    λ = 10^6,
)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 -h1/2
        2 L1 -h1/2
        3 L1 -(h2)/2
        4 L1+L2 -(h2)/2
        5 L1+L2 h2/2
        6 L1 h2/2
        7 L1 h1/2
        8 0 h1/2
    ]

    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 5 0
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 1 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]
    CCSeg = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 1 0 0 0
        5 0 0 0 0
        6 0 0 0 0
        7 0 0 0 0
        8 0 1 0 0
    ]

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function facing_step(
    ne = 15,
    tipo = 2,
    L = 5,
    h = 1,
    y_step = 0.5,
    x_step = 1,
    Re = 1,
    λ = 10^6,
)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 y_step
        2 x_step y_step
        3 x_step 0
        4 L 0
        5 L h
        6 0 h
    ]

    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 5 0
        5 5 6 0
        6 6 1 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
    ]
    CCSeg = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 1 0 0 0
        5 0 0 0 0
        6 0 1 0 0
    ]

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function channel(ne = 15, tipo = 2, h = 1, L = 1, Re = 50, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 L 0
        3 L h
        4 0 h
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 20*ne tipo
        2 ne tipo
        3 20*ne tipo
        4 ne tipo
    ]
    CCSeg = [
        0 0 0 0 0
        2 1 0 0 0
        3 0 0 0 0
        4 0 1 0 0
    ]
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function channel_simetria(ne = 15, tipo = 2, h = 0.01, L = 0.04, Re = 50, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 L 0
        3 L h
        4 0 h
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        0 0 0 0 0
        2 1 0 0 0
        3 1 0 0 0
        4 0 1 0 0
    ]
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function gradiente_pressao(ne = 15, tipo = 2, L = 1, Re = 50, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 L 0
        3 L L
        4 0 L
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        0 0 0 0 0
        2 1 0 0 0
        3 0 0 0 0
        4 1 0 0 0
    ]
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function gravidade(ne = 15, tipo = 2, L = 0.1, h = 0.01, g = 9.81, Re = 1, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 L 0
        3 L h
        4 0 h
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        0 0 0 0 0
        2 1 0 0 0
        3 1 0 0 0
        4 1 0 0 0
    ]
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function placa_plana_2(ne = 15, tipo = 2, L = 1, Ls = 0.4, h = 10, Re = 0.1, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 Ls 0
        3 L+Ls 0
        4 L+Ls h
        5 0 h
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 5 0
        5 5 1 0
    ]
    MALHA = [
        1 Int(ne / 2) tipo
        2 2*ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
    ]
    CCSeg = [
        0 1 0 0 0
        2 0 0 0 0
        3 1 0 1 0
        4 0 1 0 0
        5 0 1 0 0
    ]
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function placa_plana(ne = 15, tipo = 2, L = 1, h = 10, Re = 0.1, λ = 10^6)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 0 0
        2 L 0
        3 L h
        4 0 h
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        0 0 0 0 0
        2 1 0 0 0
        3 0 1 0 0
        4 0 1 0 0
    ]
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end

function von_karman(ne = 15, tipo = 2, L = 0.1, h = 0.05, r = 0.01, xc = 0.04, Re = 50)

    μ = 1 / Re

    v = λ / (2 * (μ + λ))
    E = 2 * (1 + v) * μ

    PONTOS = [
        1 r+xc 0
        2 xc -r
        3 -r+xc 0
        4 xc r
        5 L -h
        6 L h
        7 0 h
        8 0 -h
    ]

    SEGMENTOS = [
        1 1 2 -r
        2 2 3 -r
        3 3 4 -r
        4 4 1 -r
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 5 0
    ]

    MALHA = [
        1 nelem_circle tipo
        2 nelem_circle tipo
        3 nelem_circle tipo
        4 nelem_circle tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]

    CCSeg = [
        1 0 0 0 0
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0
        5 1 0 1 0
        6 0 0 0 0
        7 0 1 0 0
        8 0 0 0 0
    ]
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)
end


# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno

function gravidade_poisson(ne = 15, tipo = 2, h = 0.1, L = 1)

    PONTOS = [
        1 0 0
        2 L 0
        3 L h
        4 0 h
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    CCSeg = [
        1 0 0
        2 1 0
        3 1 0
        4 1 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function Poiseullie_poisson(ne = 15, tipo = 2, L = 1)

    PONTOS = [
        1 0 0
        2 L 0
        3 L L
        4 0 L
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        0 0 0
        2 1 0
        3 0 0
        4 1 0
    ]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function Taylor_Couette(ne = 15, tipo = 2, ra = 1, rb = 2, omega_a = -1, omega_b = 2)

    PONTOS = [
        1 ra 0
        2 rb 0
        3 0 rb
        4 0 ra
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 rb
        3 3 4 0
        4 4 1 -ra
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    CCSeg = [
        1 1 0
        2 0 omega_b*rb
        3 1 0
        4 0 omega_a*ra
    ]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function Poiseullie_Couette_poisson(ne = 15, tipo = 2, L = 1, u₀ = 1)

    PONTOS = [
        1 0 0
        2 L 0
        3 L L
        4 0 L
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    CCSeg = [
        1 0 0
        2 1 0
        3 0 u₀
        4 1 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function Couette_poisson(ne = 15, tipo = 2, L = 1, u₀ = 1)

    PONTOS = [
        1 0 0
        2 L 0
        3 L L
        4 0 L
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    CCSeg = [
        1 0 0
        2 1 0
        3 0 u₀
        4 1 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function von_karman_p(ne = 15, tipo = 2, L = 0.1, h = 0.05, r = 0.01, xc = 0.04)


    PONTOS = [
        1 r+xc 0
        2 xc -r
        3 -r+xc 0
        4 xc r
        5 L -h
        6 L h
        7 0 h
        8 0 -h
    ]

    SEGMENTOS = [
        1 1 2 -r
        2 2 3 -r
        3 3 4 -r
        4 4 1 -r
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 5 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]

    CCSeg = [
        1 1 0
        2 1 0
        3 1 0
        4 1 0
        5 1 0
        6 1 0
        7 0 0
        8 1 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function cavity_poisson(ne, tipo, Re, variavel)

    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]

    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    if variavel == "u"
        CCSeg = [
            1 0 0
            2 0 0
            3 0 1
            4 0 0
        ]
    elseif variavel == "v"
        CCSeg = [
            1 0 0
            2 0 0
            3 0 0
            4 0 0
        ]
    else
        CCSeg = [
            1 1 0
            2 1 0
            3 0 0
            4 1 0
        ]
    end


    # Condutividade Térmica do material
    k = 1 / Re
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function cavity_ψ(ne = 15, tipo = 2)
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]

    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]


    CCSeg = [
        1 0 0
        2 0 0
        3 0 0
        4 0 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function cavity_ω(ne = 15, tipo = 2)
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]

    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    CCSeg = [
        1 0 0
        2 0 0
        3 0 0
        4 0 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function von_karman_ω(ne = 15, tipo = 2, L = 0.1, h = 0.05, r = 0.01, xc = 0.04)

    PONTOS = [
        1 r+xc 0
        2 xc -r
        3 -r+xc 0
        4 xc r
        5 L -h
        6 L h
        7 0 h
        8 0 -h
    ]

    SEGMENTOS = [
        1 1 2 -r
        2 2 3 -r
        3 3 4 -r
        4 4 1 -r
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 5 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]

    CCSeg = [
        1 0 0
        2 0 0
        3 0 0
        4 0 0
        5 0 0
        6 0 0
        7 0 0
        8 0 0
    ]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function channel_poisson(ne = 15, tipo = 2, L = 5, h = 1, Re = 1, variavel = "u")

    PONTOS = [
        1 0 0
        2 L 0
        3 L h
        4 0 h
    ]

    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    if variavel == "u"
        CCSeg = [
            1 0 0
            2 1 0
            3 0 0
            4 0 1
        ]

    elseif variavel == "v"
        CCSeg = [
            1 0 0
            2 0 0
            3 0 0
            4 0 0
        ]
    else
        CCSeg = [
            1 1 0
            2 0 0
            3 1 0
            4 1 0
        ]
    end

    # Condutividade Térmica do material
    k = 1 / Re
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function NACA_2(ne = 3, tipo = 2)

    L = 20

    PONTOS = [
        NACA
        n+1 L/2 -L/2
        n+2 L/2 L/2
        n+3 -L/2 L/2
        n+4 -L/2 -L/2
    ]

    SEGMENTOS = [
        colum1 colum1 [colum2[begin:n-1]; 1] zeros(n)
        n+1 n+1 n+2 0
        n+2 n+2 n+3 0
        n+3 n+3 n+4 0
        n+4 n+4 n+1 0
    ]

    MALHA = [
        Int.(colum1) fill(ne, n) fill(tipo, n)
        Int.(n + 1) 4*ne tipo
        Int.(n + 2) 4*ne tipo
        Int.(n + 3) 4*ne tipo
        Int.(n + 4) 4*ne tipo
    ]

    CCSeg = [
        Int.(colum1) Int.(ones(n)) Int.(zeros(n))
        n+1 0 L/2
        n+2 0 0
        n+3 0 -L/2
        n+4 0 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function circular_cylinder(ne = 15, tipo = 2)
    r = 1
    L = 10
    PONTOS = [
        1 r 0
        2 0 -r
        3 -r 0
        4 0 r
        5 L/2 -L/2
        6 L/2 L/2
        7 -L/2 L/2
        8 -L/2 -L/2
    ]

    SEGMENTOS = [
        1 1 2 -r
        2 2 3 -r
        3 3 4 -r
        4 4 1 -r
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 5 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]

    CCSeg = [
        1 1 0
        2 1 0
        3 1 0
        4 1 0
        5 0 L/2
        6 0 0
        7 0 -L/2
        8 0 0
    ]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function cilindro_chorin(
    ne = 15,
    tipo = 2,
    r = 0.1,
    L = 3,
    h = 1,
    xc = 0.5,
    Re = 1,
    variavel = "u",
)

    yc = h / 2

    PONTOS = [
        1 r+xc +yc
        2 xc -r+yc
        3 -r+xc yc
        4 xc r+yc
        5 0 0
        6 L 0
        7 L h
        8 0 h
    ]

    SEGMENTOS = [
        1 1 2 -r
        2 2 3 -r
        3 3 4 -r
        4 4 1 -r
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 5 0
    ]

    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
    ]

    if variavel == "u"
        CCSeg = [
            1 0 0
            2 0 0
            3 0 0
            4 0 0
            5 0 0
            6 1 0
            7 0 0
            8 0 1
        ]
    elseif variavel == "v"
        CCSeg = [
            1 0 0
            2 0 0
            3 0 0
            4 0 0
            5 0 0
            6 1 0
            7 0 0
            8 0 0
        ]
    else
        CCSeg = [
            1 1 0
            2 1 0
            3 1 0
            4 1 0
            5 1 0
            6 0 0
            7 1 0
            8 1 0
        ]
    end
    # Condutividade Térmica do material
    k = 1 / Re
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function NACA_potencial(ne = 3, tipo = 2)

    PONTOS = [
        NACA
        n+1 10 0.1
        n+2 -10 0
        n+3 10 -0.1
    ]

    SEGMENTOS = [
        colum1 colum1 colum2 zeros(n)
        n+1 n+1 n+2 10
        n+2 n+2 n+3 10
        n+3 n+3 1 0
    ]

    MALHA = [
        Int.(colum1) [fill(ne, n - 1); 4 * ne] fill(tipo, n)
        Int.(n + 1) 4*ne tipo
        Int.(n + 2) 4*ne tipo
        Int.(n + 3) 4*ne tipo
    ]

    CCSeg = [
        Int.(colum1) Int.([ones(n - 1); 0]) Int.([zeros(n - 1); 1])
        n+1 1 1
        n+2 1 1
        n+3 1 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function NACA_artigo_v(ne = 3, tipo = 2)

    L = 20

    PONTOS = [
        NACA
        n+1 L/2 -L/2
        n+2 L/2 L/2
        n+3 -L/2 L/2
        n+4 -L/2 -L/2
    ]

    SEGMENTOS = [
        colum1 colum1 [colum2[begin:n-1]; 1] zeros(n)
        n+1 n+1 n+2 0
        n+2 n+2 n+3 0
        n+3 n+3 n+4 0
        n+4 n+4 n+1 0
    ]

    MALHA = [
        Int.(colum1) fill(ne, n) fill(tipo, n)
        Int.(n + 1) 4*ne tipo
        Int.(n + 2) 4*ne tipo
        Int.(n + 3) 4*ne tipo
        Int.(n + 4) 4*ne tipo
    ]

    CCSeg = [
        Int.(colum1) Int.(ones(n)) Int.(zeros(n))
        n+1 1 0
        n+2 1 0
        n+3 0 0
        n+4 0 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function NACA_artigo_w(ne = 3, tipo = 2)

    L = 20

    PONTOS = [
        NACA
        n+1 L/2 -L/2
        n+2 L/2 L/2
        n+3 -L/2 L/2
        n+4 -L/2 -L/2
    ]

    SEGMENTOS = [
        colum1 colum1 [colum2[begin:n-1]; 1] zeros(n)
        n+1 n+1 n+2 0
        n+2 n+2 n+3 0
        n+3 n+3 n+4 0
        n+4 n+4 n+1 0
    ]

    MALHA = [
        Int.(colum1) fill(ne, n) fill(tipo, n)
        Int.(n + 1) 4*ne tipo
        Int.(n + 2) 4*ne tipo
        Int.(n + 3) 4*ne tipo
        Int.(n + 4) 4*ne tipo
    ]

    CCSeg = [
        Int.(colum1) Int.(ones(n)) Int.(zeros(n))
        n+1 1 1
        n+2 1 0
        n+3 0 0
        n+4 0 0
    ]

    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end