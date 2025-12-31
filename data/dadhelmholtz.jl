# Entrada de dados para análise do problemade helmholtz pelo
# método dos elementos de contorno
function helm1d(ne = 15, tipo = 2, FR = 1)
    PONTOS = [1 0 0; 2 1 0; 3 1 1; 4 0 1]
    SEGMENTOS = [1 1 2 0; 2 2 3 0; 3 3 4 0; 4 4 1 0]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    CCSeg = [
        1 1 0 0
        2 1 1 0
        3 1 0 0
        4 0 0 0
    ]


    # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
    # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
    # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

    GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
    RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
    CW = sqrt(GE / RO) # Velocidade de propagação de onda
    # FR = 0.1

    # Malha de pontos internos
    return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR = FR, CW = CW, GE = GE)
end

function ANA_helm1d(dad, array)
    sin.(dad.k.FR * array[:, 1]) / dad.k.FR / cos(dad.k.FR)
end


# Entrada de dados para análise do problemade helmholtz pelo
# método dos elementos de contorno
function helmdirechlet(ne = 15, tipo = 2, FR = 1)
    PONTOS = [1 0 0; 2 1 0; 3 1 1; 4 0 1]
    SEGMENTOS = [1 1 2 0; 2 2 3 0; 3 3 4 0; 4 4 1 0]
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
    ]    # Condutividade Térmica do material


    # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
    # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
    # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

    GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
    RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
    CW = sqrt(GE / RO) # Velocidade de propagação de onda
    #FR = 0.1

    # Malha de pontos internos
    return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR = FR, CW = CW, GE = GE)
end


function corrigeCDC_helmdirechlet(dad)
    NELEM = length(dad.ELEM[:, 1])
    N = floor(NELEM / 4)

    for j = 1:NELEM
        if j >= (N + 1) && j <= (2 * N)
            k = length(dad.ELEM[j].indices)
            for f = 1:k
                NPonto = dad.ELEM[j].indices[f]
                x, y = dad.NOS[NPonto, :]
                u = sin(π * y)
                dad.ELEM[j].valorCDC[f] = u
            end
        end
    end
end
function ANA_helmdirechlet(dad, array)
    sin.(array[:, 1] * √(dad.k.FR^2 - π^2)) .* sin.(π * array[:, 2]) /
    sin(√(dad.k.FR^2 - π^2))
end

# Entrada de dados para análise do problemade helmholtz pelo
# método dos elementos de contorno

#A função NR calcula o numero de divisoes radiais para que tenhamos N pontos internos totais
function calcula_NR(NT, R, N; Δt_interno = 0.05, Δt_externo = 0.05)
    r_NR = R - Δt_externo
    r_1 = Δt_interno

    f(x) =
        x / 2 * (x + 1) * (r_NR - r_1) / (x - 1) +
        x * (r_1 - (r_NR - r_1) / (x - 1)) +
        -N * R / NT

    sol = find_zero(f, 0.5)  # ponto inicial 0.5
    return sol
end

function pontos_internos_circulo(NR, NT, R; Δr_interno = 0.05, Δr = 0.05)
    @assert NR ≥ 1 "NR deve ser ≥ 1"
    @assert NT ≥ 3 "NT deve ser ≥ 3"
    @assert Δr_interno < R "Δr_interno deve ser menor que R"

    # Divisões radiais
    r_values = range(Δr_interno, stop = R - Δr, length = NR)[1:end]

    pontos = Float64[]

    for r in r_values
        # número de pontos no círculo de raio r (mínimo 3)
        Ni = max(3, round(Int, NT * r / R))

        θ_values = range(0, stop = 2π, length = Ni + 1)[1:end-1]

        # acumular coordenadas (x,y)
        for θ in θ_values
            append!(pontos, (r * cos(θ), r * sin(θ)))
        end
    end

    # Converter para matriz N×2
    npoints = length(pontos) ÷ 2
    pontos_internos = reshape(pontos, 2, npoints)'  # transposto -> (N×2)

    return pontos_internos
end


# function helmcirculo(ne = 15, tipo = 2, FR = 1)
#     PONTOS = [1 -1 0; 2 0 -1; 3 1 0; 4 0 1]
#     SEGMENTOS = [1 1 2 1; 2 2 3 1; 3 3 4 1; 4 4 1 1]
#     MALHA = [
#         1 ne tipo
#         2 ne tipo
#         3 ne tipo
#         4 ne tipo
#     ]

#     CCSeg = [
#         1 0 1
#         2 0 1
#         3 0 1
#         4 0 1
#     ]


#     # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
#     # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
#     # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

#     GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
#     RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
#     CW = sqrt(GE / RO) # Velocidade de propagação de onda

#     # Malha de pontos internos
#     return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR = FR, CW = CW, GE = GE)
# end

function helmcirculo(ne = 15, tipo = 2, FR = 1)
    PONTOS = [1 0 0; 2 1 0; 3 0 1]
    SEGMENTOS = [1 1 2 0; 2 2 3 1; 3 3 1 0]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
    ]

    CCSeg = [
        1 1 0
        2 0 1
        3 1 0
    ]


    # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
    # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
    # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

    GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
    RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
    CW = sqrt(GE / RO) # Velocidade de propagação de onda

    # Malha de pontos internos
    return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR = FR, CW = CW, GE = GE)
end

function ANA_helmcirculo(dad, array)
    r = sqrt.(array[:, 1] .^ 2 + array[:, 2] .^ 2)
    besselj0.(dad.k.FR * r) / besselj0(dad.k.FR)
end

# Entrada de dados para análise do problemade helmholtz pelo
# método dos elementos de contorno
function helmcirculoinfinito(ne = 15, tipo = 2, FR = 1)
    PONTOS = [1 -1 0; 2 0 -1; 3 1 0; 4 0 1]
    SEGMENTOS = [1 1 4 -1; 2 4 3 -1; 3 3 2 -1; 4 2 1 -1]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]

    CCSeg = [
        1 1 1 0
        2 1 1 0
        3 1 1 0
        4 1 1 0
    ]    # Condutividade Térmica do material


    # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
    # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
    # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

    GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
    RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
    CW = sqrt(GE / RO) # Velocidade de propagação de onda

    # Malha de pontos internos
    return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR = FR, CW = CW, GE = GE)
end
