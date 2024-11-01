# Entrada de dados para análise do problemade helmholtz pelo
# método dos elementos de contorno
function helm1d(ne=15, tipo=2)
    PONTOS = [1 0 0; 2 1 0; 3 1 1; 4 0 1]
    SEGMENTOS = [1 1 2 0; 2 2 3 0; 3 3 4 0; 4 4 1 0]
    MALHA = [1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo]

    CCSeg = [1 0 0 0
        2 1 0 0
        3 1 1 0
        4 1 0 0]


    # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
    # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
    # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

    GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
    RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
    CW = sqrt(GE / RO) # Velocidade de propagação de onda
    FR = 0.1

    # Malha de pontos internos
    return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR=FR, CW=CW, GE=GE)
end

function ANA_helm1d(dad)
    sin.(dad.k.FR * dad.NOS[:, 2]) / dad.k.FR / cos(dad.k.FR)
end


# Entrada de dados para análise do problemade helmholtz pelo
# método dos elementos de contorno
function helmdirechlet(ne=15, tipo=2)
    PONTOS = [1 0 0; 2 1 0; 3 1 1; 4 0 1]
    SEGMENTOS = [1 1 2 0; 2 2 3 0; 3 3 4 0; 4 4 1 0]
    MALHA = [1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo]
    CCSeg = [1 0 0
        2 0 0
        3 0 1
        4 0 0]    # Condutividade Térmica do material


    # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
    # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
    # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

    GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
    RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
    CW = sqrt(GE / RO) # Velocidade de propagação de onda
    FR = 0.1

    # Malha de pontos internos
    return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR=FR, CW=CW, GE=GE)
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
function ANA_helmdirechlet(dad)
    sin.(dad.NOS[:, 1] * √(dad.k.FR^2 - π^2)) .* sin(π * dad.NOS[:, 2]) / sin(√(dad.k.FR^2 - π^2))
end

# Entrada de dados para análise do problemade helmholtz pelo
# método dos elementos de contorno
function helmcirculo(ne=15, tipo=2, FR=1)
    PONTOS = [1 -1 0; 2 0 -1; 3 1 0; 4 0 1]
    SEGMENTOS = [1 1 2 1; 2 2 3 1; 3 3 4 1; 4 4 1 1]
    MALHA = [1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo]

    CCSeg = [1 1 0 0
        2 1 1 0
        3 1 1 0
        4 1 1 0]


    # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
    # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
    # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

    GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
    RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
    CW = sqrt(GE / RO) # Velocidade de propagação de onda

    # Malha de pontos internos
    return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR=FR, CW=CW, GE=GE)
end

function ANA_helmcirculo(dad)
    x = dad.NOS[:, 1]
    y = dad.NOS[:, 2]

    r = sqrt.(x .^ 2 + y .^ 2)
    besselj0.(dad.k.CW * r) / besselj0(dad.k.CW)
end

# Entrada de dados para análise do problemade helmholtz pelo
# método dos elementos de contorno
function helmcirculoinfinito(ne=15, tipo=2, FR=1)
    PONTOS = [1 -1 0; 2 0 -1; 3 1 0; 4 0 1]
    SEGMENTOS = [1 1 4 -1; 2 4 3 -1; 3 3 2 -1; 4 2 1 -1]
    MALHA = [1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo]

    CCSeg = [1 1 1 0
        2 1 1 0
        3 1 1 0
        4 1 1 0]    # Condutividade Térmica do material


    # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
    # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
    # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

    GE = 1.E0    # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
    RO = 1       # Densidade do meio  (RO = 1 no caso de ac�stica)
    CW = sqrt(GE / RO) # Velocidade de propagação de onda

    # Malha de pontos internos
    return helmholtz, PONTOS, SEGMENTOS, MALHA, CCSeg, (FR=FR, CW=CW, GE=GE)
end