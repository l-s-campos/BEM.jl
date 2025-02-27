# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno
function potencial1d(ne = 15, tipo = 2)
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    # CCSeg = [1 1 0%correto
    #     2 1 -1
    #     3 1 0
    #     4 0 0]
    CCSeg = [
        1 1 0
        2 0 1
        3 1 0
        4 0 0
    ]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function membrana_triangular(ne = 15, tipo = 2)
    PONTOS = [
        1 2.5 -5*sqrt(3)/6
        2 0 sqrt(25 - 2.5^2)-5*sqrt(3)/6
        3 -2.5 -5*sqrt(3)/6
    ]
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 1 0
    ]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
    ]
    CCSeg = [
        1 0 0
        2 0 0
        3 0 0
    ]
    k = 1
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function placacomfuro(ne = 15, tipo = 2)
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
        5 0.25 0.5
        6 0.75 0.5
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
        5 5 6 -0.25
        6 6 5 -0.25
    ]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
    ]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [
        1 1 0
        2 1 1
        3 1 0
        4 0 0
        5 1 0
        6 1 0
    ]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return PONTOS, SEGMENTOS, MALHA, CCSeg, k
end


# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno
function potencial1diso(ne = 15, tipo = 2)
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [
        1 1 0
        2 1 -1
        3 1 0
        4 0 0
    ]
    # CCSeg = [1 0 1
    # 2 0 1
    # 3 0 1
    # 4 0 1]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial_iga, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end
function T_ana_1d(NOS)
    NOS[:, 1]
end
function placa_moulton(ne = 8, tipo = 2)
    PONTOS = [
        1 -1 0
        2 0 0
        3 1 0
        4 1 1
        5 -1 1
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 5 0
        5 5 1 0
    ]

    MALHA = [
        1 floor(Int, ne / 2) tipo
        2 floor(Int, ne / 2) tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
    ]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido

    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [
        1 0 0
        2 1 0
        3 1 0
        4 1 0
        5 1 0
    ]
    # CCSeg = [1 0 1
    # 2 0 1
    # 3 0 1
    # 4 0 1]
    # Condutividade Térmica do material
    k = 1
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function corrigeCDC_moulton(dad)
    # floor(Int,ne/2)
    NELEM = length(dad.ELEM[:, 1])
    NN = length(dad.NOS[:, 1])
    N = floor(NELEM / 4)


    for j = 1:NELEM
        # if j>=(1) && j<(N)
        # 	k=length(dad.ELEM[j].indices)
        # 	for f=1:k
        # 		#  @show 1,i,j
        # 		# dad.ELEM[j].valorCDC[f]=dad.ELEM[j].valorCDC[f]
        # 	end
        # end
        if j >= (N + 1) && j <= (2 * N)
            k = length(dad.ELEM[j].indices)
            for f = 1:k
                NPonto = dad.ELEM[j].indices[f]
                x, y = dad.NOS[NPonto, :]
                r = sqrt(x^2 + y^2)
                theta = atan(y, x)
                # @show theta/pi*180
                q = -cos(theta / 2) / (2 * sqrt(r))
                dad.ELEM[j].valorCDC[f] = q

            end
        elseif j >= (2 * N + 1) && j <= (3 * N)
            k = length(dad.ELEM[j].indices)
            for f = 1:k
                NPonto = dad.ELEM[j].indices[f]
                x, y = dad.NOS[NPonto, :]
                r = sqrt(x^2 + y^2)
                theta = atan(y, x)
                # @show theta/pi*180
                q = -sin(theta / 2) / (2 * sqrt(r))
                dad.ELEM[j].valorCDC[f] = q
                # @show NPonto
                # @show 3,i,j
            end
        elseif j >= (3 * N + 1) && j <= (4 * N)
            k = length(dad.ELEM[j].indices)
            for f = 1:k
                NPonto = dad.ELEM[j].indices[f]
                x, y = dad.NOS[NPonto, :]
                r = sqrt(x^2 + y^2)
                theta = atan(y, x)
                # @show theta/pi*180
                q = cos(theta / 2) / (2 * sqrt(r))
                dad.ELEM[j].valorCDC[f] = q
                # @show NPonto
                # @show 4,i,j
            end
        end
    end
    dad
end


function T_ana_moulton(NOS)
    NN = length(NOS[:, 1])
    T = zeros(length(NOS[:, 1]))
    for i = 1:NN
        x, y = NOS[i, :]
        r = sqrt(x^2 + y^2)
        xc = 0
        yc = 0
        x1 = 1
        y1 = 0
        teta1, teta = BEM.calcula_arco(x1, y1, x, y, xc, yc, r)
        T[i] = sqrt(r) * cos(teta / 2)
    end
    return T
end



function dad_laquini1(ne = 5, tipo = 2)
    # Matriz para definição de pontos que definem a geometria
    # PONTOS = [número do ponto, coord. x do ponto, coord. y do ponto]
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS=[N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    #= Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
                              inicial para o ponto final)
                         < 0 -> O centro é à direita do segmento (do ponto
                              inicial para o ponto final)
                         = 0 -> O segmento é uma linha reta        =#
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    # Matriz para definição da malha
    # MALHA = [n° do segmento, n° de elementos no segmento]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    #= Condições de contorno nos segmentos
       CCSeg=[no do segmento, tipo da CDC, valor da CDC]
       tipo da CDC = 0 => a temperatura é conhecida
       tipo da CDC = 1 => O fluxo é conhecido         =#
    CCSeg = [
        1 0 0
        2 0 0
        3 1 -1
        4 0 0
    ]
    k = 1  # Condutividade Térmica do material

    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function dad_laquini2(ne = 15, tipo = 2)
    # Matriz para definição de pontos que definem a geometria
    # PONTOS = [número do ponto, coord. x do ponto, coord. y do ponto]
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS=[N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    #= Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
                              inicial para o ponto final)
                         < 0 -> O centro é à direita do segmento (do ponto
                              inicial para o ponto final)
                         = 0 -> O segmento é uma linha reta        =#
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    # Matriz para definição da malha
    # MALHA = [n° do segmento, n° de elementos no segmento]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    #= Condições de contorno nos segmentos
       CCSeg=[no do segmento, tipo da CDC, valor da CDC]
       tipo da CDC = 0 => a temperatura é conhecida
       tipo da CDC = 1 => O fluxo é conhecido         =#
    CCSeg = [
        1 0 0
        2 1 -1
        3 1 -1
        4 0 0
    ]
    k = 1  # Condutividade Térmica do material

    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function dad_laquini3(ne = 5, tipo = 2)
    # Matriz para definição de pontos que definem a geometria
    # PONTOS = [número do ponto, coord. x do ponto, coord. y do ponto]
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS=[N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    #= Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
                              inicial para o ponto final)
                         < 0 -> O centro é à direita do segmento (do ponto
                              inicial para o ponto final)
                         = 0 -> O segmento é uma linha reta        =#
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    # Matriz para definição da malha
    # MALHA = [n° do segmento, n° de elementos no segmento]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    #= Condições de contorno nos segmentos
       CCSeg=[no do segmento, tipo da CDC, valor da CDC]
       tipo da CDC = 0 => a temperatura é conhecida
       tipo da CDC = 1 => O fluxo é conhecido         =#
    CCSeg = [
        1 0 0
        2 0 0
        3 0 1
        4 0 0
    ]
    k = 1  # Condutividade Térmica do material

    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end
#________________________________________________________________________________________
"Funcao para calcular a temperatura do primeiro problema da tese do Laquini
Entrada: NOS
Saída: Temperatura analítica"
function fa1(NOS, ns = 500)
    nc = size(NOS, 1)
    Ts = zeros(nc)
    for i = 1:nc
        x = NOS[i, 1]
        y = NOS[i, 2]
        T_ana = 0
        for n = 1:ns
            # inc = ((((-1)^(n + 1)) + 1) / n) * sin(n * π * x) * (sinh(n * π * y) / (n * π * cosh(n * π)))
            inc =
                ((((-1)^(n + 1)) + 1) / n) *
                sin(n * π * x) *
                ((1 - exp(-2n * π * y)) / (n * π * (1 + exp(-2n * π)))) *
                exp(n * π * (y - 1))
            if isnan(inc)
                break
            end
            T_ana += inc
        end
        T_ana *= 2 / π
        Ts[i] = T_ana
    end
    Ts
end
"Funcao para calcular a temperatura do segundo problema da tese do Laquini
Entrada: NOS
Saída: Temperatura analítica"
function fa2(NOS, ns = 500)
    nc = size(NOS, 1)
    Ts = zeros(nc)
    for i = 1:nc
        x = NOS[i, 1]
        y = NOS[i, 2]
        T_ana = 0.0
        for n = 1:ns
            inc = (
                (sin(n * π * x / 2)) * (
                    (
                        (2 / (n * π)) *
                        (((-1)^(n + 1) + 1) / (((n * π) / 2) * cosh(n * π / 2))) +
                        (8 / ((n^2) * (π^2))) * sin(n * π / 2) * tanh(n * π / 2)
                    ) * sinh(n * π * y / 2) -
                    ((8 / ((n^2) * (π^2))) * sin(n * π / 2) * cosh(n * π * y / 2))
                )
            )
            if isnan(inc)
                # @show n, inc
                break
            end
            T_ana += inc
        end
        T_ana += x
        Ts[i] = T_ana
    end
    Ts
end

#________________________________________________________________________________________
"Funcao para calcular a temperatura do terceiro problema da tese do Laquini
Entrada: NOS
Saída: Temperatura analítica"
function fa3(NOS, ns = 1000)
    nc = size(NOS, 1)
    Ts = zeros(nc)
    for i = 1:nc
        x = NOS[i, 1]
        y = NOS[i, 2]
        T_ana = 0.0
        for n = 1:ns
            # if y == 1
            #     inc = ((((-1)^(n + 1)) + 1) / n) * sin(n * π * x)
            # else
            #     inc = ((((-1)^(n + 1)) + 1) / n) * sin(n * π * x) * (sinh(n * π * y) / (sinh(n * π)))
            # end
            inc =
                ((((-1)^(n + 1)) + 1) / n) *
                sin(n * π * x) *
                ((1 - exp(-2n * π * y)) / (1 - exp(-2n * π))) *
                exp(n * π * (y - 1))
            if isnan(inc)
                #  @show i, n, inc
                break
            end
            T_ana += inc
        end
        T_ana *= 2 / π
        Ts[i] = T_ana
    end
    Ts
end

# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno
function quarto_circ(ne = 15, tipo = 2)
    PONTOS = [
        1 1 0
        2 2 0
        3 0 2
        4 0 1
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [
        1 1 2 0
        2 2 3 2
        3 3 4 0
        4 4 1 -1
    ]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    # CCSeg = [1 1 0%correto
    #     2 1 -1
    #     3 1 0
    #     4 0 0]
    CCSeg = [
        1 1 0
        2 1 -200
        3 1 0
        4 0 100
    ]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function ana_quarto_circ(NOS, ri = 1, re = 2, ti = 100, qe = -200)
    NN = length(NOS[:, 1])
    T = zeros(length(NOS[:, 1]))
    for i = 1:NN
        x, y = NOS[i, :]
        r = sqrt(x^2 + y^2)
        T[i] = ti - qe * re * log(r / ri)

    end
    return T

end



function difusao_placa(ne = 15, tipo = 2)
    #Exercício 2 aula 7
    PONTOS = [
        1 0 0
        2 1 0
        3 1 1
        4 0 1
    ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [
        1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [
        1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
    ]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [
        1 1 0
        2 1 0
        3 0 1
        4 1 0
    ]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial, PONTOS, SEGMENTOS, MALHA, CCSeg, k
end

function ana_difusao_placa(tempo, ys)
    serie = 0
    for k = 0:1000
        serie +=
            ((-1)^k) / (2 * k + 1) *
            exp(-((2 * k + 1)^2 * pi^2 * 1 * tempo) / (4 * 1^2)) *
            cos((2 * k + 1) * pi * ys / (2 * 1))
    end
    T_analit = 1 - 4 / pi * serie
    return T_analit
end