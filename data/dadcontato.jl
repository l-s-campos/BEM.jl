function sapata(ref=11, tipo=2)
    # Matriz para definição de pontos
    w = 6.5       # largura da sapata
    R = 70        # Raio da sapata
    theta = asin(w / R)
    x = R * cos(theta)
    y = R - x

    # Definição dos pontos
    PONTOS = [
        1 -w y
        2 w y
        3 w 2*w
        4 -w 2*w
    ]

    # Segmentos que definem a geometria
    # SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final, Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é à esquerda do segmento
    #                   < 0 -> O centro é à direita do segmento
    #                   = 0 -> O segmento é uma linha reta

    SEGMENTOS = [
        1 1 2 R
        2 2 3 0
        3 3 4 0
        4 4 1 0
    ]

    # Matriz para definição da malha
    MALHA = [
        1 3*ref+1 tipo
        2 1*ref tipo
        3 2*ref+1 tipo
        4 1*ref tipo
    ]

    cargav = 136.8551  # Carga vertical

    # Condições de contorno nos segmentos
    # CCSeg=[Segmento,tipo da CDC em n, valor da CDC em n , tipo da CDC em t, valor da CDC em t]
    # tipo da CDC = 0 => o deslocamento é conhecido
    # tipo da CDC = 1 => a força de superfície é conhecida
    # tipo da CDC = 2 => tudo desconhecido

    CCSeg = [1 2 0 2 0
        2 1 0 1 0
        3 1 0 1 -cargav
        4 1 0 1 0]

    E = 73.4e3  # Módulo de Elasticidade
    v = 0.33   # Coeficiente de Poisson
    μ = 0.1  # Coeficiente de atrito

    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v, μ=μ)
end

function sapata_conforme(ref=11, tipo=2)

    # Matriz para definição de pontos
    w = 6.5  # largura da sapata
    PONTOS = [1 -w 0; 2 w 0; 3 w 2w; 4 -w 2w]

    # Segmentos que definem a geometria
    #  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
    #            Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é à direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 0; 2 2 3 0; 3 3 4 0; 4 4 1 0]

    # Matriz para definição da malha
    MALHA = [1 ref tipo; 2 ref tipo; 3 ref tipo; 4 ref tipo]

    cargav = 136.8551  # Carga vertical

    # Condições de contorno nos segmentos
    # CCSeg=[Segmento,tipo da CDC em n, valor da CDC em n , tipo da CDC em t,
    #        valor da CDC em t]
    # tipo da CDC = 0 => o deslocamento é conhecido
    # tipo da CDC = 1 => a força de superfície é conhecida
    CCSeg = [
        1 2 0 2 0
        2 1 0 1 0
        3 1 0 1 -cargav
        4 1 0 1 0
    ]

    E = 73.4e3  # Módulo de Elasticidade
    v = 0.33  # Coeficiente de Poisson
    μ = 0.1  # Coeficiente de atrito


    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v, μ=μ)

end

function sapata_conforme_multicorpo(ref=11, tipo=2)

    # Matriz para definição de pontos
    w = 6.5  # largura da sapata
    PONTOS = [1 -w 0; 2 w 0; 3 w 2w; 4 -w 2w; 5 -w -2w; 6 w -2w; 7 w 0; 8 -w 0]

    # Segmentos que definem a geometria
    #  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
    #            Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é à direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 0; 2 2 3 0; 3 3 4 0; 4 4 1 0; 5 5 6 0; 6 6 7 0; 7 7 8 0; 8 8 5 0]

    # Matriz para definição da malha
    MALHA = [1 ref tipo; 2 ref tipo; 3 ref tipo; 4 ref tipo; 5 ref tipo; 6 ref tipo; 7 ref tipo; 8 ref tipo]

    cargav = 136.8551  # Carga vertical

    # Condições de contorno nos segmentos
    # CCSeg=[Segmento,tipo da CDC em n, valor da CDC em n , tipo da CDC em t,
    #        valor da CDC em t]
    # tipo da CDC = 0 => o deslocamento é conhecido
    # tipo da CDC = 1 => a força de superfície é conhecida
    # tipo da CDC = 2 contato -> 1=superfície correlata 2 = master ou slave
    CCSeg = [
        1 2 8 2 0
        2 1 0 1 0
        3 1 0 1 -cargav
        4 1 0 1 0
        5 0 0 0 0
        6 1 0 1 0
        7 2 1 2 1
        8 1 0 1 0
    ]

    E = 73.4e3  # Módulo de Elasticidade
    v = 0.33  # Coeficiente de Poisson
    μ = 0.1  # Coeficiente de atrito
    subregioes = BEM.define_SubRegioes_contato(MALHA, CCSeg, [4, 8])

    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=[E, 1.1E], nu=[v, v], μ=μ), subregioes
end