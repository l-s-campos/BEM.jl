# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno
function potencial1d(ne=15,tipo=2)
    PONTOS  = [1 0 0
        2 1 0
        3 1 1
        4 0 1 ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [1 ne tipo 
            2 ne tipo 
            3 ne tipo 
            4 ne tipo]
        # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [1 1 0
        2 1 -1
        3 1 0
        4 0 0]
        # CCSeg = [1 0 1
        # 2 0 1
        # 3 0 1
        # 4 0 1]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial,PONTOS,SEGMENTOS,MALHA,CCSeg,k
end


function placacomfuro(ne=15,tipo=2)
       PONTOS  = [1 0 0
        2 1 0
        3 1 1
        4 0 1 
        5 0.25 0.5
        6 0.75 0.5]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0
        5 5 6 -0.25
        6 6 5 -0.25]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [1 ne tipo 
            2 ne tipo 
            3 ne tipo 
            4 ne tipo
            5 ne tipo
            6 ne tipo]
        # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [1 1 0
        2 1 1
        3 1 0
        4 0 0
        5 1 0
        6 1 0]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return PONTOS,SEGMENTOS,MALHA,CCSeg,k
end


# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno
function potencial1diso(ne=15,tipo=2)
    PONTOS  = [1 0 0
        2 1 0
        3 1 1
        4 0 1 ]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 0
        2 2 3 0
        3 3 4 0
        4 4 1 0]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [1 ne tipo 
            2 ne tipo 
            3 ne tipo 
            4 ne tipo]
        # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [1 1 0
        2 1 -1
        3 1 0
        4 0 0]
        # CCSeg = [1 0 1
        # 2 0 1
        # 3 0 1
        # 4 0 1]
    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return potencial_iso,PONTOS,SEGMENTOS,MALHA,CCSeg,k
end