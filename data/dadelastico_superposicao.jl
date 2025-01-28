function elastico1d_inclusao(ne = 15, tipo = 2)
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
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida
    CCSeg = [
        1 1 0 1 0
        2 1 1 1 0
        3 1 0 1 0
        4 0 0 0 0
    ]
    # CCSeg = [1 0 1
    # 2 0 1
    # 3 0 1
    # 4 0 1]
    # Condutividade Térmica do material
    E = 1
    v = 0.0
    # Malha de pontos internos
    PONTOSi = [
        1 0.25 0.25
        2 0.75 0.25
        3 0.75 0.75
        4 0.25 0.75
    ]
    CCSegi = [
        1 1 0 1 0
        2 1 0 1 0
        3 1 0 1 0
        4 1 0 1 0
    ]
    mult = 10
    return [
        [elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E = E, nu = v)],
        [
            elastico,
            PONTOSi,
            SEGMENTOS,
            MALHA,
            CCSegi,
            (E = mult * E - E, nu = mult * v - v),
        ],
    ]
end
