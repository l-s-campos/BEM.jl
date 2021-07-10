function elastico1d(ne=15,tipo=2)
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
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida
    CCSeg = [1 1 0 1 0
        2 1 1 1 0
        3 1 0 0 0
        4 0 0 1 0]
        # CCSeg = [1 0 1
        # 2 0 1
        # 3 0 1
        # 4 0 1]
    # Condutividade Térmica do material
    E = 1
    v = 0.0
    # Malha de pontos internos
    return elastico,PONTOS,SEGMENTOS,MALHA,CCSeg,[E,v]
end

function elastico1diso(ne=15,tipo=2)
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
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida
    # CCSeg = [1 1 0 1 0
    #     2 1 1 1 0
    #     3 1 0 0 0
    #     4 0 0 1 0]
        CCSeg=[1 0 0 0 0
        2 1 0 1 0
        3 1 0 1 1
        4 1 0 1 0];
    # Condutividade Térmica do material
    E = 1
    v = 0.0
    # Malha de pontos internos
    return elastico_iso,PONTOS,SEGMENTOS,MALHA,CCSeg,[E,v]
end


function elastico_aniso_1d(ne=15,tipo=2)
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
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida
    CCSeg = [1 1 0 1 0
        2 1 1 1 0
        3 1 0 0 0
        4 0 0 1 0]
    # Condutividade Térmica do material
    Material = [1 2.200000    4.4   0.7692      0.4286  15 1;
    2 2.200000    4.4   0.7692      0.4286 -25 1];
    const_material = Compute_Material(Material)
    mi,A,q,g = const_material
    k=zeros(Complex,7,2)
    k[1,:]=mi
    k[2:3,:]=A
    k[4:5,:]=q
    k[6:7,:]=g
    # Malha de pontos internos
    return elastico_aniso,PONTOS,SEGMENTOS,MALHA,CCSeg,k
end