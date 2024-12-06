function elastico1d(ne=15, tipo=2; nt=false)
    PONTOS = [1 0 0
        2 1 0
        3 1 1
        4 0 1]
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
        2 2ne tipo
        3 ne tipo
        4 2ne tipo]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida
    if nt
        CCSeg = [1 1 0 0 0
            2 1 1 1 0
            3 1 0 1 0
            4 0 0 1 0]
    else
        CCSeg = [1 1 0 0 0
            2 1 1 1 0
            3 1 0 1 0
            4 0 0 1 0]
    end

    E = 1
    v = 0.3
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v)
end

function elastico1diga(ne=15, tipo=2)
    PONTOS = [1 0 0
        2 1 0
        3 1 1
        4 0 1]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 1
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
    CCSeg = [1 0 0 0 0
        2 1 0 1 0
        3 1 0 1 1
        4 1 0 1 0]
    # Condutividade Térmica do material
    E = 1
    v = 0.0
    # Malha de pontos internos
    return elastico_iga, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v)
end


function elastico_aniso_1d(ne=15, tipo=2)
    PONTOS = [1 0 0
        2 1 0
        3 1 1
        4 0 1]
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
    Material = [1 1.001 1.0 0.5 0.0 0 1.0]
    # const_material = Compute_Material(Material)
    # mi,A,q,g = const_material
    # k=zeros(Complex,7,2)
    # k[1,:]=mi
    # k[2:3,:]=A
    # k[4:5,:]=q
    # k[6:7,:]=g
    # Malha de pontos internos
    return elastico_aniso, PONTOS, SEGMENTOS, MALHA, CCSeg, Compute_Material(Material)
end


function elastico_aniso_1d_iga(ne=15, tipo=2)
    PONTOS = [1 0 0
        2 1 0
        3 1 1
        4 0 1]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 1
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
    Material = [1 1.000001 1 1 1 0 1]
    # const_material = Compute_Material(Material)
    # mi,A,q,g = const_material
    # k=zeros(Complex,7,2)
    # k[1,:]=mi
    # k[2:3,:]=A
    # k[4:5,:]=q
    # k[6:7,:]=g
    # Malha de pontos internos
    return elastico_aniso_iga, PONTOS, SEGMENTOS, MALHA, CCSeg, Compute_Material(Material)
end
function Subregiao1d(ne=15, tipo=2)
    a = 1

    PONTOS = [1 0 0
        2 a 0
        3 a a/2
        4 0 a/2
        5 a a
        6 0 a]

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
        5 4 3 0
        6 3 5 0
        7 5 6 0
        8 6 4 0]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]

    MALHA = [1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida

    carga = 1



    CCSeg = [1 1 0 0 0
        2 1 carga 1 0
        3 1 0 1 0
        4 0 0 1 0
        5 1 0 1 0
        6 1 carga 1 0
        7 1 0 1 0
        8 0 0 1 0]


    E = 1
    v = 0.3

    # @show vcola
    # Malha de pontos internos
    subregioes = define_SubRegioes(SEGMENTOS, MALHA)


    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=[E, E], nu=[v, v]), subregioes
end
function dad_3PBT(ne=6, tipo=2)

    hs = 1.08
    L = 50
    l = 25
    e = 23
    f = 22
    ha = 4
    d = 1

    PONTOS = [1 -L/2 0 0
        2 -e 0 0
        3 -f 0 0
        4 -l/2 0 0
        5 l/2 0 0
        6 f 0 0
        7 e 0 0
        8 L/2 0 0
        9 L/2 hs 0
        10 d/2 hs 0
        11 -d/2 hs 0
        12 -L/2 hs 0
        13 -l/2 -(ha) 0
        14 l/2 -(ha) 0]

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
        4 4 5 0# interface
        5 5 6 0
        6 6 7 0
        7 7 8 0
        8 8 9 0
        9 9 10 0
        10 10 11 0
        11 11 12 0
        12 12 1 0
        13 13 14 0
        14 14 5 0
        15 5 4 0# interface
        16 4 13 0]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [1 2*ne tipo
        2 ne tipo
        3 4*ne tipo
        4 20*ne tipo
        5 4*ne tipo
        6 ne tipo# interface de baixo
        7 2*ne tipo
        8 ne tipo
        9 8*ne tipo# interface de cima
        10 ne tipo
        11 8*ne tipo# interface de baixo
        12 ne tipo
        13 8*ne tipo
        14 4*ne tipo
        15 20*ne tipo
        16 4*ne tipo] # interface de cima
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida

    carga = 1000 / 32


    CCSeg = [1 1 0 1 0 0
        2 0 0 0 0 0
        3 1 0 1 0 0
        4 1 0 1 0 0
        5 1 0 1 0 0
        6 0 0 0 0 0
        7 1 0 1 0 0
        8 0 0 0 0 0
        9 1 0 1 0 0
        10 1 0 1 -carga 0
        11 1 0 1 0 0
        12 1 0 1 0 0
        13 1 0 1 0 0
        14 1 0 1 0 0
        15 1 0 1 0 0
        16 1 0 1 0 0]


    # CCSeg = [1 0 1
    # 2 0 1
    # 3 0 1
    # 4 0 1]
    # Condutividade Térmica do material
    E = 73.1e3
    v = 0.3
    # Ecola = E
    # vcola = 0.3
    Ecola = 1.12e3
    vcola = 0.34
    # @show vcola
    # Malha de pontos internos
    subregioes = define_SubRegioes(SEGMENTOS, MALHA)


    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=[E, Ecola], nu=[v, vcola]), subregioes
end
function Subregiao(ne=15, tipo=2)
    hh = 0.2
    HH = 2
    L = 85
    # % L = 0.85;
    l = 15
    GR = 7.5
    thickness = 25

    PONTOS = [1 -l/2 -(HH + hh / 2) 0
        2 L-l/2 -(HH + hh / 2) 0
        3 (L-l/2) -hh/2 0
        4 l/2 -hh/2 0
        5 -l/2 -hh/2 0
        6 l/2 hh/2 0
        7 l/2 HH+hh/2 0
        8 -L+l/2 HH+hh/2 0
        9 -L+l/2 hh/2 0
        10 -l/2 hh/2 0
        11 L-l/2-GR -(HH + hh / 2) 0
        12 (L-l/2)-GR -hh/2 0
        13 -L+l/2+GR hh/2 0
        14 -L+l/2+GR HH+hh/2 0]

    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 11 0
        2 11 2 0
        3 2 3 0
        4 3 12 0
        5 12 4 0
        6 4 5 0
        7 5 1 0
        8 4 6 0
        9 6 10 0
        10 10 5 0
        11 5 4 0
        12 6 7 0
        13 7 14 0
        14 14 8 0
        15 8 9 0
        16 9 13 0
        17 13 10 0
        18 10 6 0]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    ne = 6
    nead = 4 * ne #% ne da interface do adesivo
    MALHA = [1 2*ne tipo
        2 ne tipo
        3 2 tipo
        4 ne tipo
        5 2*ne tipo
        6 nead tipo
        7 2 tipo
        8 1 tipo
        9 nead tipo
        10 1 tipo
        11 nead tipo
        12 2 tipo
        13 2*ne tipo
        14 ne tipo
        15 2 tipo
        16 ne tipo
        17 2*ne tipo
        18 nead tipo]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida

    carga = 1000

    Area = 2 * thickness * GR

    CCSeg = [1 1 0 1 0 0
        2 1 carga/Area 0 0 0
        3 1 0 1 0 0
        4 1 carga/Area 0 0 0
        5 1 0 1 0 0
        6 1 0 1 0 0
        7 1 0 1 0 0
        8 1 0 1 0 0
        9 1 0 1 0 0
        10 1 0 1 0 0
        11 1 0 1 0 0
        12 1 0 1 0 0
        13 1 0 1 0 0
        14 1 0 0 0 0
        15 0 0 1 0 0
        16 1 0 0 0 0
        17 1 0 1 0 0
        18 1 0 1 0 0]


    # CCSeg = [1 0 1
    # 2 0 1
    # 3 0 1
    # 4 0 1]
    # Condutividade Térmica do material
    E = 73.1e3
    v = 0.3
    # Ecola = E
    # vcola = 0.3
    Ecola = 1.12e3
    vcola = 0.34
    # @show vcola
    # Malha de pontos internos
    subregioes = define_SubRegioes(SEGMENTOS, MALHA)


    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=[E, Ecola, E], nu=[v, vcola, v]), subregioes
end

function telastico1d(ne=15, tipo=2)
    x = 1
    y = 1
    PONTOS = [1 0 0
        2 x 0
        3 x y
        4 0 y]
    # PONTOS = [1 -x/2 -y/2
    #     2 x/2 -y/2
    #     3 x/2 y/2
    #     4 -x/2 y/2]
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
    # CCSeg = [1 1 0 0 0#preso em x
    #     2 0 0 1 0
    #     3 1 0 1 0
    #     4 0 0 1 0]
    # CCSeg = [1 1 0 0 0#preso em x e y
    #     2 0 0 1 0
    #     3 1 0 0 0
    #     4 0 0 1 0]
    CCSeg = [1 0 0 0 0#engastado x e y 
        2 0 0 0 0
        3 0 0 0 0
        4 0 0 0 0]
    # propriedades do material
    E = 1
    v = 0.3
    E = (1 - v^2) * E
    v = v / (1 + v)
    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v)
end

function telasticogao(ne=15, tipo=2; planestress=false)
    H = 1

    PONTOS = [1 -H 0.5
        2 -H -0.5
        3 H -0.5
        4 H 0.5]
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
        2 2ne tipo
        3 ne tipo
        4 2ne tipo]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida
    CCSeg = [1 0 0 1 0
        2 1 0 0 0
        3 0 0 1 0
        4 1 0 1 0]
    # CCSeg = [1 0 1
    # 2 0 1
    # 3 0 1
    # 4 0 1]
    # Condutividade Térmica do material
    E = 10000
    v = 0.3
    k = 1e-5
    if planestress
        E = E * (1 + 2v) / (1 + v)^2
        v = v / (1 + v)
        k = k / (1 + v)
    end

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v, k=k)
end

function telasticocilindro(ne=15, tipo=2; planestress=false)
    ri = 30
    re = 80
    PONTOS = [1 ri 0.0
        2 re 0.0
        3 0 re
        4 0 ri]
    # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS = [1 1 2 0
        2 2 3 re
        3 3 4 0
        4 4 1 -ri]
    # Matriz para definição da malha
    # MALHA = [número do segmento, número de elementos no segmento]
    MALHA = [1 ne tipo
        2 2ne tipo
        3 ne tipo
        4 ne tipo]
    # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDCx, valor da CDCx, tipo da CDCy, valor da CDCy]
    # tipo da CDC = 0 => deslocamento é conhecido
    # tipo da CDC = 1 => força é conhecida
    CCSeg = [1 1 0 0 0
        2 1 0 1 0
        3 0 0 1 0
        4 1 0 1 0]
    Ti = 100
    Te = 0
    CCSegterm = [1 1 0
        2 0 Te
        3 1 0
        4 0 Ti]
    # CCSeg = [1 0 1
    # 2 0 1
    # 3 0 1
    # 4 0 1]
    # Condutividade Térmica do material
    # E = 210e9
    # v = 0.0
    # α = 1.1e-5
    # k = 52.3
    E = 2.1e5
    v = 0.33
    α = 1.1e-5
    k = 60
    if planestress
        E = E * (1 + 2v) / (1 + v)^2
        v = v / (1 + v)
        k = k / (1 + v)
    end
    # Malha de pontos internos
    return [[elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v, k=α)], [potencial, PONTOS, SEGMENTOS, MALHA, CCSegterm, k]]
end


function telasticosquare(ne=15, tipo=2; planestress=false)

    H = 0.5

    PONTOS = [1 -H 0.5
        2 -H -0.5
        3 H -0.5
        4 H 0.5]
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
    CCSeg = [1 0 0 1 0
        2 1 0 0 0
        3 0 0 1 0
        4 1 0 0 0]
    # CCSeg = [1 0 1
    # 2 0 1
    # 3 0 1
    # 4 0 1]
    # Condutividade Térmica do material
    E = 210e9
    v = 0.3
    k = 1.2e-5
    if planestress
        E = E * (1 + 2v) / (1 + v)^2
        v = v / (1 + v)
        k = k / (1 + v)
    end

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v, k=k)
end

function telasticobar(ne=15, tipo=2; planestress=true)

    H = 0.5
    L = 4
    PONTOS = [1 0 0
        2 L 0
        3 L H
        4 0 H]
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
    CCSeg = [1 1 0 0 0
        2 1 0 1 0
        3 1 0 1 0
        4 0 0 1 0]

    # Condutividade Térmica do material
    E = 210e9
    v = 0.3
    k = 1.2e-5
    if planestress
        E = E * (1 + 2v) / (1 + v)^2
        v = v / (1 + v)
        k = k / (1 + v)
    end

    # Malha de pontos internos
    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=E, nu=v, k=k)
end


function Turbina(ne=5, tipo=2)

    PONTOS = [1.0000 0.0000 -23.2303 0
        2.0000 0 0 0
        3.0000 -8.1100 -0.0000 0
        4.0000 -8.7042 -17.9864 0
        5.0000 -11.7050 -22.0073 0
        6.0000 -14.6650 -22.0058 0
        7.0000 -14.6650 -23.2231 0
        8.0000 -12.8085 -6.7101 0
        9.0000 -22.1970 1.7552 0
        10.0000 -14.6970 1.7552 0
        11.0000 -13.4925 -4.8307 0
        12.0000 -9.3883 -16.1070 0
        13.0000 -13.1470 -21.4751 0
        14.0000 -22.1970 -21.4751 0]

    SEGMENTOS = [1 1 2 0
        2 2 3 0
        3 3 8 5
        4 8 12 0
        5 12 4 0
        6 4 5 -3
        7 5 6 0
        8 6 7 0
        9 7 1 0
        10 9 14 0
        11 14 13 0
        12 13 12 4
        13 12 8 0
        14 8 11 0
        15 11 10 -20
        16 10 9 0]

    CCSeg = [1 0 0 1 0
        2 1 0 1 0
        3 1 0 1 0
        4 1 0 1 0
        5 1 0 1 0
        6 1 0 1 0
        7 1 0 1 0
        8 1 0 1 0
        9 1 0 0 0
        10 0 0 1 0
        11 1 0 1 0
        12 1 0 1 0
        13 1 0 1 0
        14 1 0 1 0
        15 1 0 1 0
        16 1 0 1 1]


    MALHA = [1 ne tipo
        2 ne tipo
        3 ne tipo
        4 ne tipo # interface de baixo
        5 ne tipo
        6 ne tipo
        7 ne tipo
        8 ne tipo
        9 ne tipo
        10 ne tipo
        11 ne tipo
        12 ne tipo
        13 ne tipo # interface de baixo
        14 ne tipo
        15 ne tipo
        16 ne tipo
    ]

    NPI = [1 1]

    # Propriedades geométricas e do material
    # Material = [E E G G G ni k h]
    E = 73.1e3 # Módulo de elasticidade
    ni = 0.3
    # Malha de pontos internos
    subregioes = define_SubRegioes(SEGMENTOS, MALHA)


    return elastico, PONTOS, SEGMENTOS, MALHA, CCSeg, (E=[E, E], nu=[ni, ni]), subregioes
end