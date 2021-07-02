# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno
function helm1d(ne=15,tipo=2)
    PONTOS    = [1 0 0 ;2 1 0 ;3 1 1 ;4 1 0];
	SEGMENTOS = [1 1 2 0 ;2 2 3 0 ;3 3 4 0 ;4 4 1 0]
    MALHA = [1 ne tipo 
            2 ne tipo 
            3 ne tipo 
            4 ne tipo]
    CCSeg=[1 0 0
            2 1 0
            3 0 1
            4 1 0];    # Condutividade Térmica do material
    k = 1
    # Malha de pontos internos
    return helmholtz,PONTOS,SEGMENTOS,MALHA,CCSeg,k
end
