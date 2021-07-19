# Entrada de dados para análise de temperatura pelo
# método dos elementos de contorno
function helm1d(ne=15,tipo=2)
    PONTOS    = [1 0 0 ;2 1 0 ;3 1 1 ;4 0 1];
	SEGMENTOS = [1 1 2 0 ;2 2 3 0 ;3 3 4 0 ;4 4 1 0]
    MALHA = [1 ne tipo 
            2 ne tipo 
            3 ne tipo 
            4 ne tipo]
    CCSeg=[1 0 0
            2 1 0
            3 0 1
            4 1 0];    # Condutividade Térmica do material
            CCSeg=[1 0 0 0
            2 1 0 0
            3 1 1 0
            4 1 0 0];
        
        
        # RGE=1.E6; # Módulo elástico do meio (RGE=1 no caso de acústica)
        # DAM=0.05; # Amortecimento  (um número bem pequeno no caso de acústica)
        # RO=100; # Densidade do meio  (RO = 1 no caso de acústica)

        RGE=1.E0; 				   # M�dulo el�stico do meio (RGE=1 no caso de ac�stica)
        DAM=0.05; 			   	   # Amortecimento  (um n�mero bem pequeno no caso de ac�stica)
        RO=1; 				   	   # Densidade do meio  (RO = 1 no caso de ac�stica)

        GE=complex(RGE,RGE*2*DAM); # Módulo complexo
        CW=sqrt(GE/RO); # Velocidade de propagação de onda
FR=0.1

        # Malha de pontos internos
    return helmholtz,PONTOS,SEGMENTOS,MALHA,CCSeg,(FR=FR,CW=CW,GE=GE)
end
