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
    return potencial_iga,PONTOS,SEGMENTOS,MALHA,CCSeg,k
end

function placa482_iga(ne=8,tipo=2)
    PONTOS=[1 -1 0
    2 0 0
    3 1 0
    4 1 1
    5 -1 1]
 # Segmentos que definem a geometria
    # SEGMENTOS = [N° do segmento, N° do ponto inicial, N° do ponto final
    #                                                  Raio, tipo do elemento]
    # Raio do segmento: > 0 -> O centro é a esquerda do segmento (do ponto
    #                          inicial para o ponto final)
    #                   < 0 -> O centro é a direita do segmento (do ponto
    #                          inicial para o ponto final)
    #                   = 0 -> O segmento é uma linha reta
    SEGMENTOS=[1 1 2 0
    2 2 3 0
    3 3 4 0
    4 4 5 0
    5 5 1 0]

    MALHA = [1 floor(Int,ne/2) tipo 
    2 floor(Int,ne/2) tipo 
    3 ne tipo 
    4 ne tipo
    5 ne tipo]
# Condições de contorno nos segmentos
# CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
# tipo da CDC = 0 => A temperatura é conhecida
# tipo da CDC = 1 => O fluxo é conhecido

        # Condições de contorno nos segmentos
    # CCSeg = [N° do segmento, tipo da CDC, valor da CDC]
    # tipo da CDC = 0 => A temperatura é conhecida
    # tipo da CDC = 1 => O fluxo é conhecido
    CCSeg = [1 0 0
        2 1 0
        3 1 0
        4 1 0
        5 1 0]
        # CCSeg = [1 0 1
        # 2 0 1
        # 3 0 1
        # 4 0 1]
    # Condutividade Térmica do material
    k = 1
    return potencial_iga,PONTOS,SEGMENTOS,MALHA,CCSeg,k
end

function aplicaCDCplaca_iga(dad)
	# floor(Int,ne/2)
	NELEM=length(dad.ELEM[:,1])
	NN=length(dad.NOS[:,1])
	N=floor(length(dad.NOS[:,1])/4)
	np=length(dad.valorCDC)
		
		for j=1:np
			if j>=(1) && j<(N)
				# k=length(dad.ELEM[i].indices)
				# for j=1:k
				# 	# @show 12
					dad.valorCDC[j]=dad.valorCDC[j]
				# end
			end
			if j>=(N+1) && j<(2*N)
				# k=length(dad.ELEM[i].indices)
				# for f=1:k
					# NPonto=dad.ELEM[i].indices[f]
					x,y=dad.NOS[j,:]
					r=sqrt(x^2+y^2)
					xc=0;
					yc=0;
					x1=1;
					y1=0;
					teta=atan(y,x);
					q=-(0.5*sqrt(r))*(cos(0.5*teta)*cos(teta)+sin(0.5*teta)*sin(teta));
					dad.valorCDC[j]=q
					
					#  @show NPonto,q
				# end
			end
			if j>=(2*N+1) && j<(3*N)
				# k=length(dad.ELEM[i].tipoCDC)
				# for f=1:k
					# NPonto=dad.ELEM[i].indices[f]
					x,y=dad.NOS[j,:]
					r=sqrt(x^2+y^2)
					xc=0;
					yc=0;
					x1=1;
					y1=0;
					teta=atan(y,x);
					q=-(0.5*sqrt(r))*(cos(0.5*teta)*sin(teta)+sin(0.5*teta)*cos(teta));
					dad.valorCDC[j]=q
					# @show NPonto
				# end
			end
			if j>=(3*N+1) && j<(4*N)
				# k=length(dad.ELEM[i].tipoCDC)
				# for f=1:k
				# 	NPonto=dad.ELEM[i].indices[f]
					x,y=dad.NOS[j,:]
					r=sqrt(x^2+y^2)
					xc=0;
					yc=0;
					x1=1;
					y1=0;
					teta=atan(y,x);
					q=-(0.5*sqrt(r))*(cos(0.5*teta)*cos(teta)+sin(0.5*teta)*sin(teta));
					dad.valorCDC[j]=q
					# @show NPonto
				# end
			end
		
	end
	dad
end

function placa482treal(dad)
	NELEM=length(dad.ELEM[:,1])
	NN=length(dad.NOS[:,1])
	T=zeros(length(dad.NOS[:,1]))
	for i=1:NN
		x,y=dad.NOS[i,:]
		r=sqrt(x^2+y^2)
		xc=0;
		yc=0;
		x1=1;
		y1=0;
		teta1,teta=BEM.calcula_arco(x1,y1,x,y,xc,yc,r);
		T[i]=sqrt(r)*cos(teta/2)
	end
	return T
end