function testa_ponto(xpi,ypi,xl1,yl1,xl2,yl2,lx,raio)
# Rotina para verifica��o da condi��o de um ponto com rela��o
# � geometria (interno ou externo).
# O algoritmo utilizado � baseado no n�mero de intersec��es entre
# o contorno e uma linha vertical que parte do ponto interno.
# N. �mpar de intersec��es - ponto interno
# N. par de intersec��es   - ponto externo
#
#   Autor: Frederico Lourenao
#   Data de cria��o setembro de 1999
#   Revis�o 0.0


pc = 0;					# Contador de pontos do contorno
npc = length(xl1);	# N. de pontos do contorno
sai = false

while (sai==false && pc < npc)
  pc = pc + 1
  if xpi == xl1[pc]
    xpi = xpi + lx*1e-2
    sai = true
  end
end

interv = 0;		# N. de intersec��es entre o contorno e a linha
					# vertical abaixo do ponto interno.
l=1
while l<= npc; # for over the lines that defines the boundary
   x1 = xl1[l]; y1 = yl1[l];
   x2 = xl2[l]; y2 = yl2[l]
   if(raio[l]==0) # The segment is a straight line
      if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1) # Check only if the x coordinate of the point is
                                        # between the x coordinate of the intial and final points of the recta
         m = (y2-y1)/(x2-x1);		# Angular coeficient of the recta
         yi = m*(xpi-x1)+y1;			# y coordinate in the intersection
         if ypi  >yi # If this is true; the point is above the line
            interv = interv + 1;	# Counter the intersection
         end
      end
   else   # The segment is an arc
      xc,yc=calcula_centro(x1,y1,x2,y2,raio[l]); # compute the center of the arc
      if(xpi<xc+abs(raio[l]) && xpi>xc-abs(raio[l])) # check only if the x coordinate of the point is between
                 # the values xc-radius and xc+radius
         teta_i,teta_f=calcula_arco(x1,y1,x2,y2,xc,yc,raio[l]); # compute the arc between the line that defines the arc and
         # the horizontal direction (-π<teta<π)
         tetac=zeros(2,1)
         y=zeros(2,1)
         tetac[1]=acos((xpi-xc)/abs(raio[l])); # first intersection of the horizontal line that cross the
         # point with the circunference (angle between the radius that cross the intersection point and
         # the horizontal direction)
         tetac[2]=-tetac[1];  # angle of the second intersection
         y[1]=abs(raio[l])*sin(tetac[1])+yc; # y coordinate of the first intersection point
         y[2]=abs(raio[l])*sin(tetac[2])+yc; # y coordinate of the second intersection point
         for k=1:2 # check if the angles of the two intersection points are between the angle of the initial and the
            # final angles of the arc. If yes so the vertical line from the the point to the botton direction
            # intercept the arc (the counter should be increased).
            if(raio[l]>0) # the center is on the left of the arc (from the initial point to the final point)
               if(teta_f>teta_i)
                  if(tetac[k]>teta_i && tetac[k]<teta_f)
                     if(y[k]<ypi)
                        interv=interv+1
                     end
                  end
               else # teta_f<teta_i
                  if(tetac[k]>teta_i || tetac[k]<teta_f)
                     if(y[k]<ypi)
                        interv=interv+1
                     end
                  end
               end
            else # raio[l] < 0 the center is on the right of the arc (from the initial point to the final point)
               if(teta_i > teta_f)
                  if(tetac[k]>teta_f && tetac[k]<teta_i)
                     if(y[k]<ypi)
                        interv=interv+1
                     end
                  end
               else # teta_i < teta_f
                  if(tetac[k]>teta_f || tetac[k]<teta_i)
                     if(y[k]<ypi)
                        interv=interv+1
                     end
                  end
               end
            end
         end
      end
   end
   l=l+1
end


if rem(interv,2) != 0	# Resto da divis�o de interv por 2
  ponto = "interno"
else
  ponto = "externo"
end
return  xpi,ponto
end

function  aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio)
  # Fun??o para testar a proximidade de um ponto interno ao contorno
  # Se o ponto estiver mais proximo que o desejado (d<d_min) o ponto
  # ? reprovado (aceita = nao). Caso contrr?rio (aceita = sim)


pc = 0;					# Contador de pontos do contorno
npc = length(xl1);	# N. de pontos do contorno
aceita = "sim"

# Verificando a proximidade do ponto interno aos pontos do contorno
while (aceita=="sim"  && pc < npc)
  pc = pc + 1
  x1 = xl1[pc]
  y1 = yl1[pc]
  d = √((xpi-x1)^2+(ypi-y1)^2)
  if d < d_min
    aceita = "nao"
  end
end

# Verificando a proximidade ?s linhas do contorno
l = 0;	# Contador das linhas do contorno

while (aceita=="sim" && l < npc)
   l = l + 1
   x1 = xl1[l]; y1 = yl1[l]
   x2 = xl2[l]; y2 = yl2[l]
   if(raio[l]==0) # The segment is a straight line

      if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1)
         m = (y2-y1)/(x2-x1);		# Angular coeficient of the recta
         yi = m*(xpi-x1)+y1;			# y coordinate of the intersection of a recta normal to the current boundary
         # recta that cross the point
         dy = ypi - yi
         d = abs(dy*cos(atan(m))); # distance from the point to the recta
         if d < d_min
            aceita = "nao"
         end
      end

      if (ypi > y1 && ypi < y2) || (ypi > y2 && ypi < y1)
         if x1 == x2
            d = abs(xpi-x1)
         else
            m = (y2-y1)/(x2-x1);			# Angular coeficient of the recta
            xi = 1/m*(ypi-y1)+x1;	# x coordinate of the intersection of a recta normal to the current boundary
    # recta that cross the point
            dx = xpi - xi
            d = abs(dx*sin(atan(m))); # distance from the point to the recta
         end
         if d < d_min
            aceita = "nao"
         end
      end
   else    # The segment is an arc
      xc,yc=calcula_centro(x1,y1,x2,y2,raio[l]); # Center of the arc
      teta_i,teta_f=calcula_arco(x1,y1,x2,y2,xc,yc,raio[l]); # angle of the lines that defines the arc with the horizontal direction
      teta_i,teta_p=calcula_arco(x1,y1,xpi,ypi,xc,yc,raio[l]);  # teta_p angle of the line that cross the center point and the
      # internal point with the horizontal direction
      if(raio[l]>0) # The center is in the left side of the arc (from the initial to the end point)
         if(teta_f>teta_i)
            if(teta_p>teta_i && teta_p<teta_f)
               d=abs(raio[l]-√((xpi-xc)^2+(ypi-yc)^2)); # distance from the point to the arc
               if d < d_min
                  aceita = "nao"
               end
            end
         else # teta_f<teta_i
            if(teta_p>teta_i || teta_p<teta_f)
               d=abs(raio[l]-√((xpi-xc)^2+(ypi-yc)^2)); # distance from the point to the arc
               if d < d_min
                  aceita = "nao"
               end
            end
         end
      else # raio[l] < 0 # The center is in the right side of the arc (from the initial to the end point)
         if(teta_i > teta_f)
           if(convert(Float64,teta_p)>convert(Float64,teta_f) && convert(Float64,teta_p)<convert(Float64,teta_i))
               d=abs(abs(raio[l])-√((xpi-xc)^2+(ypi-yc)^2)); # distance from the point to the arc
               if d < d_min
                  aceita = "nao"
               end
            end
         else # teta_i < teta_f
            if(teta_p>teta_f || teta_p<teta_i)
               d=abs(abs(raio[l])-√((xpi-xc)^2+(ypi-yc)^2)); # distance from the point to the arc
               if d < d_min
                  aceita = "nao"
               end
            end
         end
      end
   end
end
# Retorna a aprova??o ou nao do ponto interno
return aceita
end

function gera_p_in(NPX,NPY,PONTO,SEGMENTOS)
  #PONTOS_INT
# Programa para cria��o de pontos internos a uma  geometria
# gen�rica formada por retas
#
#   Autor: Frederico Lourenao
#   Data de cria��o setembro de 1999
#   Revis�o 0.0

# Defini��o da �rea m�xima para cria��o de pontos internos
xmin = minimum(PONTO[:,2])
xmax = maximum(PONTO[:,2])
ymin = minimum(PONTO[:,3])
ymax = maximum(PONTO[:,3])
lx = xmax - xmin;		# Largura do ret�ngulo que cont�m a geometria
ly = ymax - ymin;		# Altura do ret�ngulo que cont�m a geometria
n_SEGMENTOSs = length(SEGMENTOS[:,1])

# Defini��o da maior SEGMENTOS do problema
max_dl = 0
for lin = 1 : length(SEGMENTOS[:,1])
    p1 = round(UInt64,SEGMENTOS[lin,2])
    p2 = round(UInt64,SEGMENTOS[lin,3])
    xp1 = PONTO[p1,2]
    yp1 = PONTO[p1,3]
    xp2 = PONTO[p2,2]
    yp2 = PONTO[p2,3]
    dl = √((xp1-xp2)^2+(yp1-yp2)^2)
    if dl > max_dl
        max_dl = dl
    end
end


d_min = 0.003*max_dl;	# Dist�ncia m�nima dos pontos internos ao contorno
npx = NPX+1;				# N. de pontos na horizontal
npy = NPY+1;				# N. de pontos na vertical

PONTOS_INT =zeros(NPX*NPY,2)
# Atribui��o dos pontos finais e iniciais das SEGMENTOSs aos
# vetores xl1; xl2; yl1 e yl2
xl1=zeros(n_SEGMENTOSs)
xl2=zeros(n_SEGMENTOSs)
yl1=zeros(n_SEGMENTOSs)
yl2=zeros(n_SEGMENTOSs)
raio=zeros(n_SEGMENTOSs)

for t = 1 : n_SEGMENTOSs		# Percorre todas as SEGMENTOSs
    xl1[t] = PONTO[round(UInt64,SEGMENTOS[t,2]),2]
    xl2[t] = PONTO[round(UInt64,SEGMENTOS[t,3]),2]
    yl1[t] = PONTO[round(UInt64,SEGMENTOS[t,2]),3]
    yl2[t] = PONTO[round(UInt64,SEGMENTOS[t,3]),3]
    raio[t]= SEGMENTOS[t,4]
end

npi = 0;	# Inicializa��o no n. de pontos internos
for i = 1 : NPY
    # Cria��o do candidato a ponto interno (xpi,ypi)
    ypi = ymin + (ly/npy)*i;	# y dentro do ret�ngulo
    for j = 1 : NPX
        xpi = xmin + (lx/npx)*j;	# x dentro do ret�ngulo

        # In�cio dos testes para valida��o do ponto interno

        # 1. Verificando se o ponto est� dentro da geometria
        xpi,ponto = testa_ponto(xpi,ypi,xl1,yl1,xl2,yl2,lx,raio)
        # 2. Verificando se o ponto est� muito pr�ximo do contorno
        if (ponto=="interno")
            aceita = aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio)
        else
            aceita = "nao"
        end
        # Armazenando os dados do ponto interno
        if (aceita=="sim")	# O ponto est� dentro da geometria e
            npi = npi + 1;		# e respeita a dist�ncia ao contorno
            PONTOS_INT[npi,:] = [xpi ypi]

        end
    end
end
PONTOS_INT=PONTOS_INT[1:npi,:]
return  PONTOS_INT
end
# Nesse ponto est�o calculados os pontos internos
