function  Compute_Material(Material)
   n_laminae = length(Material[:,1]);
   stiffness = zeros(Float64,3,3);
   aux=zeros(Complex, 2)
   mi=zeros(Complex, 2)
   A=zeros(Complex,2,2)
   total_layers = 0;  # Número de lâminas

   for i = 1:n_laminae
      theta = Material[i,6];
      n_layers = Material[i,7];
      total_layers = total_layers+n_layers;
      E1 = Material[i,2];
      E2 = Material[i,3];
      G12 = Material[i,4];
      nu12 = Material[i,5];
      nu21 = nu12*E2/E1;
      #   if [[1-nu12*nu21]<=0]
      #       error('Error in at least one of laminae constants')
      #   end;  # Verificar a função erro se ela é escrtia dessa forma
      Q = Assembly_Q(E1,E2,G12,nu12);
      Qbar = Compute_Qbar(Q,theta);
      stiffness = stiffness + Qbar*2*n_layers;
   end;
   stiffness = stiffness./[2*total_layers];

   compliance = inv(stiffness);
   a11 = compliance[1,1];
   a12 = compliance[1,2];
   a16 = compliance[1,3];
   a22 = compliance[2,2];
   a26 = compliance[2,3];
   a66 = compliance[3,3];
#   charac_poly = [a11, -2*a16, 2*a12+a66 ,-2*a26, a22]; # Characteristic polynomial
   charac_poly = [a22, -2*a26, 2*a12+a66 ,-2*a16, a11]; # Characteristic polynomial
   roots_poly = roots(charac_poly); # Verificar essa função  
   j = 1;
   for i = 1:4
      if(imag(roots_poly[i])>0)
         aux[j] = roots_poly[i]
         j = j+1;
      end;
   end;

   if(real(aux[1])>real(aux[2]))
      mi[1] = aux[2];
      mi[2] = aux[1];
   else
      mi[1] = aux[1];
      mi[2] = aux[2];
   end;

   q = [a11*mi[1]^2+a12-a16*mi[1] a11*mi[2]^2+a12-a16*mi[2]
   a12*mi[1]+a22/mi[1]-a26     a12*mi[2]+a22/mi[2]-a26];
   Matriz = [ 1           -1            1            -1
      mi[1]       -conj(mi[1])      mi[2]        -conj(mi[2])
      q[1,1]  -conj(q[1,1])   q[1,2]    -conj(q[1,2])
      q[2,1]   -conj(q[2,1])   q[2,2]    -conj(q[2,2])];
   b1 = [0; -1/(2*pi*sqrt(Complex(-1))); 0; 0];
   b2 = [1/(2*pi*sqrt(Complex(-1))); 0; 0; 0];
   A1 = inv(Matriz)*b1;
   A2 = inv(Matriz)*b2;
   A[1,1] = A1[1];
   A[1,2] = A1[3];
   A[2,1] = A2[1];
   A[2,2] = A2[3];
   g = [mi[1] mi[2]
      -1     -1];
return (mi=mi,A=A,q=q,g=g)
end

function  Compute_T(theta)
   theta = theta*pi/180;
   m = cos(theta);
   n = sin(theta);
   T = [  m^2  n^2   2*m*n
         n^2  m^2  -2*m*n
         -m*n  m*n  m^2-n^2];
return T
end

function  Compute_Qbar(Q,theta)
   T = Compute_T(theta);
   Qbar = inv(T)*Q*inv(T');
return Qbar
end

function Assembly_Q(El,Et,Glt,nult)
   nutl = nult*Et/El;
   Q11 = El/[1-nult*nutl];
   Q22 = Et/[1-nult*nutl];
   Q66 = Glt;
   Q16 = 0;
   Q26 = 0;
   Q12 = nutl*El/[1-nult*nutl];
   Q = [Q11 Q12 Q16
      Q12 Q22 Q26
      Q16 Q26 Q66];
return  Q
end