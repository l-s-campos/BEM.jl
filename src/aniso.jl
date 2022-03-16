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




function integrabeziersing(pf,cf,we,eta,w,elem::bezier,dad::elastico_aniso_iga,eet)
   h = zeros(Float64,2,2*size(elem))
   A=dad.k.A
   q=dad.k.q
   g=dad.k.g
   mi=dad.k.mi   


   basisrc,dN=calc_fforma(eet,elem,we)       
   dxdqsi = cf*dN   # dx/dξ & dy/dξ
   fonte = cf*basisrc    # Ponto de gauss interpolador

   dgamadqsif = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano    
   sx=dxdqsi[1]/dgamadqsif # vetor tangente dx/dΓ
   sy=dxdqsi[2]/dgamadqsif # vetor tangente dy/dΓ

   nfonte=[sy,-sx]

   
   
   Nsing,dNsing = calc_fforma(eet,elem,we)
 
   z0 =[fonte[1]+fonte[2]*dad.k.mi[1]*im,fonte[1]+fonte[2]*dad.k.mi[2]*im]
   mi_n=[(dad.k.mi[1]*nfonte[1]-nfonte[2]),(dad.k.mi[2]*nfonte[1]-nfonte[2])]
 

 
   Nm = zeros(Float64,2,2*size(elem))
   Nmsing = zeros(Float64,2,2*size(elem))
   Nmsing[1,1:2:end]=Nsing
   Nmsing[2,2:2:end]=Nsing
   Hsing  = zeros(Float64,2,2*size(elem)) 
   integral=log(1-eet)-log(1+eet)
 
   for k = 1:size(w,1)
     N,dN = calc_fforma(eta[k],elem,we)

     delta=eta[k]-eet;
     kernelsing=[(1/delta) 0
             0 (1/delta)]

     pg = cf*N    # Ponto de gauss interpolador
     dxdqsi = cf*dN   # dx/dξ & dy/dξ
     dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
     sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
     sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
     nx=sy
     ny=-sx
     Nm[1,1:2:end]=N
     Nm[2,2:2:end]=N
     z1 = pg[1] - pf[1] + mi[1] * (pg[2] - pf[2]);
     z2 = pg[1] - pf[1] + mi[2] * (pg[2] - pf[2]); 
     mi_n_z     = [(mi[1]*nx-ny)/z1         0;
                                 0         (mi[2]*nx-ny)/z2];    
     kernelreg = 2*real(A*(mi_n_z)*conj(g)')*dgamadqsi*Nm/2
	 
     kernelexp = 2*real(A*(kernelsing)*conj(g)')*Nmsing
     Hsing=Hsing+(kernelreg-kernelexp)*w[k]
   # @show (kernelreg[2]-kernelexp[2]),eet-eta[k]
   #   @infiltrate

   end

if abs(eet)==1
    beta_m=2/dgamadqsif
    Hsing+=2*real(A*conj(g)')*Nmsing*log(abs(2/beta_m))*sign(-eet)
    #        println("h = $(htermMatrix*log(abs(2/beta_m))*sign(-eet))")
   else
    Hsing+=2*real(A*conj(g)')*Nmsing*log(abs((1-eet)/(1+eet)));
    2*real(A*conj(g)')*Nmsing*log(abs((1-eet)/(1+eet)))
   end   
   Hsing
 end
 
function calsolfund(pg,pf,n,dad::Union{elastico_aniso, elastico_aniso_iga})
   # @infiltrate
   mi=dad.k.mi
   A=dad.k.A
   q=dad.k.q
   g=dad.k.g
 
   xcampo=pg[1]
   ycampo=pg[2]
   xf=pf[1]
   yf=pf[2]
 
 
     #Cálculo da distância do ponto fonte (xf,yf) ao ponto campo
     z1 = xcampo - xf+mi[1]*(ycampo - yf);
     z2 = xcampo - xf+mi[2]*(ycampo - yf);
     
     # Solução fundamental de deslocamento
     
     lns=[log(z1)     0
             0  log(z2)];
     
     
     uast = 2*real(A*lns*conj(q)');
     
     # Solução fundamental de forcas de superficies
     
     mi_n_z=[(mi[1]*n[1]-n[2])/z1         0
                    0          (mi[2]*n[1]-n[2])/z2];
     
     tast = 2*real(A*mi_n_z*conj(g)');
   
   uast,tast
   end