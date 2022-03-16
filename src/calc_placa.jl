function Compute_Material_Placa(Material)

n_laminae=size(Material,1);
D11=0;
D22=0;
D66=0;
D12=0;
D16=0;
D26=0;
hant=0;

for i=1:n_laminae
theta=Material[i,7];
theta=theta*pi/180;
E1=Material[i,2];
E2=Material[i,3];
G12=Material[i,4];
ni12=Material[i,5];
h = Material[i,6]; 
ni21=ni12*E2/E1;



Q11 = E1 / (1 - ni12 * ni21);
Q22 = E2 / (1 - ni12 * ni21);
Q66 = G12;
Q16 = 0;
Q26 = 0;
Q12 = ni21 * E1 / (1 - ni12 * ni21);

Q = [Q11 Q12 Q16
Q12 Q22 Q26
Q16 Q26 Q66];

m = cos(theta);
n = sin(theta);

T = [  m^2  n^2   2*m*n
n^2  m^2  -2*m*n
-m*n  m*n  m^2-n^2];

Q = inv(T) * Q * inv(T');
compliance=inv(Q);

a11=compliance[1,1];
a12=compliance[1,2];
a16=compliance[1,3];
a22=compliance[2,2];
a26=compliance[2,3];
a66=compliance[3,3];

A = [	a11	a12	a16
a12	a22	a26
a16	a26	a66];

del = 1/det(A);

B11 = del * (a22*a66 - a26^2);
B22 = del * (a11*a66 - a16^2);
B66 = del * (a11*a22 - a12^2);
B12 = del * (a16*a26 - a12*a66);
B16 = del * (a12*a26 - a22*a16);
B26 = del * (a12*a16 - a11*a26);

D11 = D11+2*((h+hant)^3-hant^3) / 3 * B11;
D22 = D22+2*((h+hant)^3-hant^3) / 3 * B22;
D66 = D66+2*((h+hant)^3-hant^3) / 3 * B66;
D12 = D12+2*((h+hant)^3-hant^3) / 3 * B12;
D16 = D16+2*((h+hant)^3-hant^3) / 3 * B16;
D26 = D26+2*((h+hant)^3-hant^3) / 3 * B26;
hant=hant+h;

end;
# using PolynomialRoots

charac_poly=[(D11),(4 * D16),(2 * D12 + 4 * D66),(4 * D26),(D22)]; # Characteristic polynomial
# charac_poly=[(D22),(4 * D26),(2 * D12 + 4 * D66),(4 * D16),(D11)]; # Characteristic polynomial

roots_poly=roots(charac_poly);
# @show roots_poly
# Conforme Lekhnitskii, p�gina 28,
# parauma placa isotr�pica:mi1 = mi2 =iemi1c = mi2c = -i

# Pega somente as ra�zes que tem a parte imagin�ria positiva.
# Exclui o conjugado complexo.
aux=zeros(Complex,0)
for i = 1 : 4
 if (imag(roots_poly[i]) > 0)
push!(aux,roots_poly[i])
 end;
end;
mi=zeros(Complex,2)
# Coloca os n�meros em ordem crescente de acordo com a parte real. 
if (real(aux[1]) > real(aux[2]))
 mi[1]=aux[2];
 mi[2]=aux[1];
else
 mi[1]=aux[1];
 mi[2]=aux[2];
end;

d=zeros(2)
e=zeros(2)

d[1] = real(mi[1])
e[1] = imag(mi[1])
d[2] = real(mi[2])
e[2] = imag(mi[2])
(d=d,e=e,D11=D11,D22=D22,D66=D66,D12=D12,D16=D16,D26=D26)

end



function calc_HeG(dad::placa_fina,npg=8)
nelem = size(dad.ELEM,1)# Quantidade de elementos discretizados no contorno
n = size(dad.NOS,1)
H=zeros(2*n,2*n)
G=zeros(2*n,2*n)
qsi,w = gausslegendre(npg)# Quadratura de gauss
for i=1:n
pf = dad.NOS[i,:] # Coordenada (x,y)dos pontos fonte
for elem_j in dad.ELEM#Laço dos elementos
x = dad.NOS[elem_j.indices,:] # Coordenada (x,y) dos nós geométricos
Δelem=x[end,:]-x[1,:] # Δx e Δy entre o primeiro e ultimo nó geometrico
eet=(elem_j.ξs[end]-elem_j.ξs[1])*dot(Δelem,pf.-x[1,:])/norm(Δelem)^2+elem_j.ξs[1]
N_geo,dN=calc_fforma(eet,elem_j)
ps=N_geo'*x
dxdqsi = dN'*x
dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
b=norm(ps'-pf)/norm(Δelem)
eta,Jt=sinhtrans(qsi,eet,b)
# eta,Jt=telles(qsi,eet)
h,g= integraelem(pf,[sy,-sx],x,eta,w.*Jt,elem_j,dad)
nosing =elem_j.indices .== i
if sum(nosing) == 1
no_pf=findfirst(nosing)
xi0=elem_j.ξs[no_pf]
# Hesing=calc_Hsing(xi0,x,dad,elem_j,qsi,w)
#  h[:,2*no_pf-1:2*no_pf]=Hesing;
 h.=0;
end
cols=[2elem_j.indices.-1 2elem_j.indices]'[:]
H[2i-1:2i,cols]=h
G[2i-1:2i,cols]=g
end
end

 # for i = 1:n#i=1:size(dad.NOS,1) #Laço dos pontos fontes
 # H[2i-1:2i,2i-1:2i].=0
 # H[2i-1:2i,2i-1:2i]=-[sum(H[2i-1:2i,1:2:end],dims=2) sum(H[2i-1:2i,2:2:end],dims=2)]
 # end
H,G
end

function calsolfund(pg,pf,n,nf,dad::Union{placa_fina})

rs=pg-pf
nx=n[1]
ny=n[2]
m1=nf[1]
m2=nf[2]
# Distance of source and field points
r1 =rs[1]
r2 = rs[2]
r = norm(rs)
# Thin plate fundamental solutions

theta=atan(r2,r1);


G = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] + dad.k.e[2])^2;

H = (dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1] - dad.k.e[2])^2;

C1 = ((dad.k.d[1] - dad.k.d[2])^2 - (dad.k.e[1]^2 - dad.k.e[1]^2)) / (G * H * dad.k.e[1]);

C2 = ((dad.k.d[1] - dad.k.d[2])^2 + (dad.k.e[1]^2 - dad.k.e[1]^2))/ (G * H * dad.k.e[2]);

C3 = 4 * (dad.k.d[1] - dad.k.d[2]) / (G * H);

a=1;
# a=constFS;

d2Rdx2=zeros(2)
d2Rdxdy=zeros(2)
d2Rdy2=zeros(2)
d2Sdx2=zeros(2)
d2Sdxdy=zeros(2)
d2Sdy2=zeros(2)
d3Rdx2dy=zeros(2)
d3Rdx3=zeros(2)
d3Rdxdy2=zeros(2)
d3Rdy3=zeros(2)
d3Sdx2dy=zeros(2)
d3Sdx3=zeros(2)
d3Sdxdy2=zeros(2)
d3Sdy3=zeros(2)
d4Rdx2dy2=zeros(2)
d4Rdx3dy=zeros(2)
d4Rdx4=zeros(2)
d4Rdxdy3=zeros(2)
d4Rdy4=zeros(2)
d4Sdx2dy2=zeros(2)
d4Sdx3dy=zeros(2)
d4Sdx4=zeros(2)
d4Sdxdy3=zeros(2)
d4Sdy4=zeros(2)
dRdx=zeros(2)
dRdy=zeros(2)
dSdx=zeros(2)
dSdy=zeros(2)
R=zeros(2)
S=zeros(2)
for i=1:2
 R[i]=r^2*((cos(theta)+dad.k.d[i]*sin(theta))^2-dad.k.e[i]^2*sin(theta)^2)*(log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2))-3)-4*r^2*dad.k.e[i]*sin(theta)*(cos(theta)+dad.k.d[i]*sin(theta))*atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));
 
 S[i]=r^2*dad.k.e[i]*sin(theta)*(cos(theta)+dad.k.d[i]*sin(theta))*(log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2))-3)+r^2*((cos(theta)+dad.k.d[i]*sin(theta))^2-dad.k.e[i]^2*sin(theta)^2)*atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));
 
 dRdx[i]=2*r*(cos(theta)+dad.k.d[i]*sin(theta))*(log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+ dad.k.e[i]^2*sin(theta)^2))-2)-4*r*dad.k.e[i]*sin(theta)*atan(dad.k.e[i]*sin(theta), (cos(theta)+dad.k.d[i]*sin(theta)));

 dRdy[i]=2*r*(dad.k.d[i]*(cos(theta)+dad.k.d[i]*sin(theta))-dad.k.e[i]^2*sin(theta))*(log(r^2/a^2*((cos(theta)+ dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2))-2)-4*r*dad.k.e[i]*(cos(theta)+2*dad.k.d[i]*sin(theta))* atan((dad.k.e[i]*sin(theta)),(cos(theta)+dad.k.d[i]*sin(theta)));

 d2Rdx2[i]=2*log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));
 
 d2Rdxdy[i]=2*dad.k.d[i]*log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)) -4*dad.k.e[i]*atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));

 d2Rdy2[i]=2*(dad.k.d[i]^2-dad.k.e[i]^2)*log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+ dad.k.e[i]^2*sin(theta)^2))-8*dad.k.d[i]*dad.k.e[i]*atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));

 d3Rdx3[i]=4*(cos(theta)+dad.k.d[i]*sin(theta))/ (r*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));
 
 d3Rdx2dy[i]=4*(dad.k.d[i]*(cos(theta)+dad.k.d[i]*sin(theta))+dad.k.e[i]^2*sin(theta))/ (r*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));

 d3Rdxdy2[i]=4*((dad.k.d[i]^2-dad.k.e[i]^2)*cos(theta)+(dad.k.d[i]^2+dad.k.e[i]^2)*dad.k.d[i]*sin(theta))/ (r*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));
 
 d3Rdy3[i]=4*(dad.k.d[i]*(dad.k.d[i]^2-3*dad.k.e[i]^2)*cos(theta)+(dad.k.d[i]^4-dad.k.e[i]^4)*sin(theta))/ (r*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));
 
 d4Rdx4[i]=-4*((cos(theta)+dad.k.d[i]*sin(theta))^2-dad.k.e[i]^2*sin(theta)^2)/ (r^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)^2);
 
 d4Rdx3dy[i]=-4/r^2*(dad.k.d[i]/((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2) + 2*dad.k.e[i]^2*sin(theta)*cos(theta)/((cos(theta)+dad.k.d[i]*sin(theta))^2+ dad.k.e[i]^2*sin(theta)^2)^2);

 d4Rdx2dy2[i]=-4/r^2*((dad.k.d[i]^2+dad.k.e[i]^2)/((cos(theta)+dad.k.d[i]*sin(theta))^2+ dad.k.e[i]^2*sin(theta)^2)-2*dad.k.e[i]^2*cos(theta)^2/ ((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)^2);
 
 d4Rdxdy3[i]=-4/r^2*(dad.k.d[i]*(dad.k.d[i]^2+dad.k.e[i]^2)/((cos(theta)+dad.k.d[i]*sin(theta))^2+ dad.k.e[i]^2*sin(theta)^2)-2*dad.k.e[i]^2*cos(theta)*(2*dad.k.d[i]*cos(theta)+ (dad.k.d[i]^2+dad.k.e[i]^2)*sin(theta))/((cos(theta)+dad.k.d[i]*sin(theta))^2+ dad.k.e[i]^2*sin(theta)^2)^2);

d4Rdy4[i]=-4/r^2*((dad.k.d[i]^4-dad.k.e[i]^4)/((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2) - 2*dad.k.e[i]^2*cos(theta)*((3*dad.k.d[i]^2-dad.k.e[i]^2)*cos(theta)+2*dad.k.d[i]*(dad.k.d[i]^2+dad.k.e[i]^2)*sin(theta))/ ((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)^2);
 
 
 
 dSdx[i]=r*dad.k.e[i]*sin(theta)*(log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+ dad.k.e[i]^2*sin(theta)^2))-2)+2*r*(cos(theta)+dad.k.d[i]*sin(theta))* atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));
 
 dSdy[i]=r*dad.k.e[i]*(cos(theta)+2*dad.k.d[i]*sin(theta))*(log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2))-2)+ 2*r*(dad.k.d[i]*(cos(theta)+dad.k.d[i]*sin(theta))-dad.k.e[i]^2*sin(theta))*atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));
 
 d2Sdx2[i]=2*atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));
 
 d2Sdxdy[i]=dad.k.e[i]*(log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)))+ 2*dad.k.d[i]*atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));
 
 
 d2Sdy2[i]=2*dad.k.d[i]*dad.k.e[i]*log(r^2/a^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2))+ 2*(dad.k.d[i]^2-dad.k.e[i]^2)*atan(dad.k.e[i]*sin(theta),(cos(theta)+dad.k.d[i]*sin(theta)));

 d3Sdx3[i]=-2*dad.k.e[i]*sin(theta)/(r*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));


 d3Sdx2dy[i]=2*dad.k.e[i]*cos(theta)/(r*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));

 d3Sdxdy2[i]=2*dad.k.e[i]*(2*dad.k.d[i]*(cos(theta)+dad.k.d[i]*sin(theta))-(dad.k.d[i]^2-dad.k.e[i]^2)*sin(theta))/ (r*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));

 d3Sdy3[i]=2*dad.k.e[i]*((3*dad.k.d[i]^2-dad.k.e[i]^2)*cos(theta)+2*dad.k.d[i]*(dad.k.d[i]^2+dad.k.e[i]^2)*sin(theta))/ (r*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2));

 d4Sdx4[i]=4*dad.k.e[i]*sin(theta)*(cos(theta)+dad.k.d[i]*sin(theta))/ (r^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)^2);

 d4Sdx3dy[i]=2*dad.k.e[i]/r^2*(1/((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)- 2*cos(theta)*(cos(theta)+dad.k.d[i]*sin(theta))/ ((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)^2);
 
 d4Sdx2dy2[i]=-4*dad.k.e[i]*cos(theta)*(dad.k.d[i]*(cos(theta)+dad.k.d[i]*sin(theta))+dad.k.e[i]^2*sin(theta))/ (r^2*((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)^2);

 d4Sdxdy3[i]=-2*dad.k.e[i]/r^2*((dad.k.d[i]^2+dad.k.e[i]^2)/((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)+ (2*(dad.k.d[i]^2+dad.k.e[i]^2)*cos(theta)*(cos(theta)+dad.k.d[i]*sin(theta))-4*dad.k.e[i]^2*cos(theta)^2)/ ((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)^2);
 

 d4Sdy4[i]=-4*dad.k.e[i]/r^2*(dad.k.d[i]*(dad.k.d[i]^2+dad.k.e[i]^2)/((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)+ cos(theta)*(dad.k.d[i]*(dad.k.d[i]^2-3*dad.k.e[i]^2)*cos(theta)+(dad.k.d[i]^4-dad.k.e[i]^4)*sin(theta))/ ((cos(theta)+dad.k.d[i]*sin(theta))^2+dad.k.e[i]^2*sin(theta)^2)^2); 
end


w=1/(8*pi*dad.k.D22)*(C1*R[1]+C2*R[2]+C3*(S[1]-S[2]));

dwdx=1/(8*pi)*(C1*dRdx[1]+C2*dRdx[2]+C3*(dSdx[1]-dSdx[2]));
dwdy=1/(8*pi)*(C1*dRdy[1]+C2*dRdy[2]+C3*(dSdy[1]-dSdy[2]));
 
d2wdx2=1/(8*pi)*(C1*d2Rdx2[1]+C2*d2Rdx2[2]+C3*(d2Sdx2[1]-d2Sdx2[2]));
d2wdxdy=1/(8*pi)*(C1*d2Rdxdy[1]+C2*d2Rdxdy[2]+C3*(d2Sdxdy[1]-d2Sdxdy[2]));
d2wdy2=1/(8*pi)*(C1*d2Rdy2[1]+C2*d2Rdy2[2]+C3*(d2Sdy2[1]-d2Sdy2[2]));

d3wdx3=1/(8*pi)*(C1*d3Rdx3[1]+C2*d3Rdx3[2]+C3*(d3Sdx3[1]-d3Sdx3[2]));
d3wdx2dy=1/(8*pi)*(C1*d3Rdx2dy[1]+C2*d3Rdx2dy[2]+C3*(d3Sdx2dy[1]-d3Sdx2dy[2]));
d3wdxdy2=1/(8*pi)*(C1*d3Rdxdy2[1]+C2*d3Rdxdy2[2]+C3*(d3Sdxdy2[1]-d3Sdxdy2[2]));
d3wdy3=1/(8*pi)*(C1*d3Rdy3[1]+C2*d3Rdy3[2]+C3*(d3Sdy3[1]-d3Sdy3[2]));
 
d4wdx4=1/(8*pi)*(C1*d4Rdx4[1]+C2*d4Rdx4[2]+C3*(d4Sdx4[1]-d4Sdx4[2]));
d4wdx3dy=1/(8*pi)*(C1*d4Rdx3dy[1]+C2*d4Rdx3dy[2]+C3*(d4Sdx3dy[1]-d4Sdx3dy[2]));
d4wdx2dy2=1/(8*pi)*(C1*d4Rdx2dy2[1]+C2*d4Rdx2dy2[2]+C3*(d4Sdx2dy2[1]-d4Sdx2dy2[2]));
d4wdxdy3=1/(8*pi)*(C1*d4Rdxdy3[1]+C2*d4Rdxdy3[2]+C3*(d4Sdxdy3[1]-d4Sdxdy3[2]));
d4wdy4=1/(8*pi)*(C1*d4Rdy4[1]+C2*d4Rdy4[2]+C3*(d4Sdy4[1]-d4Sdy4[2]));

f1=dad.k.D11*nx^2+2*dad.k.D16*nx*ny+dad.k.D12*ny^2;
f2=2*(dad.k.D16*nx^2+2*dad.k.D66*nx*ny+dad.k.D26*ny^2);
f3=dad.k.D12*nx^2+2*dad.k.D26*nx*ny+dad.k.D22*ny^2;

 
 
h1=dad.k.D11*nx*(1+ny^2)+2*dad.k.D16*ny^3-dad.k.D12*nx*ny^2;
h2=4*dad.k.D16*nx+dad.k.D12*ny*(1+nx^2)+4*dad.k.D66*ny^3-dad.k.D11*nx^2*ny-2*dad.k.D26*nx*ny^2;
h3=4*dad.k.D26*ny+dad.k.D12*nx*(1+ny^2)+4*dad.k.D66*nx^3-dad.k.D22*nx*ny^2-2*dad.k.D16*nx^2*ny;
h4=dad.k.D22*ny*(1+nx^2)+2*dad.k.D26*nx^3-dad.k.D12*nx^2*ny;



dwdn=(dwdx*nx+dwdy*ny)/dad.k.D22;
mn=-(f1*d2wdx2+f2*d2wdxdy+f3*d2wdy2)/dad.k.D22;
vn=-(h1*d3wdx3+h2*d3wdx2dy+h3*d3wdxdy2+h4*d3wdy3)/dad.k.D22;


dmndx=-(f1*d3wdx3+f2*d3wdx2dy+f3*d3wdxdy2);
dmndy=-(f1*d3wdx2dy+f2*d3wdxdy2+f3*d3wdy3);

dvndx=-(h1*d4wdx4+h2*d4wdx3dy+h3*d4wdx2dy2+h4*d4wdxdy3);
dvndy=-(h1*d4wdx3dy+h2*d4wdx2dy2+h3*d4wdxdy3+h4*d4wdy4);

dwdm=-(dwdx*m1+dwdy*m2)/dad.k.D22;
d2wdndm=-(d2wdx2*nx*m1+d2wdxdy*(nx*m2+ny*m1)+d2wdy2*ny*m2)/dad.k.D22;
dmndm=-(dmndx*m1+dmndy*m2)/dad.k.D22;
dvndm=-(dvndx*m1+dvndy*m2)/dad.k.D22; 

# Assembly of matrices that contain fundamental solutions.
u = [w -dwdn
 dwdm -d2wdndm];

 p = [ vn -mn
 dvndm -dmndm];
u,p
end


function integraelem(pf,nf,x,eta,w,elem,dad::placa_fina)
    h = zeros(Float64,2,2*size(elem))
    g = zeros(Float64,2,2*size(elem))
    Nm = zeros(Float64,2,2*size(elem))
    for k = 1:size(w,1)
      N,dN=calc_fforma(eta[k],elem)
        pg = N'*x    # Ponto de gauss interpolador
        dxdqsi = dN'*x   # dx/dξ & dy/dξ
        dgamadqsi = norm(dxdqsi)  # dΓ/dξ = J(ξ) Jacobiano
        sx=dxdqsi[1]/dgamadqsi # vetor tangente dx/dΓ
        sy=dxdqsi[2]/dgamadqsi # vetor tangente dy/dΓ
        uast,tast=calsolfund(pg',pf,[sy,-sx],nf,dad)
        Nm[1,1:2:end]=N
        Nm[2,2:2:end]=N
        h+=tast*Nm*dgamadqsi*w[k]
        g+=uast*Nm*dgamadqsi*w[k]
        # @infiltrate
  
  end
  h,g
  end