function sladek03(ne=15,tipo=2)

POINTS  = [ 1   0     0 
        2   .254  0 
        3   .254  .254
        4   0     .254];
          
# Lines (how the points are joined)
# 
SEGMENTS = [1  1 2  0
	   	    2 2 3 0
         3 3 4 0
	     4 4 1 0];

     
# Discretization (elements per line)

MESH = [1 ne tipo
	    2 ne tipo
        3 ne tipo
        4 ne tipo];

#  Boundary conditions of lines
# BC = [type value type value type value type value]
# 1 line of the matrix per line of the geometry
# Se a CDC imposta for borda livre -> 1
#1,0,1,0
# Se a CDC imposta for borda apoiada -> 2
# 0,0,1,0
#  Se a CDC imposta for borda esgastada -> 3
# 0,0,0,0

# #
# BC_Segments = [1  2
#                2  2
#                3  2
#                4  2];   #simplesmente apoiado
BC_Segments = [1  0 0 0 0
               2  0 0 0 0
               3  0 0 0 0
               4  0 0 0 0
                'c' 0 0 0 0]
                 #engastado
# Boundary condition of corners (thin plate)
# cornerbc = 'f': Corners free
# cornerbc = 'c': Corners clamped

# type_cornerbc = 'c';
aprfun=3; # Type of approximation function
dA=32/2;

# Domain load
A=0;
B=0;
C=2.07e6;

# n_dt = 300;
# dt =10e-3/n_dt;
# thickness   = .0127;
# rho = 7.166e3;
#
E2=.6895e10;
E1=2*E2;
nu12=.3;
G12=E2/(2*(1+nu12));
thickness = .0127;
h=thickness;
nu21=nu12*E2/E1;
rho = 7.166e3;


D22=E2*h^3/(12*(1-nu12*nu21));
a=.254;
#-------------------------------------------------------------------------#
# Termo normalizador
#-------------------------------------------------------------------------#
# to=a^2*sqrt(rho*h/D22)/4;
to=0.02;
#-------------------------------------------------------------------------#
# N�mero(n_dt) de intervalos e Passos(dt) de Tempo
#-------------------------------------------------------------------------#
n_dt =200;

# dt = 0.9*to/n_dt; #simplesmente apoiado
# dt =0.5*to/n_dt; #engastada
dt =to/n_dt; #

#-------------------------------------------------------------------------#
# Matriz de propriedades do material
#-------------------------------------------------------------------------#
Material =[1 E1  E2 G12  nu12  thickness/2  0]


placa_fina,POINTS,SEGMENTS,MESH,BC_Segments,Compute_Material_Placa(Material)
end
