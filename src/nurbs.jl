mutable struct Curve
    knots
    coefs
    number
    order
    C
    conn
    range
    fontes
end
mutable struct Surf
    uknots
    vknots
    coefs
    number
    order
    C
    conn
    range
    fontes
end
mutable struct Fontes
    pts
    coords
    basis
    normal
end

###############################################################################

function nrbmak(coefs,knots)
    np = size(coefs);
    dim = np[1];
    number = np[2];
    if (dim < 3)
        # nurbs.coefs = repmat([0.0 0.0 1.0]',[1 np[2]]);
        # nurbs.coefs[1:dim,:] = coefs;
        coefs=[coefs;ones(1,np[2])]
    end
    order = length(knots)-np[2]; # m=n+p+1  <=> m-n=p+1
#    println(knots)
#    knots = sort(knots,dims=1);  # sort(A,dims=1 ou 2 resp colunas linhas) 
    knots = (knots.-knots[1])/(knots[end]-knots[1]);  
    C,nb,conn,uk=bezierExtraction(knots,order-1)

    range=[uk[1:end-1] uk[2:end]]
#    x1,wx1=gausslegendre(number)   # -1 a 1  
collocPts=zeros(number,1);
z=0;
p=order-1;

collocPts=greville(knots,p,0.25)

    fontes=Array{Fontes}(undef,size(collocPts,1))
    nurbs=Curve(knots,coefs,number,order,C,conn,range,0)
    for i=1:length(collocPts)
        collocCoord,B,normal= evalcrv(nurbs,collocPts[i])  # B são as funções de forma Nurbs
        fontes[i]=Fontes(collocPts[i],collocCoord,B,normal)
    end
    nurbs=Curve(knots,coefs,number,order,C,conn,range,fontes)
end

function   bspderiv(d,c,k) 
    #   bspdeval - Evaluate the control points and knot sequence of the derivative 
    #              of a univariate B-Spline.
    # Calling Sequence:
    #   dc,dk = bspderiv(d,c,k)
    # Parameters:
    #   d	: Degree of the B-Spline.
    #   c	: Control Points, matrix of size (dim,nc).
    #   k	: Knot sequence, row vector of size nk.
    #   dc	: Control points of the derivative
    #   dk	: Knot sequence of the derivative
    # Description:
    #   Evaluate the derivative of univariate B-Spline, which is itself a B-Spline.

    mc,nc = size(c)
    nk = length(k) 
    dc = zeros(mc,nc-1)
    for i=0:nc-2
        tmp = d / (k[i+d+2] - k[i+2])
        for j=0:mc-1
            dc[j+1,i+1] = tmp*(c[j+1,i+2] - c[j+1,i+1])
        end
    end
    dk = zeros(1,nk-2)
    for i=1:nk-2
        dk[i] = k[i+1]
    end
    return dc,dk
end


function nrbeval(nurbs,tt,basis=0)
    # Function Name:
    #
    #   nrbeval - Evaluate a NURBS at parameteric points
    #
    # Calling Sequence:
    #
    #   p = nrbeval(crv,ut)
    #   [p,w] = nrbeval(srf,{ut,vt})
    #
    # Parameters:
    #
    #   crv		: NURBS curve, see nrbmak.
    #
    #   srf		: NURBS surface, see nrbmak.
    #
    #   ut		: Parametric evaluation points along U direction.
    #
    #   vt		: Parametric evaluation points along V direction.
    #
    #   p		: Evaluated points on the NURBS curve or surface as cartesian
    # 		coordinates (x,y,z). If w is included on the lhs argument list
    # 		the points are returned as homogeneous coordinates (wx,wy,wz).
    #
    #   w		: Weights of the homogeneous coordinates of the evaluated
    # 		points. Note inclusion of this argument changes the type
    # 		of coordinates returned in p (see above).

    # NURBS structure represents a curve
    #  tt represent a vector of parametric points in the u direction

    val = bspeval(nurbs.order-1,nurbs.coefs,nurbs.knots,tt);
#     println(" val = $val ")
    w = val[4,:];
    #println(" w = $w ")
    p = val[1:3,:];   # EmerEdited 1:2  para 1:3
    #println(" p = $p ")
    p = p./repeat(w',3,1); # EmerEdited 2,1  para 3,1

    p
end

###############################################################################

function  bspeval(d,c,k,u)
    #
    # Function Name:
    #
    #   bspeval - Evaluate a univariate B-Spline.
    #
    # Calling Sequence:
    #
    #   p = bspeval(d,c,k,u)
    #
    # Parameters:
    #
    #   d	: Degree of the B-Spline.
    #
    #   c	: Control Points, matrix of size (dim,nc).
    #
    #   k	: Knot sequence, row vector of size nk.
    #
    #   u	: Parametric evaluation points, row vector of size nu.
    #
    #   p	: Evaluated points, matrix of size (dim,nu)
    #
    # Description:
    #
    #   Evaluate a univariate B-Spline. This function provides an interface to
    #   a toolbox 'C' routine.
    
    nu = length(u);
    mc,nc = size(c);
    #   int bspeval(int d, double *c, int mc, int nc, double *k, int nk, double *u,int nu, double *p){
    #   int ierr = 0;
    #   int i, s, tmp1, row, col;
    #   double tmp2;
    #
    #   // Construct the control points
    #   double **ctrl = vec2mat(c,mc,nc);
    #
    #   // Contruct the evaluated points
    p = zeros(mc,nu);                               #   double **pnt = vec2mat(p,mc,nu);
    #   // space for the basis functions
    N = zeros(d+1,1);                               #   double *N = (double*) mxMalloc((d+1)*sizeof(double));
    #
    #   // for each parametric point i
    for col=1:nu                                    #   for (col = 0; col < nu; col++) {
        s = findspan(nc-1, d, u[col], k);           #     s = findspan(nc-1, d, u[col], k);
        N = basisfun(s,u[col],d,k);                 #     basisfun(s, u[col], d, k, N);
                                                    #
        tmp1 = s - d + 1;                           #     tmp1 = s - d;
        for row=1:mc                                #     for (row = 0; row < mc; row++)  {
            tmp2 = 0;                               #       tmp2 = 0.0;
            for i=0:d                               #       for (i = 0; i <= d; i++)
                tmp2 = tmp2 + N[i+1]*c[row,tmp1+i];  # 	tmp2 += N[i] * ctrl[tmp1+i][row];
            end                                     #
            p[row,col] = tmp2;                      #       pnt[col][row] = tmp2;
        end                                         #     }
    end                                             #   }
    #
    #   mxFree(N);
    #   freevec2mat(pnt);
    #   freevec2mat(ctrl);
    #
    #   return ierr;
    p                                                #   }
end

###############################################################################

function  findspan(n,p,u,U)
    # FINDSPAN  Find the span of a B-Spline knot vector at a parametric point
    # -------------------------------------------------------------------------
    # ADAPTATION of FINDSPAN from C
    # -------------------------------------------------------------------------
    #
    # Calling Sequence:
    #
    #   s = findspan(n,p,u,U)
    #
    #  INPUT:
    #
    #    n - number of control points - 1
    #    p - spline degree
    #    u - parametric point
    #    U - knot sequence
    #
    #  RETURN:
    #
    #    s - knot span
    #
    #  Algorithm A2.1 from 'The NURBS BOOK' pg68

    # int findspan(int n, int p, double u, double *U) {

    #   int low, high, mid;
    #   // special case
    if u==U[n+2]
        s=n
        return s
    end               #   if (u == U[n+1]) return(n);
    #
    #   // do binary search
    low = p;                                        #   low = p;
    high = n + 1;                                   #   high = n + 1;
    mid = floor(Int,(low + high) / 2);                  #   mid = (low + high) / 2;
    while (u < U[mid+1] || u >= U[mid+2])           #   while (u < U[mid] || u >= U[mid+1])  {
        if (u < U[mid+1])                           #     if (u < U[mid])
            high = mid;                             #       high = mid;
        else                                        #     else
            low = mid;                              #       low = mid;
        end
        mid = floor(Int,(low + high) / 2);              #     mid = (low + high) / 2;
    end                                             #   }
    #
    s = mid;                                        #   return(mid);
    s                                           #   }
end

###############################################################################

function basisfun(i,u,p,U)
    # BASISFUN  Basis function for B-Spline
    # -------------------------------------------------------------------------
    # ADAPTATION of BASISFUN from C Routine
    # -------------------------------------------------------------------------
    #
    # Calling Sequence:
    #
    #   N = basisfun(i,u,p,U)
    #
    #    INPUT:
    #
    #      i - knot span  ( from FindSpan() )
    #      u - parametric point
    #      p - spline degree
    #      U - knot sequence
    #
    #    OUTPUT:
    #
    #      N - Basis functions vector[p+1]
    #
    #    Algorithm A2.2 from 'The NURBS BOOK' pg70.

    #   void basisfun(int i, double u, int p, double *U, double *N) {
    #   int j,r;
    #   double saved, temp;
    i = i + 1;
    #   // work space
    left = zeros(p+1,1);                              #   double *left  = (double*) mxMalloc((p+1)*sizeof(double));
    right = zeros(p+1,1);                             #   double *right = (double*) mxMalloc((p+1)*sizeof(double));
    N = zeros(p+1,1);
    N[1] = 1;                                         #   N[0] = 1.0;
    for j=1:p                                         #   for (j = 1; j <= p; j++) {
        left[j+1] = u - U[i+1-j];                     #   left[j]  = u - U[i+1-j];
        right[j+1] = U[i+j] - u;                      #   right[j] = U[i+j] - u;
        saved = 0;                                    #   saved = 0.0;

        for r=0:j-1                                   #   for (r = 0; r < j; r++) {
            temp = N[r+1]/(right[r+2] + left[j-r+1]); #   temp = N[r] / (right[r+1] + left[j-r]);
            N[r+1] = saved + right[r+2]*temp;         #   N[r] = saved + right[r+1] * temp;
            saved = left[j-r+1]*temp;                 #   saved = left[j-r] * temp;
        end                                           #   }

        N[j+1] = saved;                               #   N[j] = saved;
    end                                               #   }

    #   mxFree(left);
    #   mxFree(right);
    #   }
    N
end

###############################################################################

function  bspline_deriv(order, knots, ctrl) # prever o caso de ctrl de duas colunas <-- EmerEdited
    # Knots and control points associated with the derivative of B-spline curve.
    #
    # Input arguments:
    # order:
    #    B-spline order (2 for linear, 3 for quadratic, etc.)
    # knots:
    #    knot vector
    # ctrl:
    #    control points, typically 2-by-m, 3-by-m, or 4-by-m (for weights)
    #
    # Output arguments:
    # dctrl:
    #    control points of the derivative of the input B-spline curve
    # dknots:
    #    the new knot vector associated with the derivative B-spline curve
    p = order - 1;
    tmp = size(ctrl);
    n = tmp[2]-1;
    dim = tmp[1];

    # derivative knots
    dknots = knots[2:end-1];

    # derivative control points
    dctrl = zeros(dim,n);
    for i = 1 : n
        dctrl[:,i] = (p / (knots[i+p+1] - knots[i+1])) * (ctrl[:,i+1] - ctrl[:,i]);
    end
    if n==1
        dctrl[end]=1 #[dctrl;ones(1,n)]
    end
    dctrl,dknots
end

###############################################################################

function  nrbderiv(nurbs)
    #
    # Function Name:
    #
    #   nrbderiv - Construct the first derivative representation of a
    #              NURBS curve or surface.

    # NURBS structure represents a curve
    dcoefs,dknots = bspline_deriv(nurbs.order,nurbs.knots,nurbs.coefs);
    dnurbs = nrbmak(dcoefs, dknots);
end

###############################################################################

function nrbdeval(nurbs, dnurbs, tt)
    val = bspeval(nurbs.order-1,nurbs.coefs,nurbs.knots,tt);
    cw = val[4,:];
    cp = val[1:3,:];

    # NURBS is a curve
    temp =repeat(cw, outer = [3, 1]);
    pnt = cp./temp;

    # first derivative
    val = bspeval(dnurbs.order-1,dnurbs.coefs,dnurbs.knots,tt);
    cuw = val[4,:];
    cup = val[1:3,:];
    temp1 = repeat(cuw, outer = [3, 1]);
    jac = (cup-temp1.*pnt)./temp;
    pnt,jac
end

###############################################################################

function nrbbasisfun(nrb,points)
    n    = size(nrb.coefs, 2) -1;
    p    = nrb.order -1;
    U    = nrb.knots;
    w    = nrb.coefs[4,:];
    npoints=length(points)
#     println("npoints = $npoints")
    B=zeros(nrb.order,npoints)
    nbfu=zeros(Int,nrb.order,npoints) 
    for i=1:npoints
        spu  =  findspan(n, p, points[i], U); 
#         println("spu =$spu")
#         println("p =$p")
        nbfu[:,i] =  (spu-p+1) .+ collect(0:p)
#         println("nbfu =$nbfu")
        N     = w[nbfu[:,i]] .* basisfun(spu, points[i], p, U);
        B[:,i]     = N./sum(N);    # funções bases NURBS
#         @show(B)
    end
    B, nbfu
end

###############################################################################

function nrbplot(crv,divs=100)
#     p=nrbeval(crv,linspace(0,1,divs))
    p=nrbeval(crv,range(0,length=divs,stop=1))
    f=plot(p[1,:],p[2,:],c=:black)
    f=scatter(crv.coefs[1,:]./crv.coefs[4,:],crv.coefs[2,:]./crv.coefs[4,:],c=:blue)   # < ----- emeredited
#     f=scatter(crv.coefs[1,:]./crv.coefs[3,:],crv.coefs[2,:]./crv.coefs[3,:],c=:black)
    
end

function nrbplot2(crv,divs=100)
    #     p=nrbeval(crv,linspace(0,1,divs))
        p=nrbeval(crv,range(0,length=divs,stop=1))
        f=plot(p[1,:],p[2,:],c=:black)
        #f=scatter(crv.coefs[1,:]./crv.coefs[4,:],crv.coefs[2,:]./crv.coefs[4,:],c=:black)   # < ----- emeredited
    #     f=scatter(crv.coefs[1,:]./crv.coefs[3,:],crv.coefs[2,:]./crv.coefs[3,:],c=:black)
        
    end

###############################################################################
   
###############################################################################

function bincoeff(n,k)
#  Computes the binomial coefficient.
#
#      ( n )      n!
#      (   ) = --------
#      ( k )   k!(n-k)!
#
#  b = bincoeff(n,k)
#
#  Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.

                                                          # double bincoeff(int n, int k)
                                                          # {
b = floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));      #   return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
                                                          # }
return b
end

###############################################################################

function factln(n)
# computes ln(n!)
if n <= 1
    f = 0;
else
    f = log(factorial(n)) #log(factorial(n));
end
return f
end

###############################################################################

###############################################################################

function vectrans(v)
    # v = [vector[:],0,0];
    [1 0 0 v[1]; 0 1 0 v[2]; 0 0 1 v[3]; 0 0 0 1];

end

###############################################################################

function vecmag(v)
    # v = [vector[:],0,0];
    sqrt.(sum(v.^2,1));
end

###############################################################################

function vecangle(num,dem)
    ang = atan2.(num,dem)
    if length(num)>1
        ang[ang .< 0.0] += 2*pi
    elseif ang[1] <0.0
        ang+=2*pi
    end
    ang
end

###############################################################################

function vecrotx(angle)
    sn = sin(angle);
    cn = cos(angle);
    [1 0 0 0; 0 cn -sn 0; 0 sn cn 0; 0 0 0 1]
end

###############################################################################

function vecroty(angle)
    sn = sin(angle);
    cn = cos(angle);
    [cn 0 sn 0; 0 1 0 0; -sn 0 cn 0; 0 0 0 1]
end

###############################################################################

function vecrotz(angle)
    sn = sin(angle);
    cn = cos(angle);
    [cn -sn 0 0; sn cn 0 0; 0 0 1 0; 0 0 0 1]
end

###############################################################################

function nrbtform!(nurbs::Curve,tmat)
    dim,nu = size(nurbs.coefs);
    nurbs.coefs = reshape(tmat*reshape(nurbs.coefs,dim,nu),dim, nu);
end

###############################################################################

function nrbtform!(nurbs::Surf,tmat)
    dim,nu,nv = size(nurbs.coefs);
    nurbs.coefs = reshape(tmat*reshape(nurbs.coefs,dim,nu*nv),dim, nu,nv);
end

###############################################################################

function vecscale(vector)
    s = [vector[:];0;0];
    [s[1] 0 0 0; 0 s[2] 0 0; 0 0 s[3] 0; 0 0 0 1]
end

###############################################################################

function nrbcirc(radius,center,sang,eang)  # levar em conta a rotação no contorno interno EmerEdited
#     @show eang, sang                         # <------- EmerEdited
    sweep = eang - sang;       # sweep angle of arc
#     println("sweep = $sweep")                # <------- EmerEdited
    if abs(sweep) <= pi/2
        narcs = 1;                # number of arc segments
        knots = [0, 0, 0, 1, 1, 1];
    elseif abs(sweep) <= pi
        narcs = 2;
        knots = [0, 0, 0, 0.5, 0.5, 1, 1, 1];
    elseif abs(sweep) <= 3*pi/2
        narcs = 3;
        knots = [0, 0, 0, 1/3, 1/3, 2/3, 2/3, 1, 1, 1];
    else
        narcs = 4;
        knots = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1];
    end

    dsweep = sweep/(2*narcs);     # arc segment sweep angle/2
    # determine middle control point and weight
    wm = cos(dsweep);  # would not be sin(dsweep)  ??  EmerEdited
    x  = radius*wm;
    y  = radius*sin(dsweep);
    xm = x+y*tan(dsweep);

    # arc segment control points
    ctrlpt = [ x wm*xm   x;    # w*x - coordinate
              -y   0     y;    # w*y - coordinate
               0   0     0;    # w*z - coordinate
               1   wm    1];   # w   - coordinate

    # build up complete arc from rotated segments
    coefs = zeros(4,2*narcs+1);   # nurbs control points of arc
    xx = vecrotz(sang + dsweep);
    coefs[:,1:3] = xx*ctrlpt;     # rotate to start angle
    if narcs >= 2
        xx = vecrotz(2*dsweep);
#     println(xx)
        for n = 2:narcs
            m = 2*n .+ [0,1];                  # Emer EDITED
    #        println(coefs[:,m .- 2])
            coefs[:,m] = xx*coefs[:,m .- 2]
        end
    end

    # vectrans arc if necessary
    xx = vectrans(center);
    coefs = xx*coefs;
    nrbmak(coefs,knots)
end

###############################################################################

function nrbrevolve(curve,pnt,vec,theta = 2.0*pi)
    # Translate and rotate the curve into alignment with the z-axis
    T  = vectrans(-pnt);
    angx = vecangle(vec[1],vec[3]);
    RY = vecroty(-angx);
    vectmp = RY*[vec/norm(vec);1.0];
    angy = vecangle(vectmp[2],vectmp[3]);
    RX = vecrotx(angy);
    nrbtform!(curve,RX*RY*T);
    arc = nrbcirc(1.0,[0,0,0],0.0,theta);
    coefs = zeros(4,arc.number,curve.number);
    angle = vecangle(curve.coefs[2,:],curve.coefs[1,:])
    radius = vecmag(curve.coefs[1:2,:])
    for i = 1:curve.number
        coefs[:,:,i] = vecrotz(angle[i])*vectrans([0.0 0.0 curve.coefs[3,i]])*vecscale([radius[i] radius[i]])*arc.coefs
        coefs[4,:,i] = coefs[4,:,i]*curve.coefs[4,i]
    end
    surf = nrbmak(coefs,arc.knots, curve.knots);
    T = vectrans(pnt);
    RX = vecrotx(-angy);
    RY = vecroty(angx);
    nrbtform!(surf,T*RY*RX)
    surf
end

###############################################################################

function nrbline(p1,p2)
#
# NRBLINE: Construct a straight line.
#
# Calling Sequence:
#
#   crv = nrbline()
#   crv = nrbline(p1,p2)
#
# Parameters:
#
# p1        : 2D or 3D cartesian coordinate of the start point.
#
# p2            : 2D or 3D cartesian coordinate of the end point.
#
# crv       : NURBS curve for a straight line.
#
# Description:
#
#   Constructs NURBS data structure for a straight line. If no rhs
#   coordinates are included the function returns a unit straight
#   line along the x-axis.

#  D.M. Spink
#  Copyright (c) 2000.

coefs = [zeros(3,2); ones(1,2)];
coefs[1:length(p1),1] = p1[:];
coefs[1:length(p2),2] = p2[:];
knots=[0,0,1,1]
curve = nrbmak(coefs, knots);
return curve
end
















###############################################################################
###############################################################################
###########################   Backup   ####################################################
###############################################################################
###############################################################################
###############################################################################

#   #  < - ---- - - -- - - -- - emer   removed
function  bspkntins(d,c,k,u)

# BSPKNTINS:  Insert knots into a B-Spline
#
# Calling Sequence:
#
#   [ic,ik] = bspkntins(d,c,k,u)
#
#  INPUT:
#
#    d - spline degree             integer
#    c - control points            double  matrix(mc,nc)
#    k - knot sequence             double  vector(nk)
#    u - new knots                 double  vector(nu)
#
#  OUTPUT:
#
#    ic - new control points double  matrix(mc,nc+nu)
#    ik - new knot sequence  double  vector(nk+nu)
#
#  Modified version of Algorithm A5.4 from 'The NURBS BOOK' pg164.
#

mc,nc = size(c);
nu = length(u);
nk = length(k);
                                                     #
                                                     # int bspkntins(int d, double *c, int mc, int nc, double *k, int nk,
                                                     #               double *u, int nu, double *ic, double *ik)
                                                     # {
                                                     #   int ierr = 0;
                                                     #   int a, b, r, l, i, j, m, n, s, q, ind;
                                                     #   double alfa;
                                                     #
                                                     #   double **ctrl  = vec2mat(c, mc, nc);
ic = zeros(mc,nc+nu);                                #   double **ictrl = vec2mat(ic, mc, nc+nu);
ik = zeros(nk+nu);
                                                     #
n = size(c,2) - 1;                                   #   n = nc - 1;
r = length(u) - 1;                                   #   r = nu - 1;
                                                     #
m = n + d + 1;                                       #   m = n + d + 1;
a = findspan(n, d, u[1], k);                         #   a = findspan(n, d, u[0], k);
b = findspan(n, d, u[r+1], k);                       #   b = findspan(n, d, u[r], k);
b = b+1;                                             #   ++b;
                                                     #
for q=0:mc-1                                         #   for (q = 0; q < mc; q++)  {
   for j=0:a-d
       ic[q+1,j+1] = c[q+1,j+1]
   end        #     for (j = 0; j <= a-d; j++) ictrl[j][q] = ctrl[j][q];
   for j=b-1:n
       ic[q+1,j+r+2] = c[q+1,j+1]
   end      #     for (j = b-1; j <= n; j++) ictrl[j+r+1][q] = ctrl[j][q];
end                                                  #   }

for j=0:a
   ik[j+1] = k[j+1]
end                     #   for (j = 0; j <= a; j++)   ik[j] = k[j];
for j=b+d:m
   ik[j+r+2] = k[j+1]
end                 #   for (j = b+d; j <= m; j++) ik[j+r+1] = k[j];
                                                     #
i = b + d - 1;                                       #   i = b + d - 1;
s = b + d + r;                                       #   s = b + d + r;

for j=r:-1:0                                         #   for (j = r; j >= 0; j--) {
   while u[j+1] <= k[i+1] && i > a                   #     while (u[j] <= k[i] && i > a) {
       for q=0:mc-1                                  #       for (q = 0; q < mc; q++)
           ic[q+1,s-d] = c[q+1,i-d]                 #         ictrl[s-d-1][q] = ctrl[i-d-1][q];
       end
       ik[s+1] = k[i+1]                             #       ik[s] = k[i];
       s = s - 1;                                    #       --s;
       i = i - 1;                                    #       --i;
   end                                               #     }

   for q=0:mc-1                                      #     for (q = 0; q < mc; q++)
       ic[q+1,s-d] = ic[q+1,s-d+1]                  #       ictrl[s-d-1][q] = ictrl[s-d][q];
   end

   for l=1:d                                         #     for (l = 1; l <= d; l++)  {
       ind = s - d + l;                              #       ind = s - d + l;
       alfa = ik[s+l+1] - u[j+1]                    #       alfa = ik[s+l] - u[j];
       if abs(alfa) == 0                             #       if (fabs(alfa) == 0.0)
           for q=0:mc-1                              #         for (q = 0; q < mc; q++)
               ic[q+1,ind] = ic[q+1,ind+1]          #           ictrl[ind-1][q] = ictrl[ind][q];
           end
       else                                          #       else  {
           alfa = alfa/(ik[s+l+1] - k[i-d+l+1]);     #         alfa /= (ik[s+l] - k[i-d+l]);
           for q=0:mc-1                              #         for (q = 0; q < mc; q++)
               tmp = (1-alfa)*ic[q+1,ind+1];
               ic[q+1,ind] = alfa*ic[q+1,ind] + tmp; #           ictrl[ind-1][q] = alfa*ictrl[ind-1][q]+(1.0-alfa)*ictrl[ind][q];
           end
       end                                           #       }
   end                                               #     }
   #
   ik[s+1] = u[j+1]                                 #     ik[s] = u[j];
   s = s - 1;                                        #     --s;
end
return ic,ik
end
  #   emer removed
###############################################################################
###############################################################################
###############################################################################
   #  < - ---- - - -- - - -- - emer   removed
function bspdegelev(d,c,k,t)

# BSPDEGELEV:  Degree elevate a univariate B-Spline.
#
# Calling Sequence:
#
#   [ic,ik] = bspdegelev(d,c,k,t)
#
# Parameters:
#
#   d - Degree of the B-Spline.
#   c - Control points, matrix of size (dim,nc).
#   k - Knot sequence, row vector of size nk.
#   t - Raise the B-Spline degree t times.
#
#   ic - Control points of the new B-Spline.
#   ik - Knot vector of the new B-Spline.
#
mc,nc = size(c);
# nc = number of control points
# nknot = number of knots
ic = zeros(mc,nc);                                  #   double **ictrl = vec2mat(ic, mc, nc*(t+1));
n = nc - 1;                                               #   n = nc - 1;
                                                          #
bezalfs =  zeros(d+1,d+t+1);                              #   bezalfs = matrix(d+1,d+t+1);
bpts = zeros(mc,d+1);                                     #   bpts = matrix(mc,d+1);
ebpts = zeros(mc,d+t+1);                                  #   ebpts = matrix(mc,d+t+1);
Nextbpts = zeros(mc,d+1);                                 #   Nextbpts = matrix(mc,d+1);
alfs = zeros(d,1);                                        #   alfs = (double *) mxMalloc(d*sizeof(double));
                                                          #
m = n + d + 1;                                            #   m = n + d + 1;
ph = d + t;                                               #   ph = d + t;
ph2 = Int32(floor(ph / 2));                                      #   ph2 = ph / 2;
ik=zeros(ph+1)  # Duvida !!!!!!!!!!!!

                                                          #
                                                          #   // compute bezier degree elevation coefficeients
bezalfs[1,1] = 1;                                         #   bezalfs[0][0] = bezalfs[ph][d] = 1.0;
bezalfs[d+1,ph+1] = 1;                                    #

for i=1:1:ph2
   inv = 1/bincoeff(ph,i)
                                                            #     inv = 1.0 / bincoeff(ph,i);
   mpi = minimum([d,i]);                                        #     mpi = min(d,i);
                                                          #
   for j=maximum([0,i-t]):mpi                                   #     for (j = max(0,i-t); j <= mpi; j++)
       bezalfs[j+1,i+1] = inv*bincoeff(d,j)*bincoeff(t,i-j)  #       bezalfs[i][j] = inv * bincoeff(d,j) * bincoeff(t,i-j);
   end
end                                                       #   }
                                                          #
for i=ph2+1:ph-1                                          #   for (i = ph2+1; i <= ph-1; i++) {
   mpi = minimum([d,i]);                                        #     mpi = min(d, i);
   for j=maximum([0,i-t]):mpi                                   #     for (j = max(0,i-t); j <= mpi; j++)
       bezalfs[j+1,i+1] = bezalfs[d-j+1,ph-i+1];          #       bezalfs[i][j] = bezalfs[ph-i][d-j];
   end
end                                                       #   }
                                                          #
mh = ph;                                                  #   mh = ph;
kind = ph+1;                                              #   kind = ph+1;
r = -1;                                                   #   r = -1;
a = d;                                                    #   a = d;
b = d+1;                                                  #   b = d+1;
cind = 1;                                                 #   cind = 1;
ua = k[1];                                                #   ua = k[0];
                                                          #
for ii=0:mc-1                                             #   for (ii = 0; ii < mc; ii++)
   ic[ii+1,1] = c[ii+1,1]                                #     ictrl[0][ii] = ctrl[0][ii];
end                                                       #
for i=0:ph                                                #   for (i = 0; i <= ph; i++)
   ik[i+1] = ua;                                          #     ik[i] = ua;
end                                                       #
                                                          #   // initialise first bezier seg
for i=0:d                                                 #   for (i = 0; i <= d; i++)
   for ii=0:mc-1                                          #     for (ii = 0; ii < mc; ii++)
      bpts[ii+1,i+1] = c[ii+1,i+1];                       #       bpts[i][ii] = ctrl[i][ii];
   end
end                                                       #
                                                          #   // big loop thru knot vector
while b < m                                               #   while (b < m)  {
   i = b;                                                 #     i = b;
   while b < m && k[b+1] == k[b+2]                        #     while (b < m && k[b] == k[b+1])
      b = b + 1;                                          #       b++;
   end                                                    #
   mul = b - i + 1;                                       #     mul = b - i + 1;
   mh = mh + mul + t;                                     #     mh += mul + t;
   ub = k[b+1];                                           #     ub = k[b];
   oldr = r;                                              #     oldr = r;
   r = d - mul;                                           #     r = d - mul;
                                                          #
                                                          #     // insert knot u(b) r times
   if oldr > 0                                            #     if (oldr > 0)
      lbz = floor((oldr+2)/2);                            #       lbz = (oldr+2) / 2;
   else                                                   #     else
      lbz = 1;                                            #       lbz = 1;
   end                                                    #

   if r > 0                                               #     if (r > 0)
      rbz = ph - floor((r+1)/2);                          #       rbz = ph - (r+1)/2;
   else                                                   #     else
      rbz = ph;                                           #       rbz = ph;
   end                                                    #

   if r > 0                                               #     if (r > 0) {
                                                          #       // insert knot to get bezier segment
      numer = ub - ua;                                    #       numer = ub - ua;
      for q=d:-1:mul+1                                    #       for (q = d; q > mul; q--)
         alfs[q-mul] = numer / (k[a+q+1]-ua);             #         alfs[q-mul-1] = numer / (k[a+q]-ua);
      end

      for j=1:r                                           #       for (j = 1; j <= r; j++)  {
         save = r - j;                                    #         save = r - j;
         s = mul + j;                                     #         s = mul + j;
                                                          #
         for q=d:-1:s                                     #         for (q = d; q >= s; q--)
            for ii=0:mc-1                                 #           for (ii = 0; ii < mc; ii++)
               tmp1 = alfs[q-s+1]*bpts[ii+1,q+1]
               tmp2 = (1-alfs[q-s+1])*bpts[ii+1,q]
               bpts[ii+1,q+1] = tmp1 + tmp2;              #             bpts[q][ii] = alfs[q-s]*bpts[q][ii]+(1.0-alfs[q-s])*bpts[q-1][ii];
            end
         end                                              #

         for ii=0:mc-1                                    #         for (ii = 0; ii < mc; ii++)
            Nextbpts[ii+1,save+1] = bpts[ii+1,d+1]       #           Nextbpts[save][ii] = bpts[d][ii];
         end
      end                                                 #       }
   end                                                    #     }
                                                          #     // end of insert knot
                                                          #
                                                          #     // degree elevate bezier
   for i=lbz:ph                                           #     for (i = lbz; i <= ph; i++)  {
      for ii=0:mc-1                                       #       for (ii = 0; ii < mc; ii++)
         ebpts[ii+1,i+1] = 0;                             #         ebpts[i][ii] = 0.0;
      end
      mpi = minimum([d, i]);                                    #       mpi = min(d, i);
      for j=maximum([0,i-t]):mpi                                #       for (j = max(0,i-t); j <= mpi; j++)
         for ii=0:mc-1                                    #         for (ii = 0; ii < mc; ii++)
            tmp1 = ebpts[ii+1,i+1]
            tmp2 = bezalfs[j+1,i+1]*bpts[ii+1,j+1]
            ebpts[ii+1,i+1] = tmp1 + tmp2;                #           ebpts[i][ii] = ebpts[i][ii] + bezalfs[i][j]*bpts[j][ii];
         end
      end
   end                                                    #     }
                                                          #     // end of degree elevating bezier
                                                          #
   if oldr > 1                                            #     if (oldr > 1)  {
                                                          #       // must remove knot u=k[a] oldr times
      first = kind - 2;                                                    #       first = kind - 2;
      last = kind;                                        #       last = kind;
      den = ub - ua;                                      #       den = ub - ua;
      bet = floor((ub-ik[kind]) / den);                   #       bet = (ub-ik[kind-1]) / den;
                                                          #
                                                          #       // knot removal loop
      for tr=1:oldr-1                                     #       for (tr = 1; tr < oldr; tr++)  {
         i = first;                                       #         i = first;
         j = last;                                        #         j = last;
         kj = j - kind + 1;                               #         kj = j - kind + 1;
         while j-i > tr                                   #         while (j - i > tr)  {
                                                          #           // loop and compute the new control points
                                                          #           // for one removal step
            if i < cind                                   #           if (i < cind)  {
               alf = (ub-ik[i+1])/(ua-ik[i+1]);           #             alf = (ub-ik[i])/(ua-ik[i]);
                  nic=size(ic,2)                          #             for (ii = 0; ii < mc; ii++)
                  tmp1 = alf*ic[:,i+1]
                  tmp2 = (1-alf)*ic[:,i]
                  if(nic>=i+1)
                       ic[:,i+1] = tmp1 + tmp2;
                  else
                       ic=[ic (tmp1+tmp2)]                     #               ictrl[i][ii] = alf * ictrl[i][ii] + (1.0-alf) * ictrl[i-1][ii];
                  end
            end                                           #           }
            if j >= lbz                                   #           if (j >= lbz)  {
               if j-tr <= kind-ph+oldr                    #             if (j-tr <= kind-ph+oldr) {
                  gam = (ub-ik[j-tr+1]) / den;            #               gam = (ub-ik[j-tr]) / den;
                  for ii=0:mc-1                           #               for (ii = 0; ii < mc; ii++)
                     tmp1 = gam*ebpts[ii+1,kj+1]
                     tmp2 = (1-gam)*ebpts[ii+1,kj+2]
                     ebpts[ii+1,kj+1] = tmp1 + tmp2;      #                 ebpts[kj][ii] = gam*ebpts[kj][ii] + (1.0-gam)*ebpts[kj+1][ii];
                  end                                     #             }
               else                                       #             else  {
                  for ii=0:mc-1                           #               for (ii = 0; ii < mc; ii++)
                     tmp1 = bet*ebpts[ii+1,kj+1];
                     tmp2 = (1-bet)*ebpts[ii+1,kj+2];
                     ebpts[ii+1,kj+1] = tmp1 + tmp2;      #                 ebpts[kj][ii] = bet*ebpts[kj][ii] + (1.0-bet)*ebpts[kj+1][ii];
                  end
               end                                        #             }
            end                                           #           }
            i = i + 1;                                    #           i++;
            j = j - 1;                                    #           j--;
            kj = kj - 1;                                  #           kj--;
         end                                              #         }
                                                          #
         first = first - 1;                               #         first--;
         last = last + 1;                                 #         last++;
      end                                                 #       }
   end                                                    #     }
                                                          #     // end of removing knot n=k[a]
                                                          #
                                                          #     // load the knot ua
   if a != d                                              #     if (a != d)
      for i=0:ph-oldr-1                                   #       for (i = 0; i < ph-oldr; i++)  {
         nik=length(ik)
         if(nik>=kind+1)
             ik[kind+1] = ua;                                 #         ik[kind] = ua;
         else
             append!(ik,ua)
         end
         kind = kind + 1;                                 #         kind++;
      end
   end                                                    #       }
                                                          #
                                                          #     // load ctrl pts int
      for j=lbz:rbz                                       #     for (j = lbz; j <= rbz; j++)  {
            nic=size(ic,2)
            if(nic>=cind+1)
              ic[:,cind+1] = ebpts[:,j+1]            #         ictrl[cind][ii] = ebpts[j][ii];
            else
              ic=[ic ebpts[:,j+1]]
            end
         cind = cind + 1;                                 #       cind++;
      end                                                 #     }
                                                          #
      if b < m                                            #     if (b < m)  {
                                                          #       // setup for next pass thru loop
         for j=0:r-1                                      #       for (j = 0; j < r; j++)
            for ii=0:mc-1                                 #         for (ii = 0; ii < mc; ii++)
               bpts[ii+1,j+1] = Nextbpts[ii+1,j+1]       #           bpts[j][ii] = Nextbpts[j][ii];
            end
         end
         for j=r:d                                        #       for (j = r; j <= d; j++)
            for ii=0:mc-1                                 #         for (ii = 0; ii < mc; ii++)
               bpts[ii+1,j+1] = c[ii+1,b-d+j+1]          #           bpts[j][ii] = ctrl[b-d+j][ii];
            end
         end
         a = b;                                           #       a = b;
         b = b+1;                                         #       b++;
         ua = ub;                                         #       ua = ub;
                                                          #     }
      else                                                #     else
                                                          #       // end knot
         for i=0:ph
            nik=length(ik)
            if(nik>=kind+i+1)
                ik[kind+i+1] = ub;                            #         ik[kind+i] = ub;
            else
                append!(ik,ub)
            end
         end
      end
end
return ic,ik
end

function basisfundecomp(Be,dBedxi,Ce,w)
    # Calcula R e dRdxi => Eq. (65) e (66) do artigo do Borden
    Wb=Ce'*w  # Bézier weights (nodal values)
    # @infiltrate
    wb = dot(Be,Wb);   #  Bézier weight function
    dwbdxi = dot(dBedxi,Wb);  # derivada de w com relaçao a xi
    # Shape function and derivatives
    R = diagm(w)*Ce*Be/wb; # Equação (65) do artigo do Borden R tem dim [(p+1) x 1]
    dRdxi = diagm(w)*Ce*(dBedxi[:,1]/wb-dwbdxi*Be/(wb*wb)); # Eq. (66)
    return R,dRdxi
  end

  
function evalcrv(crv, u)
    uu = unique(crv.knots)
    ind = sum(u .> uu)
    unorm = 2 * (u - uu[ind]) / (uu[ind + 1] - uu[ind]) - 1
    shapes=Bernsteins(crv.order - 1, unorm)
    derivs=dBernsteins(crv.order - 1, unorm)
  # control= reshape(sup.coefs[:,sup.conn[ind][1],sup.conn[ind][2]],4,length(sup.conn[ind][1])*length(sup.conn[ind][2]))
    control = crv.coefs[:,crv.conn[ind]]
  
  # @show u
  # @show ind
    we = control[4,:]
    ptse = control[1:3,:] ./ we'
    R, dR = basisfundecomp(shapes, derivs, crv.C[:,:,ind],we)
    x = ptse * R; # Calcula a coordenada x do ponto de integra��o
  #@show norm(x)
    dxdqsi = ptse * dR
    dgamadqsi = norm(dxdqsi);

  
    nx = dxdqsi[2] / dgamadqsi; # Componente x do vetor normal unit�rio
    ny = -dxdqsi[1] / dgamadqsi; # Componente y do vetor normal unit�rio
  # Rcheio=zeros(sup.number[1]*sup.number[2])
  # Rcheio[(ones(Int64,length(sup.conn[ind][2]))*(sup.conn[ind][1]-1)'*(sup.number[2]).+sup.conn[ind][2])[:]]=R
  # x,Rcheio,n   
  # x, R ./ we, [nx,ny]  
  x, R, [nx,ny]  
end