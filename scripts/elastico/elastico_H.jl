## Início da análise
# using Pkg
# Pkg.activate(pwd())
# Pkg.instantiate()
using DataFrames
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
# nelem = 10  #Numero de elementos
# order = 3
NPX = 10 #pontos internos na direção x
NPY = 10 #pontos internos na direção y
npg = 16    #apenas números pares
# ## Formatação dos dados ________________________________________________

prob = cilindro
ana = ana_cilindro

p = 2
i = 7
j = 2
nelem = 2 * 2^i  #Numero de elementos
order = j

dad = format_dad(prob(nelem, order), NPX, NPY) # dados
# include(datadir("dadpotencial.jl"))

corrige_CDC_cilindro(dad)

# dad.pontos_internos = testeinterno
# function compara_potencial(dad, npg)
reset_timer!()
tHeG = @timed H, G = calc_HeG(dad, npg, interno = true)  #importante
A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
tsolve = @timed x = A \ b
u, t = separa(dad, x) #importante
M = BEM.Monta_M_RIMd(dad, npg)


elements = [
    [Point2D(dad.NOS[i, 1], dad.NOS[i, 2]) for i = 1:nc(dad)]
    [Point2D(dad.pontos_internos[i, 1], dad.pontos_internos[i, 2]) for i = 1:ni(dad)]
];
# elements = [elements elements]'[:]
struct elastMatrix <: AbstractMatrix{SMatrix{2,2,Float64,4}}
    A::Matrix{Float64}
end
function Base.getindex(K::elastMatrix, i::Int, j::Int)
    # return K.A[i, j];
     return SMatrix{2,2,Float64,4}(K.A[2i-1:2i, 2j-1:2j])
end
Base.size(K::elastMatrix) = round.(Int,size(K.A)./2)

# Xclt = Yclt = ClusterTree(repeat(elements, inner = 2), BEM.PrincipalComponentSplitter())
Xclt = Yclt = ClusterTree(elements, BEM.PrincipalComponentSplitter())
adm = StrongAdmissibilityStd()
comp = PartialACA(atol = 1e-6, rank = 20, rtol = 1e-6)
# comp = TSVD(atol = 1e-6, rank = 20, rtol = 1e-6)
KA = elastMatrix(A)
KM = elastMatrix(M)
@time HA = assemble_hmatrix(KA, Xclt, Yclt; adm, comp, threads = false, distributed = false)
@time HM = assemble_hmatrix(KM, Xclt, Yclt; adm, comp, threads = false, distributed = false)
Re = 100;
@time FH = lu(HA; atol = 1e-6);
@time F = lu(A);
@time FH \ x;
@time F \ x;

@time HM * x;
y = similar(x)
@time x1 = F \ (b - Re * M * b);
@time begin
    BEM.mul!(y, HM, b, Re, 0)
    x2 = FH \ (b - y)
end;
@time begin
    x3 = FH \ (b - Re * HM * b)
end;



