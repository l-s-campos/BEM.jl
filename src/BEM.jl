module BEM
using Base: NamedTuple
using ForwardDiff: length
using LinearAlgebra,Statistics, FastGaussQuadrature, ForwardDiff,SparseArrays,Optim
using TimerOutputs#, FFTW
using GLMakie,Triangulate#, WriteVTK
using Infiltrator,Distances,ParallelKMeans,LowRankApprox
using PolynomialRoots, SpecialFunctions

include("estruturas.jl")
export potencial,helmholtz,potencial_iga,elastico,elastico_iga,elemento,elastico_aniso,elastico_aniso_iga

include("calc_HeG_elastico.jl")
include("calc_HeG_potencial.jl")
include("calc_HeG_helmholtz.jl")
export calc_HeG,format_dad,separa,aplicaCDC,calc_Ti,calc_Aeb,nc,ni

include("format.jl")
include("gera_p_in.jl")
include("calc_fforma.jl")
include("transforma.jl")
include("mostra_problema.jl")
export mostra_geometria,mostra_resultado
include("rim.jl")
export calc_HeG_interp,Hinterp,cluster,Ainterp,matvec,Akmeans,Akmeans2
include("arvore.jl")
export Monta_M_RIMd,Finterp
include("nurbs.jl")
include("format_iga.jl")
export format_dad_iga
include("aniso.jl")
export Compute_Material,Compute_T,Compute_Qbar,Assembly_Q
end

