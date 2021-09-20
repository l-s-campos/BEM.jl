module BEM
using Base: NamedTuple
using ForwardDiff: length
using LinearAlgebra,Statistics, FastGaussQuadrature, ForwardDiff,SparseArrays
using TimerOutputs#, FFTW
using GLMakie,Triangulate#, WriteVTK
using Infiltrator,Distances,ParallelKMeans,LowRankApprox
using PolynomialRoots, SpecialFunctions

export potencial,helmholtz,potencial_iga,elastico,elastico_iga,elemento,elastico_aniso,elastico_aniso_iga
export mostra_geometria,mostra_resultado
export calc_HeG,format_dad,separa,aplicaCDC,calc_Ti,calc_Aeb,nc,ni
export calc_HeG_interp,Hinterp,cluster,Ainterp,matvec,Akmeans,Akmeans2
export Monta_M_RIMd,Finterp
export Compute_Material,Compute_T,Compute_Qbar,Assembly_Q, calc_Hsing


abstract type DadosBEM end
abstract type escalar     <: DadosBEM end
abstract type vetorial     <: DadosBEM end

struct elemento
    indices ::Array{Int64,1}
    tipoCDC ::Int64
    valorCDC ::Array{Float64,1}
    ξs ::Array{Float64,1}
end
struct elementov
    indices ::Array{Int64,1}
    tipoCDC ::Array{Int64,1}
    valorCDC ::Array{Float64,2}
    ξs ::Array{Float64,1}
end
struct bezier
    indices ::Array{Int64,1}
    C ::Array{Float64,2}
    p ::Int64
    limites ::Array{Float64,2}
    Wb :: Array{Float64,1}
    sing ::Array{Int64,1}
end
struct elastico <: vetorial 
    NOS ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    ELEM ::Array{elementov,1}
    Ev ::Array{Float64,1}
end

struct elastico_aniso <: vetorial
    NOS ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    ELEM ::Array{elementov,1}
    k ::NamedTuple
end
struct potencial  <: escalar 
    NOS ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    ELEM ::Array{elemento,1}
    k ::Float64    
end
struct helmholtz  <: escalar 
    NOS ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    ELEM ::Array{elemento,1}
    k ::NamedTuple    
end

struct potencial_iga  <: escalar 
    NOS ::Array{Float64,2}
    pontos_controle ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    tipoCDC  ::Array{Int64,2}
    valorCDC  ::Array{Float64,2}
    ELEM ::Array{bezier,1}
    k ::Float64  
    E ::SparseMatrixCSC{Float64,Int64} 
end

struct elastico_iga <: vetorial 
    NOS ::Array{Float64,2}
    pontos_controle ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    tipoCDC  ::Array{Int64,2}
    valorCDC  ::Array{Float64,2}
    ELEM ::Array{bezier,1}
    Ev ::Array{Float64,1}
    E ::SparseMatrixCSC{Float64,Int64} 
end

struct elastico_aniso_iga <: vetorial 
    NOS ::Array{Float64,2}
    pontos_controle ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    tipoCDC  ::Array{Int64,2}
    valorCDC  ::Array{Float64,2}
    ELEM ::Array{bezier,1}
    k ::NamedTuple
    E ::SparseMatrixCSC{Float64,Int64} 
end


nc(dad::DadosBEM) = size(dad.NOS,1)
ni(dad::DadosBEM) = size(dad.pontos_internos,1)
include("format.jl")
include("gera_p_in.jl")
include("calc_HeG_elastico.jl")
include("calc_HeG_potencial.jl")
include("calc_HeG_helmholtz.jl")
include("calc_fforma.jl")
include("transforma.jl")
include("mostra_problema.jl")
include("rim.jl")
include("arvore.jl")
include("nurbs.jl")
export format_dad_iga
include("format_iga.jl")
include("aniso.jl")

import Base.size
size(e::elemento)=size(e.indices,1)
size(e::bezier)=size(e.indices,1)
size(e::elementov)=size(e.indices,1)
end

