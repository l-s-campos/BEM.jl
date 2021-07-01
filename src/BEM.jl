module BEM
using LinearAlgebra,Statistics, FastGaussQuadrature, ForwardDiff,SparseArrays
using TimerOutputs#, FFTW
using GLMakie,Triangulate#, WriteVTK
using Infiltrator,Distances,ParallelKMeans,LowRankApprox

export potencial,potencial_iso,elastico,elemento
export mostra_geometria,mostra_resultado
export calc_HeG,format_dad,separa,aplicaCDC,calc_Ti,calc_Aeb,nc,ni
export calc_HeG_interp,Hinterp,cluster,Ainterp,matvec,Akmeans,Akmeans2
export Monta_M_RIMd,Finterp


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
struct elastico
    NOS ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    ELEM ::Array{elementov,1}
    Ev ::Array{Float64,1}
end
struct potencial
    NOS ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    ELEM ::Array{elemento,1}
    k ::Float64    
end
struct helmholtz
    NOS ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    ELEM ::Array{elemento,1}
    k ::Array{Float64,1}    
end

struct potencial_iso
    NOS ::Array{Float64,2}
    pontos_controle ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    tipoCDC  ::Array{Int64,1}
    valorCDC  ::Array{Float64,2}
    ELEM ::Array{bezier,1}
    k ::Float64  
    E ::SparseMatrixCSC{Float64,Int64} 
end

struct elastico_iso
    NOS ::Array{Float64,2}
    pontos_controle ::Array{Float64,2}
    pontos_internos ::Array{Float64,2}
    tipoCDC  ::Array{Int64,1}
    valorCDC  ::Array{Float64,2}
    ELEM ::Array{bezier,1}
    Ev ::Array{Float64,1}
    E ::SparseMatrixCSC{Float64,Int64} 
end

abstract type DadosBEM end
abstract type escalar     <: DadosBEM end
abstract type vetorial     <: DadosBEM end

 potencial_iso     <: escalar 
 potencial     <: escalar 
 helmholtz     <: escalar 

elastico     <: vetorial 
elastico_iso     <: vetorial 

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
export format_dad_iso
include("formatiso.jl")

import Base.size
size(e::elemento)=size(e.indices,1)
size(e::bezier)=size(e.indices,1)
size(e::elementov)=size(e.indices,1)
end

