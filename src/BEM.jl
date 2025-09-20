module BEM
using Reexport
@reexport using DrWatson
@reexport using LinearAlgebra
@reexport using Statistics
@reexport using FastGaussQuadrature
@reexport using SparseArrays
@reexport using StaticArrays
@reexport using Krylov
@reexport using TimerOutputs#, FFTW,ForwardDiff, Optim,Printf
@reexport using CairoMakie
@reexport using DelaunayTriangulation#, WriteVTK,Distances
# export lines, lines!, scatter, scatter!, tricontourf
@reexport using Infiltrator
@reexport using ParallelKMeans
@reexport using LowRankApprox
@reexport using PolynomialRoots
@reexport using SpecialFunctions
@reexport using Krylov
@reexport using NonlinearSolve
using ProgressMeter
using Distributed
using Printf
include("estruturas.jl")

include("./Hmat/HMatrices.jl")

include("integra.jl")
include("calc_HeG_elastico.jl")
include("calc_HeG_elastico_direto.jl")
export Contato_sem_atrito_NL
include("calc_HeG_elastico_superposicao.jl")
include("calc_HeG_potencial.jl")
include("calc_HeG_potencial_direto.jl")
include("calc_HeG_potencial_tudoemu.jl")
include("calc_HeG_helmholtz.jl")
include("casca.jl")

export calc_HeG,
    format_dad,
    separa,
    aplicaCDC,
    calc_Ti,
    calc_Aeb,
    nc,
    ni,
    calc_tens_int,
    calc_tens_cont,
    calc_SeD,
    potencial_correlato,
    corrigediag!,
    calc_gap,
    verifica_contato_sem_atrito,
    muda_nt!,
    muda_nt,
    transforma_tensao
include("termoelasticidade.jl")
export termoelasticidade

include("CFD.jl")
export passo_não_linear, passo_não_linear_Hmat

include("format.jl")
include("gera_p_in.jl")
include("trinca.jl")
include("calc_fforma.jl")
include("transforma.jl")
include("mostra_problema.jl")
export mostra_geometria, mostra_resultado, mostra_deformação
include("rim.jl")
import Base: *
export MatrizH, matvec
include("arvore.jl")
export Monta_M_RIMd, Finterp
include("nurbs.jl")
include("format_iga.jl")
export format_dad_iga
include("aniso.jl")
export Compute_Material, Compute_T, Compute_Qbar, Assembly_Q, calc_HeGeIt
include("calc_placa.jl")
include("calc_placa_espessa.jl")
export Compute_Material_Placa,
    simula_placa_flambagem, simula_placa_tempo, Monta_M_RIMd, Monta_M_RIM, DRM
include("placa_tensao.jl")
include("radial.jl")
export define_SubRegioes, subregioes
include("SubRegioes.jl")




export nrmse, nme
nrmse(y_true::Array, y_pred) =
    sqrt(sum((y_true .- y_pred) .^ 2) / length(y_true)) /
    (maximum(y_true) - minimum(y_true))
nrmse(y_true::Number, y_pred) = sqrt(sum((y_true .- y_pred) .^ 2) / length(y_pred)) / y_true
nme(y_true::Number, y_pred) = sum(abs.(y_true .- y_pred)) / length(y_pred) / y_true
nme(y_true::Array, y_pred) =
    sum(abs.(y_true .- y_pred)) / length(y_pred) / (maximum(y_true) - minimum(y_true))
end
