module BEM
using Reexport, DrWatson
using Base: NamedTuple
using ForwardDiff: length
using LinearAlgebra, Statistics, FastGaussQuadrature, ForwardDiff, SparseArrays, Optim
using TimerOutputs#, FFTW
using GLMakie, DelaunayTriangulation#, WriteVTK
export lines, lines!, scatter, scatter!,tricontourf
using Infiltrator, Distances, ParallelKMeans, LowRankApprox
using PolynomialRoots, SpecialFunctions
using LinearSolve, IterativeSolvers
using ProgressMeter

include("estruturas.jl")
include("integra.jl")

include("calc_HeG_elastico.jl")
include("calc_HeG_elastico_superposicao.jl")
include("calc_HeG_potencial.jl")
include("calc_HeG_helmholtz.jl")
include("casca.jl")
export calc_HeG, format_dad, separa, aplicaCDC, calc_Ti, calc_Aeb, nc, ni, calc_tens_int, calc_tens_cont, calc_SeD
include("termoelasticidade.jl")
export termoelasticidade

include("format.jl")
include("gera_p_in.jl")
include("trinca.jl")
include("calc_fforma.jl")
include("transforma.jl")
include("mostra_problema.jl")
export mostra_geometria, mostra_resultado,mostra_deformação
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
export Compute_Material_Placa, simula_placa_flambagem, simula_placa_tempo, Monta_M_RIMd, Monta_M_RIM, DRM
include("placa_tensao.jl")
include("radial.jl")
export define_SubRegioes, subregioes
include("SubRegioes.jl")
export nrmse, nme
nrmse(y_true::Array, y_pred) = sqrt(sum((y_true .- y_pred) .^ 2) / length(y_true)) / (maximum(y_true) - minimum(y_true))
nrmse(y_true::Number, y_pred) = sqrt(sum((y_true .- y_pred) .^ 2) / length(y_pred)) / y_true
nme(y_true::Number, y_pred) = sum(abs.(y_true .- y_pred)) / length(y_pred) / y_true
nme(y_true::Array, y_pred) = sum(abs.(y_true .- y_pred)) / length(y_pred) / (maximum(y_true) - minimum(y_true))
end

