using BEM
using Revise, Infiltrator
using TimerOutputs
# using Plots
# gaston()
# using GLMakie
# Makie.inline!(false)
using FileIO, LinearAlgebra, Statistics
using Krylov, SparseArrays, NonlinearSolve

includet(datadir("dadpotencial.jl"))
includet(datadir("dadelastico.jl"))
includet(datadir("dadcontato.jl"))
includet(datadir("dadelastico_superposicao.jl"))
includet(datadir("dadhelmholtz.jl"))
includet(datadir("dadplacafina.jl"))
includet(datadir("dadplacaespessa.jl"))
