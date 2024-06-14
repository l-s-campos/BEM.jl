using BEM
using Revise, Infiltrator
using TimerOutputs
# using Plots
# gaston()
# using GLMakie
# Makie.inline!(false)
using FileIO, JLD2, LinearAlgebra, Statistics
using Krylov, LinearMaps, SparseArrays

includet(datadir("dadpotencial.jl"))
includet(datadir("dadelastico.jl"))
includet(datadir("dadelastico_superposicao.jl"))
includet(datadir("dadhelmholtz.jl"))
includet(datadir("dadplacafina.jl"))
includet(datadir("dadplacaespessa.jl"))
