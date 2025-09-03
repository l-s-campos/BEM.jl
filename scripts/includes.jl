using BEM
using Revise, Infiltrator
using TimerOutputs
using FileIO, LinearAlgebra, Statistics
using SparseArrays, NonlinearSolve

includet(datadir("dadpotencial.jl"))
includet(datadir("dadelastico.jl"))
includet(datadir("dadcontato.jl"))
includet(datadir("dadelastico_superposicao.jl"))
includet(datadir("dadhelmholtz.jl"))
includet(datadir("dadplacafina.jl"))
includet(datadir("dadplacaespessa.jl"))
includet(datadir("dad_CFD.jl"))
