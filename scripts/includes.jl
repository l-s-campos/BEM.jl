using BEM
using Revise,Infiltrator
using TimerOutputs
using GLMakie
using FileIO, JLD2,LinearAlgebra
using KrylovMethods,SparseArrays

includet(datadir("dadpotencial.jl"))
includet(datadir("dadelastico.jl"))
includet(datadir("dadhelmholtz.jl"))
