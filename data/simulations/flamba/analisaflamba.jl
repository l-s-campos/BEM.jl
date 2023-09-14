using DataFrames
using DrWatson
# using CairoMakie, MakiePublication
# CairoMakie.activate!()
# dados = collect_results(datadir("simulations/flamba"); rinclude=[r"problema=ex1"])
dados = collect_results(datadir("simulations/flamba"); rinclude=[r"problema=placa_furo_retangular"])
vscodedisplay(dados)

