## Início da análise
using DrWatson,
    LinearMaps, DataFrames, CSV, Plots, DelimitedFiles, HMatrices, BenchmarkTools
@quickactivate "BEM"
include(scriptsdir("includes.jl"))

#Escreve resultados
res = DataFrame(
    prob = String[],
    order = Int[],
    n_contorno = Int[],
    n_interno = Int[],
    n_gdl = Int[],
    t_BEM = Array{Float64,1}[],
    t_DIBEM = Array{Float64,1}[],
    t_DIHBEM = Array{Float64,1}[],
    ttotal_BEM = Float64[],
    ttotal_DIBEM = Float64[],
    ttotal_DIHBEM = Float64[],
    e = Float64[],
    ed = Float64[],
    eH = Float64[],
    em = Float64[],
    edm = Float64[],
    eHm = Float64[],
)

prob = [helm1d, helmdirechlet, helmcirculo]
prob_str = ["helm1d", "helmdirechlet", "helmcirculo"]
ana = [ANA_helm1d, ANA_helmdirechlet, ANA_helmcirculo]

npg = 10    #apenas números pares
FR = 4 #frequencia de excitação
atol = 1e-7
rtol = 1e-6


# function residuo!(x)

#     TH[idx_1] = x[idx_1]
#     qH[idx_0] = x[idx_0]

#     THi = x[size(TH, 1)+1:end]
#     return HH * TH - H_extend * THi + MH * [TH; THi] * ((ω / c)^2) - GH * qH + bH
# end

# function residuo2!(x)

#     TH[idx_1] = x[idx_1]
#     qH[idx_0] = x[idx_0]

#     THi = x[size(TH, 1)+1:end]
#     return HH' * TH - H_extend' * THi + MH' * [TH; THi] * ((ω / c)^2) - GH' * qH + bH
# end

for i = 2:7
    base = 1.6
    nelem = Int(round(2 * base^i)) #Numero de elementos
    order = 2

    NPX = Int(round((3 / 2) * nelem)) #pontos internos na direção x
    NPY = NPX #pontos internos na direção y    

    n_contorno = 4 * nelem * order
    n_interno = NPX * NPY
    n_total = n_contorno + n_interno

    println(
        n_contorno,
        " pontos de contorno e ",
        n_interno,
        " pontos internos",
        " resultando em ",
        n_total,
        " totais",
    )
end

for p = 1:1, i = 2:7, j = 2:2

    # nelem = 2 * 2^(i-1) #Numero de elementos
    # order = 2

    # NPX = 10*i #pontos internos na direção x
    # NPY = NPX #pontos internos na direção y 

    base = 1.6
    nelem = Int(round(2 * base^i)) #Numero de elementos
    order = 2

    NPX = Int(round((3 / 2) * nelem)) #pontos internos na direção x
    NPY = NPX #pontos internos na direção y

    dad = format_dad(prob[p](nelem, order, FR), NPX, NPY) # dados

    if p == 3
        # println("Gera pontos internos dentro do circulo")
        # NT = 100 # número de pontos na circunferência maior
        # NR = Int(round(calcula_NR(NT,1,NPX*NPY)))
        # #NR = 10  # número de divisões radiais
        # dad.pontos_internos = pontos_internos_circulo(NR, NT, 1; Δr=0.05)
    elseif p == 2
        println("Corrigindo CDC de Dirichlet")
        corrigeCDC_helmdirechlet(dad)
        if FR <= π
            println("No caso 2, FR deve ser maior que pi")
            FR = 4
        end
    end

    n_contorno = nc(dad)
    n_interno = ni(dad)

    println(
        "Caso: ",
        prob[p],
        " com ",
        n_contorno,
        " pontos de contorno e ",
        n_interno,
        " pontos internos",
    )

    dadpot, ω, c = potencial_correlato(dad::helmholtz)

    #Matrizes HeG BEM
    tHG = @elapsed H, G = BEM.calc_HeG(dadpot, npg, interno = true)

    #Matrizes HeG DIBEM
    tHGd = @elapsed Hd, Gd = BEM.calc_HeGd(dadpot, interno = true)

    tM = @elapsed M = BEM.Monta_M_RIMdd(dadpot, npg, tiporadial = "tps")

    #A e b
    tb = @elapsed A, b = aplicaCDC(H + ((ω / c)^2) * M, G, dadpot)
    tbd = @elapsed Ad, bd = aplicaCDC(Hd + ((ω / c)^2) * M, Gd, dadpot)

    tsolve = @elapsed x, f = BEM.gmres(A, b, rtol = rtol, itmax = 100)
    T, q = separa(dadpot, x)
    Ti = x[nc(dad)+1:end]

    tsolved = @elapsed xd, f = BEM.gmres(Ad, bd, rtol = rtol, itmax = 100)
    Td, qd = separa(dadpot, xd)
    Tid = xd[nc(dad)+1:end]

    #Matrizes DIHBEM

    # tHG_H = @elapsed HH, GH = BEM.calc_HeG_Hd(dadpot, atol = atol, eta = 3, nmax = 10)
    # compression_ratio_HH = BEM.compression_ratio(HH)
    # println("compression_ratio_HH = $compression_ratio_HH")

    tM_H = @elapsed MH = BEM.Monta_M_HddH(
        dadpot,
        npg,
        tiporadial = "tps",
        atol = atol,
        eta = 3,
        nmax = 10,
    )

    tHG_H = @elapsed HH, GH = BEM.preencheH(MH, dadpot, atol = atol, rtol = rtol)


    tb_H = @elapsed begin
        AH = HH + rmul!(MH, (ω / c)^2)
        bH, AH, _ = aplicaCDC(AH, GH, dadpot)
    end


    tsolveH = @elapsed xH, fH = BEM.gmres(AH, bH, rtol = rtol, itmax = 100)
    TH, qH = separa(dadpot, xH)
    TiH = xH[nc(dad)+1:end]

    # tb_H = @elapsed begin
    #     AH = HH + rmul!(MH,(ω / c)^2)
    #     bH = aplicaCDC(AH, GH, dadpot)
    # end

    # compression_ratio_AH = BEM.compression_ratio(AH)
    # println("compression_ratio_AH = $compression_ratio_AH")

    #por linear map

    # tb_H = @elapsed begin
    #     bH = aplicaCDC(HH, GH, dadpot)[1]
    #     aux = [BEM.valorCDC(dadpot); zeros(ni(dadpot))]
    #     aux[BEM.tipoCDC(dadpot)] .= 0.0

    #     bH = bH - MH * aux * ((ω / c)^2)
    # end


    # function residuo!(x)

    #     TH[idx_1] = x[idx_1]
    #     qH[idx_0] = x[idx_0]

    #     THi = x[size(TH,1)+1:end]
    #     return HH * TH - H_extend*THi + MH * [TH;THi]* ((ω / c)^2) - GH*qH  + bH
    # end

    # function residuo2!(x)

    #     TH[idx_1] = x[idx_1]
    #     qH[idx_0] = x[idx_0]

    #     THi = x[size(TH,1)+1:end]
    #     return HH' * TH - H_extend'*THi + MH' * [TH;THi] * ((ω / c)^2) - GH'*qH  + bH
    # end

    # tsolve_H = @elapsed begin

    #     idx_0 = []
    #     idx_1 = []
    #     TH = zeros(nc(dadpot))
    #     qH = zeros(nc(dadpot))
    #     val_CDC = []

    #     for elem_i in dadpot.ELEM  # Laço dos pontos fontes
    #         ind_elem = elem_i.indices
    #         if elem_i.tipoCDC == 0
    #             append!(idx_0, ind_elem)
    #             append!(val_CDC, elem_i.valorCDC)
    #         elseif elem_i.tipoCDC == 1
    #             append!(idx_1, ind_elem)
    #             append!(val_CDC, elem_i.valorCDC)
    #         end
    #     end

    #     TH[idx_0] = val_CDC[idx_0]
    #     qH[idx_1] = val_CDC[idx_1]

    #     m, n = size(MH)

    #     H_extend = sparse([zeros(nc(dad), ni(dad));I])
    #     AH = LinearMap(residuo!, residuo2!, m, n)

    #     xH, fH = BEM.gmres(AH, bH, rtol = rtol, itmax = 100) 
    #     TH, qH = separa(dadpot, xH) 
    #     TiH = xH[nc(dad)+1:end]
    # end

    #Calcula erros
    if p == 1
        Ta = -ANA_helm1d(dad, dad.pontos_internos)
    elseif p == 2
        Ta = ANA_helmdirechlet(dad, dad.pontos_internos)
    else
        Ta = ANA_helmcirculo(dad, dad.pontos_internos)
    end

    e = norm(Ti - Ta) / norm(Ta)
    ed = norm(Tid - Ta) / norm(Ta)
    eH = norm(TiH - Ta) / norm(Ta)

    em = sum(abs, Ti - Ta) / sum(abs, Ta)
    edm = sum(abs, Tid - Ta) / sum(abs, Ta)
    eHm = sum(abs, TiH - Ta) / sum(abs, Ta)

    #Tempos

    t_BEM = [tHG, tM, tb, tsolve]

    t_DIBEM = [tHGd, tM, tbd, tsolved]

    t_DIHBEM = [tHG_H, tM_H, tb_H, tsolveH]

    ttotal_BEM = sum(t_BEM)
    ttotal_DIBEM = sum(t_DIBEM)
    ttotal_DIHBEM = sum(t_DIHBEM)

    push!(
        res,
        [
            string(prob[p]),
            order,
            n_contorno,
            n_interno,
            n_contorno + n_interno,
            t_BEM,
            t_DIBEM,
            t_DIHBEM,
            ttotal_BEM,
            ttotal_DIBEM,
            ttotal_DIHBEM,
            e,
            ed,
            eH,
            em,
            edm,
            eHm,
        ],
    )
end

#Plots results

#Plota tempo especifico
t_BEM = [res.t_BEM[i][1] for i = 1:size(res, 1)]
t_DIBEM = [res.t_DIBEM[i][1] for i = 1:size(res, 1)]
t_DIHBEM = [res.t_DIHBEM[i][1] for i = 1:size(res, 1)]

Plots.plot(
    res[res.prob.==prob_str[1], :n_gdl],
    t_BEM[:],
    label = "BEM",
    marker = :o,
    yscale = :log10,
    xscale = :log10,
    xlabel = "dof",
    ylabel = "Tempo (s)",
    legend = :bottomright,
)
Plots.plot!(res[res.prob.==prob_str[1], :n_gdl], t_DIBEM[:], label = "DIBEM", marker = :o)
Plots.plot!(res[res.prob.==prob_str[1], :n_gdl], t_DIHBEM[:], label = "DIHBEM", marker = :o)

#Plota erros ou tempo total (basta mudar e ou ttotal)

variables = [:ttotal_BEM, :ttotal_DIBEM, :ttotal_DIHBEM]
# variables = [:e, :ed, :eH]

Plots.plot(
    res[(res.prob.==prob_str[1]).&(res.order.==2), :n_gdl],
    res[(res.prob.==prob_str[1]).&(res.order.==2), variables[1]],
    label = "BEM",
    marker = :o,
    yscale = :log10,
    xscale = :log10,
    xlabel = "dof",
    ylabel = "y",
    legend = :bottomright,
)

Plots.plot!(
    res[(res.prob.==prob_str[1]).&(res.order.==2), :n_gdl],
    res[(res.prob.==prob_str[1]).&(res.order.==2), variables[2]],
    label = "DIBEM",
    marker = :o,
)

Plots.plot!(
    res[(res.prob.==prob_str[1]).&(res.order.==2), :n_gdl],
    res[(res.prob.==prob_str[1]).&(res.order.==2), variables[3]],
    label = "DIHBEM",
    marker = :o,
)


# #CSV dos results

# #salva como CSV
# CSV.write("helm_res.csv", res)

# #le arquivo csv
# res = CSV.read("helm_res.csv", DataFrame)

# #escreve arquivo .dat pra fazer o pdf no overleaf

# for i = 1:3, order = 2:3
#     println("Escrevendo arquivos .dat do caso ", prob_str[i])
#     prob = prob_str[i]

#     DelimitedFiles.writedlm(
#         "$prob-BEM-$order-time.dat",
#         [res[res.prob.==prob_str[i], :n_gdl] res[res.prob.==prob_str[i], :ttotal_BEM]],
#         '\t',
#     )
#     DelimitedFiles.writedlm(
#         "$prob-DIBEM-$order-time.dat",
#         [res[res.prob.==prob_str[i], :n_gdl] res[res.prob.==prob_str[i], :ttotal_DIBEM]],
#         '\t',
#     )
#     DelimitedFiles.writedlm(
#         "$prob-DIHBEM-$order-time.dat",
#         [res[res.prob.==prob_str[i], :n_gdl] res[res.prob.==prob_str[i], :ttotal_DIHBEM]],
#         '\t',
#     )

#     DelimitedFiles.writedlm(
#         "$prob-BEM-$order-error.dat",
#         [res[res.prob.==prob_str[i], :n_gdl] res[res.prob.==prob_str[i], :e]],
#         '\t',
#     )
#     DelimitedFiles.writedlm(
#         "$prob-DIBEM-$order-error.dat",
#         [res[res.prob.==prob_str[i], :n_gdl] res[res.prob.==prob_str[i], :ed]],
#         '\t',
#     )
#     DelimitedFiles.writedlm(
#         "$prob-DIHBEM-$order-error.dat",
#         [res[res.prob.==prob_str[i], :n_gdl] res[res.prob.==prob_str[i], :eH]],
#         '\t',
#     )
# end

# #Para caso específico

# DelimitedFiles.writedlm(
#     "helm1d_BEM_order2_time.dat",
#     [res[res.prob.==prob_str[1], :n_gdl] res[res.prob.==prob_str[1], :ttotal_BEM]],
#     '\t',
# )
# DelimitedFiles.writedlm(
#     "helm1d_DIBEM_order2_time.dat",
#     [res[res.prob.==prob_str[1], :n_gdl] res[res.prob.==prob_str[1], :ttotal_DIBEM]],
#     '\t',
# )
# DelimitedFiles.writedlm(
#     "helm1d_DIHBEM_order2_time.dat",
#     [res[res.prob.==prob_str[1], :n_gdl] res[res.prob.==prob_str[1], :ttotal_DIHBEM]],
#     '\t',
# )

# DelimitedFiles.writedlm(
#     "helm1d_BEM_order2_error.dat",
#     [res[res.prob.==prob_str[1], :n_gdl] res[res.prob.==prob_str[1], :e]],
#     '\t',
# )
# DelimitedFiles.writedlm(
#     "helm1d_DIBEM_order2_error.dat",
#     [res[res.prob.==prob_str[1], :n_gdl] res[res.prob.==prob_str[1], :ed]],
#     '\t',
# )
# DelimitedFiles.writedlm(
#     "helm1d_DIHBEM_order2_error.dat",
#     [res[res.prob.==prob_str[1], :n_gdl] res[res.prob.==prob_str[1], :eH]],
#     '\t',
# )
