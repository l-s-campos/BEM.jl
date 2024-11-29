## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 6  #Numero de elementos
# NPX = 2 #pontos internos na direção x
# NPY = 2 #pontos internos na direção y
npi = 80
di = 1e-1

pontosi = [range(-7.5 + di, 7.5 - di, length=npi) zeros(npi)]
npg = 10    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
# dad = format_dad(Subregiao(nelem, 3), NPX, NPY) # dados
dad = format_dad(Subregiao(nelem, 3), pontosi) # dados
# dad = format_dad(Subregiao1d(nelem, 2), 1, 2) # dados
# dad = format_dad(elastico1d(nelem, 2), 1, 2) # dados
# dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados
function solucaoAnalitica_adesivo(dad, P, t, ta, c)
    # Goland Reissner
    E, nu, Ecola, nucola = dad.k.E[1], dad.k.nu[1], dad.k.E[2], dad.k.nu[2]
    Gcola = Ecola / (2 * (1 + nucola))
    # P is the applied load per unit width
    # c is half the length of the bonded zone
    # t is the thickness of the adhesive, nu is the Poisson's ratio
    Ga = Gcola
    Ea = Ecola
    x = dad.pontos_internos[:, 1]
    # x = range(-7.5, 7.5, length=100)

    beta = sqrt(8 * Ga * t / E / ta)
    u2 = 1 / t * sqrt(3 * (1 - nu^2) * P / 2 / t / E)

    # Flector moment factor k
    k = cosh(u2 * c) / (cosh(u2 * c) + 2 * sqrt(2) * sinh(u2 * c))

    tau = -P / (8 * c) * (beta * c / t * (1 + 3 * k) * cosh.(beta * c * x / t / c) / sinh.(beta * c / t) .+ 3 * (1 - k))

    gamma = (6 * Ea * t / E / ta)^(1 / 4)
    lambda = gamma * c / t

    delta = 1 / 2 * (sin(2 * lambda) + sinh(2 * lambda))
    kappal = k * c / t * sqrt(3 * (1 - nu^2) * P / t / E)
    R1 = cosh(lambda) * sin(lambda) + sinh(lambda) * cos(lambda)
    R2 = -cosh(lambda) * sin(lambda) + sinh(lambda) * cos(lambda)

    sigma = P * t / delta / c^2 * ((R2 * lambda^2 * k / 2 + lambda * kappal * cosh(lambda) * cos(lambda)) * cosh.(lambda * x / c) .* cos.(lambda * x / c) + (R1 * lambda^2 * k / 2 + lambda * kappal * sinh(lambda) * sin(lambda)) * sinh.(lambda * x / c) .* sin.(lambda * x / c))

    return sigma, tau
end
println("2. Montando a matriz A e o vetor b")
H, G = calc_HeG(dad, npg, subcheia=true)  #importante


A, b = BEM.aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
# A, b = aplicaCDC(H, G, dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A \ b
println("4. Separando fluxo e temperatura")
u, t = separa(dad, x) #importante
tens_int = calc_tens_int(dad, u, t, fill(2, size(dad.pontos_internos)), 30)

sigma_ana, tau_ana = solucaoAnalitica_adesivo(dad, 40, 2, 0.2, 7.5)
p = lines(dad.pontos_internos[:, 1], sigma_ana)
lines!(dad.pontos_internos[:, 1], tens_int[:, 2])
p
p2 = lines(dad.pontos_internos[:, 1], tau_ana)
lines!(dad.pontos_internos[:, 1], tens_int[:, 3])
p2

# A1,b1 = calc_Aeb(dad,npg)  #importante
# x1=A1\b1
# println("5. Calculando nos pontos internos")
# Ti = calc_Ti(dad,T,q,npg);

# println("6. Gerando mapa de cor")
# erro=norm([T-dad.NOS[:,1];Ti-dad.pontos_internos[:,1]])

# # geo=mostra_geometria(dad)
# # mapa=mostra_resultado(dad,[T;Ti])

