
abstract type DadosBEM end
abstract type escalar <: DadosBEM end
abstract type vetorial <: DadosBEM end
mutable struct subregioes
    regiao::Array{Int64,1}
    equivale::Array{Int64,2}
    Hc::Array{Float64,2}
end
mutable struct elemento
    indices::Array{Int64,1}
    tipoCDC::Int64
    valorCDC::Array{Float64,1}
    ξs::Array{Float64,1}
    comprimento::Float64
    regiao::Int64
end
mutable struct elementov
    indices::Array{Int64,1}
    tipoCDC::Array{Int64,1}
    valorCDC::Array{Float64,2}
    ξs::Array{Float64,1}
    comprimento::Float64
    regiao::Int64
end
mutable struct bezier
    indices::Array{Int64,1}
    C::Array{Float64,2}
    p::Int64
    limites::Array{Float64,2}
    Wb::Array{Float64,1}
    sing::Array{Int64,1}
    eets::Array{Float64,1}
end
mutable struct elastico <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    k::NamedTuple
end

mutable struct elastico_aniso <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    k::NamedTuple
end
mutable struct potencial <: escalar
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elemento,1}
    k::Float64
end
mutable struct helmholtz <: escalar
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elemento,1}
    k::NamedTuple
end

mutable struct potencial_iga <: escalar
    NOS::Array{Float64,2}
    pontos_controle::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    tipoCDC::Array{Int64,2}
    valorCDC::Array{Float64,2}
    ELEM::Array{bezier,1}
    k::Float64
    E::SparseMatrixCSC{Float64,Int64}
end
mutable struct potencial_iga_3d <: escalar
    NOS::Array{Float64,2}
    pontos_controle::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    tipoCDC::Array{Int64,2}
    valorCDC::Array{Float64,2}
    ELEM::Array{bezier,1}
    k::Float64
    E::SparseMatrixCSC{Float64,Int64}
end

mutable struct elastico_iga <: vetorial
    NOS::Array{Float64,2}
    pontos_controle::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    tipoCDC::Array{Int64,2}
    valorCDC::Array{Float64,2}
    ELEM::Array{bezier,1}
    Ev::Array{Float64,1}
    E::SparseMatrixCSC{Float64,Int64}
end

mutable struct elastico_aniso_iga <: vetorial
    NOS::Array{Float64,2}
    pontos_controle::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    tipoCDC::Array{Int64,2}
    valorCDC::Array{Float64,2}
    ELEM::Array{bezier,1}
    k::NamedTuple
    E::SparseMatrixCSC{Float64,Int64}
end
mutable struct placa_fina <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    k::NamedTuple
end
mutable struct placa_fina_isotropica <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    k::NamedTuple
end

mutable struct placa_espessa <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    k::NamedTuple
end
mutable struct placa_espessa_isotropica <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    k::NamedTuple
end

mutable struct casca <: vetorial
    dadplaca::placa_fina_isotropica
    dadpe::elastico
    R11::Float64
    R22::Float64
end
mutable struct casca_aniso <: vetorial
    dadplaca::placa_fina
    dadpe::elastico_aniso
    R11::Float64
    R22::Float64
end

mutable struct hmat
    A::Matrix{Matrix{Float64}}
    block::Matrix{Int64}
    Tree1::Vector{Vector{Int64}}
    Tree2::Vector{Vector{Int64}}
    dad::DadosBEM
    cols::Vector{Vector{Int64}}
    n::Int
end
nc(dad::DadosBEM) = size(dad.NOS, 1)
ni(dad::DadosBEM) = size(dad.pontos_internos, 1)
const Point2D = BEM.SVector{2,Float64}


import Base.size
size(e::elemento) = size(e.indices, 1)
size(e::bezier) = size(e.indices, 1)
size(e::elementov) = size(e.indices, 1)
"""
KH=kernelH(dad,BEM.calc_normais(dad))

"""
mutable struct kernelH <: AbstractMatrix{Float64}
    dad::DadosBEM
    normalfonte::Matrix{Float64}
    integralelem::Vector{Float64}
end
mutable struct kernelG <: AbstractMatrix{Float64}
    dad::DadosBEM
    integralelem::Vector{Float64}
end
function Base.getindex(K::kernelH, i::Int, j::Int)
    if i > nc(K.dad)
        xi = K.dad.pontos_internos[i-nc(K.dad), 1]
        yi = K.dad.pontos_internos[i-nc(K.dad), 2]
    else
        xi = K.dad.NOS[i, 1]
        yi = K.dad.NOS[i, 2]
    end

    if i == j
        return 0.0
    end

    if j > nc(K.dad)
        xj = K.dad.pontos_internos[j-nc(K.dad), 1]
        yj = K.dad.pontos_internos[j-nc(K.dad), 2]
    else
        xj = K.dad.NOS[j, 1]
        yj = K.dad.NOS[j, 2]
    end
    if j > nc(K.dad)
        return 0.0
    end
    Qast, Tast = calsolfund([xj - xi, yj - yi], K.normalfonte[j, :], K.dad)
    return Qast * K.integralelem[j]
end
function Base.getindex(K::kernelG, i::Int, j::Int)
    if i > nc(K.dad)
        xi = K.dad.pontos_internos[i-nc(K.dad), 1]
        yi = K.dad.pontos_internos[i-nc(K.dad), 2]
    else
        xi = K.dad.NOS[i, 1]
        yi = K.dad.NOS[i, 2]
    end

    if i == j
        return 0.0
    end

    if j > nc(K.dad)
        xj = K.dad.pontos_internos[j-nc(K.dad), 1]
        yj = K.dad.pontos_internos[j-nc(K.dad), 2]
    else
        xj = K.dad.NOS[j, 1]
        yj = K.dad.NOS[j, 2]
    end
    if j > nc(K.dad)
        return 0.0
    end
    Qast, Tast = calsolfund([xj - xi, yj - yi], [0, 0], K.dad)
    return Tast * K.integralelem[j]
end
Base.size(K::kernelH) = nc(K.dad), nc(K.dad)
Base.size(K::kernelG) = nc(K.dad), nc(K.dad)


export potencial,
    helmholtz,
    potencial_iga,
    elastico,
    elastico_iga,
    elemento,
    elastico_aniso,
    elastico_aniso_iga,
    hmat,
    placa_fina,
    placa_fina_isotropica,
    casca,
    casca_aniso,
    placa_espessa,
    placa_espessa_isotropica,
    matrizA
