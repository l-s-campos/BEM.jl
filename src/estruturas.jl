
abstract type DadosBEM end
abstract type escalar <: DadosBEM end
abstract type vetorial <: DadosBEM end
@kwdef mutable struct subregioes
    regiao::Array{Int64,1} = zeros(Int64, 0)
    equivale::Array{Int64,2} = zeros(Int64, 0, 0)
    Hc::Array{Float64,2} = zeros(Int64, 0, 0)
end
@kwdef mutable struct elemento
    indices::Vector{Int64}
    tipoCDC::Int64
    valorCDC::Vector{Float64}
    ξs::Vector{Float64}
    comprimento::Float64
    regiao::Int64
end
@kwdef mutable struct elementov
    indices::Vector{Int64}
    tipoCDC::Vector{Int64}
    valorCDC::Matrix{Float64}
    ξs::Vector{Float64}
    comprimento::Float64
    regiao::Int64
end

@kwdef mutable struct bezier
    indices::Vector{Int64}
    C::Array{Float64,2}
    p::Int64
    limites::Array{Float64,2}
    Wb::Array{Float64,1}
    sing::Array{Int64,1}
    eets::Array{Float64,1}
end
@kwdef mutable struct elastico <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    normal::Array{Float64,2}
    E::Array{Float64,1}
    ν::Array{Float64,1}
    subregioes::subregioes = subregioes()
    cantos::Array{Float64,2} = zeros(Int64, 0, 0)
end

@kwdef mutable struct elastico_aniso <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    normal::Array{Float64,2}
    k::NamedTuple
end
@kwdef mutable struct potencial <: escalar
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elemento,1}
    normal::Array{Float64,2}
    k::Float64
end
@kwdef mutable struct helmholtz <: escalar
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elemento,1}
    normal::Array{Float64,2}
    k::NamedTuple
end

@kwdef mutable struct potencial_iga <: escalar
    NOS::Array{Float64,2}
    pontos_controle::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    tipoCDC::Array{Int64,2}
    valorCDC::Array{Float64,2}
    ELEM::Array{bezier,1}
    k::Float64
    E::SparseMatrixCSC{Float64,Int64}
end
@kwdef mutable struct potencial_iga_3d <: escalar
    NOS::Array{Float64,2}
    pontos_controle::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    tipoCDC::Array{Int64,2}
    valorCDC::Array{Float64,2}
    ELEM::Array{bezier,1}
    k::Float64
    E::SparseMatrixCSC{Float64,Int64}
end

@kwdef mutable struct elastico_iga <: vetorial
    NOS::Array{Float64,2}
    pontos_controle::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    tipoCDC::Array{Int64,2}
    valorCDC::Array{Float64,2}
    ELEM::Array{bezier,1}
    Ev::Array{Float64,1}
    E::SparseMatrixCSC{Float64,Int64}
end

@kwdef mutable struct elastico_aniso_iga <: vetorial
    NOS::Array{Float64,2}
    pontos_controle::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    tipoCDC::Array{Int64,2}
    valorCDC::Array{Float64,2}
    ELEM::Array{bezier,1}
    k::NamedTuple
    E::SparseMatrixCSC{Float64,Int64}
end
@kwdef mutable struct placa_fina <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    normal::Array{Float64,2}
    k::NamedTuple
end
@kwdef mutable struct placa_fina_isotropica <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    normal::Array{Float64,2}
    k::NamedTuple
end

@kwdef mutable struct placa_espessa <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    normal::Array{Float64,2}
    k::NamedTuple
end
@kwdef mutable struct placa_espessa_isotropica <: vetorial
    NOS::Array{Float64,2}
    pontos_internos::Array{Float64,2}
    ELEM::Array{elementov,1}
    normal::Array{Float64,2}
    k::NamedTuple
end

@kwdef mutable struct casca <: vetorial
    dadplaca::placa_fina_isotropica
    dadpe::elastico
    R11::Float64
    R22::Float64
end
@kwdef mutable struct casca_aniso <: vetorial
    dadplaca::placa_fina
    dadpe::elastico_aniso
    R11::Float64
    R22::Float64
end

@kwdef mutable struct hmat
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
@kwdef mutable struct kernelH <: AbstractMatrix{Float64}
    dad::DadosBEM
    integralelem::Vector{Float64}
end
@kwdef mutable struct kernelG <: AbstractMatrix{Float64}
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
    Qast, Tast = calsolfund([xj - xi, yj - yi], K.dad.normal[j, :], K.dad)
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


@kwdef mutable struct kernelHv <: AbstractMatrix{Float64}
    dad::DadosBEM
    integralelem::Vector{Float64}
    Hp::SparseMatrixCSC
end
@kwdef mutable struct kernelGv <: AbstractMatrix{Float64}
    dad::DadosBEM
    integralelem::Vector{Float64}
    Gp::SparseMatrixCSC
end
function Base.getindex(K::kernelHv, ii::Int, jj::Int)
   i= ceil(Int,ii/2)
   j= ceil(Int,jj/2)
    if i > nc(K.dad)
        xi = K.dad.pontos_internos[i-nc(K.dad), :]
    else
        xi = K.dad.NOS[i, :]
    end

    if K.Hp[ii,jj] ==0.0
        if j > nc(K.dad)
            xj = K.dad.pontos_internos[j-nc(K.dad), :]
        else
            xj = K.dad.NOS[j, :]
        end
        if j > nc(K.dad)
            return 0.0
        end
        uast, tast = calsolfund(xj,xi, K.dad.normal[j, :], K.dad)
        # @show ii jj i j ii-2i+2 jj-2j+2
        return tast[ii-2i+2,jj-2j+2] * K.integralelem[j]
    else
        return K.Hp[ii,jj]
    end

end
function Base.getindex(K::kernelGv, ii::Int, jj::Int)
   i= ceil(Int,ii/2)
   j= ceil(Int,jj/2)
    if i > nc(K.dad)
        xi = K.dad.pontos_internos[i-nc(K.dad), :]
    else
        xi = K.dad.NOS[i, :]
    end

    if i == j
        return 0.0
    end

    if K.Gp[ii,jj] ==0.0
        if j > nc(K.dad)
            xj = K.dad.pontos_internos[j-nc(K.dad), :]
        else
            xj = K.dad.NOS[j, :]
        end
        if j > nc(K.dad)
            return 0.0
        end
        uast, tast = calsolfund(xj,xi, K.dad.normal[j, :], K.dad)
        return uast[ii-2i+2,jj-2j+2] * K.integralelem[j]
    else
        return K.Gp[ii,jj]
    end
end
Base.size(K::kernelHv) = 2*nc(K.dad), 2*nc(K.dad)
Base.size(K::kernelGv) = 2*nc(K.dad), 2*nc(K.dad)


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
    Point2D
