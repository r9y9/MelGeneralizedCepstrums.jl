function mgcep{T<:FloatingPoint,N}(x::Array{T,N}, order::Int, α::Float64,
                                   γ::Float64; kargs...)
    raw = SPTK.mgcep(x, order, α, γ; kargs...)
    MelGeneralizedCepstrum{Mel,GeneralizedLog,T,N}(α, γ, raw)
end

function mcep{T<:FloatingPoint,N}(x::Array{T,N}, order::Int, α::Float64;
                                  kargs...)
    raw = SPTK.mcep(x, order, α; kargs...)
    MelGeneralizedCepstrum{Mel,StandardLog,T,N}(α, zero(T), raw)
end

function gcep{T<:FloatingPoint,N}(x::Array{T,N}, order::Int, γ::Float64;
                                  kargs...)
    raw = SPTK.gcep(x, order, γ; kargs...)
    MelGeneralizedCepstrum{Linear,GeneralizedLog,T,N}(zero(T), γ, raw)
end
