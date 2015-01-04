function mgcep{T<:FloatingPoint,N}(x::Array{T,N}, order::Int, α::FloatingPoint,
                                   γ::FloatingPoint; kargs...)
    raw = SPTK.mgcep(x, order, α, γ; kargs...)
    MelGeneralizedCepstrum{Mel,GeneralizedLog,T,N}(α, γ, raw)
end

function gcep{T<:FloatingPoint,N}(x::Array{T,N}, order::Int, γ::FloatingPoint;
                                  kargs...)
    raw = SPTK.gcep(x, order, γ; kargs...)
    MelGeneralizedCepstrum{Linear,GeneralizedLog,T,N}(zero(T), γ, raw)
end
