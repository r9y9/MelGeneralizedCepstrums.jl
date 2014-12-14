# Mel-generalized cepstrum

immutable MelGeneralizedCepstrum{FS<:FrequencyScale, L<:LogFunc} <: AbstractMelGeneralizedCepstrum{FS,L}
    α::Float64
    γ::Float64
    data::Vector{Float64}
    
    function MelGeneralizedCepstrum(α::Float64, γ::Float64, 
                                    data::Vector{Float64})
        abs(α) < 1 || error("|α| < 1 is supported")
        (-1 <= γ <= 0) || error("-1 <= γ <= 0 is supported")
        @assert length(data) > 1
        new(α, γ, data)
    end
end

function MelGeneralizedCepstrum(α::Float64, γ::Float64, data::Vector{Float64})
    MelGeneralizedCepstrum{Mel, GeneralizedLog}(α, γ, data)
end

Base.eltype(c::MelGeneralizedCepstrum) = eltype(c.data)
Base.length(c::MelGeneralizedCepstrum) = length(c.data)
Base.size(c::MelGeneralizedCepstrum) = size(c.data)
Base.getindex(c::MelGeneralizedCepstrum, i::Real) = getindex(c.data)
rawdata(c::MelGeneralizedCepstrum) = c.data

function mgcep(x::Vector{Float64}, order::Int, α::Float64, γ::Float64)
    raw = SPTK.mgcep(x, order, α, γ)
    MelGeneralizedCepstrum(α, γ, raw)
end
