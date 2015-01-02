# Common types and functions

# Mel-Generalized Cepstrum is parametrized by type of frequency scale (F) and
# log function (L). It is designed to be a subtype of AbstractArray{T,N}.
abstract AbstractMelGeneralizedCepstrumArray{F,L,T,N} <: AbstractArray{T,N}

abstract Frequency
type Mel <: Frequency
end
type Linear <: Frequency
end

abstract Log
type GeneralizedLog <: Log
end
type StandardLog <: Log
end
type AllPoleLog <: Log
end
type AllZeroLog <: Log
end

frequency_scale{F<:Frequency,L<:Log,T,N}(::Type{AbstractMelGeneralizedCepstrumArray{F,L,T,N}}) = F
frequency_scale{T<:AbstractMelGeneralizedCepstrumArray}(::Type{T}) = frequency_scale(super(T))

log_func{F<:Frequency,L<:Log,T,N}(::Type{AbstractMelGeneralizedCepstrumArray{F,L,T,N}}) = L
log_func{T<:AbstractMelGeneralizedCepstrumArray}(::Type{T}) = log_func(super(T))

immutable MelGeneralizedCepstrum{F<:Frequency,L<:Log,T<:FloatingPoint,N} <: AbstractMelGeneralizedCepstrumArray{F,L,T,N}
    α::T
    γ::T
    data::Array{T,N}

    function MelGeneralizedCepstrum(α::T, γ::T, data::Array{T,N})
        abs(α) < 1 || error("|α| < 1 is supported")
        (-1 <= γ <= 0) || error("-1 <= γ <= 0 is supported")
        @assert length(data) > 1
        new(α, γ, data)
    end
end

function MelGeneralizedCepstrum{T,N}(α::T, γ::T, data::Array{T,N})
    F = (α == zero(T)) ? Linear : Mel

    L = GeneralizedLog
    if γ == zero(T)
        L = StandardLog
    elseif γ == -one(T)
        L = AllPoleLog
    elseif γ == one(T)
        L = AllZeroLog
    end

    MelGeneralizedCepstrum{F,L,T,N}(α, γ, data)
end

typealias MelFrequencyCepstrum{L,T,N} MelGeneralizedCepstrum{Mel,L,T,N}
typealias LinearFrequencyCepstrum{L,T,N} MelGeneralizedCepstrum{Linear,L,T,N}
typealias GeneralizedLogCepstrum{F,T,N} MelGeneralizedCepstrum{F,GeneralizedLog,T,N}
typealias StandardLogCepstrum{F,T,N} MelGeneralizedCepstrum{F,StandardLog,T,N}
typealias AllPoleCepstrum{F,T,N} MelGeneralizedCepstrum{F,AllPoleLog,T,N}
typealias AllZeroCepstrum{F,T,N} MelGeneralizedCepstrum{F,AllZeroLog,T,N}

typealias MelCepstrum{T,N} MelGeneralizedCepstrum{Mel,StandardLog,T,N}
typealias GeneralizedCepstrum{T,N} MelGeneralizedCepstrum{Linear,GeneralizedLog,T,N}

## AbstractArray implementation ##

Base.eltype(c::MelGeneralizedCepstrum) = eltype(c.data)
Base.length(c::MelGeneralizedCepstrum) = length(c.data)
Base.size(c::MelGeneralizedCepstrum) = size(c.data)
Base.getindex(c::MelGeneralizedCepstrum, i::Real) = getindex(c.data, i)
Base.getindex(c::MelGeneralizedCepstrum, i::Real...) = getindex(c.data, i...)
Base.getindex(c::MelGeneralizedCepstrum, i::Range,j::Real) = getindex(c.data, i, j)
Base.getindex(c::MelGeneralizedCepstrum, i::Real, j::Range) = getindex(c.data, i, j)
Base.eltype(c::MelGeneralizedCepstrum) = Base.eltype(g.data)

## Mel-generalized cepstrum related functions ##

rawdata(c::MelGeneralizedCepstrum) = c.data
order(c::MelGeneralizedCepstrum) = size(c.data, 1) - 1
allpass_alpha(c::MelGeneralizedCepstrum) = c.α
glog_gamma(c::MelGeneralizedCepstrum) = c.γ
powercoef{F,L,T}(c::MelGeneralizedCepstrum{F,L,T,1}) = first(rawdata(c))
powercoef{F,L,T}(c::MelGeneralizedCepstrum{F,L,T,2}) = rawdata(c)[1,:]
