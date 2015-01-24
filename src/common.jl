# Common types and functions

import Base: eltype, length, size, getindex

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

frequency_scale{F<:Frequency,L<:Log,T,N}(::Type{AbstractMelGeneralizedCepstrumArray{F,L,T,N}}) = F
frequency_scale{T<:AbstractMelGeneralizedCepstrumArray}(::Type{T}) = frequency_scale(super(T))

log_func{F<:Frequency,L<:Log,T,N}(::Type{AbstractMelGeneralizedCepstrumArray{F,L,T,N}}) = L
log_func{T<:AbstractMelGeneralizedCepstrumArray}(::Type{T}) = log_func(super(T))

immutable MelGeneralizedCepstrum{F<:Frequency,L<:Log,T<:FloatingPoint,N} <: AbstractMelGeneralizedCepstrumArray{F,L,T,N}
    α::T
    γ::T
    data::Array{T,N}

    function MelGeneralizedCepstrum(α::FloatingPoint, γ::FloatingPoint,
                                    data::Array{T,N})
        abs(α) < 1 || throw(ArgumentError("|α| < 1 is supported"))
        (-1 <= γ <= 0) || throw(ArgumentError("-1 <= γ <= 0 is supported"))
        @assert size(data, 1) > 1
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
    end

    MelGeneralizedCepstrum{F,L,T,N}(α, γ, data)
end

function MelGeneralizedCepstrum(mgc::MelGeneralizedCepstrum)
    α = allpass_alpha(mgc)
    γ = glog_gamma(mgc)
    MelGeneralizedCepstrum(α, γ, rawdata(mgc))
end

typealias MelFrequencyCepstrum{L,T,N} MelGeneralizedCepstrum{Mel,L,T,N}
typealias LinearFrequencyCepstrum{L,T,N} MelGeneralizedCepstrum{Linear,L,T,N}
typealias GeneralizedLogCepstrum{F,T,N} MelGeneralizedCepstrum{F,GeneralizedLog,T,N}
typealias StandardLogCepstrum{F,T,N} MelGeneralizedCepstrum{F,StandardLog,T,N}
typealias AllPoleLogCepstrum{F,T,N} MelGeneralizedCepstrum{F,AllPoleLog,T,N}

typealias LinearCepstrum{T,N} MelGeneralizedCepstrum{Linear,StandardLog,T,N}
typealias MelCepstrum{T,N} MelGeneralizedCepstrum{Mel,StandardLog,T,N}
typealias GeneralizedCepstrum{T,N} MelGeneralizedCepstrum{Linear,GeneralizedLog,T,N}
typealias AllPoleCepstrum{T,N} MelGeneralizedCepstrum{Linear,AllPoleLog,T,N}
typealias MelAllPoleCepstrum{T,N} MelGeneralizedCepstrum{Mel,AllPoleLog,T,N}


immutable MelGeneralizedCepstrumFilterCoef{F<:Frequency,L<:Log,T<:FloatingPoint,N} <: AbstractMelGeneralizedCepstrumArray{F,L,T,N}
    α::T
    γ::T
    data::Array{T,N}

    function MelGeneralizedCepstrumFilterCoef(α::FloatingPoint,
                                              γ::FloatingPoint,
                                              data::Array{T,N})
        abs(α) < 1 || throw(ArgumentError("|α| < 1 is supported"))
        (-1 <= γ <= 0) || throw(ArgumentError("-1 <= γ <= 0 is supported"))
        @assert size(data, 1) > 1
        new(α, γ, data)
    end
end

function MelGeneralizedCepstrumFilterCoef{T,N}(α::T, γ::T, data::Array{T,N})
    F = (α == zero(T)) ? Linear : Mel

    L = GeneralizedLog
    if γ == zero(T)
        L = StandardLog
    elseif γ == -one(T)
        L = AllPoleLog
    end

    MelGeneralizedCepstrumFilterCoef{F,L,T,N}(α, γ, data)
end

function MelGeneralizedCepstrumFilterCoef(mgc::MelGeneralizedCepstrumFilterCoef)
    α = allpass_alpha(mgc)
    γ = glog_gamma(mgc)
    MelGeneralizedCepstrumFilterCoef(α, γ, rawdata(mgc))
end

typealias LMADFCoef{T,N} MelGeneralizedCepstrumFilterCoef{Linear,StandardLog,T,N}
typealias MLSADFCoef{T,N} MelGeneralizedCepstrumFilterCoef{Mel,StandardLog,T,N}
typealias MGLSADFCoef{T,N} MelGeneralizedCepstrumFilterCoef{Linear,GeneralizedLog,T,N}

## AbstractArray implementation ##

eltype(c::AbstractMelGeneralizedCepstrumArray) = eltype(c.data)
length(c::AbstractMelGeneralizedCepstrumArray) = length(c.data)
size(c::AbstractMelGeneralizedCepstrumArray) = size(c.data)
getindex(c::AbstractMelGeneralizedCepstrumArray, i::Real) = getindex(c.data, i)
getindex(c::AbstractMelGeneralizedCepstrumArray, i::Real...) = getindex(c.data, i...)
function getindex(c::MelGeneralizedCepstrum, i::Range, j::Real)
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    MelGeneralizedCepstrum(α, γ, getindex(c.data, i, j))
end
function getindex(c::MelGeneralizedCepstrumFilterCoef, i::Range, j::Real)
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    MelGeneralizedCepstrumFilterCoef(α, γ, getindex(c.data, i, j))
end
getindex(c::AbstractMelGeneralizedCepstrumArray, i::Real, j::Range) = getindex(c.data, i, j)
eltype(c::AbstractMelGeneralizedCepstrumArray) = Base.eltype(g.data)

## Mel-generalized cepstrum related functions ##

rawdata(c::AbstractMelGeneralizedCepstrumArray) = c.data
order(c::AbstractMelGeneralizedCepstrumArray) = size(c.data, 1) - 1
allpass_alpha(c::AbstractMelGeneralizedCepstrumArray) = c.α
glog_gamma(c::AbstractMelGeneralizedCepstrumArray) = c.γ
powercoef{F,L,T}(c::AbstractMelGeneralizedCepstrumArray{F,L,T,1}) = first(rawdata(c))
powercoef{F,L,T}(c::AbstractMelGeneralizedCepstrumArray{F,L,T,2}) = rawdata(c)[1,:]
