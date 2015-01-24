# Linear Prediction Coefficients (LPC)

abstract AbstractLinearPredictionCoef{F,T,N} <: AbstractMelGeneralizedCepstrumArray{F,AllPoleLog,T,N}

immutable MelLinearPredictionCoef{F<:Frequency,T<:FloatingPoint,N} <: AbstractLinearPredictionCoef{F,T,N}
    α::T
    γ::T
    data::Array{T,N}
    loggain::Bool

    function MelLinearPredictionCoef(α::FloatingPoint,
                                     data::Array{T,N},
                                     loggain::Bool)
        abs(α) < 1 || throw(ArgumentError("|α| < 1 is supported"))
        @assert size(data, 1) > 1
        new(α, -1.0, data, loggain)
    end
end

function getindex(c::MelLinearPredictionCoef, i::Range, j::Real)
    α = allpass_alpha(c)
    MelLinearPredictionCoef(α, getindex(c.data, i, j), c.loggain)
end

typealias LinearPredictionCoef{T,N} MelLinearPredictionCoef{Linear,T,N}

function MelLinearPredictionCoef{T,N}(α::T, data::Array{T,N}, loggain::Bool)
    F = ifelse(α == zero(T), Linear, Mel)
    MelLinearPredictionCoef{F,T,N}(α, data, loggain)
end

function MelLinearPredictionCoef(mlpc::MelLinearPredictionCoef)
    α = allpass_alpha(mlpc)
    MelLinearPredictionCoef(α, rawdata(mlpc), mlpc.loggain)
end

function lpc{T<:FloatingPoint,N}(x::AbstractArray{T,N},
                                 order::Int=40,
                                 kargs...)
    raw = _mgcep(x, order, 0.0, -1.0; kargs..., otype=5) # force output type
    r = MelLinearPredictionCoef{Linear,T,N}(0.0, raw, false)
end
