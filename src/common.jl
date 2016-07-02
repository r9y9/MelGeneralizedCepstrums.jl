import Base: eltype, size, length, getindex, setindex!, similar, copy

import SPTK

### Generic interface ###

abstract SpectralParam

function estimate(s::SpectralParam, x::AbstractArray)
    error("should provide an estimator")
end

function params(s::SpectralParam)
    error("should provide get access to parameters")
end

# subtypes of SpectralParam should have field `order`
param_order(s::SpectralParam) = s.order

### Types that characterize spectral parameters ###

abstract FrequencyForm
type MelFrequency <: FrequencyForm
end
type LinearFrequency <: FrequencyForm
end

abstract LogForm
type GeneralizedLog <: LogForm
end
type StandardLog <: LogForm
end
type AllPoleLog <: LogForm
end

### Definitions of spectral envelope parameters ###

immutable MelGeneralizedCepstrum{F<:FrequencyForm,L<:LogForm} <: SpectralParam
    order::Int
    α::Real
    γ::Real
    function MelGeneralizedCepstrum(order::Int, α::Real, γ::Real)
        order > 0 || throw(ArgumentError("order should be larger than 0"))
        abs(α) < 1 || throw(ArgumentError("|α| < 1 is supported"))
        (-1 <= γ <= 0) || throw(ArgumentError("-1 <= γ <= 0 is supported"))
        new(order, α, γ)
    end
end

typealias MelFrequencyCepstrum{L} MelGeneralizedCepstrum{MelFrequency,L}
typealias LinearFrequencyCepstrum{L} MelGeneralizedCepstrum{LinearFrequency,L}
typealias GeneralizedLogCepstrum{F} MelGeneralizedCepstrum{F,GeneralizedLog}
typealias StandardLogCepstrum{F} MelGeneralizedCepstrum{F,StandardLog}
typealias AllPoleLogCepstrum{F} MelGeneralizedCepstrum{F,AllPoleLog}

typealias GeneralizedCepstrum MelGeneralizedCepstrum{LinearFrequency,GeneralizedLog}
typealias MelCepstrum MelGeneralizedCepstrum{MelFrequency,StandardLog}
typealias LinearCepstrum MelGeneralizedCepstrum{LinearFrequency,StandardLog}
typealias AllPoleCepstrum MelGeneralizedCepstrum{LinearFrequency,AllPoleLog}
typealias MelAllPoleCepstrum MelGeneralizedCepstrum{MelFrequency,AllPoleLog}

freq_form{F<:FrequencyForm,L<:LogForm}(::Type{MelGeneralizedCepstrum{F,L}}) = F
freq_form{T<:MelGeneralizedCepstrum}(::Type{T}) = freq_form(super(T))

log_form{F<:FrequencyForm,L<:LogForm}(::Type{MelGeneralizedCepstrum{F,L}}) = L
log_form{T<:MelGeneralizedCepstrum}(::Type{T}) = freq_form(super(T))

# Generic fallback constructor, which determines type parmaters from α and γ
function MelGeneralizedCepstrum(order::Int, α::Real, γ::Real)
    F = (α == zero(α)) ? LinearFrequency : MelFrequency

    L = GeneralizedLog
    if γ == zero(γ)
        L = StandardLog
    elseif γ == -one(γ)
        L = AllPoleLog
    end

    MelGeneralizedCepstrum{F,L}(order, α, γ)
end

function MelGeneralizedCepstrum(mgc::MelGeneralizedCepstrum)
    order = param_order(mgc)
    α = allpass_alpha(mgc)
    γ = glog_gamma(mgc)
    MelGeneralizedCepstrum(order, α, γ)
end

@compat function (::Type{GeneralizedCepstrum})(order::Int, γ::Real)
    MelGeneralizedCepstrum{LinearFrequency,GeneralizedLog}(order, 0.0, γ)
end

@compat function (::Type{MelCepstrum})(order::Int, α::Real)
    MelGeneralizedCepstrum{MelFrequency,StandardLog}(order, α, 0.0)
end

@compat function (::Type{LinearCepstrum})(order::Int)
    MelGeneralizedCepstrum{LinearFrequency,StandardLog}(order, 0.0, 0.0)
end

@compat function (::Type{AllPoleCepstrum})(order::Int)
    MelGeneralizedCepstrum{LinearFrequency,AllPoleLog}(order, 0.0, -1.0)
end

@compat function (::Type{MelAllPoleCepstrum})(order::Int, α::Real)
    MelGeneralizedCepstrum{MelFrequency,AllPoleLog}(order, α, -1.0)
end

params(c::MelGeneralizedCepstrum) = (c.order, c.α, c.γ)
params(c::GeneralizedCepstrum) = (c.order, c.γ)
params(c::MelCepstrum) = (c.order, c.α)
params(c::LinearCepstrum) = (c.order,)
params(c::AllPoleCepstrum) = (c.order,)
params(c::MelAllPoleCepstrum) = (c.order, c.α)
allpass_alpha(c::MelGeneralizedCepstrum) = c.α
glog_gamma(c::MelGeneralizedCepstrum) = c.γ

### Definitions of LPC variants ###

abstract LinearPredictionCoefVariants <: SpectralParam

# LPC, LSP and PARCOR
immutable LinearPredictionCoef <: LinearPredictionCoefVariants
    order
end

immutable LineSpectralPair <: LinearPredictionCoefVariants
    order
end

immutable PartialAutoCorrelation <: LinearPredictionCoefVariants
    order
end

params(s::LinearPredictionCoefVariants) = (s.order,)

### State, which keeps actual computation results in arrays ###

abstract AbstractParamState{T,N} <: AbstractArray{T,N}

function paramdef(s::AbstractParamState)
    error("should provide get access for parameter definition")
end

function rawdata(s::AbstractParamState)
    error("should provide get access for raw data")
end

type SpectralParamState{S<:SpectralParam,T,N} <: AbstractParamState{T,N}
    def::S
    data::Array{T,N}

    has_loggain::Bool
    gain_normalized::Bool
    ready_to_filt::Bool  # state is explicitly converted to filter coef or not.

    function SpectralParamState(s::S, data::Array{T,N},
                                has_loggain::Bool,
                                gain_normalized::Bool,
                                ready_to_filt::Bool)
        new(s, data, has_loggain, gain_normalized, ready_to_filt)
    end
end

typealias SpectralParamStateVector{S,T} SpectralParamState{S,T,1}
typealias SpectralParamStateMatrix{S,T} SpectralParamState{S,T,2}

function SpectralParamState{S<:SpectralParam,T,N}(s::S, data::Array{T,N},
                                                  has_loggain::Bool,
                                                  gain_normalized::Bool;
                                                  ready_to_filt::Bool=false)
    SpectralParamState{S,T,N}(s, data, has_loggain, gain_normalized,
                              ready_to_filt)
end

eltype(s::SpectralParamState) = eltype(s.data)
size(s::SpectralParamState) = size(s.data)
length(s::SpectralParamState) = length(s.data)
getindex(s::SpectralParamState, i::Real) = getindex(s.data, i)
getindex(s::SpectralParamState, i::Real...) = getindex(s.data, i...)

# colon indexing like s[:,1]
function getindex(s::SpectralParamState, i::Colon, j...)
    SpectralParamState(paramdef(s), getindex(s.data, i, j...),
                       has_loggain(s), gain_normalized(s),
                       ready_to_filt=ready_to_filt(s))
end

setindex!(s::SpectralParamState, i::Real) = setindex!(s.data, i)
setindex!(s::SpectralParamState, i::Real...) = setindex!(s.data, i...)
setindex!(s::SpectralParamState, i::Colon, j...) = setindex!(s.data, i, j...)

paramdef(s::SpectralParamState) = s.def
rawdata(s::SpectralParamState) = s.data
has_loggain(s::SpectralParamState) = s.has_loggain
gain_normalized(s::SpectralParamState) = s.gain_normalized
ready_to_filt(s::SpectralParamState) = s.ready_to_filt

function loggain!(s::SpectralParamStateVector)
    has_loggain(s) && return s
    data = rawdata(s)
    data[1] = log(data[1])
    s.has_loggain = true
    s
end

function loggain!(s::SpectralParamStateMatrix)
    has_loggain(s) && return s
    data = rawdata(s)
    for i in 1:size(s, 2)
        data[1,i] = log(data[1,i])
    end
    s.has_loggain = true
    s
end

function unloggain!(s::SpectralParamStateVector)
    !has_loggain(s) && return s
    data = rawdata(s)
    data[1] = exp(data[1])
    s.has_loggain = false
    s
end

function unloggain!(s::SpectralParamStateMatrix)
    !has_loggain(s) && return s
    data = rawdata(s)
    for i in 1:size(s, 2)
        data[1,i] = exp(data[1,i])
    end
    s.has_loggain = false
    s
end

function similar(s::SpectralParamState)
    def = paramdef(s)
    data = rawdata(s)
    SpectralParamState(def, similar(data), has_loggain(s), gain_normalized(s),
                       ready_to_filt=ready_to_filt(s))
end

function copy(s::SpectralParamState)
    def = paramdef(s)
    data = rawdata(s)
    SpectralParamState(def, copy(data), has_loggain(s), gain_normalized(s),
                       ready_to_filt=ready_to_filt(s))
end

# Force re-type mel-generalized cepstrums by its parameter values
# needed?
function retype{T<:MelGeneralizedCepstrum}(s::SpectralParamState{T})
    def = paramdef(s)
    newdef = MelGeneralizedCepstrum(def)
    SpectralParamState(newdef, rawdata(s), has_loggain(s), gain_normalized(s),
                       ready_to_filt=ready_to_filt(s))
end


### Asserts ###

function assert_ready_to_filt(state::SpectralParamState)
    ready_to_filt(state) && return nothing
    throw(ArgumentError("""ready_to_filt = $(ready_to_filt(state))
                        filter coefficient input is expected"""))
end

function assert_not_ready_to_filt(state::SpectralParamState)
    !ready_to_filt(state) && return nothing
    throw(ArgumentError("""ready_to_filt = $(ready_to_filt(state))
                        filter coefficient input (can be obtain by `mgc2b` or `mc2b`) is not allowed"""))
end

function assert_gain_normalized(state::SpectralParamState)
    gain_normalized(state) && return nothing
    throw(ArgumentError("""gain_normalized: $(gain_normalized(state))
                        gain is assumed to be normalized"""))
end

function assert_gain_unnormalized(state::SpectralParamState)
    gain_normalized(state) || return nothing
    throw(ArgumentError("""gain_normalized: $(gain_normalized(state))
                        gain is assumed to be un-normalized"""))
end
