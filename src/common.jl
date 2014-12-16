# Common types and functions

abstract FrequencyScale
type Mel <: FrequencyScale end
type Linear <: FrequencyScale end

abstract LogFunc
type GeneralizedLog <: LogFunc end
type StandardLog <: LogFunc end

# Generic type
abstract AbstractMelGeneralizedCepstrum{FS<:FrequencyScale,L<:LogFunc} <: AbstractVector{Float64}

frequency_scale{FS<:FrequencyScale,L<:LogFunc}(::Type{AbstractMelGeneralizedCepstrum{FS,L}}) = FS
frequency_scale{T<:AbstractMelGeneralizedCepstrum}(::Type{T}) = frequency_scale(super(T))

log_func{FS<:FrequencyScale,L<:LogFunc}(::Type{AbstractMelGeneralizedCepstrum{FS,L}}) = L
log_func{T<:AbstractMelGeneralizedCepstrum}(::Type{T}) = log_func(super(T))

rawdata(c::AbstractMelGeneralizedCepstrum) = error("not implemented")

# order of mel-generalized cepstrum
order(c::AbstractMelGeneralizedCepstrum) = error("not implemented")

# all-pass constant for mel-cepstrum analysis
allpass_alpha(c::AbstractMelGeneralizedCepstrum) = 0.0

# paramter of generalized log function
glog_gamma(c::AbstractMelGeneralizedCepstrum) = 0.0

powercoef(c::AbstractMelGeneralizedCepstrum) = error("not implemented")
