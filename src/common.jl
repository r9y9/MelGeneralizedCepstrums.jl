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

function rawdata(c::AbstractMelGeneralizedCepstrum)
    error("All concrete types must implement rawdata for accessing raw data.")
end

# order of mel-generalized cepstrum
order(c::AbstractMelGeneralizedCepstrum) = length(rawdata(c))-1

# all-pass constant for mel-cepstrum analysis
alpha(c::AbstractMelGeneralizedCepstrum) = c.α

# paramter of generalized log function
gamma(c::AbstractMelGeneralizedCepstrum) = c.γ

powercoef(c::AbstractMelGeneralizedCepstrum) = first(rawdata(c))
