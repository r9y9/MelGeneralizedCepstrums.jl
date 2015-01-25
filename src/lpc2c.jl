function lpc2c!{T<:FloatingPoint}(c::AbstractVector{T}, a::AbstractVector{T})
    @assert length(a) == length(c)
    m = length(a) - 1
    fill!(c, zero(T))
    c[1] = log(a[1])
    c[2] = -a[2]
    for k=3:m+1
        d = 0.0
        upperlimit = (k >= m+1) ? m+1 : k
        for i=2:upperlimit
            d += (i-1) *  c[i] * a[k-i+1]
        end
        c[k] = -d / (k-1) - a[k]
    end

    c
end

function lpc2c{T<:FloatingPoint}(a::AbstractVector{T})
    c = similar(a)
    lpc2c!(c, a)
end

function lpc2c{F,T,N}(a::MelLinearPredictionCoef{F,T,N})
    isa(a, LinearPredictionCoef) || throw(ArgumentError("unexpected lpc form"))
    α = allpass_alpha(a)
    raw = lpc2c(rawdata(a))
    MelGeneralizedCepstrum{F,StandardLog,T,N}(α, -1.0, raw)
end
