function gnorm!{T<:FloatingPoint}(c::AbstractVector{T}, γ::FloatingPoint)
    if γ == zero(T)
        c[1] = exp(c[1])
        return c
    end

    gain = one(T) + γ*c[1]
    c[2:end] = c[2:end]/gain
    c[1] = gain^(one(T)/γ)

    c
end

function gnorm{T<:FloatingPoint}(c::AbstractVector{T}, γ::FloatingPoint)
    normalizedc = copy(c)
    gnorm!(normalizedc, γ)
end

function gnorm(c::MelGeneralizedCepstrum)
    γ = glog_gamma(c)
    raw = gnorm(rawdata(c), γ)
    MelGeneralizedCepstrum(allpass_alpha(c), γ, raw)
end
