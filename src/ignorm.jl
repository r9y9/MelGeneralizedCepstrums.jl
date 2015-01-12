function ignorm!{T<:FloatingPoint}(c::AbstractVector{T},
                                   γ::FloatingPoint)
    if γ == zero(T)
        c[1] = log(c[1])
        return c
    end

    gain = c[1]^γ
    for i=2:length(c)
        @inbounds c[i] *= gain
    end
    c[1] = (gain - one(T))/γ

    c
end

function ignorm{T<:FloatingPoint}(normalizedc::AbstractVector{T},
                                  γ::FloatingPoint)
    c = copy(normalizedc)
    ignorm!(c, γ)
end

function ignorm(c::MelGeneralizedCepstrum)
    γ = glog_gamma(c)
    raw = ignorm(rawdata(c), γ)
    MelGeneralizedCepstrum(allpass_alpha(c), γ, raw)
end
