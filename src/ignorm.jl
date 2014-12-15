function ignorm!(c::Vector{Float64}, γ::Float64)
    if γ == 0.0
        c[1] = log(c[1])
        return c
    end

    gain = c[1]^γ
    c[1:end] = c[1:end]*gain
    c[1] = (gain - 1.0)/γ
    
    c
end

function ignorm(normalizedc::Vector{Float64}, γ::Float64)
    c = copy(normalizedc)
    ignorm!(c, γ)
end

function ignorm{FS,L}(c::MelGeneralizedCepstrum{FS,L})
    γ = glog_gamma(c)
    raw = ignorm(rawdata(c), γ)
    MelGeneralizedCepstrum{FS,L}(allpass_alpha(c), γ, raw)
end
