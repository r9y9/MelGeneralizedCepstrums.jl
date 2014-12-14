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
    raw = ignorm(rawdata(c), gamma(c))
    MelGeneralizedCepstrum{FS,L}(alpha(c), gamma(c), raw)
end
