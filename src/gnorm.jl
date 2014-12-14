function gnorm!(c::Vector{Float64}, γ::Float64)
    if γ == 0.0
        c[1] = exp(c[1])
        return normalizedc
    end

    gain = 1.0 + γ*c[1]
    c[2:end] = c[2:end]/gain
    c[1] = gain^(1.0/γ)
    
    c
end

function gnorm(c::Vector{Float64}, γ::Float64)
    normalizedc = copy(c)
    gnorm!(normalizedc, γ)
end

function gnorm{FS,L}(c::MelGeneralizedCepstrum{FS,L})
    raw = gnorm(rawdata(c), gamma(c))
    MelGeneralizedCepstrum{FS,L}(alpha(c), gamma(c), raw)
end
