function gnorm(c::Vector{Float64}, γ::Float64)
    normalizedc = Array(Float64, length(c))

    if γ == 0.0
        copy!(normalizedc, c)
        normalizedc[1] = exp(c[1])
        return normalizedc
    end

    gain = 1.0 + γ*c[1]
    normalizedc[2:end] = c[2:end]/gain
    normalizedc[1] = gain^(1.0/γ)
    
    normalizedc
end

function gnorm{FS,L}(c::MelGeneralizedCepstrum{FS,L})
    raw = gnorm(rawdata(c), gamma(c))
    MelGeneralizedCepstrum{FS,L}(alpha(c), gamma(c), raw)
end
