function ignorm(normalizedc::Vector{Float64}, γ::Float64)
    c = Array(Float64, length(normalizedc))

    if γ == 0.0
        copy!(c, normalizedc)
        c[1] = log(normalizedc[1])
        return c
    end

    gain = normalizedc[1]^γ
    c[1:end] = normalizedc[1:end]*gain
    c[1] = (gain - 1.0)/γ
    
    c
end

function ignorm{FS,L}(c::MelGeneralizedCepstrum{FS,L})
    raw = ignorm(rawdata(c), gamma(c))
    MelGeneralizedCepstrum{FS,L}(alpha(c), gamma(c), raw)
end
