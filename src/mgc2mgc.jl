function mgc2mgc(c1::Vector{Float64}, α¹::Float64, γ¹::Float64,
                 m2::Int, α²::Float64, γ²::Float64)
    c2 = zeros(Float64, m2+1)

    α = (α²-α¹) / (1.0-α²*α¹)

    if α == 0.0
        c2 = gnorm(c1, γ¹)
        c2 = gc2gc(c2, γ¹, m2, γ²)
    else
        c2 = freqt(c1, m2, α)
        c2 = gnorm(c2, γ¹)
        c2 = gc2gc(c2, γ¹, m2, γ²)
        c2 = ignorm(c2, γ²)
    end

    c2
end

function mgc2mgc(c::MelGeneralizedCepstrum, m2::Int, α²::Float64,
                 γ²::Float64)
    α¹ = allpass_alpha(c)
    γ¹ = glog_gamma(c)
    raw = mgc2mgc(rawdata(c), α¹, γ¹, m2, α², γ²)
    MelGeneralizedCepstrum(α², γ², raw)
end
