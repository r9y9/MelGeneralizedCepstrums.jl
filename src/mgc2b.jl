# mgc2b converts mel generalized cesptrum to MGLSADF filter coefficients.
function mgc2b!(mgc::Vector{Float64}, α::Float64, γ::Float64)
    b = mgc
    mc2b!(b, α)

    # when gamma = 0, mel-generalized cespstrum corresponds to mel cepstrum
    if γ == 0.0
        return b
    end

    gnorm!(b, γ)

    # scale by gamma
    b[1] = log(b[1])
    b[2:end] *= γ

    b
end

function mgc2b(mgc::Vector{Float64}, α::Float64, γ::Float64)
    b = copy(mgc)
    mgc2b!(b, α, γ)
end

function mgc2b(c::MelGeneralizedCepstrum)
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    raw = mgc2b(rawdata(c), α, γ)
    MelGeneralizedCepstrum(α, γ, raw)
end
