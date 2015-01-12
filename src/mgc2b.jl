# mgc2b converts mel generalized cesptrum to MGLSADF filter coefficients.
function mgc2b!{T<:FloatingPoint}(mgc::AbstractVector{T}, α::FloatingPoint,
                                  γ::FloatingPoint)
    b = mgc
    mc2b!(b, α)

    # when gamma = 0, mel-generalized cespstrum corresponds to mel cepstrum
    if γ == zero(T)
        return b
    end

    gnorm!(b, γ)

    # scale by gamma
    b[1] = log(b[1])
    for i=2:length(b)
        @inbounds b[i] *= γ
    end

    b
end

function mgc2b{T<:FloatingPoint}(mgc::AbstractVector{T}, α::FloatingPoint,
                                 γ::FloatingPoint)
    b = copy(mgc)
    mgc2b!(b, α, γ)
end

function mgc2b(c::MelGeneralizedCepstrum)
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    raw = mgc2b(rawdata(c), α, γ)
    MelGeneralizedCepstrum(α, γ, raw)
end
