function mgc2sp(mgc::AbstractVector, α::FloatingPoint, γ::FloatingPoint,
                fftlen::Int)
    # transform to cesptrum
    c = mgc2mgc(mgc, α, γ, fftlen>>1, 0.0, 0.0)

    # zero padding
    buf = zeros(eltype(c), fftlen)
    copy!(buf, 1, c, 1, length(c))

    # FFT
    real(rfft(buf))
end

function mgc2sp(mgc::MelGeneralizedCepstrum, fftlen::Int)
    α = allpass_alpha(mgc)
    γ = glog_gamma(mgc)
    mgc2sp(rawdata(mgc), α, γ, fftlen)
end
