function mgc2sp(mgc::AbstractVector, α=0.35, γ=0.0, fftlen=256)
    # transform to cesptrum
    c = mgc2mgc(mgc, α, γ, fftlen>>1, 0.0, 0.0)

    # zero padding
    buf = zeros(eltype(c), fftlen)
    copy!(buf, 1, c, 1, length(c))

    # FFT
    rfft(buf)
end

function mgc2sp{T<:AbstractFloat}(mgc::AbstractMatrix{T}, α=0.0, γ=0.0, fftlen=256)
    sp = Array{Complex{T}}(fftlen>>1 + 1, size(mgc, 2))
    for i in 1:size(mgc, 2)
        @inbounds sp[:,i] = mgc2sp(view(mgc, :, i), α, γ, fftlen)
    end
    sp
end

function mgc2sp{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T}, fftlen=256)
    assert_not_ready_to_filt(state)

    def = paramdef(state)
    α = allpass_alpha(def)
    γ = glog_gamma(def)
    data = rawdata(state)
    mgc2sp(data, α, γ, fftlen)
end
