# mc2e computes energy from mel-cepstrum.
# TODO: theoretical background
function mc2e(mc::AbstractVector, α=0.35, len=256)
    # back to linear frequency domain
    c = SPTK.freqt(mc, len-1, -α)

    # compute impule response from cepsturm
    ir = c2ir(c, len)

    sum(abs2, ir)
end

function mc2e(mc::AbstractMatrix, α=0.35, len=256)
    r = Array{eltype(mc)}(size(mc, 2))
    for i in 1:length(r)
        r[i] = mc2e(view(mc, :, i), α, len)
    end
    r
end

function mc2e(state::SpectralParamState, len=256)
    error("not implemented")
end

function mc2e{T<:Union{MelCepstrum,LinearCepstrum}}(state::SpectralParamState{T},
                                                   len=256)
    assert_not_ready_to_filt(state)

    def = paramdef(state)
    α = allpass_alpha(def)
    data = rawdata(state)
    mc2e(data, α, len)
end
