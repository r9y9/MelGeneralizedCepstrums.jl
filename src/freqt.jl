function freqt!(wc::AbstractVector, c::AbstractVector, α;
                prev=Array(eltype(wc), length(wc)))
    fill!(wc, zero(eltype(wc)))
    dst_order = length(wc) - 1

    m1 = length(c)-1
    for i in -m1:0
        copy!(prev, wc)
        if dst_order >= 0
            @inbounds wc[1] = c[-i+1] + α*prev[1]
        end
        if dst_order >= 1
            wc[2] = (1.0-α*α)*prev[1] + α*prev[2]
        end
        for m=3:dst_order+1
            @inbounds wc[m] = prev[m-1] + α*(prev[m] - wc[m-1])
        end
    end

    wc
end

function freqt(c::AbstractVector, order=25, α=0.35)
    wc = Array(eltype(c), order+1)
    freqt!(wc, c, α)
end

function freqt{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T},
                                          order=param_order(paramdef(state)),
                                          α²=allpass_alpha(paramdef(state)))
    assert_not_ready_to_filt(state)

    def = paramdef(state)
    α¹ = allpass_alpha(def)
    α = (α²-α¹) / (1.0-α²*α¹)
    γ = glog_gamma(def)
    data = rawdata(state)
    newdef = MelGeneralizedCepstrum(order, α², γ)
    SpectralParamState(newdef, freqt(data, order, α),
                       has_loggain(state), gain_normalized(state))
end
