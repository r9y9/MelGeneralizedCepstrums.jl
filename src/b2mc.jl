function b2mc!(mc::AbstractVector, α)
    m = length(mc)
    o = zero(eltype(mc))

    mc[m] = mc[m]
    d = mc[m]
    for i=m-1:-1:1
        o = mc[i] + α*d
        d = mc[i]
        mc[i] = o
    end

    mc
end

b2mc(b::AbstractVector, α=0.35) = b2mc!(copy(b), α)

function b2mc!{T<:MelCepstrum}(state::SpectralParamState{T})
    assert_ready_to_filt(state)

    def = paramdef(state)
    α = allpass_alpha(def)
    data = rawdata(state)
    state.data = b2mc!(data, α)
    state.ready_to_filt = false
    state
end

function b2mc!{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T})
    throw(ArgumentError("""
                        filter coefficients based on mel-generalized cesptrum cannot be converted to mel-cepstrum"""))

end

function b2mc{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T})
    b2mc!(copy(state))
end
