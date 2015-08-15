function mc2b!(mc::AbstractVector, α)
    b = mc

    for i=length(mc)-1:-1:1
        @inbounds b[i] = b[i] - α*b[i+1]
    end
    b
end

mc2b(mc::AbstractVector, α) = mc2b!(copy(mc), α)

function mc2b!{T<:MelCepstrum}(state::SpectralParamState{T})
    assert_not_ready_to_filt(state)

    α = allpass_alpha(paramdef(state))
    state.ready_to_filt = true
    α == zero(α) && return state

    data = rawdata(state)
    state.data = mc2b!(data, α)
    state
end

function mc2b!{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T})
    throw(ArgumentError("""
                        use `mgc2b` instead of `mc2b` for mel-generalized cepstrums"""))
end

function mc2b{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T})
    mc2b!(copy(state))
end
