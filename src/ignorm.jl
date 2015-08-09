function ignorm!(c::AbstractVector, γ)
    T = eltype(c)
    if γ == zero(T)
        c[1] = log(c[1])
        return c
    end

    gain = c[1]^γ
    for i=2:length(c)
        @inbounds c[i] *= gain
    end
    c[1] = (gain - one(T))/γ

    c
end

ignorm(normalizedc::AbstractVector, γ=0.0) = ignorm!(copy(normalizedc), γ)

function ignorm!{T<:GeneralizedLogCepstrum}(state::SpectralParamState{T})
    assert_not_ready_to_filt(state)
    assert_gain_normalized(state)

    γ = glog_gamma(paramdef(state))
    state.data = ignorm!(rawdata(state), γ)
    state.gain_normalized = false
    state
end

function ignorm!{T<:Union{AllPoleLogCepstrum,StandardLogCepstrum}}(state::SpectralParamState{T})
    assert_not_ready_to_filt(state)
    assert_gain_normalized(state)
    state
end

function ignorm{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T})
    ignorm!(copy(state))
end
