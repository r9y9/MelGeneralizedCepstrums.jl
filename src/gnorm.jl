
function gnorm!(c::AbstractVector, γ)
    T = eltype(c)
    if γ == zero(T)
        c[1] = exp(c[1])
        return c
    end

    gain = one(T) + γ*c[1]
    for i=2:length(c)
        @inbounds c[i] /= gain
    end
    c[1] = gain^(one(T)/γ)

    c
end

gnorm(c::AbstractVector, γ=0.0) = gnorm!(copy(c), γ)

function gnorm!{T<:GeneralizedLogCepstrum}(state::SpectralParamState{T})
    assert_not_ready_to_filt(state)
    assert_gain_unnormalized(state)

    γ = glog_gamma(paramdef(state))
    state.data = gnorm!(rawdata(state), γ)
    state.gain_normalized = true
    state
end

function gnorm!{T<:Union{AllPoleLogCepstrum,StandardLogCepstrum}}(state::SpectralParamState{T})
    assert_not_ready_to_filt(state)
    assert_gain_unnormalized(state)
    state
end

function gnorm{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T})
    gnorm!(copy(state))
end
