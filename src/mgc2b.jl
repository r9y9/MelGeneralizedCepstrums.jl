# mgc2b converts mel generalized cesptrum to MGLSADF filter coefficients.
# TODO: consider `otype` optional argument?
function mgc2b!(mgc::AbstractVector, α, γ)
    b = mgc
    mc2b!(b, α)

    # when gamma = 0, mel-generalized cespstrum corresponds to mel cepstrum
    γ == zero(γ) && return mgc

    gnorm!(b, γ)

    # scale by gamma
    b[1] = log(b[1])
    for i=2:length(b)
        @inbounds b[i] *= γ
    end

    b
end

mgc2b(mgc::AbstractVector, α=0.35, γ=0.0) = mgc2b!(copy(mgc), α, γ)

function mgc2b!{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T})
    assert_not_ready_to_filt(state)

    @assert has_loggain(state)
    γ = glog_gamma(paramdef(state))
    if γ == zero(γ)
        @assert gain_normalized(state)
    else
        @assert !gain_normalized(state)
    end

    def = paramdef(state)
    α = allpass_alpha(def)
    γ = glog_gamma(def)
    data = rawdata(state)
    state.data = mgc2b!(data, α, γ)
    state.has_loggain = true
    state.gain_normalized = true
    state.ready_to_filt = true
    state
end

# mgc is assumed not to be gamma scaled (TODO: handle this explicitly?)
function mgc2b{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T})
    mgc2b!(copy(state))
end
