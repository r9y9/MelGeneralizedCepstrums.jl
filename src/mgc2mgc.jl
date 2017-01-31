function mgc2mgc!(dst_ceps::AbstractVector, dst_α, dst_γ,
                  src_ceps::AbstractVector, src_α, src_γ)
    dst_order = length(dst_ceps) - 1
    α = (dst_α-src_α) / (1.0-dst_α*src_α)

    if α == zero(α)
        dst_ceps = gnorm(src_ceps, src_γ)
        dst_ceps = gc2gc(dst_ceps, src_γ, dst_order, dst_γ)
        ignorm!(dst_ceps, dst_γ)
    else
        freqt!(dst_ceps, src_ceps, α)
        gnorm!(dst_ceps, src_γ)
        dst_ceps = gc2gc(dst_ceps, src_γ, dst_order, dst_γ)
        ignorm!(dst_ceps, dst_γ)
    end

    dst_ceps
end

function mgc2mgc(src_ceps::AbstractVector,
                 src_α=0.0,
                 src_γ=0.0,
                 dst_order=length(src_ceps) - 1,
                 dst_α=0.0,
                 dst_γ=0.0)
    mgc2mgc!(Array{eltype(src_ceps),1}(dst_order+1), dst_α, dst_γ,src_ceps, src_α, src_γ)
end

function mgc2mgc{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T},
                                            dst_order=param_order(paramdef(state)),
                                            dst_α=allpass_alpha(paramdef(state)),
                                            dst_γ=glog_gamma(paramdef(state)))
    assert_not_ready_to_filt(state)

    def = paramdef(state)
    src_α = allpass_alpha(def)
    src_γ = glog_gamma(def)
    normalized = (dst_γ == zero(dst_γ))
    data = rawdata(state)
    newdef = MelGeneralizedCepstrum(dst_order, dst_α, dst_γ)
    SpectralParamState(newdef, mgc2mgc(data, src_α, src_γ, dst_order, dst_α, dst_γ),
                       has_loggain(state), normalized)
end
