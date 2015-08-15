function gc2gc!(dst_ceps::AbstractVector, dst_γ,
                src_ceps::AbstractVector, src_γ)
    fill!(dst_ceps, zero(eltype(dst_ceps)))
    dst_order = length(dst_ceps)-1
    m1 = length(src_ceps)-1
    dst_ceps[1] = src_ceps[1]

    for m=2:dst_order+1
        ss1, ss2 = 0.0, 0.0
        min = (m1 < m-1) ? m1 : m - 2

        for k=2:min+1
            @inbounds cc = src_ceps[k] * dst_ceps[m-k+1]
            ss2 += (k-1) * cc
            ss1 += (m-k) * cc
        end

        if m <= m1+1
            @inbounds dst_ceps[m] = src_ceps[m] + (dst_γ*ss2 - src_γ*ss1)/(m-1)
        else
            @inbounds dst_ceps[m] = (dst_γ*ss2 - src_γ*ss1)/(m-1)
        end
    end

    dst_ceps
end

function gc2gc(src_ceps::AbstractVector, src_γ=0.0,
               dst_order=length(src_ceps)-1, dst_γ=0.0)
    gc2gc!(Array(eltype(src_ceps), dst_order+1), dst_γ, src_ceps, src_γ)
end

function gc2gc{T<:MelGeneralizedCepstrum}(state::SpectralParamState{T},
                                          dst_order=param_order(paramdef(state)),
                                          dst_γ=glog_gamma(paramdef(state)))
    assert_not_ready_to_filt(state)

    def = paramdef(state)
    order = param_order(def)
    src_α = allpass_alpha(def)
    src_γ = glog_gamma(def)
    data = rawdata(state)
    normalized = (dst_γ == zero(dst_γ))
    newdef = MelGeneralizedCepstrum(order, src_α, dst_γ)
    SpectralParamState(newdef, gc2gc(data, src_γ, dst_order, dst_γ),
                       has_loggain(state), normalized)
end
