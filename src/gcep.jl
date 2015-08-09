function estimate(mgc::GeneralizedCepstrum, x::AbstractArray;
                  norm=false,
                  kargs...)
    order = param_order(mgc)
    γ = glog_gamma(mgc)
    data = SPTK.gcep(x, order, γ; norm=norm, kargs...)
    SpectralParamState(mgc, data, true, norm)
end

function gcep(x::AbstractArray, order=25, γ=0.0; kargs...)
    estimate(MelGeneralizedCepstrum(order, 0.0, γ), x; kargs...)
end
