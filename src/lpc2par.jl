function lpc2par(state::SpectralParamState{LinearPredictionCoef};
                 kargs...)
    @assert !has_loggain(state)
    lpcdef = paramdef(state)
    data = rawdata(state)
    data = SPTK.lpc2par(data; kargs...)
    newdef = PartialAutoCorrelation(param_order(lpcdef))
    SpectralParamState(newdef, data,
                       has_loggain(state), gain_normalized(state))
end
