function par2lpc(state::SpectralParamState{PartialAutoCorrelation};
                 kargs...)
    @assert !has_loggain(state)
    pardef = paramdef(state)
    data = rawdata(state)
    data = SPTK.par2lpc(data; kargs...)
    newdef = LinearPredictionCoef(param_order(pardef))
    SpectralParamState(newdef, data,
                       has_loggain(state), gain_normalized(state))
end
