function lpc2lsp(state::SpectralParamState{LinearPredictionCoef};
                 loggain=false, kargs...)
    @assert !has_loggain(state)
    lpcdef = paramdef(state)
    data = rawdata(state)
    data = SPTK.lpc2lsp(data; loggain=loggain, kargs...)
    newdef = LineSpectralPair(param_order(lpcdef))
    SpectralParamState(newdef, data,
                       has_loggain(state), gain_normalized(state))
end
