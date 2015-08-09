# Linear Prediction Coefficients (LPC)

function estimate(def::LinearPredictionCoef, x::AbstractArray;
                  use_mgcep::Bool=false, kargs...)
    order = param_order(def)
    data = if use_mgcep
        _mgcep(x, order, 0.0, -1.0; otype=5, kargs...)
    else
        SPTK.lpc(x, order; kargs...)
    end
    SpectralParamState(def, data, false, true)
end

function lpc(x::AbstractArray, order=25; kargs...)
    def = LinearPredictionCoef(order)
    estimate(def, x; kargs...)
end
