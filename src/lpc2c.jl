function lpc2c!(c::AbstractVector, a::AbstractVector, loggain::Bool=false)
    @assert length(a) == length(c)
    m = length(a) - 1
    fill!(c, zero(eltype(c)))
    if loggain
        c[1] = a[1]
    else
        c[1] = log(a[1])
    end
    c[2] = -a[2]
    for k=3:m+1
        d = 0.0
        upperlimit = (k >= m+1) ? m+1 : k
        for i=2:upperlimit
            d += (i-1) *  c[i] * a[k-i+1]
        end
        c[k] = -d / (k-1) - a[k]
    end

    c
end

lpc2c(a::AbstractVector, loggain::Bool=false) = lpc2c!(similar(a), a, loggain)

function lpc2c(state::SpectralParamState{LinearPredictionCoef})
    assert_not_ready_to_filt(state)

    def = paramdef(state)
    newdef = LinearCepstrum(param_order(def))
    data = rawdata(state)
    SpectralParamState(newdef, lpc2c(data, has_loggain(state)),
                       true, gain_normalized(state))
end
