# mc2e computes energy from mel-cepstrum.
function mc2e{T<:FLoatingPoint}(mc::AbstractVector{T}, alpha::FloatingPoint,
                                len::Int)
    # back to linear frequency domain
    c = freqt(mc, len-1, -alpha)

    # compute impule response from cepsturm
    ir = c2ir(c, len)

    sumabs2(ir)
end

function mc2e(c::MelGeneralizedCepstrum, len::Int)
    α = allpass_alpha(c)
    mc2e(rawdata(c), α, len)
end
