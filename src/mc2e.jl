# mc2e computes energy from mel-cepstrum.
function mc2e{T<:FloatingPoint}(mc::AbstractVector{T}, α::FloatingPoint,
                                len::Int)
    # back to linear frequency domain
    c = freqt(mc, len-1, -α)

    # compute impule response from cepsturm
    ir = c2ir(c, len)

    sumabs2(ir)
end

function mc2e{T<:FloatingPoint}(mc::AbstractMatrix{T}, α::FloatingPoint,
                                len::Int)
    [mc2e(sub(mc, 1:size(mc, 1), i), α, len) for i=1:size(mc, 2)]
end

function mc2e(c::MelGeneralizedCepstrum, len::Int)
    α = allpass_alpha(c)
    mc2e(rawdata(c), α, len)
end
