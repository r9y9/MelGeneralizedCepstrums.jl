# mc2e computes energy from mel-cepstrum.
function mc2e(mc::Vector{Float64}, alpha::Float64, len::Int)
    # back to linear frequency domain
    c = freqt(mc, len-1, -alpha)

    # compute impule response from cepsturm
    ir = c2ir(c, len)

    sumabs2(ir)
end
