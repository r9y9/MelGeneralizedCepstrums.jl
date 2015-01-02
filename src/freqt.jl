function freqt{T<:FloatingPoint}(c::Vector{T}, order::Int, α::Float64)
    wc = zeros(T, order+1)
    prev = zeros(T, order+1)

    m1 = length(c)-1
    for i=-m1:0
        copy!(prev, wc)
        if order >= 0
            wc[1] = c[-i+1] + α*prev[1]
        end
        if order >= 1
            wc[2] = (1.0-α*α)*prev[1] + α*prev[2]
        end
        for m=3:order+1
            @inbounds wc[m] = prev[m-1] + α*(prev[m] - wc[m-1])
        end
    end

    wc
end

function freqt(c::MelGeneralizedCepstrum, order::Int, α::Float64)
    raw = rawdata(c)
    cc = freqt(raw, order, α)
    MelGeneralizedCepstrum(allpass_alpha(c)+α, glog_gamma(c), cc)
end
