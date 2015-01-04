function gc2gc{T<:FloatingPoint}(c1::AbstractVector{T}, γ¹::FloatingPoint,
                                 m2::Int, γ²::FloatingPoint)
    m1 = length(c1)
    c2 = zeros(T, m2+1)
    c2[1] = c1[1]

    for m=2:m2+1
        ss1, ss2 = 0.0, 0.0
        min = m1

        if m1 >= m-1
            min = m-1
        end

        for k=2:min
            @inbounds cc = c1[k] * c2[m-k+1]
            ss2 += (k-1) * cc
            ss1 += (m-k) * cc
        end

        if m <= m1
            @inbounds c2[m] = c1[m] + (γ²*ss2 - γ¹*ss1)/(m-1)
        else
            c2[m] = (γ²*ss2 - γ¹*ss1)/(m-1)
        end
    end

    c2
end

function gc2gc(c::MelGeneralizedCepstrum, m2::Int, γ²::FloatingPoint)
    γ¹ = glog_gamma(c)
    raw = gc2gc(rawdata(c), γ¹, m2, γ²)
    MelGeneralizedCepstrum(allpass_alpha(c), γ², raw)
end
