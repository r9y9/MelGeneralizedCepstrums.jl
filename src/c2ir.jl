function c2ir(c::Vector{Float64}, len::Int)
    h = Array(Float64, len)
    m = length(c) - 1

    h[1] = exp(c[1])
    for n=2:len
        d = 0.0
        upperlimit = n
        if n >= m+1
            upperlimit = m
        end
        for k=2:upperlimit+1
            @inbounds d += (k-1) * c[k] * h[n-k+1]
        end
        @inbounds h[n] = d / (n-1)
    end

    h
end

c2ir(c::MelGeneralizedCepstrum, len::Int) = c2ir(rawdata(c), len)
