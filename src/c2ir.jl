function c2ir(c::AbstractVector, len=256)
    h = Array{eltype(c),1}(len)
    m = length(c) - 1

    h[1] = exp(c[1])
    for n=2:len
        d = 0.0
        upperlimit = (n >= m+1) ? m+1 : n
        for k=2:upperlimit
            @inbounds d += (k-1) * c[k] * h[n-k+1]
        end
        @inbounds h[n] = d / (n-1)
    end

    h
end
