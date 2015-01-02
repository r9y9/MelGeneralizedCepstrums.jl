function b2mc{T<:FloatingPoint}(b::Vector{T}, α::Float64)
    mc = Array(T, length(b))
    m = length(b)

    mc[m] = b[m]
    d = mc[m]

    for i=m-1:-1:1
        o = b[i] + α*d
        d = b[i]
        mc[i] = o
    end

    mc
end

# TODO(ryuichi): add b2mc using MelGeneralizedCepstrum type
