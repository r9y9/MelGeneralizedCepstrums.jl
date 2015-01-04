function b2mc{T<:FloatingPoint}(b::AbstractVector{T}, α::FloatingPoint)
    mc = Array(T, length(b))
    m = length(b)

    mc[m] = b[m]

    for i=m-1:-1:1
        @inbounds mc[i] = b[i] + α*b[i+1]
    end

    mc
end

# TODO(ryuichi): add b2mc using MelGeneralizedCepstrum type
