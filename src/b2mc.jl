function b2mc!{T<:FloatingPoint}(mc::AbstractVector{T}, α::FloatingPoint)
    m = length(mc)-1
    o = zero(T)

    mc[m] = mc[m]
    d = mc[m]
    for i=m-1:-1:1
        o = mc[i] + α*d
        d = mc[i]
        mc[i] = o
    end

    mc
end

function b2mc{T<:FloatingPoint}(b::AbstractVector{T}, α::FloatingPoint)
    mc = copy(b)
    b2mc!(mc, α)
end

# TODO(ryuichi): add b2mc using MelGeneralizedCepstrum type
