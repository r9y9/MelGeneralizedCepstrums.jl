function b2mc!{T<:FloatingPoint}(mc::AbstractVector{T}, α::FloatingPoint)
    m = length(mc)
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

function b2mc{F,L,T,N}(c::MelGeneralizedCepstrumFilterCoef{F,L,T,N})
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    raw = b2mc(rawdata(c), α)
    MelGeneralizedCepstrumFilter{F,L,T,N}(α, γ, raw)
end
