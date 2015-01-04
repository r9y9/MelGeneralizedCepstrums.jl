function mc2b!{T<:FloatingPoint}(mc::AbstractVector{T}, α::FloatingPoint)
    for i=length(mc)-1:-1:1
        @inbounds mc[i] = mc[i] - α*mc[i+1]
    end
    mc
end

function mc2b{T<:FloatingPoint}(mc::AbstractVector{T}, α::FloatingPoint)
    b = copy(mc)
    mc2b!(b, α)
end

function mc2b(c::MelGeneralizedCepstrum)
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    raw = mc2b(rawdata(c), α)
    MelGeneralizedCepstrum(α, γ, raw)
end
