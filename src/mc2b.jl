function mc2b!{T<:FloatingPoint}(mc::AbstractVector{T}, α::FloatingPoint)
    b = mc

    for i=length(mc)-1:-1:1
        @inbounds b[i] = b[i] - α*b[i+1]
    end
    b
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
