function mc2b!(mc::Vector{Float64}, α::Float64)
    for i=length(mc)-1:-1:1
        @inbounds mc[i] = mc[i] - α*mc[i+1]
    end
    mc
end

function mc2b(mc::Vector{Float64}, α::Float64)
    b = copy(mc)
    mc2b!(b, α)
end

function mc2b(c::MelGeneralizedCepstrum)
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    raw = mc2b(rawdata(c), α)
    MelGeneralizedCepstrum(α, γ, raw)
end
