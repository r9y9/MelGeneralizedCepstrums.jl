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

function mc2b{FS,L}(c::MelGeneralizedCepstrum{FS,L})
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    raw = mc2b(rawdata(c), α)
    MelGeneralizedCepstrum{FS,L}(α, γ, raw)
end
