function mc2b!(coef::Vector{Float64}, mc::Vector{Float64}, α::Float64)
    @assert length(coef) == length(mc)
    coef[end] = mc[end]
    for i=length(mc)-1:-1:1
        @inbounds coef[i] = mc[i] - α*coef[i+1]
    end    
    coef
end

function mc2b(mc::Vector{Float64}, α::Float64)
    coef = Array(Float64, length(mc))
    mc2b!(coef, mc, α)
end

function mc2b{FS,L}(c::MelGeneralizedCepstrum{FS,L})
    α = allpass_alpha(c)
    γ = glog_gamma(c)
    raw = mc2b(rawdata(c), α)
    MelGeneralizedCepstrum{FS,L}(α, γ, raw)
end
