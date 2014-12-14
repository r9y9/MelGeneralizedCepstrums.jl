function mc2b!(coef::Vector{Float64}, mc::Vector{Float64}, α::Float64)
    @assert length(coef) == length(mc)
    coef[end] = mc[end]
    for i=length(mc)-1:-1:1
        coef[i] = mc[i] - α*coef[i+1]
    end    
    coef
end

function mc2b(mc::Vector{Float64}, α::Float64)
    coef = Array(Float64, length(mc))
    mc2b!(coef, mc, α)
end

function mc2b{FS,L}(c::MelGeneralizedCepstrum{FS,L})
    raw = mc2b(rawdata(c), alpha(c))
    MelGeneralizedCepstrum{FS,L}(alpha(c), gamma(c), raw)
end
