function b2mc(b::Vector{Float64}, α::Float64)
    mc = Array(Float64, length(b))
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
