typealias MelCepstrum MelGeneralizedCepstrum{Mel,StandardLog}

function MelCepstrum(α::Float64, data::Vector{Float64})
    MelGeneralizedCepstrum{Mel,StandardLog}(α, 0.0, data)    
end

function mcep(x::Vector{Float64}, order::Int, α::Float64)
    raw = SPTK.mcep(x, order, α)
    MelCepstrum(α, raw)
end
