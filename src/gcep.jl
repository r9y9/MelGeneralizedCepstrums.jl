typealias GeneralizedCepstrum MelGeneralizedCepstrum{Linear,GeneralizedLog}

function GeneralizedCepstrum(γ::Float64, data::Vector{Float64})
    MelGeneralizedCepstrum{Linear, GeneralizedLog}(0.0, γ, data)
end

function gcep(x::Vector{Float64}, order::Int, γ::Float64)
    SPTK.gcep(x, order, γ)
end
