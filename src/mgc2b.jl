# mgc2b converts mel generalized cesptrum to MGLSADF filter coefficients.
function mgc2b(mgc::Vector{Float64}, α::Float64, γ::Float64)
    b = mc2b(mgc, α)

    # when gamma = 0, mel-generalized cespstrum corresponds to mel cepstrum
    if γ == 0.0
        return b
    end

    b = gnorm(b, γ) # TODO replace with inplace version
    
    # scale by gamma
    b[1] = log(b[1])
    b[2:end] *= γ
    
    b
end
