VERSION >= v"0.4.0-dev+6521" && __precompile__()

module MelGeneralizedCepstrums

using Compat

export
    # Spectral parameter types
    SpectralParam,
    MelGeneralizedCepstrum,
    GeneralizedCepstrum,
    MelCepstrum,
    LinearCepstrum,
    AllPoleCepstrum,
    MelAllPoleCepstrum,
    LinearPredictionCoefVariants,
    LinearPredictionCoef,
    LineSpectralPair,
    PartialAutoCorrelation,

    # State, which keeps actual spectral parameter values
    SpectralParamState, # subtype of AbstractArray

    # typealiases for SpectralParamState
    SpectralParamStateVector,
    SpectralParamStateMatrix,

    # Generic interfaces
    estimate,        # estimate spectral parameter values
    params,          # returns parameters of a spectral parameter
    param_order,     # returns order of a spectral parameter

    # Mel-generalized cepstrum properties
    allpass_alpha,   # all-pass constant (α)
    glog_gamma,      # parameter of generalized log function (γ)

    # Spectral state properties
    paramdef,        # returns parameter definition from state
    rawdata,         # returns rawdata in state
    has_loggain,
    gain_normalized,
    ready_to_filt,

    loggain!,
    unloggain!,

    # Spectral parameter estimation
    gcep,            # generalized cepstrum analysis
    mcep,            # mel-cepstrum analysis
    mgcep,           # mel-generalized cepstrum analysis
    lpc,             # linear prediction (pipe to mgcep for now)
    periodogram2mcep,# periodogram to mel-cepstrum

    mcepalpha,       # determine appropriate all pass constant in mel-cesptum

    # LPC, LSP and PARCOR conversions
    lpc2c,           # LPC to cepstrum
    lpc2lsp,         # LPC -> LSP
    lpc2par,         # LPC -> PARCOR
    par2lpc,         # PARCOR -> LPC

    # Mel-generalized cepstrum conversions
    mgc2b,           # mel-generalized cepstrum -> MGLSADF coef.
    mgc2b!,          #
    mc2b,            # mel-cepstrum -> MLSADF coef.
    mc2b!,           #
    b2mc,            # MLSADF coef. -> mel-cepstrum
    b2mc!,           #
    c2ir,            # cepstrum -> impulse response
    gc2gc,           # generalized cepstrum -> generalized cepstrum
    gnorm,           # gain normalization
    gnorm!,          #
    ignorm,          # inverse gain normalization
    ignorm!,         #
    freqt,           # frequency transform
    mgc2mgc,         # mel-generalized cepstrum -> mel-generalized cepstrum
    mgc2sp           # mel-generalized cepstrum -> spectrum envelope

for fname in [
              "common",
              "gcep",
              "mcep",
              "mgcep",
              "lpc",
              "lpc2c",
              "lpc2lsp",
              "lpc2par",
              "par2lpc",
              "mcepalpha",
              "mgc2b",
              "b2c",
              "mc2e",
              "frqtr",
              "mc2b",
              "b2mc",
              "c2ir",
              "gc2gc",
              "gnorm",
              "ignorm",
              "freqt",
              "mgc2mgc",
              "mgc2sp",
              "extend"
    ]
    include(string(fname, ".jl"))
end

end # module
