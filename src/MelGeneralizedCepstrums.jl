module MelGeneralizedCepstrums

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
    SpectralParamState,

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

    # Spectral parameter estimation
    gcep,            # generalized cepstrum analysis
    mcep,            # mel-cepstrum analysis
    mgcep,           # mel-generalized cepstrum analysis
    lpc,             # linear prediction (pipe to mgcep for now)
    periodogram2mcep,# periodogram to mel-cepstrum

    mcepalpha,       # determine appropriate all pass constant in mel-cesptum

    # Conversions
    lpc2c,           # LPC to cepstrum
    lpc2lsp,         # LPC -> LSP
    lpc2par,         # LPC -> PARCOR
    mc2b,            # mel-cepstrum -> mlsadf filter coef.
    mc2b!,           #
    mc2e,            # mel-cesptrum to energy (TODO: need more tests)
    mgc2b,           # mel-generalized cepstrum -> mglsadf coef.
    mgc2b!,          #
    mgc2sp,          # mel-generalized cepstrum -> spectrum envelope
    b2mc,            # mlsadf filter coef. -> mel-cepstrum
    b2mc!,           #
    gnorm,           # gain normalization for mel-generalized cepstrum
    gnorm!,          #
    ignorm,          # inverse gain normalization for mel-generalized cepstrum
    ignorm!,         #
    c2ir,            # cepstrum -> impulse response
    freqt,           # frequency transform
    gc2gc,           # generalized cepstrum -> genralized cepstrum
    mgc2mgc          # mel-generalized cepstrum -> mel-generalized cepstrum

for fname in [
              "common",
              "gcep",
              "mcep",
              "mgcep",
              "lpc",
              "lpc2c",
              "lpc2lsp",
              "lpc2par",
              "mcepalpha",
              "mc2b",
              "mc2e",
              "mgc2b",
              "mgc2sp",
              "b2mc",
              "b2c",
              "gnorm",
              "ignorm",
              "c2ir",
              "freqt",
              "frqtr",
              "gc2gc",
              "mgc2mgc",
              "extend"
    ]
    include(string(fname, ".jl"))
end

end # module
