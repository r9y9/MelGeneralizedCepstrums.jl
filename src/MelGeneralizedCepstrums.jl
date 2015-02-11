module MelGeneralizedCepstrums

export
    # Types
    MelGeneralizedCepstrum,
    MelFrequencyCepstrum,
    LinearFrequencyCepstrum,
    GeneralizedLogCepstrum,
    StandardLogCepstrum,
    AllPoleCepstrum,
    LinearCepstrum,
    MelCepstrum,
    GeneralizedCepstrum,
    MelAllPoleCepstrum,

    MelLinearPredictionCoef,
    LinearPredictionCoef,

    MelGeneralizedCepstrumFilterCoef,
    LMADFCoef,
    MLSADFCoef,
    MGLSADFCoef,

    # Basic property of Mel-generalized cepstrum
    order,           # order of cepstrum (not including 0-th coef.)
    allpass_alpha,   # all-pass constant (alpha)
    glog_gamma,      # parameter of generalized log function
    powercoef,       # power coef. (0-th order of mgcep)

    # Feature extraction
    mgcep,           # mel-generalized cepstrum analysis
    mcep,            # mel-cepstrum analysis
    lpc,             # linear prediction (pipe to mgcep for now)
    periodogram2mcep,# periodogram to mel-cepstrum

    mcepalpha,       # determine appropriate all pass constant in mel-cesptum

    # Conversions
    lpc2c,           # LPC to cepstrum
    lpc2c!,          #
    mc2b,            # mel-cepstrum -> mlsadf filter coef.
    mc2b!,           #
    mc2e,            # mel-cesptrum to energy
    mgc2b,           # mel-generalized cepstrum -> mglsadf coef.
    mgc2b!,          #
    mgc2sp,          # mel-generalized cepstrum -> spectrum envelope
    b2mc,            # mlsadf filter coef. -> mel-cepstrum
    b2mc!,           #
    b2c,             # b2c in _mgcep.c
    b2c!,            #
    gnorm,           # gain normalization for mel-generalized cepstrum
    gnorm!,          #
    ignorm,          # inverse gain normalization for mel-generalized cepstrum
    ignorm!,         #
    c2ir,            # cepstrum -> impulse response
    freqt,           # frequency transform
    freqt!,          #
    frqtr,           # frqtr in _mcep.c
    frqtr!,          #
    gc2gc,           # generalized cepstrum -> genralized cepstrum
    gc2gc!,          #
    mgc2mgc          # mel-generalized cepstrum -> mel-generalized cepstrum

for fname in [
              "common",
              "mgcep",
              "mcep",
              "lpc",
              "lpc2c",
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
