module MelGeneralizedCepstrums

# TODO: remove this dependency
import SPTK

export
    # Types
    MelGeneralizedCepstrum,

    MelFrequencyCepstrum,
    LinearFrequencyCepstrum,
    GeneralizedLogCepstrum,
    StandardLogCepstrum,
    AllPoleCepstrum,
    AllZeroCepstrum,
    MelCepstrum,
    GeneralizedCepstrum,

    # Basic property of Mel-generalized cepstrum
    order,
    allpass_alpha,   # all-pass constant (alpha)
    glog_gamma,      # parameter of generalized log function
    powercoef,

    # Feature extraction
    mgcep,
    mcep,
    gcep,

    # Conversions
    mc2b,
    mc2b!,
    mc2e,
    mgc2b,
    mgc2b!,
    b2mc,
    gnorm!,
    gnorm,
    ignorm!,
    ignorm,
    c2ir,
    freqt,
    freqt!,
    gc2gc,
    mgc2mgc

for fname in [
              "common",
              "mgcep",
              "mc2b",
              "mgc2b",
              "b2mc",
              "gnorm",
              "ignorm",
              "c2ir",
              "freqt",
              "gc2gc",
              "mgc2mgc",
              "extend"
    ]
    include(string(fname, ".jl"))
end

end # module
