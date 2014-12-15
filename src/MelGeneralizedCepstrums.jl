module MelGeneralizedCepstrums

# TODO: remove this dependency
import SPTK

export
  # Types
  AbstractMelGeneralizedCepstrum,
  MelGeneralizedCepstrum,
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
  mgc2b,
  b2mc,
  gnorm!,
  gnorm,
  ignorm!,
  ignorm,
  c2ir,
  freqt,
  gc2gc,
  mgc2mgc

for fname in [
              "common",
              "mgcep",
              "mcep",
              "gcep",
              "mc2b",
              "mgc2b",
              "b2mc",
              "gnorm",
              "ignorm",
              "c2ir",
              "freqt",
              "gc2gc",
              "mgc2mgc"
    ]
    include(string(fname, ".jl"))
end

end # module
