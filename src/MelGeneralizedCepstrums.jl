module MelGeneralizedCepstrums

# TODO: remove this dependency
import SPTK

export
  # Types
  AbstractMelGeneralizedCepstrum,
  MelGeneralizedCepstrum,
  MelCepstrum,
  GeneralizedCepstrum,

  # Frequency scales
  Mel,
  Linear,

  # Log type
  StandardLog,
  GeneralizedLog,

  # Basic property of Mel-generalized cepstrum
  order,
  alpha,
  gamma,
  powercoef,

  # Feature extraction
  mgcep,
  mcep,
  gcep,

  # Conversions
  mc2b,
  mgc2b,
  b2mc,
  gnorm,
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
