# MelGeneralizedCepstrums

[![Build Status](https://travis-ci.org/r9y9/MelGeneralizedCepstrums.jl.svg?branch=master)](https://travis-ci.org/r9y9/MelGeneralizedCepstrums.jl)

MelGeneralizedCepstrums.jl provides a mel generalized-log cepstrum anlysis for spectrum envelope estimation, which includes Linear Predicition, Mel-Cepstrum analysis and Generalized Cepstrum analysis for Julia. The core is re-writen by pure Julia language based on [Speech Signal Processing Toolkit (SPTK)](http://sp-tk.sourceforge.net/).

## Type design

MelGeneralizedCepstrums.jl provides a general type `MelGeneralizedCepstrum{F,L,T,N}` which is a subtype of `AbstractArray{T,N}`, where `F` represents a type of frequency scale and `L` represents a type of log function. Linear Prediction-based cepstrum, mel-cepstrum, generalized cepstrum and other kind of cepstrum are represented as a special type with specific type pareamters of `MelGeneralizedCepstrum{F,L,T,N}`.

The actual type definition is as follows:

```julia
immutable MelGeneralizedCepstrum{F<:Frequency,L<:Log,T<:FloatingPoint,N} <: AbstractMelGeneralizedCepstrumArray{F,L,T,N}
    α::T
    γ::T
    data::Array{T,N}

    function MelGeneralizedCepstrum(α::T, γ::T, data::Array{T,N})
        abs(α) < 1 || error("|α| < 1 is supported")
        (-1 <= γ <= 0) || error("-1 <= γ <= 0 is supported")
        @assert size(data, 1) > 1
        new(α, γ, data)
    end
end
```

where 

```julia
abstract Frequency
type Mel <: Frequency
end
type Linear <: Frequency
end

abstract Log
type GeneralizedLog <: Log
end
type StandardLog <: Log
end
type AllPoleLog <: Log
end
```

Following the above definition, for example, mel-cepstrum is represented as 

```julia
typealias MelCepstrum{T,N} MelGeneralizedCepstrum{Mel,StandardLog,T,N}
```

and generalized cepstrum is represented as:

```julia
typealias GeneralizedCepstrum{T,N} MelGeneralizedCepstrum{Linear,GeneralizedLog,T,N}
```

For more information about types, please check [common.jl]](src/common.jl). 
