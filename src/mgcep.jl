# Mel generalized cepstrum analysis
# re-coded from SPTK

# Overview:
# Tokuda, Keiichi, et al. "Mel-generalized cepstral analysis-a unified approach
# to speech spectral estimation." ICSLP. 1994.
# Details:
# 徳田恵一, et al. "メル一般化ケプストラム分析による音声のスペクトル推定."
# 電子情報通信学会論文誌 A 75.7 (1992): 1124-1134.

# memo(ryuichi):
# Definition of generalized log function:
# sᵧ(ω) =
# (ωᵞ - 1)/γ, 0 < |γ| ≤ 1,
# log(ω), γ = 0.
#
# Definition of inverse generalized log function:
# sᵧ⁻¹(ω) =
# (1 + γω)¹/ᵞ, 0 < |γ| ≤ 1,
# exp(ω), γ = 0.
#
# All-pass filter model:
# z̃⁻¹ = Ψ(z) = (z⁻¹ - α)/(1 - αz⁻¹)
#
# A spectral model based on a mel generalized cepstrum:
# H(z) = sᵧ⁻¹(∑ₘ cₐ,ᵧ(m)z̃⁻ᵐ) =
# (1 + γ∑ₘ cₐ,ᵧ(m)z̃⁻ᵐ)¹/ᵞ, 0 < |γ| ≤ 1,
# exp (∑ₘ cₐ,ᵧ(m) z̃⁻ᵐ), γ = 0.
#
# Another representation with explicit gain `K`:
# H(z) = sᵧ⁻¹(∑ₘ bₐ,ᵧ(m)Φₘ(z)) = K⋅D(z)
# where
# K = sᵧ⁻¹(bₐ,ᵧ(0)),
# D(z) =
#  (1 + γ∑ₘ bₐ,ᵧ'(m)z̃⁻ᵐ)¹/ᵞ, 0 < |γ| ≤ 1,
#  exp (∑ₘ bₐ,ᵧ'(m) z̃⁻ᵐ), γ = 0,
#
# Gain normalization (e.g. b -> b'):
# bₐ,ᵧ'(m) = bₐ,ᵧ(m) / (1 + γbₐ,ᵧ(0))

function ptrans!{T}(p::AbstractVector{T},  m::Int, α::FloatingPoint)
    d, o = zero(T), zero(T)

    d = p[m+1]
    @inbounds for i=m-1:-1:1
        o = p[i+1] + α * d
        d = p[i+1]
        p[i+1] = o
    end

    o = α * d
    p[1] = (one(T) - α*α) * p[1] + 2o

    p
end

function qtrans!{T}(q::AbstractVector{T},  m::Int, α::FloatingPoint)
    d = q[2]
    @inbounds for i=2:2m
        o = q[i+1] + α*d
        d = q[i+1]
        q[i+1] = o
    end

    q
end

function gain{T}(er::AbstractVector{T}, c::AbstractVector{T}, m::Int,
                 g::FloatingPoint)
    t = zero(T)
    if g != zero(T)
        for i=2:m+1
            @inbounds t += er[i] * c[i]
        end
        return er[1] + g*t
    else
        return er[1]
    end
end

function newton!{T}(c::AbstractVector{T}, # mel-generalized cepstrum stored
                    x::AbstractVector{T}, # modified periodogram
                    order::Int,           # order of cepstrum
                    α::FloatingPoint,     # allpass constant
                    γ::FloatingPoint,     # parameter of generalized log function
                    n::Int,               # the order of recursion
                    iter::Int,            # current iter #
                    y_fft::Array{Complex{T}}, # *length must be equal to length(x)*
                    z_fft::Array{T},      # *length must be equal to length(x)*
                    bplan::FFTW.Plan,     # FFTW.Plan used in backward fft
                    cr::Vector{T} = zeros(T, length(x)),
                    pr::Vector{T} = zeros(T, length(x)),
                    rr::Vector{T} = zeros(T, length(x)),
                    ri::Vector{T} = zeros(T, length(x)),
                    qr::Vector{T} = zeros(T, length(x)),
                    qi::Vector{T} = zeros(T, length(x)),
                    Tm::Matrix{T} = Array(T, order, order),
                    Hm::Matrix{T} = Array(T, order, order)
    )
    @assert length(x) > length(c)
    @assert n < length(x)
    @assert length(y_fft) == length(x)
    @assert length(z_fft) == length(x)

    copy!(cr, 2, c, 2, order)

    if α != zero(T)
        b2c!(sub(cr, 1:n+1), cr[1:order+1], -α)
    end

    y = fft(cr)

    γ⁻¹ = one(T)/γ
    if γ == -one(T)
        copy!(pr, 1, x, 1, length(x))
    elseif γ == zero(T)
        for i=1:length(x)
            @inbounds pr[i] = x[i] / exp(2real(y[i]))
        end
    else
        @inbounds for i=1:length(x)
            tr = one(T) + γ*real(y[i])
            ti = γ*imag(y[i])
            trr, tii = tr*tr, ti*ti
            s = trr + tii
            t = x[i] * s^(-γ⁻¹)
            t /= s
            pr[i] = t
            rr[i] = tr * t
            ri[i] = ti * t
            t /= s
            qr[i] = (trr - tii) * t
            s = tr * ti * t
            qi[i] = s + s
        end
    end

    copy!(y_fft, pr)
    FFTW.execute(bplan.plan, y_fft, pr)
    scale!(pr, FFTW.normalization(pr))

    if α != zero(T)
        b2c!(sub(pr, 1:2order+1), pr[1:n+1], α)
    end

    if γ == zero(T) || γ == -one(T)
        copy!(qr, 1, pr, 1, 2order+1)
        copy!(rr, 1, pr, 1, order+1)
    else
        for i=1:length(qr)
            @inbounds y_fft[i] = Complex(qr[i], qi[i])
        end
        FFTW.execute(bplan.plan, y_fft, qr)
        scale!(qr, FFTW.normalization(qr))
        for i=1:length(rr)
            @inbounds y_fft[i] = Complex(rr[i], ri[i])
        end
        FFTW.execute(bplan.plan, y_fft, rr)
        scale!(rr, FFTW.normalization(rr))

        if α != zero(T)
            b2c!(sub(qr, 1:n+1), qr[1:n+1], α)
            b2c!(sub(rr, 1:order+1), rr[1:n+1], α)
        end
    end

    if α != zero(T)
        ptrans!(pr, order, α)
        qtrans!(qr, order, α)
    end

    ϵ = zero(T)
    if γ != -one(T)
        ϵ = gain(rr, c, order, γ)
        c[1] = √ϵ
    end

    if γ == -one(T)
        fill!(qr, zero(T))
    elseif γ != zero(T)
        for i=3:2order+1
            @inbounds qr[i] *= one(T) + γ
        end
    end

    te = sub(pr, 1:order)
    fill_toeplitz!(Tm, te)
    he = sub(qr, 3:2order+1)
    fill_hankel!(Hm, he)

    # solve linear equation
    b = (Tm + Hm) \ sub(rr, 2:order+1)

    for i=2:order+1
        @inbounds c[i] += b[i-1]
    end

    if γ == -one(T)
        ϵ = gain(rr, c, order, γ)
        c[1] = √ϵ
    end

    log(ϵ)
end


# mgcepnorm! changes form of mel-generalized cepstrums
# input:
# K, bₐ,ᵧ'(1), ..., bₐ,ᵧ'(m)
#
# output(otype):
# 0: cₐ,ᵧ(0), cₐ,ᵧ(1), ..., cₐ,ᵧ(m)
# 1: bₐ,ᵧ(0), bₐ,ᵧ(1), ..., bₐ,ᵧ(m)
# 2: Kₐ, cₐ,ᵧ'(1), ..., cₐ,ᵧ'(m)
# 3: K, bₐ,ᵧ'(1), ..., bₐ,ᵧ'(m)
# 4: Kₐ, γcₐ,ᵧ'(1), ..., γcₐ,ᵧ'(m)
# 5: K, γbₐ,ᵧ'(1), ..., γbₐ,ᵧ'(m)
#
# For simplicity, we represent bₐ,ᵧ as bᵧ.
function mgcepnorm!{T<:FloatingPoint}(bᵧ′::AbstractVector{T},
                                      α::FloatingPoint,
                                      γ::FloatingPoint,
                                      otype::Int)
    0<=otype<=5 || throw(ArgumentError("0 ≤ otype ≤ 5 are supported"))

    mgc = bᵧ′

    if otype == 0 || otype == 1 || otype == 2 || otype == 4
        # K, bᵧ' -> bᵧ
        ignorm!(mgc, γ)
    end

    if otype == 0 || otype == 2 || otype == 4
        # bᵞ -> cᵞ
        b2mc!(mgc, α)
    end

    if otype == 2 || otype == 4
        # cᵧ -> cᵧ'
        gnorm!(mgc, γ)
    end

    if otype == 4 || otype == 5
        # cᵧ' -> γcᵧ'
        for i=2:length(mgc)
            @inbounds mgc[i] *= γ
        end
    end

    mgc
end

function _mgcep{T<:FloatingPoint}(x::AbstractVector{T},          # a *windowed* signal
                                  order::Int=40,                 # order of mgcep
                                  α::FloatingPoint=0.41,         # all-pass constant
                                  γ::FloatingPoint=0.0;          # parameter of generalized log
                                  n::Int=length(x)-1,            # order of recursion
                                  miniter::Int=2,
                                  maxiter::Int=30,
                                  criteria::FloatingPoint=0.001, # stopping criteria
                                  e::T=zero(T),                  # floor of
                                  otype::Int=0,                  # output type
                                  verbose::Bool=false
    )
    @assert n < length(x)

    # Periodogram
    periodogram = abs2(fft(x)) + e

    # Allocate memory
    cr = zeros(T, length(x))
    pr = zeros(T, length(x))
    rr = zeros(T, length(x))
    ri = zeros(T, length(x))
    qr = zeros(T, length(x))
    qi = zeros(T, length(x))
    Tm = Array(T, order, order)
    Hm = Array(T, order, order)

    # FFT workspace
    y = Array(Complex{T}, length(x))
    z = Array(T, length(x))
    bplan = FFTW.Plan(y, z, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)

    bᵧ′ = zeros(T, order+1)
    ϵ⁰ = newton!(bᵧ′, periodogram, order, α, -one(T), n, 1, y, z, bplan,
                 cr, pr, rr, ri, qr, qi, Tm, Hm)

    if γ != -one(T)
        d = Array(T, order+1)
        if α != zero(T)
            ignorm!(bᵧ′, -1.0)
            b2mc!(bᵧ′, α)
            copy!(d, bᵧ′)
            gnorm!(d, -1.0)
        else
            copy!(d, bᵧ′)
        end
        bᵧ′ = gc2gc(d, -1.0, order, γ)

        if α != zero(T)
            ignorm!(bᵧ′, γ)
            mc2b!(bᵧ′, α)
            gnorm!(bᵧ′, γ)
        end
    end

    if γ != -one(T)
        ϵᵗ = ϵ⁰
        for i=1:maxiter
            ϵ = newton!(bᵧ′, periodogram, order, α, γ, n, i, y, z, bplan,
                        cr, pr, rr, ri, qr, qi, Tm, Hm)
            if i >= miniter
                err = abs((ϵᵗ - ϵ)/ϵ)
                verbose && println("nmse: $err")
                if err < criteria
                    break
                end
            end
            ϵᵗ = ϵ
        end
    end

    mgcepnorm!(bᵧ′, α, γ, otype)
end

function mgcep{T<:FloatingPoint,N}(x::AbstractArray{T,N},
                                   order::Int=40,
                                   α::FloatingPoint=0.41,
                                   γ::FloatingPoint=0.0;
                                   kargs...)
    raw = _mgcep(x, order, α, γ; kargs...)
    MelGeneralizedCepstrum(α, γ, raw)
end
