# Mel generalized cepstrum analysis
# re-coded from SPTK

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
                    iter::Int             # current iter #
    )
    @assert length(x) > length(c)
    @assert n < length(x)

    # Allocate memory
    cr = zeros(T, length(x))
    pr = zeros(T, length(x))
    rr = zeros(T, length(x))
    ri = zeros(T, length(x))
    qr = zeros(T, length(x))
    qi = zeros(T, length(x))
    copy!(cr, 2, c, 2, order)

    if α != zero(T)
        b2c!(sub(cr, 1:n+1), cr[1:order+1], -α)
    end

    y = fft(cr)

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
            t = x[i] * s^(-one(T)/γ)
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

    pr = real(ifft(pr))

    if α != zero(T)
        b2c!(sub(pr, 1:2order+1), pr[1:n+1], α)
    end

    if γ == zero(T) || γ == -one(T)
        copy!(qr, 1, pr, 1, 2order+1)
        copy!(rr, 1, pr, 1, order+1)
    else
        q1 = ifft(complex(qr, qi))
        qr, qi = real(q1), imag(qi)
        q2 = ifft(complex(rr, ri))
        rr, ri = real(q2), imag(q2)

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
        c[1] = sqrt(ϵ)
    end

    if γ == -one(T)
        fill!(qr, zero(T))
    elseif γ != zero(T)
        for i=3:2order+1
            @inbounds qr[i] *= one(T) + γ
        end
    end

    Tm = Array(T, order, order)
    Hm = Array(T, order, order)
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
        c[1] = sqrt(ϵ)
    end

    log(ϵ)
end

function mgcepnorm!{T<:FloatingPoint}(mgc::AbstractVector{T},
                                      α::FloatingPoint,
                                      γ::FloatingPoint,
                                      otype::Int)
    order = length(mgc)-1

    if otype == 0 || otype == 1 || otype == 2 || otype == 4
        ignorm!(mgc, γ)
    end

    if otype == 0 || otype == 2 || otype == 4
        b2mc!(mgc, α)
    end

    if otype == 2 || otype == 4
        gnorm!(mgc, γ)
    end

    if otype == 4 || otype == 5
        mgc[2:end] *= γ
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

    b = zeros(order+1)
    ϵ⁰ = newton!(b, periodogram, order, α, -one(T), n, 1)

    if γ != -one(T)
        d = Array(T, order+1)
        if α != zero(T)
            ignorm!(b, -1.0)
            b2mc!(b, α)
            copy!(d, b)
            gnorm!(d, -1.0)
        else
            copy!(d, b)
        end
        b = gc2gc(d, -1.0, order, γ)

        if α != zero(T)
            ignorm!(b, γ)
            mc2b!(b, α)
            gnorm!(b, γ)
        end
    end

    if γ != -one(T)
        ϵᵗ = ϵ⁰
        for i=1:maxiter
            ϵ = newton!(b, periodogram, order, α, γ, n, i)
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

    mgcepnorm!(b, α, γ, otype)
end

function mgcep{T<:FloatingPoint,N}(x::AbstractArray{T,N},
                                   order::Int=40,
                                   α::FloatingPoint=0.41,
                                   γ::FloatingPoint=0.0;
                                   kargs...)
    raw = _mgcep(x, order, α, γ; kargs...)
    MelGeneralizedCepstrum(α, γ, raw)
end
