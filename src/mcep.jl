function fill_al!{T<:FloatingPoint}(al::Vector{T}, α::Float64)
    al[1] = one(T)
    for i=2:length(al)
        @inbounds al[i] = -α*al[i-1]
    end
    al
end

function fill_toeplitz!{T}(A::Matrix{T}, t::Vector{T})
    n = length(t)
    for i=1:n, j=1:n
        if i-j+1 >= 1
            @inbounds A[i,j] = t[i-j+1]
        else
            @inbounds A[i,j] = t[j-i+1]
        end
    end
    A
end

function fill_hankel!{T}(A::Matrix{T}, h::Vector{T})
    n = div(length(h),2) + 1
    for i=1:n, j=1:n
        @inbounds A[i,j] = h[i+j-1]
    end
    A
end

function fill_real_part!{T}(y::AbstractVector{Complex{T}}, v::AbstractVector{T})
    for i=1:length(v)
        @inbounds y[i] = Complex(v[i], 0.0)
    end
end

function mcep2{T<:FloatingPoint}(x::Vector{T}, order::Int, α::Float64;
                                 miniter::Int=2,
                                 maxiter::Int=30,
                                 threshold::Float64=0.001,
                                 e::Float64=0.0)
    const xh = div(length(x),2)

    y = Array(Complex{T}, xh+1)
    c = Array(T, length(x))

    # create FFT plan
    fplan = FFTW.Plan(x, y, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
    iplan = FFTW.Plan(y, c, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)

    # Periodgram
    FFTW.execute(fplan.plan, x, y)
    periodgram = abs2(y)
    logperiodgram = log(periodgram + e)

    # Initial value of cepstrum
    fill_real_part!(y, logperiodgram)
    FFTW.execute(iplan.plan, y, c)
    scale!(c, FFTW.normalization(c))
    c[1] /= 2.0
    c[xh+1] /= 2.0

    # Initial value of mel-cesptrum
    mc = freqt(sub(c, 1:xh+1), order, α)
    s = c[1]

    # solving linear equation (Tm + Hm)d = b to compute derivative
    Tm = Array(T, order+1, order+1)
    Hm = Array(T, order+1, order+1)
    b = Array(T, order+1)  # elements of toeplitz matrix
    h = Array(T, 2order+1) # elements of hankel matrix

    al = Array(T, order+1)
    fill_al!(al, α)

    # Newton raphson roop
    for i=1:maxiter
        fill!(c, zero(T))
        freqt!(sub(c, 1:xh+1), mc, -α)

        FFTW.execute(fplan.plan, c, y)
        for i=1:length(y)
            @inbounds y[i] = Complex(periodgram[i] / exp(2y[i].re), 0.0)
        end
        FFTW.execute(iplan.plan, y, c)
        scale!(c, FFTW.normalization(c))

        frqtr!(sub(c, 1:2order+1), c[1:xh+1], α)

        # check convergence
        t = c[1]
        if i >= miniter
            if abs((t-s)/t) < threshold
                break
            end
            s = t
        end

        for j=1:order+1
            @inbounds b[j] = c[j] - al[j]
        end
        for j=1:2order+1
            @inbounds h[j] = c[j]
        end
        for j=1:2:2order+1
            @inbounds h[j] -= c[1]
        end
        for j=3:2:order+1
            @inbounds c[j] += c[1]
        end
        c[1] += c[1]

        fill_hankel!(Hm, h)
        fill_toeplitz!(Tm, c[1:order+1])

        # solve!
        mc += (Tm + Hm) \ b
    end

    mc
end
