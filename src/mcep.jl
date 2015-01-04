# Mel-cepstrum analysis
# re-coded from SPTK

function fill_al!{T<:FloatingPoint}(al::Vector{T}, α::Float64)
    al[1] = one(T)
    for i=2:length(al)
        @inbounds al[i] = -α*al[i-1]
    end
    al
end

function fill_toeplitz!{T}(A::AbstractMatrix{T}, t::AbstractVector{T})
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

function fill_hankel!{T}(A::AbstractMatrix{T}, h::AbstractVector{T})
    n = div(length(h),2) + 1
    for i=1:n, j=1:n
        @inbounds A[i,j] = h[i+j-1]
    end
    A
end


function fill_only_real_part!{T}(y::AbstractVector{Complex{T}},
                                 v::AbstractVector{T})
    for i=1:length(v)
        @inbounds y[i] = Complex(v[i], zero(T))
    end
    y
end

function update_hankel_elements!(he::AbstractVector, c::AbstractVector)
    for j=1:length(he)
        @inbounds he[j] = c[j]
    end
    for j=1:2:length(he)
        @inbounds he[j] -= c[1]
    end
    he
end

function update_toeplitz_elements!(te::AbstractVector, c::AbstractVector)
    for j=1:2:length(te)
        @inbounds te[j] += c[1]
    end
    te
end

function mcep{T<:FloatingPoint}(x::Vector{T}, order::Int, α::T;
                                miniter::Int=2,
                                maxiter::Int=30,
                                threshold::T=0.001,
                                e::T=zero(T))
    const xh = div(length(x),2)

    y = Array(Complex{T}, xh+1)
    c = Array(T, length(x))

    # create FFT plan
    fplan = FFTW.Plan(c, y, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
    iplan = FFTW.Plan(y, c, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)

    # Periodgram
    FFTW.execute(fplan.plan, x, y)
    periodgram = abs2(y)
    logperiodgram = log(periodgram + e)

    # Initial value of cepstrum
    fill_only_real_part!(y, logperiodgram)
    FFTW.execute(iplan.plan, y, c)
    scale!(c, FFTW.normalization(c))
    c[1] /= 2.0
    c[xh+1] /= 2.0

    # Initial value of mel-cesptrum
    mc = freqt(sub(c, 1:xh+1), order, α)
    s = c[1]

    # Allocate memory for solving linear equation (Tm + Hm)d = b
    Tm = Array(T, order+1, order+1)
    Hm = Array(T, order+1, order+1)
    he = Array(T, 2order+1) # elements of hankel matrix
    te = Array(T, order+1)  # elements of toeplitz matrix
    b = Array(T, order+1)   # right side of linear equation

    al = Array(T, order+1)
    fill_al!(al, α)

    # Newton raphson roop
    for i=1:maxiter
        fill!(c, zero(T))
        freqt!(sub(c, 1:xh+1), mc, -α)

        FFTW.execute(fplan.plan, c, y)
        for i=1:length(y)
            @inbounds y[i] = Complex(periodgram[i] / exp(2real(y[i])), zero(T))
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

        copy!(te, 1, c, 1, order+1)

        for j=1:order+1
            @inbounds b[j] = c[j] - al[j]
        end

        update_hankel_elements!(he, c)
        update_toeplitz_elements!(te, c)

        fill_hankel!(Hm, he)
        fill_toeplitz!(Tm, te)

        # solve linear equation and add derivative
        mc += (Tm + Hm) \ b
    end

    mc
end
