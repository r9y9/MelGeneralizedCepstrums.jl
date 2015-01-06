function frqtr!{T<:FloatingPoint}(wc::AbstractVector{T},
                                  c::AbstractVector{T},
                                  α::FloatingPoint;
                                  prev::Vector{T}=Array(T,length(wc)))
    fill!(wc, zero(T))
    desired_order = length(wc) - 1

    m1 = length(c)-1
    for i=-m1:0
        copy!(prev, wc)
        if desired_order >= 0
            @inbounds wc[1] = c[-i+1]
        end
        for m=2:desired_order+1
            @inbounds wc[m] = prev[m-1] + α*(prev[m] - wc[m-1])
        end
    end

    wc
end

function frqtr{T<:FloatingPoint}(c::AbstractVector{T}, order::Int,
                                 α::FloatingPoint)
    wc = Array(T, order+1)
    frqtr!(wc, c, α)
end

function frqtr(c::MelGeneralizedCepstrum, order::Int, α::FloatingPoint)
    raw = rawdata(c)
    cc = frqtr(raw, order, α)
    MelGeneralizedCepstrum(allpass_alpha(c)+α, glog_gamma(c), cc)
end
