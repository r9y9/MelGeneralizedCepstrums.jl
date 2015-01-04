# Npte that wc must be initialized.
function freqt!{T<:FloatingPoint}(wc::AbstractVector{T}, c::AbstractVector{T},
                                  α::Float64;
                                  prev::Vector{T}=Array(T,length(wc)))
    desired_order = length(wc) - 1

    m1 = length(c)-1
    for i=-m1:0
        copy!(prev, wc)
        if desired_order >= 0
            wc[1] = c[-i+1] + α*prev[1]
        end
        if desired_order >= 1
            wc[2] = (1.0-α*α)*prev[1] + α*prev[2]
        end
        for m=3:desired_order+1
            @inbounds wc[m] = prev[m-1] + α*(prev[m] - wc[m-1])
        end
    end

    wc
end

function freqt{T<:FloatingPoint}(c::AbstractVector{T}, order::Int, α::Float64)
    wc = zeros(T, order+1)
    freqt!(wc, c, α)
end

function freqt(c::MelGeneralizedCepstrum, order::Int, α::Float64)
    raw = rawdata(c)
    cc = freqt(raw, order, α)
    MelGeneralizedCepstrum(allpass_alpha(c)+α, glog_gamma(c), cc)
end

# Npte that wc must be initialized.
function frqtr!{T<:FloatingPoint}(wc::AbstractVector{T},
                                  c::AbstractVector{T},
                                  α::Float64;
                                  prev::Vector{T}=Array(T,length(wc)))
    desired_order = length(wc) - 1

    m1 = length(c)-1
    for i=-m1:0
        copy!(prev, wc)
        if desired_order >= 0
            wc[1] = c[-i+1]
        end
        for m=2:desired_order+1
            @inbounds wc[m] = prev[m-1] + α*(prev[m] - wc[m-1])
        end
    end

    wc
end

function frqtr{T<:FloatingPoint}(c::AbstractVector{T}, order::Int, α::Float64)
    wc = zeros(T, order+1)
    frqtr!(wc, c, α)
end

function frqtr(c::MelGeneralizedCepstrum, order::Int, α::Float64)
    raw = rawdata(c)
    cc = frqtr(raw, order, α)
    MelGeneralizedCepstrum(allpass_alpha(c)+α, glog_gamma(c), cc)
end
