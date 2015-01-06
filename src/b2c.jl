function b2c!{T<:FloatingPoint}(wc::AbstractVector{T}, c::AbstractVector{T},
                                α::FloatingPoint;
                                prev::Vector{T}=Array(T,length(wc)))
    fill!(wc, zero(T))
    desired_order = length(wc) - 1

    m1 = length(c)-1
    for i=-m1:0
        copy!(prev, wc)
        wc[1] = c[-i+1]

        if desired_order >= 1
            wc[2] = (1.0-α*α)*prev[1] + α*prev[2]
        end
        for m=3:desired_order+1
            @inbounds wc[m] = prev[m-1] + α*(prev[m] - wc[m-1])
        end
    end

    wc
end

function b2c{T<:FloatingPoint}(c::AbstractVector{T}, order::Int,
                               α::FloatingPoint)
    wc = Array(T, order+1)
    b2c!(wc, c, α)
end
