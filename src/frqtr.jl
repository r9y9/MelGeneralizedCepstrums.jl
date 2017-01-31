function frqtr!(wc::AbstractVector, c::AbstractVector, α,
                prev::Vector=Array{eltype(wc),1}(length(wc)))
    fill!(wc, zero(eltype(wc)))
    dst_order = length(wc) - 1

    m1 = length(c)-1
    for i in -m1:0
        copy!(prev, wc)
        if dst_order >= 0
            @inbounds wc[1] = c[-i+1]
        end
        for m=2:dst_order+1
            @inbounds wc[m] = prev[m-1] + α*(prev[m] - wc[m-1])
        end
    end

    wc
end

function frqtr(c::AbstractVector, order=25, α=0.35)
    frqtr!(Array{eltype(c),1}(order+1), c, α)
end
