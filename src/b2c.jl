function b2c!(wc::AbstractVector, c::AbstractVector, α,
              prev=Array(eltype(wc),length(wc)))
    T = eltype(wc)

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

b2c(c::AbstractVector, order, α) = b2c!(Array(eltype(c), order+1), c, α)
