## Extend functions for matrix input (col-wise) ##

for f in [
          :lpc2c!,
          :mc2b!,
          :mgc2b!,
          :b2mc!,
          :b2c!,
          :gnorm!,
          :ignorm!,
          :freqt!,
          :frqtr!,
          :gc2gc!
          ]
    @eval begin
        function ($f){T<:FloatingPoint}(x::AbstractMatrix{T}, args...; kargs...)
            for i = 1:size(x, 2)
                @inbounds x[:, i] = $f(x[:, i], args...; kargs...)
            end
            x
        end
        function ($f){T<:FloatingPoint}(x::Matrix{T}, args...; kargs...)
            for i = 1:size(x, 2)
                @inbounds $f(sub(x, 1:size(x, 1), i), args...; kargs...)
            end
            x
        end
    end
end

for f in [
          :_mcep,
          :_mgcep,
          :lpc2c!,
          :mc2b,
          :mgc2b,
          :mgc2sp,
          :gnorm,
          :ignorm,
          :c2ir,
          :freqt,
          :gc2gc,
          :mgc2mgc
          ]
    @eval begin
        function ($f){T<:FloatingPoint}(x::AbstractMatrix{T}, args...; kargs...)
            r = $f(x[:, 1], args...; kargs...)
            ret = Array(eltype(x), size(r, 1), size(x, 2))
            for i = 1:length(r)
                @inbounds ret[i, 1] = r[i]
            end
            for i = 2:size(x, 2)
                @inbounds ret[:, i] = $f(x[:, i], args...; kargs...)
            end
            ret
        end
    end
end
