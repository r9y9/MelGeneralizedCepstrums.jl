## Extend functions for matrix input (col-wise) ##

for f in [
          :mc2b!,
          :mgc2b!,
          :gnorm!,
          :ignorm!
          ]
    @eval begin
        function ($f){T<:FloatingPoint}(x::AbstractMatrix{T}, args...; kargs...)
            for i = 1:size(x, 2)
                @inbounds x[:, i] = $f(x[:, i], args...; kargs...)
            end
            x
        end
    end
end

for f in [
          :_mcep,
          :_mgcep,
          :mc2b,
          :mc2e,
          :mgc2b,
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
