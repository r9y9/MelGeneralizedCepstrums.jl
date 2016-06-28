## Extend functions for matrix input (col-wise) ##

const vec2vec_inplace = [
                         :lpc2c!,
                         :mc2b!,
                         :mgc2b!,
                         :b2mc!,
                         :b2c!,
                         :gnorm!,
                         :ignorm!,
                         :freqt!,
                         :frqtr!,
                         :gc2gc!,
                         :mgc2mgc!
                         ]

for f in vec2vec_inplace
    @eval begin
        function ($f)(x::AbstractMatrix, args...; kargs...)
            for i = 1:size(x, 2)
                @inbounds $f(view(x, :, i), args...; kargs...)
            end
            x
        end
    end
end

const vec2vec = [
                 :_mcep,
                 :_mgcep,
                 :lpc2c,
                 :mc2b,
                 :mgc2b,
                 :b2mc,
                 :b2c,
                 :gnorm,
                 :ignorm,
                 :c2ir,
                 :freqt,
                 :gc2gc,
                 :mgc2mgc
                 ]

for f in vec2vec
    @eval begin
        function ($f)(x::AbstractMatrix, args...; kargs...)
            r = $f(x[:, 1], args...; kargs...)
            ret = Array(eltype(x), size(r, 1), size(x, 2))
            for i = 1:length(r)
                @inbounds ret[i, 1] = r[i]
            end
            for i = 2:size(x, 2)
                @inbounds ret[:, i] = $f(view(x, :, i), args...; kargs...)
            end
            ret
        end
    end
end
