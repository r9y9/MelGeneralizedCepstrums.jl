using MelGeneralizedCepstrums

import MelGeneralizedCepstrums: _mcep, _mgcep
import SPTK


function perf(f1, f2, params...; name="foo", n=100, N=1024)
    println("benchmark: $name")

    srand(98765)
    x = rand(N)

    @show params
    f1(x, params...)
    f2(x, params...)

    @time begin
        elapsed = @elapsed for i=1:n
            f1(x, params...)
        end
    end

    @time begin
        elapsed_sptk = @elapsed for i=1:n
            f2(x, params...)
        end
    end

    r = elapsed/elapsed_sptk
    println("$r x slower than SPTK implementation")
end

gc_enable(false)

gc()
perf(_mgcep, SPTK.mgcep, 25, 0.41, -0.1; name="mgcep", n=200, N=1024)

gc()
perf(_mcep, SPTK.mcep, 25, 0.41; name="mcep", n=2000, N=1024)

gc()
perf(mc2b, SPTK.mc2b, 0.41; name="mb2b", n=20000, N=25)

gc()
perf(b2mc, SPTK.b2mc, 0.41; name="b2mc", n=20000, N=25)

gc()
perf(freqt, SPTK.freqt, 30, 0.41; name="freqt", n=20000, N=25)
