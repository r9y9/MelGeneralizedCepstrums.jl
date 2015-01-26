using MelGeneralizedCepstrums: _mcep, _mgcep
import SPTK

function perf_mcep()
    println("benchmark: mcep")

    srand(98765)
    x = rand(1024)

    n = 500
    order = 25
    α = 0.41

    @time begin
        elapsed = @elapsed for i=1:n
            _mcep(x, order, α)
        end
    end

    @time begin
        elapsed_sptk = @elapsed for i=1:n
            SPTK.mcep(x, order, α)
        end
    end

    r = elapsed/elapsed_sptk
    println("$r x slower than SPTK implementation")
end

function perf_mgcep()
    println("benchmark: mgcep")

    srand(98765)
    x = rand(1024)

    n = 100
    order = 25
    α = 0.41
    γ = -0.1

    @time begin
        elapsed = @elapsed for i=1:n
            _mgcep(x, order, α, γ)
        end
    end

    @time begin
        elapsed_sptk  = @elapsed for i=1:n
            SPTK.mgcep(x, order, α, γ)
        end
    end

    r = elapsed/elapsed_sptk
    println("$r x slower than SPTK implementation")
end

perf_mcep()
perf_mgcep()
