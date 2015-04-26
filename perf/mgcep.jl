using MelGeneralizedCepstrums: _mcep, _mgcep, mc2b
import SPTK

function perf_mc2b()
    println("benchmark: mc2b")

    srand(98765)
    mc = rand(21)
    α = 0.41
    n = 50000

    mc2b(mc, α)
    SPTK.mc2b(mc, α)

    @time begin
        elapsed = @elapsed for i=1:n
            mc2b(mc, α)
        end
    end

    @time begin
        elapsed_sptk = @elapsed for i=1:n
            SPTK.mc2b(mc, α)
        end
    end

    r = elapsed/elapsed_sptk
    println("$r x slower than SPTK implementation")
end

function perf_mcep()
    println("benchmark: mcep")

    srand(98765)
    x = rand(1024)

    n = 2000
    order = 25
    α = 0.41

    _mcep(x, order, α)
    SPTK.mcep(x, order, α)

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

    n = 200
    order = 25
    α = 0.41
    γ = -0.1

    _mgcep(x, order, α, γ)
    SPTK.mgcep(x, order, α, γ)

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

gc_disable()

gc()
perf_mcep()

gc()
perf_mgcep()

gc()
perf_mc2b()
