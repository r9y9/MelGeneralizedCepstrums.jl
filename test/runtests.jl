using MelGeneralizedCepstrums
using Base.Test

import SPTK
import MelGeneralizedCepstrums: frequency_scale, log_func, rawdata, Mel, Linear,
  StandardLog, GeneralizedLog, AllPoleLog

function test_mcep_type()
    srand(98765)
    x = rand(1024)
    order = 20

    # Linear frequency cepstrum
    mc_typed = mcep(x, order, 0.0)

    # For type stability, mcep always returns a Type{::MelCepstrum}-typed value
    # even if α = 0.0. To get expected type that can be estimated from α,
    # just pass it to the generic constructor as follows:
    mc_typed = MelGeneralizedCepstrum(mc_typed)

    # turn out to be linear cepstrum
    @test isa(mc_typed, LinearCepstrum)

    mc_typed = mcep(x, order, 0.41)
    @test isa(mc_typed, MelCepstrum)
end

function test_mgcep_type()
    srand(98765)
    x = rand(1024)
    order = 20

    # Linear frequency cepstrum
    mgc_typed = mgcep(x, order, 0.0, 0.0)
    @test isa(mgc_typed, LinearCepstrum)

    # Mel-cepstrum
    mgc_typed = mgcep(x, order, 0.41, 0.0)
    @test isa(mgc_typed, MelCepstrum)

    # Generalized cepstrum
    mgc_typed = mgcep(x, order, 0.0, -0.1)
    @test isa(mgc_typed, GeneralizedCepstrum)

    # All-pole cepstrum (Linear prediction)
    mgc_typed = mgcep(x, order, 0.0, -1.0)
    @test isa(mgc_typed, AllPoleCepstrum)

    # Mel all-pole cepstrum (Warped linear prediction)
    mgc_typed = mgcep(x, order, 0.41, -1.0)
    @test isa(mgc_typed, MelAllPoleCepstrum)

    # Mel-Generalized Cesptrum
    mgc_typed = mgcep(x, order, 0.41, -0.1)
    @test isa(mgc_typed, MelGeneralizedCepstrum)
end

function test_mgcep_basics()
    srand(98765)
    c = rand(21)

    mgc = MelGeneralizedCepstrum(0.41, -0.01, c)
    @test isa(mgc, MelFrequencyCepstrum)
    @test isa(mgc, GeneralizedLogCepstrum)
    @test allpass_alpha(mgc) == 0.41
    @test glog_gamma(mgc) == -0.01
    @test order(mgc) == 20
    @test powercoef(mgc) == mgc[1]
    @test size(mgc) == size(c)

    @test frequency_scale(typeof(mgc)) == Mel
    @test log_func(typeof(mgc)) == GeneralizedLog
end

function test_mcep_basics()
    srand(98765)
    c = rand(21)

    mc = MelGeneralizedCepstrum(0.41, 0.0, c)
    @test isa(mc, MelCepstrum)
    @test isa(mc, StandardLogCepstrum)
    @test allpass_alpha(mc) == 0.41
    @test glog_gamma(mc) == 0.0
    @test order(mc) == 20
    @test powercoef(mc) == mc[1]
    @test size(mc) == size(c)

    @test frequency_scale(typeof(mc)) == Mel
    @test log_func(typeof(mc)) == StandardLog
end

function test_gcep_basics()
    srand(98765)
    c = rand(21)

    gc = MelGeneralizedCepstrum(0.0, -0.01, c)
    @test isa(gc, GeneralizedCepstrum)
    @test isa(gc, LinearFrequencyCepstrum)
    @test allpass_alpha(gc) == 0.0
    @test glog_gamma(gc) == -0.01
    @test order(gc) == 20
    @test powercoef(gc) == gc[1]
    @test size(gc) == size(c)

    @test frequency_scale(typeof(gc)) == Linear
    @test log_func(typeof(gc)) == GeneralizedLog

    gc = MelGeneralizedCepstrum(0.0, -1.0, c)
    @test log_func(typeof(gc)) == AllPoleLog
end

function test_mcep(order::Int, α::Float64)
    srand(98765)
    x = rand(1024)
    mc = SPTK.mcep(x, order, α)
    mĉ = MelGeneralizedCepstrums._mcep(x, order, α)
    @test_approx_eq mc mĉ
end

function test_mgcep(order::Int, α::Float64, γ::Float64)
    srand(98765)
    x = rand(1024)
    mgc = SPTK.mgcep(x, order, α, γ; dd=0.001)
    mgĉ = MelGeneralizedCepstrums._mgcep(x, order, α, γ; threshold=0.001)
    @test_approx_eq mgc mgĉ
end

function test_gnorm(γ::Float64)
    srand(98765)
    x = rand(100)
    mc = rand(21)

    g = SPTK.gnorm(mc, γ)
    ĝ = gnorm(mc, γ)
    @test_approx_eq g ĝ
    ĝ = copy(mc)
    gnorm!(ĝ, γ)
    @test_approx_eq g ĝ

    mc = MelGeneralizedCepstrum(0.0, γ, mc)
    ĝ = gnorm(mc)
    @test_approx_eq g rawdata(ĝ)
end

function test_ignorm(γ::Float64)
    srand(98765)
    x = rand(100)
    mc = rand(21)

    g = SPTK.ignorm(mc, γ)
    ĝ = ignorm(mc, γ)
    @test_approx_eq g ĝ
    ĝ = copy(mc)
    ignorm!(ĝ, γ)
    @test_approx_eq g ĝ

    mc = MelGeneralizedCepstrum(0.0, γ, mc)
    ĝ = ignorm(mc)
    @test_approx_eq g rawdata(ĝ)
end

function test_mc2b(α::Float64)
    srand(98765)
    x = rand(100)
    mc = rand(21)

    b = SPTK.mc2b(mc, α)
    b̂ = mc2b(mc, α)
    @test_approx_eq b b̂
    b̂ = copy(mc)
    mc2b!(b̂, α)
    @test_approx_eq b b̂
end

function test_mgc2b(α::Float64, γ::Float64)
    srand(98765)
    x = rand(100)
    mgc = rand(21)

    b = mgc2b(mgc, α, γ)
    b̂ = copy(mgc)
    mgc2b!(b̂, α, γ)
    @test_approx_eq b b̂
    mgc = MelGeneralizedCepstrum(α, γ, mgc)
    b̂ = mgc2b(mgc)
    @test_approx_eq b rawdata(b̂)
end

function test_b2mc(α::Float64)
    srand(98765)
    x = rand(100)
    b = rand(21)

    mc = SPTK.b2mc(b, α)
    mĉ = b2mc(b, α)
    @test_approx_eq mc mĉ
end

function test_c2ir(len::Int=512)
    srand(98765)
    x = rand(100)
    c = rand(21)

    ir = SPTK.c2ir(c, len)
    ir̂ = c2ir(c, len)
    @test_approx_eq ir ir̂
end

function test_freqt(order::Int, α::Float64)
    srand(98765)
    x = rand(100)
    mc = rand(21)

    m = SPTK.freqt(mc, order, α)
    m̂ = freqt(mc, order, α)
    @test_approx_eq m m̂

    mc = MelGeneralizedCepstrum(0.0, 0.0, mc)
    m̂ = freqt(mc, order, α)
    @test_approx_eq m rawdata(m̂)
end

function test_b2c(order::Int, α::Float64)
    srand(98765)
    x = rand(100)
    mc = rand(21)

    m = SPTK.b2c(mc, order, α)
    m̂2 = b2c(mc, order, α)
    @test_approx_eq m m̂2
end

function test_frqtr(order::Int, α::Float64)
    srand(98765)
    x = rand(100)
    mc = rand(21)

    m = SPTK.frqtr(mc, order, α)
    m̂ = frqtr(mc, order, α)
    @test_approx_eq m m̂

    mc = MelGeneralizedCepstrum(0.0, 0.0, mc)
    m̂ = frqtr(mc, order, α)
    @test_approx_eq m rawdata(m̂)
end

function test_gc2gc(order::Int, γ::Float64)
    srand(98765)
    x = rand(100)
    gc = rand(21)

    gc2 = SPTK.gc2gc(gc, 0.0, order, γ)
    gc2̂ = gc2gc(gc, 0.0, order, γ)
    @test_approx_eq gc2 gc2̂

    gc = MelGeneralizedCepstrum(0.0, 0.0, gc)
    gc2̂ = gc2gc(gc, order, γ)
    @test_approx_eq gc2 rawdata(gc2̂)
end

function test_gc2gc(order::Int, α::Float64, γ::Float64)
    srand(98765)
    x = rand(100)
    mgc = rand(21)

    mgc2 = SPTK.mgc2gc(gc, 0.0, 0.0, order, α, γ)
    mgc2̂ = mgc2gc(gc, 0.0, 0.0, order, α, γ)
    @test_approx_eq mgc2 mgc2̂

    mgc = MelGeneralizedCepstrum(0.0, 0.0, gc)
    mgc2̂ = mgc2gc(gc, order, α, γ)
    @test_approx_eq mgc2 rawdata(mgc2̂)
    @test allpass_alpha(mgc2²) == α
    @test glog_gamma(mgc2²) == γ
end

test_mcep_type()
test_mgcep_type()

test_mgcep_basics()
test_mcep_basics()
test_gcep_basics()

for order in 10:2:30
    for α in [-0.544, -0.41, -0.35, 0.0, 0.35, 0.41, 0.544]
        println("mcep: testing with order=$order, α=$α")
        test_mcep(order, α)
    end
end

for order in 25:5:35
    for α in [-0.544, -0.41, -0.35, 0.0, 0.35, 0.41, 0.544]
        for γ in [-1.0, -0.75, -0.5, -0.25, 0.0]
            println("mgcep: testing with order=$order, α=$α, γ=$γ")
            test_mgcep(order, α, γ)
        end
    end
end

for γ in [-1.0, -0.75, -0.5, -0.25, 0.0]
    println("gnorm: testing with γ=$γ")
    test_gnorm(γ)
end

for γ in [-1.0, -0.75, -0.5, -0.25, 0.0]
    println("ignorm: testing with γ=$γ")
    test_ignorm(γ)
end

for α in [-0.544, -0.41, -0.35, 0.0, 0.35, 0.41, 0.544]
    println("mc2b: testing with α=$α")
    test_mc2b(α)
end

# TODO(ryuichi): tests don't passed when α -0.544, -0.41, -0.35, γ = -1.0
for α in [0.0, 0.35, 0.41, 0.544]
    for γ in [-1.0, -0.75, -0.5, -0.25, 0.0]
        println("mgc2b: testing with α=$α, γ=$γ")
        test_mgc2b(α, γ)
    end
end

for α in [-0.544, -0.41, -0.35, 0.0, 0.35, 0.41, 0.544]
    println("mc2b: testing with α=$α")
    test_b2mc(α)
end

for len in [128, 256, 512, 1024]
    println("c2ir: testing with len=$len")
    test_c2ir(len)
end

for order in 10:2:30
    for α in [-0.544, -0.41, -0.35, 0.0, 0.35, 0.41, 0.544]
        println("freqt: testing with order=$order, α=$α")
        test_freqt(order, α)
    end
end

for order in 10:2:30
    for α in [-0.544, -0.41, -0.35, 0.0, 0.35, 0.41, 0.544]
        println("b2c: testing with order=$order, α=$α")
        test_b2c(order, α)
    end
end

for order in 10:2:30
    for α in [-0.544, -0.41, -0.35, 0.0, 0.35, 0.41, 0.544]
        println("frqtr: testing with order=$order, α=$α")
        test_frqtr(order, α)
    end
end

for order in 15:5:35
    for γ in [-1.0, -0.75, -0.5, -0.25, 0.0]
        println("gc2gc: testing with order=$order, γ=$γ")
        test_gc2gc(order, γ)
    end
end

for order in 15:5:35
    for α in [-0.544, -0.41, -0.35, 0.35, 0.41, 0.544]
        for γ in [-1.0, -0.75, -0.5, -0.25, 0.0]
            println("mgc2mgc: testing with order=$order, α=$α, γ=$γ")
            test_gc2gc(order, γ)
        end
    end
end
