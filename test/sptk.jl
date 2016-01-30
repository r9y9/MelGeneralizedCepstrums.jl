# Check consistency with SPTK

module SPTKTest

using MelGeneralizedCepstrums
using Base.Test

import SPTK
import MelGeneralizedCepstrums: lpc2c!, b2c, frqtr, frqtr!

function test_mcep(order::Int, α::Float64)
    srand(98765)
    x = rand(512)
    mc = SPTK.mcep(x, order, α)
    mĉ = MelGeneralizedCepstrums._mcep(x, order, α)
    @test_approx_eq mc mĉ
end

function test_mgcep(order::Int, α::Float64, γ::Float64)
    srand(98765)
    x = rand(512)
    mgc = SPTK.mgcep(x, order, α, γ; threshold=0.001)
    mgĉ = MelGeneralizedCepstrums._mgcep(x, order, α, γ; criteria=0.001)
    @test_approx_eq mgc mgĉ
end

function test_lpc2c()
    srand(98765)
    a = rand(21)

    c = SPTK.lpc2c(a)
    ĉ = lpc2c(a)
    @test_approx_eq c ĉ
    ĉ = similar(a)
    lpc2c!(ĉ, a)
    @test_approx_eq c ĉ
end

function test_gnorm(γ::Float64)
    srand(98765)
    mc = rand(21)

    g = SPTK.gnorm(mc, γ)
    ĝ = gnorm(mc, γ)
    @test_approx_eq g ĝ
    ĝ = copy(mc)
    gnorm!(ĝ, γ)
    @test_approx_eq g ĝ
end

function test_ignorm(γ::Float64)
    srand(98765)
    mc = rand(21)

    g = SPTK.ignorm(mc, γ)
    ĝ = ignorm(mc, γ)
    @test_approx_eq g ĝ
    ĝ = copy(mc)
    ignorm!(ĝ, γ)
    @test_approx_eq g ĝ
end

function test_mc2b(α::Float64)
    srand(98765)
    mc = rand(21)

    b = SPTK.mc2b(mc, α)
    b̂ = mc2b(mc, α)
    @test_approx_eq b b̂
    b̂ = copy(mc)
    mc2b!(b̂, α)
    @test_approx_eq b b̂
end


function test_b2mc(α::Float64)
    srand(98765)
    b = rand(21)

    mc = SPTK.b2mc(b, α)
    mĉ = b2mc(b, α)
    @test_approx_eq mc mĉ
end

function test_c2ir(len::Int=512)
    srand(98765)
    c = rand(21)

    ir = SPTK.c2ir(c, len)
    ir̂ = c2ir(c, len)
    @test_approx_eq ir ir̂
end

function test_freqt(order::Int, α::Float64)
    srand(98765)
    mc = rand(21)

    m = SPTK.freqt(mc, order, α)
    m̂ = freqt(mc, order, α)
    @test_approx_eq m m̂
end

function test_b2c(order::Int, α::Float64)
    srand(98765)
    mc = rand(21)

    m = SPTK.b2c(mc, order, α)
    m̂2 = b2c(mc, order, α)
    @test_approx_eq m m̂2
end

function test_frqtr(order::Int, α::Float64)
    srand(98765)
    mc = rand(21)

    m = SPTK.frqtr(mc, order, α)
    m̂ = frqtr(mc, order, α)
    @test_approx_eq m m̂
end

function test_gc2gc(order::Int, γ::Float64)
    srand(98765)
    gc = rand(21)

    gc2 = SPTK.gc2gc(gc, 0.0, order, γ)
    gc2̂ = gc2gc(gc, 0.0, order, γ)
    @test_approx_eq gc2 gc2̂
end

function test_mgc2mgc(order::Int, α::Float64, γ::Float64)
    srand(98765)
    mgc = rand(21)

    mgc2 = SPTK.mgc2mgc(mgc, 0.41, 0.0, order, α, γ)
    mgc2̂ = mgc2mgc(mgc, 0.41, 0.0, order, α, γ)
    @test_approx_eq mgc2 mgc2̂
end

function test_mgc2sp(α::Float64, γ::Float64)
    srand(98765)
    mgc = rand(21)

    fftlen = 512
    sp = SPTK.mgc2sp(mgc, α, γ, fftlen)
    sp̂ = mgc2sp(mgc, α, γ, fftlen)
    @test_approx_eq sp sp̂
    @test length(sp̂) == fftlen>>1 + 1
    @test eltype(sp̂) == Complex{eltype(mgc)}
end

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
            test_mgc2mgc(order, α, γ)
        end
    end
end

for α in [-0.544, -0.41, -0.35, 0.35, 0.41, 0.544]
    for γ in [-0.75, -0.5, -0.25, 0.0] # TODO(ryuichi) -1.0
        println("mgc2sp: testing with α=$α and γ=$γ")
        test_mgc2sp(α, γ)
    end
end

println("lpc2c: testing")
test_lpc2c()

end # module SPTKTestModule
