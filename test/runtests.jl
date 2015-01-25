using MelGeneralizedCepstrums
using Base.Test

import SPTK
import MelGeneralizedCepstrums: frequency_scale, log_func, rawdata, Mel, Linear,
  StandardLog, GeneralizedLog, AllPoleLog, mgcepnorm!

function test_mcep_type()
    srand(98765)
    x = rand(512)
    order = 20

    # Linear frequency cepstrum
    mc_typed = mcep(x, order, 0.0)

    # For type stability, mcep always returns a Type{::MelCepstrum}-typed value
    # even if α = 0.0. To get expected type that can be estimated from α,
    # just pass it to the generic constructor as follows:
    @test !isa(mc_typed, LinearCepstrum)
    mc_typed = MelGeneralizedCepstrum(mc_typed)

    # turn out to be linear cepstrum
    @test isa(mc_typed, LinearCepstrum)

    mc_typed = mcep(x, order, 0.41)
    @test isa(mc_typed, MelCepstrum)
end

function test_mgcep_type()
    srand(98765)
    x = rand(512)
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

function test_lpc_type()
    srand(98765)
    x = rand(512)
    order = 20

    lpc_typed = lpc(x, order)
    @test isa(lpc_typed, LinearPredictionCoef)
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

    c = repmat(c, 1, 2)
    mgc = MelGeneralizedCepstrum(0.41, -0.01, c)
    @test size(mgc) == (21, 2)
    @test isa(mgc[:,1], MelGeneralizedCepstrum)

    @test_throws ArgumentError MelGeneralizedCepstrum(1.0, -0.01, c)
    @test_throws ArgumentError MelGeneralizedCepstrum(0.41, 0.01, c)
    @test_throws ArgumentError mgcepnorm!(rand(2), 0.41, -0.01, -1)
    @test_throws ArgumentError mgcepnorm!(rand(2), 0.41, -0.01, 6)
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

function test_lpc_basics()
    srand(98765)
    c = rand(21)

    l = MelLinearPredictionCoef(0.0, c, false)
    @test isa(l, LinearPredictionCoef)
    l = MelLinearPredictionCoef(0.41, c, false)
    @test !isa(l, LinearPredictionCoef)
    @test isa(l, MelLinearPredictionCoef)

    @test allpass_alpha(l) == 0.41
    @test glog_gamma(l) == -1.0
    @test order(l) == 20
    @test powercoef(l) == l[1]
    @test size(l) == size(c)

    @test frequency_scale(typeof(l)) == Mel
    @test log_func(typeof(l)) == AllPoleLog

    c = repmat(c, 1, 2)
    l = MelLinearPredictionCoef(0.0, c, false)
    @test size(l) == (21, 2)
    @test isa(l[:,1], MelLinearPredictionCoef)

    @test_throws ArgumentError MelLinearPredictionCoef(1.0, c, false)
end

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
    mgc = SPTK.mgcep(x, order, α, γ; dd=0.001)
    mgĉ = MelGeneralizedCepstrums._mgcep(x, order, α, γ; criteria=0.001)
    @test_approx_eq mgc mgĉ
end

function test_mgcep_otypes(α::Float64, γ::Float64)
    srand(98765)
    x = rand(512)
    order = 20

    mgc_2 = MelGeneralizedCepstrums._mgcep(x, order, α, γ; otype=2)
    mgc_3 = MelGeneralizedCepstrums._mgcep(x, order, α, γ; otype=3)
    mgc_4 = MelGeneralizedCepstrums._mgcep(x, order, α, γ; otype=4)
    mgc_5 = MelGeneralizedCepstrums._mgcep(x, order, α, γ; otype=5)

    mgc_2̂ = MelGeneralizedCepstrums._mgcep(10x, order, α, γ; otype=2)
    mgc_3̂ = MelGeneralizedCepstrums._mgcep(10x, order, α, γ; otype=3)
    mgc_4̂ = MelGeneralizedCepstrums._mgcep(10x, order, α, γ; otype=4)
    mgc_5̂ = MelGeneralizedCepstrums._mgcep(10x, order, α, γ; otype=5)

    # check if gain normalized
    @test_approx_eq_eps mgc_2[2:end] mgc_2̂[2:end] 1.0e-3
    @test_approx_eq_eps mgc_3[2:end] mgc_3̂[2:end] 1.0e-3
    @test_approx_eq_eps mgc_4[2:end] mgc_4̂[2:end] 1.0e-3
    @test_approx_eq_eps mgc_5[2:end] mgc_5̂[2:end] 1.0e-3
end

function test_extend()
    srand(98765)
    mc = rand(21)
    mc_mat = repmat(mc, 1, 2)
    mc_mat2 = copy(mc_mat)
    mc_submat = sub(mc_mat2, 1:size(mc_mat2, 1), 1:size(mc_mat2, 2))

    # case 1: Vector input
    r1 = mc2b(mc, 0.41)
    # case 2: Matrix input
    r2 = mc2b(mc_mat, 0.41)
    # case 3: SubArray input
    r3 = mc2b(mc_submat, 0.41)

    # make sure non-inplace function doesn't destoy input
    @test r1 != mc
    @test r2 != mc_mat
    @test r3 != mc_mat2

    @test r1 == r2[:,1]
    @test r1 == r3[:,1]
end

function test_inplace_extend()
    srand(98765)
    mc = rand(21)
    mc_mat = repmat(mc, 1, 2)
    mc_mat2 = copy(mc_mat)
    mc_submat = sub(mc_mat2, 1:size(mc_mat2, 1), 1:size(mc_mat2, 2))

    mc_org = copy(mc)

    # case 1: Vector input
    mc2b!(mc, 0.41)
    # case 2: Matrix input
    mc2b!(mc_mat, 0.41)
    # case 3: SubArray input
    mc2b!(mc_submat, 0.41)

    @test mc != mc_org

    # all of the cases above must get an equal result for same input
    @test mc == mc_mat[:,1]
    @test mc == mc_submat[:,1]
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

    a_typed = MelLinearPredictionCoef(0.0, a, false)
    c = lpc2c(a_typed)
    @test isa(c, LinearCepstrum)
    @test c[1] == log(a_typed[1])

    # assume a has loggain
    a_typed = MelLinearPredictionCoef(0.0, a, true)
    c = lpc2c(a_typed)
    @test c[1] == a_typed[1]

    a_typed = MelLinearPredictionCoef(0.41, a, false)
    @test_throws ArgumentError lpc2c(a_typed)
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

    mc = MelGeneralizedCepstrum(0.0, γ, mc)
    ĝ = gnorm(mc)
    @test_approx_eq g rawdata(ĝ)
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

    mc = MelGeneralizedCepstrum(0.0, γ, mc)
    ĝ = ignorm(mc)
    @test_approx_eq g rawdata(ĝ)
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

    if α != 0.0
        mc = MelGeneralizedCepstrum(α, 0.0, mc)
        b̂ = mc2b(mc)
        @test isa(b̂, MLSADFCoef)
    end
end

function test_mgc2b(α::Float64, γ::Float64)
    srand(98765)
    mgc = rand(21)

    b = mgc2b(mgc, α, γ)
    b̂ = copy(mgc)
    mgc2b!(b̂, α, γ)
    @test_approx_eq b b̂
    mgc = MelGeneralizedCepstrum(α, γ, mgc)
    b̂ = mgc2b(mgc)
    @test_approx_eq b rawdata(b̂)
    if α != 0.0 && γ != 0.0 && γ != -1.0
        @test isa(b̂, MGLSADFCoef)
    end
end

function test_b2mc(α::Float64)
    srand(98765)
    b = rand(21)

    mc = SPTK.b2mc(b, α)
    mĉ = b2mc(b, α)
    @test_approx_eq mc mĉ

    mc = MelGeneralizedCepstrum(α, 0.0, b)
    b = mc2b(mc)
    mĉ = b2mc(b)
    if α != 0.0
        @test isa(mĉ, MelCepstrum)
        @test rawdata(mc) != rawdata(b)
    end
    @test_approx_eq rawdata(mc) rawdata(mĉ)
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

    mc = MelGeneralizedCepstrum(0.0, 0.0, mc)
    m̂ = freqt(mc, order, α)
    @test_approx_eq m rawdata(m̂)
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

    mc = MelGeneralizedCepstrum(0.0, 0.0, mc)
    m̂ = frqtr(mc, order, α)
    @test_approx_eq m rawdata(m̂)
end

function test_gc2gc(order::Int, γ::Float64)
    srand(98765)
    gc = rand(21)

    gc2 = SPTK.gc2gc(gc, 0.0, order, γ)
    gc2̂ = gc2gc(gc, 0.0, order, γ)
    @test_approx_eq gc2 gc2̂

    gc = MelGeneralizedCepstrum(0.0, 0.0, gc)
    gc2̂ = gc2gc(gc, order, γ)
    @test_approx_eq gc2 rawdata(gc2̂)
end

function test_mgc2mgc(order::Int, α::Float64, γ::Float64)
    srand(98765)
    mgc = rand(21)

    mgc2 = SPTK.mgc2mgc(mgc, 0.41, 0.0, order, α, γ)
    mgc2̂ = mgc2mgc(mgc, 0.41, 0.0, order, α, γ)
    @test_approx_eq mgc2 mgc2̂

    mgc = MelGeneralizedCepstrum(0.41, 0.0, mgc)
    mgc2̂ = mgc2mgc(mgc, order, α, γ)
    @test_approx_eq mgc2 rawdata(mgc2̂)
    @test allpass_alpha(mgc2̂) == α
    @test glog_gamma(mgc2̂) == γ
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

    mgc = MelGeneralizedCepstrum(α, γ, mgc)
    sp̂ = mgc2sp(mgc, fftlen)
    @test length(sp̂) == fftlen>>1 + 1
end

test_mcep_type()
test_mgcep_type()
test_lpc_type()

test_mgcep_basics()
test_mcep_basics()
test_gcep_basics()
test_lpc_basics()

test_extend()
test_inplace_extend()

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

for α in [-0.35, 0.0, 0.35]
    for γ in [-1.0, -0.5, 0.0]
        println("mgcep_otypes: testing with α=$α and γ=$γ")
        test_mgcep_otypes(α, γ)
    end
end

println("lpc2c: testing")
test_lpc2c()

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

for f in (:uint, :int, :float)
    @eval begin
        @test_approx_eq  mcepalpha($f(8000))  0.312
        @test_approx_eq  mcepalpha($f(11025)) 0.357
        @test_approx_eq  mcepalpha($f(16000)) 0.41
        @test_approx_eq  mcepalpha($f(22050)) 0.455
        @test_approx_eq  mcepalpha($f(44100)) 0.544
        @test_approx_eq  mcepalpha($f(48000)) 0.554
    end
end
