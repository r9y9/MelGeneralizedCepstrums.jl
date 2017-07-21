using MelGeneralizedCepstrums
using Base.Test

import MelGeneralizedCepstrums: freq_form, log_form, rawdata, MelFrequency,
    LinearFrequency, StandardLog, GeneralizedLog, AllPoleLog, mgcepnorm!, retype,
    mc2e

@static is_linux() ? include("sptk.jl") : nothing

### Test functions ###

function test_spectral_param_state()
    srand(98765)

    # state that doesn't have loggain
    s = SpectralParamState(LinearPredictionCoef(20), rand(21, 2), false, true)
    s_copy= copy(s)
    loggain!(s_copy)
    @test s_copy[1,1] == log(s[1,1])
    @test s_copy[1,2] == log(s[1,2])

    unloggain!(s_copy)
    @test s_copy[1,1] == s[1,1]
    @test s_copy[1,2] == s[1,2]

    # state that has loggain
    s = SpectralParamState(LinearCepstrum(20), rand(21, 2), true, true)
    s_copy= copy(s)

    loggain!(s_copy)
    @test s_copy[1,1] == s[1,1]
    @test s_copy[1,2] == s[1,2]
end

function test_mcep_special_cases(order=20)
    # Linear frequency cepstrum
    mc = MelCepstrum(order, 0.0)
    mc_retyped = MelGeneralizedCepstrum(mc)
    # should turn out to be linear cepstrum
    @test isa(mc_retyped, LinearCepstrum)
end

function test_mcep_basics(order=20, α=0.41)
    srand(98765)
    x = rand(512)

    mc = MelCepstrum(order, α)
    @test isa(mc, MelCepstrum)
    @test params(mc) == (order, α)
    @test allpass_alpha(mc) == α
    @test glog_gamma(mc) == 0.0
    @test param_order(mc) == order

    @test freq_form(typeof(mc)) == MelFrequency
    @test log_form(typeof(mc)) == StandardLog

    state = estimate(mc, x)
    @test all(isfinite.(state))
    @test mcep(x, order, α) ≈ state

    # details
    mc1 = MelGeneralizedCepstrums._mcep(x, order, α)
    mc2 = periodogram2mcep(abs2.(rfft(x)), order, α)
    @test mc1 ≈ mc2
end

function test_gcep_special_cases(order=20)
    gc = GeneralizedCepstrum(order, 0.0)
    gc_retyped = MelGeneralizedCepstrum(gc)
    @test isa(gc_retyped, LinearCepstrum)

    gc = GeneralizedCepstrum(order, -1.0)
    gc_retyped = MelGeneralizedCepstrum(gc)
    @test isa(gc_retyped, AllPoleCepstrum)
end

function test_gcep_basics(order=20, γ=-0.01)
    srand(98765)
    x = rand(512)

    gc = GeneralizedCepstrum(order, γ)
    @test isa(gc, GeneralizedCepstrum)
    @test params(gc) == (order, γ)
    @test allpass_alpha(gc) == 0.0
    @test glog_gamma(gc) == γ
    @test param_order(gc) == order

    @test freq_form(typeof(gc)) == LinearFrequency
    @test log_form(typeof(gc)) == GeneralizedLog

    state = estimate(gc, x, norm=false)
    @test all(isfinite.(state))
    @test gcep(x, order, γ, norm=false) ≈ state
    @test !gain_normalized(state)

    # gain-unnormalized
    # should have different values for different scaled signals
    @test rawdata(estimate(gc, 10x))[2:end] != rawdata(state)[2:end]

    # gain normalized
    state = estimate(gc, x, norm=true)
    @test gain_normalized(state)
    # should have same values for different scaled signals
    @test rawdata(estimate(gc, 10x, norm=true))[2:end] ≈ rawdata(state)[2:end]
end

function test_mgcep_special_cases(order=20)
    # Linear frequency cepstrum
    mgc = MelGeneralizedCepstrum(order, 0.0, 0.0)
    mgc_retyped = MelGeneralizedCepstrum(mgc)
    # should turn out to be LinearCepstrum
    @test isa(mgc_retyped, LinearCepstrum)

    # Mel-cepstrum
    mgc = MelGeneralizedCepstrum(order, 0.41, 0.0)
    mgc_retyped = MelGeneralizedCepstrum(mgc)
    # should turn out to be MelCepstrum
    @test isa(mgc_retyped, MelCepstrum)

    # Generalized cepstrum
    mgc = MelGeneralizedCepstrum(order, 0.0, -0.1)
    mgc_retyped = MelGeneralizedCepstrum(mgc)
    @test isa(mgc_retyped, GeneralizedCepstrum)

    # All-pole cepstrum
    mgc = MelGeneralizedCepstrum(order, 0.0, -1.0)
    mgc_retyped = MelGeneralizedCepstrum(mgc)
    @test isa(mgc_retyped, AllPoleCepstrum)

    # Mel all-pole cepstrum
    mgc = MelGeneralizedCepstrum(order, 0.41, -1.0)
    mgc_retyped = MelGeneralizedCepstrum(mgc)
    @test isa(mgc_retyped, MelAllPoleCepstrum)
end

function test_mgcep_basics(order=20, α=0.41, γ=-0.01)
    srand(98765)
    x = rand(512)

    mgc = MelGeneralizedCepstrum(order, α, γ)
    @test isa(mgc, MelGeneralizedCepstrum)
    @test params(mgc) == (order, α, γ)
    @test allpass_alpha(mgc) == α
    @test glog_gamma(mgc) == γ
    @test param_order(mgc) == order

    @test freq_form(typeof(mgc)) == MelFrequency
    @test log_form(typeof(mgc)) == GeneralizedLog

    state = estimate(mgc, x)
    @test all(isfinite.(state))
    @test mgcep(x, order, α, γ) ≈ state

    # should only accept otype=0
    for otype in [1,2,3,4,5]
        @test_throws ArgumentError estimate(mgc, x, otype=otype)
    end

    @test_throws ArgumentError MelGeneralizedCepstrum(0, -.41, -0.01)
    @test_throws ArgumentError MelGeneralizedCepstrum(order, 1.0, -0.01)
    @test_throws ArgumentError MelGeneralizedCepstrum(order, 0.41, 0.01)
    @test_throws ArgumentError mgcepnorm!(rand(2), 0.41, -0.01, -1)
    @test_throws ArgumentError mgcepnorm!(rand(2), 0.41, -0.01, 6)
end

function test_lpc_basics()
    srand(98765)
    x = rand(512)
    order = 20

    lpcdef = LinearPredictionCoef(order)
    @test isa(lpcdef, LinearPredictionCoef)
    state = estimate(lpcdef, x)
    @test all(isfinite.(state))
    @test lpc(x, order) ≈ state
    @test isa(paramdef(v), LinearPredictionCoef)

    # LSP
    lspstate = lpc2lsp(state)
    @test isa(paramdef(lspstate), LineSpectralPair)
    @test !has_loggain(lpsv)
    @test !ready_to_filt(lspstate)
    @test has_loggain(lpc2lsp(state, loggain=true))

    # PARCOR
    parstate = lpc2par(state)
    @test isa(paramdef(parstate), PartialAutoCorrelation)
    @test !has_loggain(parstate)
    @test !ready_to_filt(parstate)

    # LPC cepstrum
    mgcdef = AllPoleCepstrum(order)
    @test isa(mgcdef, AllPoleCepstrum)
    @test allpass_alpha(mgcdef) == 0.0
    @test glog_gamma(mgcdef) == -1.0
    @test param_order(mgcdef) == order

    @test freq_form(typeof(mgcdef)) == LinearFrequency
    @test log_form(typeof(mgcdef)) == AllPoleLog

    state = estimate(mgcdef, x)
    @test all(isfinite.(state))
    @test has_loggain(state)
    @test gain_normalized(state)
    @test !ready_to_filt(state)
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
    @test isapprox(mgc_2[2:end], mgc_2̂[2:end], atol=1.0e-3)
    @test isapprox(mgc_3[2:end], mgc_3̂[2:end], atol=1.0e-3)
    @test isapprox(mgc_4[2:end], mgc_4̂[2:end], atol=1.0e-3)
    @test isapprox(mgc_5[2:end], mgc_5̂[2:end], atol=1.0e-3)
end

function test_lpc2c(order=20)
    srand(98765)
    x = rand(512)

    state = estimate(LinearPredictionCoef(order), x)
    @assert !has_loggain(state)
    @test gain_normalized(state)
    @test !ready_to_filt(state)

    state2 = lpc2c(state)
    @test rawdata(state2) ≈ lpc2c(rawdata(state))
    @test has_loggain(state2)
    @test gain_normalized(state2)
    @test !ready_to_filt(state2)
end

function test_lpc_and_par_conversion(order=20)
    srand(98765)
    x = rand(512)

    lpc_state1 = estimate(LinearPredictionCoef(order), x)
    @test all(isfinite.(rawdata(lpc_state1)))
    par_state = lpc2par(lpc_state1)
    @test all(isfinite.(rawdata(par_state)))
    lpc_state2 = par2lpc(par_state)

    @test rawdata(lpc_state1) ≈ rawdata(lpc_state2)
end

function test_mc2b(order=20, α=0.41)
    srand(98765)
    x = rand(512)

    mcdef = MelCepstrum(order, α)
    state = estimate(mcdef, x)
    @test has_loggain(state)
    @test gain_normalized(state)
    @test !ready_to_filt(state)

    b = mc2b(state)
    @test rawdata(b) ≈ mc2b(rawdata(state), α)
    @test has_loggain(b)
    @test gain_normalized(b)
    @test ready_to_filt(b)

    mc2b!(state)
    @test rawdata(state) ≈ rawdata(b)
end

function test_mc2e(order=20, α=0.41, len=512)
    srand(98765)
    x = rand(512)

    state = estimate(MelCepstrum(order, α), x)

    e = mc2e(state, len)
    @test all(isfinite.(e))
    @test e ≈ mc2e(rawdata(state), α, len)
end

function test_mgc2b(order=20, α=0.41, γ=-0.01)
    srand(98765)
    x = rand(512)

    state_mgc = estimate(MelGeneralizedCepstrum(order, α, γ), x)
    state_b = mgc2b(state_mgc)

    @test has_loggain(state_b)
    @test gain_normalized(state_b)
    @test ready_to_filt(state_b)
    @test rawdata(state_b) ≈ mgc2b(rawdata(state_mgc), α, γ)

    # mgc2b should not accept filter coefficients
    @test_throws ArgumentError mgc2b(state_b)
    @test_throws ArgumentError mgc2b!(state_b)

    state2 = copy(state_mgc)
    mgc2b!(state2)
    @test rawdata(state2) ≈ rawdata(state_b)
end

function test_mgc2sp(order=20, α=0.41, γ=-0.01, fftlen=512)
    srand(98765)
    x = rand(512)

    state = estimate(MelGeneralizedCepstrum(order, α, γ), x)

    logsp = mgc2sp(state, fftlen)
    @test logsp ≈ mgc2sp(rawdata(state), α, γ, fftlen)
    @test length(logsp) == fftlen>>1 + 1
    @test all(isfinite.(logsp))
end

function test_b2mc(order=20, α=0.41)
    srand(98765)
    x = rand(512)

    state = estimate(MelCepstrum(order, α), x)

    b = mc2b(state)
    @test ready_to_filt(b)

    mc = b2mc(b)
    @test rawdata(mc) ≈ b2mc(rawdata(b), α)
    @test has_loggain(mc)
    @test gain_normalized(mc)
    @test !ready_to_filt(mc)

    # check invertibility
    @test rawdata(state) ≈ rawdata(b2mc(mc2b(state)))

    b2 = copy(b)
    b2mc!(b2)
    @test rawdata(b2) ≈ rawdata(mc)
end

function test_gnorm_and_ignorm(order=20, γ=-0.01)
    srand(98765)
    x = rand(512)

    state1 = estimate( GeneralizedCepstrum(order, γ), x)
    @test !gain_normalized(state1)
    state2 = gnorm(state1)
    @test gain_normalized(state2)
    state3 = ignorm(state2)
    @test !gain_normalized(state3)

    # check invertibility
    @test rawdata(state1) ≈ rawdata(state3)

    # inplace version
    state4 = copy(state1)
    gnorm!(state4)
    @test gain_normalized(state4)
    @test rawdata(state4) ≈ rawdata(state2)

    ignorm!(state4)
    @test rawdata(state4) ≈ rawdata(state1)
end

function test_freqt(src_order, src_α, dst_order=20, dst_α=0.35)
    srand(98765)
    x = rand(512)

    state1 = estimate(MelCepstrum(src_order, src_α), x)
    state2 = freqt(state1, dst_order, dst_α)
    @test allpass_alpha(paramdef(state2)) == dst_α
    @test length(state2) == dst_order + 1
    @test all(isfinite.(state1))
    @test all(isfinite.(state2))

    state4 = estimate(MelGeneralizedCepstrum(20, 0.41, -0.02), x)
    state5 = freqt(state4, dst_order, dst_α)
    @test allpass_alpha(paramdef(state5)) == dst_α
    @test length(state5) == dst_order + 1
    @test all(isfinite.(state4))
    @test all(isfinite.(state5))
end

function test_gc2gc(src_order=20, src_α=0.35, src_γ=0.0,
                    dst_order=20, dst_γ=-0.05)
    srand(98765)
    x = rand(512)

    def = MelGeneralizedCepstrum(src_order, src_α, src_γ)

    state1 = estimate(def, x)
    state2 = gc2gc(state1, dst_order, dst_γ)
    @test glog_gamma(paramdef(state2)) == dst_γ
    if dst_γ != -1.0 && dst_γ != 0.0
        @test !gain_normalized(state2)
    end
    @test !ready_to_filt(state2)
end

function test_mgc2mgc(src_order, src_α, src_γ, dst_order, dst_α, dst_γ)
    srand(98765)
    x = rand(512)

    state = estimate(MelGeneralizedCepstrum(src_order, src_α, src_γ), x)
    state2 = mgc2mgc(state, dst_order, dst_α, dst_γ)
    if log_form(typeof(paramdef(state2))) == StandardLog
        @test gain_normalized(state2)
    else
        @test !gain_normalized(state2)
    end
    @test !ready_to_filt(state2)
end

function test_extend()
    srand(98765)
    mc = rand(21)
    mc_mat = repmat(mc, 1, 2)
    mc_mat2 = copy(mc_mat)
    mc_submat = view(mc_mat2, 1:size(mc_mat2, 1), 1:size(mc_mat2, 2))

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
    mc_submat = view(mc_mat2, 1:size(mc_mat2, 1), 1:size(mc_mat2, 2))

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

###  Perform tests ###

println("testing: spectral param state")
test_spectral_param_state()

let
    println("testing: gcep")
    test_gcep_special_cases()
    test_gcep_basics()
end

let
    println("testing: gcep")
    test_mcep_special_cases()
    test_mcep_basics()
end

let
    println("testing: mgcep")
    test_mgcep_special_cases()
    test_mgcep_basics()
end

for α in [-0.35, 0.0, 0.35]
    for γ in [-1.0, -0.5, 0.0]
        println("mgcep: testing with α=$α and γ=$γ")
        test_mgcep_otypes(α, γ)
    end
end

println("testing: mcepalpha")
for elty in (UInt, Int, Float32, Float64)
    @test mcepalpha(convert(elty, 8000)) ≈ 0.312
    @test mcepalpha(convert(elty, 11025)) ≈ 0.357
    @test mcepalpha(convert(elty, 16000)) ≈ 0.41
    @test mcepalpha(convert(elty, 22050)) ≈ 0.455
    @test mcepalpha(convert(elty, 44100)) ≈ 0.544
    @test mcepalpha(convert(elty, 48000)) ≈ 0.554
end

for order in 15:5:25
    println("lpc2c: testing with order=$order")
    test_lpc2c(order)
end

for order in 15:5:25
    println("par2lpc and lpc2par: testing with order=$order")
    test_lpc_and_par_conversion(order)
end

for order in 15:5:25
    for α in [0.0, 0.35]
        println("mc2b: testing with order=$order, α=$α")
        test_mc2b(order, α)
    end
end

let
    println("testing: mc2b exceptions")
    srand(98765)
    x = rand(512)
    b = mc2b(estimate(MelCepstrum(20, 0.41), x))
    # mc2b assumes input is not filter coef.
    @test_throws ArgumentError mc2b(b)
    @test_throws ArgumentError mc2b!(b)
    @test_throws ArgumentError mc2b(mgcep(x, 20, 0.41, -0.2))
end

for order in [20, 25, 30]
    for α in [0.0, 0.35, 0.41, 0.544]
        for fftlen in [256, 512]
            println("mc2e: testing with order=$order, α=$α, fftlen=$fftlen")
            test_mc2e(order, α, fftlen)
        end
    end
end

let
    println("testing: mc2e exceptions")
    srand(98765)
    x = rand(512)
    state = estimate(MelCepstrum(20, 0.35), x)
    mc2e(state)
    state = estimate(LinearCepstrum(20), x)
    mc2e(state)
    state = estimate(GeneralizedCepstrum(20, -0.01), x)
    @test_throws ErrorException mc2e(state)
    state = estimate(MelGeneralizedCepstrum(20, 0.35, -0.01), x)
    @test_throws ErrorException mc2e(state)
end

for order in [20, 25, 30]
    for α in [-0.554, -0.41, 0.0, 0.41, 0.544]
        for γ in [-1.0, -0.5, 0.0]
            println("mgc2b: testing with order=$order, α=$α, γ=$γ")
            test_mgc2b(order, α, γ)
        end
    end
end

let
    println("testing: mgc2b exceptions")
    srand(98765)
    x = rand(512)
    b = mgc2b(estimate(MelGeneralizedCepstrum(20, 0.41, -0.3), x))
    # mgc2b assumes input is not filter coef.
    @test_throws ArgumentError mgc2b(b)
end

for order in [20, 25, 30]
    for α in [-0.554, -0.41, 0.0, 0.41, 0.544]
        for γ in [-1.0, -0.5, 0.0]
            for fftlen in [256, 512]
                println("mgc2sp: testing with order=$order, α=$α, γ=$γ, fftlen=$fftlen")
                test_mgc2sp(order, α, γ, fftlen)
            end
        end
    end
end

# mgc2sp exceptions
let
    println("testing: mgc2sp exceptions")
    srand(98765)
    x = rand(512)

    state = mgc2b(estimate(MelGeneralizedCepstrum(20, 0.41, -0.1), x))
    @test_throws ArgumentError mgc2sp(state, 1024)
end

for order in 15:5:25
    for α in [0.0, 0.35]
        println("b2mc: testing with order=$order, α=$α")
        test_b2mc(order, α)
    end
end

# b2mc exceptions
let
    println("testing: b2mc exceptions")
    srand(98765)
    x = rand(512)

    state = estimate(MelCepstrum(20, 0.41), x)

    # b2mc should not accept mel-cepstrum
    @test_throws ArgumentError b2mc(state)
    @test_throws ArgumentError b2mc!(state)
end

for order in [20, 25, 30]
    for γ in [-1.0, -0.5, 0.0]
        println("gnorm and ignorm: testing with order=$order, γ=$γ")
        test_gnorm_and_ignorm(order, γ)
    end
end

let
    println("testing: gnorm and ignorm exceptions")
    srand(98765)
    x = rand(512)

    state1 = estimate( GeneralizedCepstrum(20, -0.1), x)
    state2 = gnorm(state1)
    # double normalization should raise error
    @test_throws ArgumentError gnorm(gnorm(state1))
    @test_throws ArgumentError ignorm(ignorm(state2))

    # should not accept filter coef.
    @test_throws ArgumentError gnorm(mgc2b(state1))
    @test_throws ArgumentError ignorm(mgc2b(state1))
end

for src_order in 15:5:25
    for src_α in [0.0, 0.35]
        for dst_order in 15:5:25
            for dst_α in [0.0, 0.35]
                println("freqt: testing with src_order=$src_order, src_α=$src_α, dst_order=$dst_order, dst_α=$dst_α")
                test_freqt(src_order, src_α, dst_order,dst_α )
            end
        end
    end
end

let
    println("testing: freqt exceptions")
    srand(98765)
    x = rand(512)

    # should not accept filter coefficients
    state1 = mc2b(estimate(MelCepstrum(20, 0.41), x))
    @test_throws ArgumentError freqt(state1, 20, 0.35)
end

for src_order in 15:5:25
    for src_α in [0.0, 0.35]
        for src_γ in [-1.0, -0.5, 0.0]
            for dst_order in 15:5:25
                for dst_γ in [-1.0, -0.5, 0.0]
                    println("gc2gc: testing with src_order=$src_order, src_α=$src_α, src_γ=$src_γ, dst_order=$dst_order, dst_γ=$dst_γ")
                    test_gc2gc(src_order,src_α, src_γ, dst_order, dst_γ)
                end
            end
        end
    end
end

for src_order in 15:5:25
    for src_α in [0.0, 0.35]
        for src_γ in [-1.0, -0.5, 0.0]
            for dst_order in 15:5:25
                for dst_α in [0.0, 0.35]
                    for dst_γ in [-1.0, -0.5, 0.0]
                        println("mgc2mgc: testing with src_order=$src_order, src_α, src_γ=$src_γ, dst_order=$dst_order, dst_α=$dst_α, dst_γ=$dst_γ")
                        test_mgc2mgc(src_order, src_α, src_γ, dst_order, dst_α, dst_γ)
                    end
                end
            end
        end
    end
end

let
    println("testing: mgc2mgc exceptions")
    srand(98765)
    x = rand(512)

    state = estimate(MelGeneralizedCepstrum(20, 0.41, -0.1), x)
    mgc2b!(state)

    # should not accept filter coefficients
    @test_throws ArgumentError mgc2mgc(state, 0.35, -1.0)
end

let
    println("testing: mat2mat functions")
    test_extend()
    test_inplace_extend()
end
