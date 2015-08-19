# This code is a julia translation of the following project:
# https://bitbucket.org/happyalu/mcep_alpha_calc

# mcepalpha computes appropriate α for a given sampling frequency.
function mcepalpha(fs::Real;
                   start::Float64=0.0,
                   stop::Float64=1.0,
                   step::Float64=0.001,
                   numpoints::Integer=1000)
    α_candidates = start:step:stop
    mel = melscale_vector(fs, numpoints)
    distances = [rms_distance(mel, warping_vector(α, numpoints)) for
                 α in α_candidates]
    return α_candidates[indmin(distances)]
end

function melscale_vector(fs::Real, len::Integer)
    step = (fs / 2.0) / len
    melscalev = 1000.0/log(2)*log(1 + step.*(1:len)./1000.0)
    return melscalev / melscalev[end]
end

function warping_vector(α::Float64, len::Integer)
    step = π / len
    ω = step .* (1:len)
    num = (1-α*α) * sin(ω)
    den = (1+α*α) * cos(ω) - 2*α
    warpfreq = atan(num./den)
    warpfreq[warpfreq .< 0] += π
    return warpfreq / warpfreq[end]
end

rms_distance(v1, v2) = sumabs2(v1 - v2) / length(v1)
