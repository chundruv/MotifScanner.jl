### Per-motif score calibration
#
# A fixed fraction of the maximum score (PrMax) is not comparable across motifs: the
# genome hit rate at PrMax >= 0.85 spans >1,000-fold across JASPAR vertebrate motifs and
# tracks information content. Here each motif gets its exact null score distribution
# under the zero-order background (dynamic programming over integer-discretised
# log-odds, as in FIMO/MOODS), giving
#   - a per-placement, single-strand p-value for any window score;
#   - an integer score threshold for a chosen p (used for calls and site counts);
#   - a probability-normalised score, CDF(score) = 1 - p, comparable to motifDiff "probNorm".

@inline baseidx(b) = b == DNA_A ? 1 : b == DNA_C ? 2 : b == DNA_G ? 3 : b == DNA_T ? 4 : 0

"""
    nullsf(ipbg, background) -> (minint, sf)

Survival function of the integer motif score for one random placement, with i.i.d.
bases drawn from `background`. `sf[k] = P(S >= minint + k - 1)`.
"""
function nullsf(ipbg::AbstractMatrix{<:Integer}, background)
    cmin = vec(minimum(ipbg, dims=1))
    cmax = vec(maximum(ipbg, dims=1))
    dist = [1.0]
    for j in axes(ipbg, 2)
        new = zeros(length(dist) + cmax[j] - cmin[j])
        for b in 1:4
            sh = ipbg[b, j] - cmin[j]
            pb = background[b]
            @inbounds for k in eachindex(dist)
                new[k + sh] += dist[k] * pb
            end
        end
        dist = new
    end
    sf = reverse(cumsum(reverse(dist)))
    sum(cmin), sf
end

"""
    calibrate(motif; background=background_human_zero_markov(), resolution=0.01, pthresh=1e-4)

Return `motif` (a NamedTuple from any loader) with calibration fields added:
`ipbg`, `ipbg_rc` (log-odds rounded to `resolution` bits), `minint`, `sf`, `resolution`,
`pthresh`, `thrint` (smallest integer score with p <= `pthresh`) and `thrscore`
(the same threshold in log2-odds units).
Use the same `background` as the log-odds were built with.
"""
function calibrate(motif; background=background_human_zero_markov(), resolution=0.01, pthresh=1e-4)
    ipbg = round.(Int, motif.pbg ./ resolution)
    minint, sf = nullsf(ipbg, background)
    k = findfirst(<=(pthresh), sf)
    thrint = isnothing(k) ? typemax(Int) : minint + k - 1
    base = (name=motif.name, id=motif.id, pwm=motif.pwm, pbg=motif.pbg,
            pbg_rc=motifpbg_rc(motif), maxscore=motifmaxscore(motif))
    merge(base, (ipbg=ipbg, ipbg_rc=rcm(ipbg), minint=minint, sf=sf, resolution=resolution,
                 pthresh=pthresh, thrint=thrint,
                 thrscore=thrint == typemax(Int) ? Inf : thrint * resolution))
end

calibrate(motifs::AbstractVector; kwargs...) = [calibrate(m; kwargs...) for m in motifs]

iscalibrated(m) = hasproperty(m, :sf) && hasproperty(m, :thrint)

"""
    motifpval(motif, s::Integer)

Single-strand, per-placement null p-value of an integer window score `s`
(1.0 for no valid placement, i.e. `s == typemin(Int)`).
"""
@inline function motifpval(m, s::Integer)
    s == typemin(Int) && return 1.0
    s < m.minint && return 1.0
    k = s - m.minint + 1
    k > length(m.sf) ? 0.0 : m.sf[k]
end

"""
    motifcalibration_table(motifs)

One row per calibrated motif: id, name, width, max score, threshold at `pthresh`,
the p-value of the best possible k-mer, and the PrMax equivalent of the threshold.
"""
function motifcalibration_table(motifs)
    DataFrame(MotifID=[m.id for m in motifs], MotifName=[m.name for m in motifs],
              Length=[size(m.pbg, 2) for m in motifs],
              MaxScore=[m.maxscore for m in motifs],
              PThresh=[m.pthresh for m in motifs],
              ThresholdScore=[m.thrscore for m in motifs],
              ThresholdPrMax=[m.thrscore / m.maxscore for m in motifs],
              BestKmerP=[m.sf[end] for m in motifs])
end
