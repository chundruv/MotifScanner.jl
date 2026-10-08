


### variant table
# chrom start stop id type = {var, ins, del} index refseq altseq

"""
    loadvcf(vcf_file; split_multiallelic=true, drop_symbolic=true)

Load variants from a VCF file into a DataFrame with columns:
`chrom`, `start`, `stop`, `ID`, `ref`, `alt`.

`start` and `stop` are 1-based, inclusive genomic coordinates for the REF allele.
When `split_multiallelic=true`, records with multiple ALT alleles are expanded into
one row per ALT allele.
"""
function loadvcf(vcf_file::AbstractString; split_multiallelic=true, drop_symbolic=true)
    chroms = String[]
    starts = Int[]
    stops = Int[]
    ids = String[]
    refs = String[]
    alts_out = String[]

    open(vcf_file, "r") do io
        for line in eachline(io)
            isempty(line) && continue
            startswith(line, "#") && continue

            fields = split(line, '\t')
            length(fields) < 5 && continue

            chrom = fields[1]
            pos = tryparse(Int, fields[2])
            isnothing(pos) && continue

            row_id = fields[3] == "." ? string(chrom, ":", pos, ":", fields[4], ">", fields[5]) : fields[3]
            ref = uppercase(fields[4])

            alts = split_multiallelic ? split(fields[5], ',') : [fields[5]]
            for alt_raw in alts
                alt = uppercase(strip(alt_raw))
                (alt == "." || isempty(alt)) && continue
                if drop_symbolic && (occursin('<', alt) || occursin('>', alt) || alt == "*")
                    continue
                end

                start = pos
                stop = pos + length(ref) - 1
                push!(chroms, chrom)
                push!(starts, start)
                push!(stops, stop)
                push!(ids, row_id)
                push!(refs, ref)
                push!(alts_out, alt)
            end
        end
    end

    DataFrame(chrom=chroms, start=starts, stop=stops, ID=ids, ref=refs, alt=alts_out)
end

function loadrefseqs(vartable, fastafile; mml=40, reffield=:ref, altfield=:alt)
    seqs = loadrefseqs(vartable.chrom, vartable.start, vartable.stop, vartable[!, reffield], vartable[!, altfield], fastafile, mml)
    if :ID in propertynames(vartable)
        seqs[!, :ID] = string.(vartable[!, :ID])
    elseif :id in propertynames(vartable)
        seqs[!, :ID] = string.(vartable[!, :id])
    end
    seqs
end

function loadrefseqs(vcf_file::AbstractString, fastafile; mml=40, split_multiallelic=true, drop_symbolic=true)
    vartable = loadvcf(vcf_file; split_multiallelic=split_multiallelic, drop_symbolic=drop_symbolic)
    loadrefseqs(vartable, fastafile; mml=mml)
end

#### currently only works for variants where Ref is single base
# function loadrefseqs(chroms, starts, stops, refs, alts, file, mml=40)
#     reader = open(FASTA.Reader, file, index=string(file, ".fai"))
    
#     cr = 0
#     ca = 0
#     seqs = DataFrame(seqstart=Int[], seqstop=Int[], Kind=String[], refind=UnitRange{Int}[], altind=UnitRange{Int}[], refseq=LongSequence{DNAAlphabet{4}}[], altseq=LongSequence{DNAAlphabet{4}}[])
#     for (c, s, e, r, a) in zip(chroms, starts, stops, refs, alts)
#         start = s - mml
#         stop  = e + mml
#         refseq = FASTA.extract(reader, DNAAlphabet{4}(), c, start:stop) 
#         refind = (s:e) .- start .+ 1
#         altind = first(refind) .+ (1:length(a)) .- 1
#         leftind = 1:(first(refind) - 1)
#         rightind = (last(refind) + 1):length(refseq)
#         seqa = LongDNASeq(a)
#         altseq = refseq[leftind]*seqa*refseq[rightind]
#         #@assert length(refseq) == length(altseq)


#         @assert seqa == altseq[altind]
#         if string(refseq[refind]) == r
#             cr += 1
#         elseif string(refseq[refind]) == a
#             ca += 1
#         else
#             @show c, s, e, r, a, refind, altind refseq
#             break
#         end

#         if length(r) == length(a) == 1
#             kind = "Var"
#         else
#             kind = "Indel"
#         end
#         push!(seqs, (start, stop, kind, refind, altind, refseq, altseq))
#     end
    
#     close(reader)
#     seqs
# end

loadrefseqs(chroms, starts, stops, refs::Vector{<:AbstractString}, alts::Vector{<:AbstractString}, file, mml=40) = loadrefseqs(chroms, starts, stops, LongDNA{4}.(refs), LongDNA{4}.(alts), file, mml)

function loadrefseqs(chroms, starts, stops, refs, alts, file, mml=40)
    reader = open(FASTA.Reader, file, index=string(file, ".fai"))

    seqstarts = Int[]
    seqstops = Int[]
    refinds = UnitRange{Int}[]
    altinds = UnitRange{Int}[]
    refseqs = LongSequence{DNAAlphabet{4}}[]
    altseqs = LongSequence{DNAAlphabet{4}}[]
    refmatch = Bool[]
    for (c, s, e, r, a) in zip(chroms, starts, stops, refs, alts)
        start = max(s - mml, 1)
        stop  = e + mml
        # refseq = FASTA.extract(reader, DNAAlphabet{4}(), c, start:stop) 
        refseq = FASTA.extract(reader, c, start:stop) |> LongDNA{4}
        refind = (s:e) .- start .+ 1
        altind = first(refind) .+ (1:length(a)) .- 1
        
        leftind = 1:(first(refind) - 1)
        rightind = (last(refind) + 1):length(refseq)
        
        altseq = refseq[leftind]*a*refseq[rightind]
    
        @assert a == altseq[altind]
        # REF must match the reference genome (build / normalisation check)
        push!(refmatch, refseq[refind] == r)

        push!(seqstarts, start)
        push!(seqstops, stop)
        push!(refinds, refind)
        push!(altinds, altind)
        push!(refseqs, refseq)
        push!(altseqs, altseq)
    end

    close(reader)
    DataFrame(seqstart=seqstarts, seqstop=seqstops, refind=refinds, altind=altinds, refseq=refseqs, altseq=altseqs, RefMatch=refmatch)
end


"""
    motifscanall(seqtable, motifs; pthresh=1e-4, maxpval=1.0, mindlog10p=0.0,
                 minprmax=-Inf, mindeltapr=0.0, background=background_human_zero_markov(),
                 resolution=0.01)

Score every variant in `seqtable` (from `loadrefseqs`) against every motif. Motifs are
calibrated once here if they are not already (see `calibrate`). Rows with
`RefMatch == false` (REF does not match the FASTA) are skipped. A row is written when
min(RefPval, AltPval) <= `maxpval`, |Log10P_RefAlt| >= `mindlog10p`, and the legacy
`minprmax`/`mindeltapr` filters pass. Use `maxpval=pthresh` to keep only rows where at
least one allele has a calibrated site.
"""
function motifscanall(seqtable, motifs; pthresh=1e-4, maxpval=1.0, mindlog10p=0.0,
                      minprmax=-Inf, mindeltapr=0.0,
                      background=background_human_zero_markov(), resolution=0.01)

    cmotifs = [iscalibrated(m) && m.pthresh == pthresh ? m :
               calibrate(m; background=background, resolution=resolution, pthresh=pthresh)
               for m in motifs]
    nrows = nrow(seqtable)
    dfs = Vector{DataFrame}(undef, nrows)
    has_id = :ID in propertynames(seqtable)
    has_match = :RefMatch in propertynames(seqtable)

    Threads.@threads for idx in 1:nrows
        row = seqtable[idx, :]
        if has_match && !row.RefMatch
            dfs[idx] = scanmots(row.refseq, row.altseq, row.refind, row.altind, cmotifs[1:0])
        else
            dfs[idx] = scanmots(row.refseq, row.altseq, row.refind, row.altind, cmotifs;
                                maxpval=maxpval, mindlog10p=mindlog10p,
                                minprmax=minprmax, mindeltapr=mindeltapr)
        end
        if has_id
            dfs[idx][!, :ID] = fill(row.ID, nrow(dfs[idx]))
        end
    end
    vcat(dfs...)
end

## for variants only

function motifscanseq(seq, ind, motif)
    n = size(motif.pwm, 2)
    start = max(first(ind) - n + 1, 1)
    stop  = min(last(ind)  + n - 1, length(seq))
    start:stop
end

"""
    scanwindow(seq, motif, ind) -> NamedTuple

Scan the placements that start in `first(ind):(last(ind) - w + 1)` on both strands,
using a calibrated motif. Placements containing a non-ACGT base are excluded.
Returns the best placement (same tie rule as `scansummax`: + strand wins only if strictly
better), its integer score, the number of valid placements, and the number and score sum
of placements at or above the calibrated threshold (`motif.thrint`).
"""
function scanwindow(seq, motif, ind::UnitRange{Int})
    pbg = motif.pbg; rot = motifpbg_rc(motif)
    ip = motif.ipbg; irc = motif.ipbg_rc
    thrint = motif.thrint
    maxscore = motifmaxscore(motif)
    m = size(pbg, 2)
    offset = first(ind)
    lastpos = last(ind) - m + 1

    nvalid = 0; nsig = 0; nmax = 0; sumsig = 0.0
    fm = -Inf; fpos = 0; fint = typemin(Int)
    rm = -Inf; rpos = 0; rint = typemin(Int)

    @inbounds for pos in offset:lastpos
        f = 0.0; r = 0.0; fi = 0; ri = 0; valid = true
        for j in 1:m
            b = baseidx(seq[pos + j - 1])
            if b == 0
                valid = false
                break
            end
            f += pbg[b, j]; r += rot[b, j]
            fi += ip[b, j]; ri += irc[b, j]
        end
        valid || continue
        nvalid += 1
        if fi >= thrint
            nsig += 1; sumsig += f
        end
        if ri >= thrint
            nsig += 1; sumsig += r
        end
        f == maxscore && (nmax += 1)
        r == maxscore && (nmax += 1)
        if f > fm
            fm = f; fpos = pos; fint = fi
        end
        if r > rm
            rm = r; rpos = pos; rint = ri
        end
    end

    if nvalid == 0
        return (score=-Inf, intscore=typemin(Int), start=0, stop=0, strand=".",
                sumscore=0.0, nsig=0, nmax=0, nvalid=0)
    end
    if fm > rm
        s, p, si, strand = fm, fpos, fint, "+"
    else
        s, p, si, strand = rm, rpos, rint, "-"
    end
    local_start = p - offset + 1
    (score=s, intscore=si, start=local_start, stop=local_start + m - 1, strand=strand,
     sumscore=sumsig / maxscore, nsig=nsig, nmax=nmax, nvalid=nvalid)
end

callclass(refhit::Bool, althit::Bool) =
    refhit ? (althit ? "retained" : "loss") : (althit ? "gain" : "none")

"""
    scanmots(refseq, altseq, refind, altind, motifs; maxpval=1.0, mindlog10p=0.0,
             minprmax=-Inf, mindeltapr=0.0)

Score one variant against calibrated motifs (see `calibrate`). Columns:

- Legacy: MotifName, MotifID, RefMotSeqStart, AltMotSeqStart, Ref/Alt MaxScore, Start,
  Stop, Strand, SumScore, TotalMotifs, TotalMax, RefPrMax, AltPrMax, LR_RefAlt,
  PR_RefAlt, RefSeq, AltSeq.
  `TotalMotifs`/`SumScore` now count placements at the calibrated threshold
  (p <= pthresh), not every placement with score > 0.
- Calibrated: RefPval/AltPval (single-strand, per-placement null p of the best placement),
  RefProbNorm/AltProbNorm (= 1 - p, motifDiff "probNorm"-style),
  ProbNorm_AltRef (alt - ref; negative = loss, motifDiff sign),
  Log10P_RefAlt (log10 AltPval - log10 RefPval; positive = loss),
  Call (loss | gain | retained | none at `pthresh`), ThresholdScore, PThresh,
  MotifLen, RefLen, AltLen, RefNValid, AltNValid.

Indel note: the longer allele is scanned over more placements, so its maximum is
biased upwards. Compare loss/gain rates only within variant classes / indel lengths.
"""
function scanmots(refseq, altseq, refind, altind, motifs; maxpval=1.0, mindlog10p=0.0,
                  minprmax=-Inf, mindeltapr=0.0)

    out = (MotifName=String[], MotifID=String[], RefMotSeqStart=Int[], AltMotSeqStart=Int[],
           RefMaxScore=Float64[], RefStart=Int[], RefStop=Int[], RefStrand=String[],
           RefSumScore=Float64[], RefTotalMotifs=Int[], RefTotalMax=Int[],
           AltMaxScore=Float64[], AltStart=Int[], AltStop=Int[], AltStrand=String[],
           AltSumScore=Float64[], AltTotalMotifs=Int[], AltTotalMax=Int[],
           RefPrMax=Float64[], AltPrMax=Float64[], LR_RefAlt=Float64[], PR_RefAlt=Float64[],
           RefPval=Float64[], AltPval=Float64[], RefProbNorm=Float64[], AltProbNorm=Float64[],
           ProbNorm_AltRef=Float64[], Log10P_RefAlt=Float64[], Call=String[],
           ThresholdScore=Float64[], PThresh=Float64[], MotifLen=Int[], RefLen=Int[],
           AltLen=Int[], RefNValid=Int[], AltNValid=Int[],
           RefSeq=LongDNA{4}[], AltSeq=LongDNA{4}[])

    for m0 in motifs
        # back-compat: uncalibrated motifs are calibrated here with defaults (slow per call;
        # motifscanall / calibrate(motifs) do this once up front)
        m = iscalibrated(m0) ? m0 : calibrate(m0)
        maxscore = motifmaxscore(m)
        refscanind = motifscanseq(refseq, refind, m)
        altscanind = motifscanseq(altseq, altind, m)

        rr = scanwindow(refseq, m, refscanind)
        ar = scanwindow(altseq, m, altscanind)

        ref_prmax = rr.score / maxscore
        alt_prmax = ar.score / maxscore
        pr_refalt = ref_prmax - alt_prmax
        refp = motifpval(m, rr.intscore)
        altp = motifpval(m, ar.intscore)
        dlog10p = log10(altp) - log10(refp)

        if min(refp, altp) > maxpval || abs(dlog10p) < mindlog10p ||
           max(ref_prmax, alt_prmax) < minprmax || !(abs(pr_refalt) >= mindeltapr)
            continue
        end

        push!(out.MotifName, m.name); push!(out.MotifID, m.id)
        push!(out.RefMotSeqStart, first(refscanind)); push!(out.AltMotSeqStart, first(altscanind))
        push!(out.RefMaxScore, rr.score); push!(out.RefStart, rr.start); push!(out.RefStop, rr.stop)
        push!(out.RefStrand, rr.strand); push!(out.RefSumScore, rr.sumscore)
        push!(out.RefTotalMotifs, rr.nsig); push!(out.RefTotalMax, rr.nmax)
        push!(out.AltMaxScore, ar.score); push!(out.AltStart, ar.start); push!(out.AltStop, ar.stop)
        push!(out.AltStrand, ar.strand); push!(out.AltSumScore, ar.sumscore)
        push!(out.AltTotalMotifs, ar.nsig); push!(out.AltTotalMax, ar.nmax)
        push!(out.RefPrMax, ref_prmax); push!(out.AltPrMax, alt_prmax)
        push!(out.LR_RefAlt, rr.score - ar.score); push!(out.PR_RefAlt, pr_refalt)
        push!(out.RefPval, refp); push!(out.AltPval, altp)
        push!(out.RefProbNorm, 1 - refp); push!(out.AltProbNorm, 1 - altp)
        push!(out.ProbNorm_AltRef, refp - altp)
        push!(out.Log10P_RefAlt, dlog10p)
        push!(out.Call, callclass(refp <= m.pthresh, altp <= m.pthresh))
        push!(out.ThresholdScore, m.thrscore); push!(out.PThresh, m.pthresh)
        push!(out.MotifLen, size(m.pbg, 2)); push!(out.RefLen, length(refind)); push!(out.AltLen, length(altind))
        push!(out.RefNValid, rr.nvalid); push!(out.AltNValid, ar.nvalid)
        push!(out.RefSeq, refseq[refscanind]); push!(out.AltSeq, altseq[altscanind])
    end
    DataFrame(out)
end


### legacy scorer (fraction-of-max threshold); kept for back-compatibility
function scansummax(seq, motif, ind::UnitRange{Int}, thr=0)
    maxscore = motifmaxscore(motif)
    thrscore = thr * maxscore
    rot = motifpbg_rc(motif)
    pbg = motif.pbg
    m = size(pbg, 2)

    offset = first(ind)
    lastpos = last(ind) - m + 1
    @assert lastpos >= offset

    totalmotifs = 0
    totalmax = 0
    sumscore_num = 0.0

    fm = -Inf
    fi = offset
    rm = -Inf
    ri = offset

    @inbounds for pos in offset:lastpos
        f = 0.0
        r = 0.0
        for j in 1:m
            base = seq[pos + j - 1]
            if base == DNA_A
                f += pbg[1, j]
                r += rot[1, j]
            elseif base == DNA_C
                f += pbg[2, j]
                r += rot[2, j]
            elseif base == DNA_G
                f += pbg[3, j]
                r += rot[3, j]
            elseif base == DNA_T
                f += pbg[4, j]
                r += rot[4, j]
            end
        end

        if f > thrscore
            sumscore_num += f
            totalmotifs += 1
        end
        if r > thrscore
            sumscore_num += r
            totalmotifs += 1
        end

        if f == maxscore
            totalmax += 1
        end
        if r == maxscore
            totalmax += 1
        end

        if f > fm
            fm = f
            fi = pos
        end
        if r > rm
            rm = r
            ri = pos
        end
    end

    sumscore = sumscore_num / maxscore

    if fm > rm
        local_start = fi - offset + 1
        local_stop = local_start + m - 1
        return fm, local_start, local_stop, "+", sumscore, totalmotifs, totalmax
    else
        local_start = ri - offset + 1
        local_stop = local_start + m - 1
        return rm, local_start, local_stop, "-", sumscore, totalmotifs, totalmax
    end
end
