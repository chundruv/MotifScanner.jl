


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
    if :ID in names(vartable)
        seqs[!, :ID] = string.(vartable[!, :ID])
    elseif :id in names(vartable)
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
    for (c, s, e, r, a) in zip(chroms, starts, stops, refs, alts)
        start = s - mml
        stop  = e + mml
        # refseq = FASTA.extract(reader, DNAAlphabet{4}(), c, start:stop) 
        refseq = FASTA.extract(reader, c, start:stop) |> LongDNA{4}
        refind = (s:e) .- start .+ 1
        altind = first(refind) .+ (1:length(a)) .- 1
        
        leftind = 1:(first(refind) - 1)
        rightind = (last(refind) + 1):length(refseq)
        
        altseq = refseq[leftind]*a*refseq[rightind]
    
        @assert a == altseq[altind]

        push!(seqstarts, start)
        push!(seqstops, stop)
        push!(refinds, refind)
        push!(altinds, altind)
        push!(refseqs, refseq)
        push!(altseqs, altseq)
    end

    close(reader)
    DataFrame(seqstart=seqstarts, seqstop=seqstops, refind=refinds, altind=altinds, refseq=refseqs, altseq=altseqs)
end


function motifscanall(seqtable, motifs; minprmax=-Inf, mindeltapr=0.0)

    nrows = nrow(seqtable)
    dfs = Vector{DataFrame}(undef, nrows)

    Threads.@threads for idx in 1:nrows
        row = seqtable[idx, :]
        df = scanmots(row.refseq, row.altseq, row.refind, row.altind, motifs; minprmax=minprmax, mindeltapr=mindeltapr)
        if nrow(df) > 0
            df[!, :ID] .= row.ID
        end
        dfs[idx] = df
    end
    vcat(dfs...)
end
## for variants only

function motifscanseq(seq, ind, motif)
    n = size(motif.pwm, 2)
    start = max(first(ind) - n + 1, 1)
    stop  = min(last(ind)  + n - 1, length(seq))
    seq[start:stop], start
end



function scanmots(refseq, altseq, refind, altind, motifs; minprmax=-Inf, mindeltapr=0.0)

    nm = length(motifs)
    motif_names = Vector{String}(undef, nm)
    motif_ids = Vector{String}(undef, nm)
    ref_mot_seq_starts = Vector{Int}(undef, nm)
    alt_mot_seq_starts = Vector{Int}(undef, nm)
    ref_max_scores = Vector{Float64}(undef, nm)
    ref_starts = Vector{Int}(undef, nm)
    ref_stops = Vector{Int}(undef, nm)
    ref_strands = Vector{String}(undef, nm)
    ref_sum_scores = Vector{Float64}(undef, nm)
    ref_total_motifs = Vector{Int}(undef, nm)
    ref_total_max = Vector{Int}(undef, nm)
    alt_max_scores = Vector{Float64}(undef, nm)
    alt_starts = Vector{Int}(undef, nm)
    alt_stops = Vector{Int}(undef, nm)
    alt_strands = Vector{String}(undef, nm)
    alt_sum_scores = Vector{Float64}(undef, nm)
    alt_total_motifs = Vector{Int}(undef, nm)
    alt_total_max = Vector{Int}(undef, nm)
    ref_prmaxs = Vector{Float64}(undef, nm)
    alt_prmaxs = Vector{Float64}(undef, nm)
    lr_refalts = Vector{Float64}(undef, nm)
    pr_refalts = Vector{Float64}(undef, nm)
    ref_seqs = Vector{LongDNA{4}}(undef, nm)
    alt_seqs = Vector{LongDNA{4}}(undef, nm)

    k = 0
    @showprogress for m in motifs
        maxscore = motifmaxscore(m)
        refmseq, refstart = MotifScanner.motifscanseq(refseq, refind, m)
        altmseq, altstart = MotifScanner.motifscanseq(altseq, altind, m)
        
        refres = scansummax(refmseq, m)
        altres = scansummax(altmseq, m)

        ref_prmax = first(refres)/maxscore
        alt_prmax = first(altres)/maxscore
        lr_refalt = first(refres) - first(altres)
        pr_refalt = ref_prmax - alt_prmax

        # Keep motifs where normalized max score and normalized delta pass filters.
        if max(ref_prmax, alt_prmax) < minprmax || abs(pr_refalt) < mindeltapr
            continue
        end

        k += 1

        motif_names[k] = m.name
        motif_ids[k] = m.id
        ref_mot_seq_starts[k] = refstart
        alt_mot_seq_starts[k] = altstart

        ref_max_scores[k] = refres[1]
        ref_starts[k] = refres[2]
        ref_stops[k] = refres[3]
        ref_strands[k] = refres[4]
        ref_sum_scores[k] = refres[5]
        ref_total_motifs[k] = refres[6]
        ref_total_max[k] = refres[7]

        alt_max_scores[k] = altres[1]
        alt_starts[k] = altres[2]
        alt_stops[k] = altres[3]
        alt_strands[k] = altres[4]
        alt_sum_scores[k] = altres[5]
        alt_total_motifs[k] = altres[6]
        alt_total_max[k] = altres[7]

        ref_prmaxs[k] = ref_prmax
        alt_prmaxs[k] = alt_prmax
        lr_refalts[k] = lr_refalt
        pr_refalts[k] = pr_refalt
        ref_seqs[k] = refmseq
        alt_seqs[k] = altmseq
    end

    resize!(motif_names, k)
    resize!(motif_ids, k)
    resize!(ref_mot_seq_starts, k)
    resize!(alt_mot_seq_starts, k)
    resize!(ref_max_scores, k)
    resize!(ref_starts, k)
    resize!(ref_stops, k)
    resize!(ref_strands, k)
    resize!(ref_sum_scores, k)
    resize!(ref_total_motifs, k)
    resize!(ref_total_max, k)
    resize!(alt_max_scores, k)
    resize!(alt_starts, k)
    resize!(alt_stops, k)
    resize!(alt_strands, k)
    resize!(alt_sum_scores, k)
    resize!(alt_total_motifs, k)
    resize!(alt_total_max, k)
    resize!(ref_prmaxs, k)
    resize!(alt_prmaxs, k)
    resize!(lr_refalts, k)
    resize!(pr_refalts, k)
    resize!(ref_seqs, k)
    resize!(alt_seqs, k)

    DataFrame(
        MotifName=motif_names,
        MotifID=motif_ids,
        RefMotSeqStart=ref_mot_seq_starts,
        AltMotSeqStart=alt_mot_seq_starts,
        RefMaxScore=ref_max_scores,
        RefStart=ref_starts,
        RefStop=ref_stops,
        RefStrand=ref_strands,
        RefSumScore=ref_sum_scores,
        RefTotalMotifs=ref_total_motifs,
        RefTotalMax=ref_total_max,
        AltMaxScore=alt_max_scores,
        AltStart=alt_starts,
        AltStop=alt_stops,
        AltStrand=alt_strands,
        AltSumScore=alt_sum_scores,
        AltTotalMotifs=alt_total_motifs,
        AltTotalMax=alt_total_max,
        RefPrMax=ref_prmaxs,
        AltPrMax=alt_prmaxs,
        LR_RefAlt=lr_refalts,
        PR_RefAlt=pr_refalts,
        RefSeq=ref_seqs,
        AltSeq=alt_seqs,
    )
    
end