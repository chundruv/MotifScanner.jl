using MotifScanner
using DataFrames

function parse_bool_arg(x::AbstractString)
    lx = lowercase(x)
    lx in ("1", "true", "t", "yes", "y")
end

function parse_vcf_chunk!(io, chunk_size::Int, split_multiallelic::Bool, drop_symbolic::Bool)
    chroms = String[]
    starts = Int[]
    stops = Int[]
    ids = String[]
    refs = String[]
    alts = String[]

    while length(starts) < chunk_size && !eof(io)
        line = readline(io)
        isempty(line) && continue
        startswith(line, "#") && continue

        fields = split(line, '\t')
        length(fields) < 5 && continue

        chrom = fields[1]
        pos = tryparse(Int, fields[2])
        isnothing(pos) && continue

        row_id = fields[3] == "." ? string(chrom, ":", pos, ":", fields[4], ">", fields[5]) : fields[3]
        ref = uppercase(fields[4])
        row_alts = split_multiallelic ? split(fields[5], ',') : [fields[5]]

        for alt_raw in row_alts
            alt = uppercase(strip(alt_raw))
            (alt == "." || isempty(alt)) && continue
            if drop_symbolic && (occursin('<', alt) || occursin('>', alt) || alt == "*")
                continue
            end

            push!(chroms, chrom)
            push!(starts, pos)
            push!(stops, pos + length(ref) - 1)
            push!(ids, row_id)
            push!(refs, ref)
            push!(alts, alt)

            if length(starts) >= chunk_size
                break
            end
        end
    end

    DataFrame(chrom=chroms, start=starts, stop=stops, ID=ids, ref=refs, alt=alts)
end

function write_scores!(io, df::DataFrame)
    for row in eachrow(df)
        fields = [
            row.MotifName,
            row.MotifID,
            string(row.RefMotSeqStart),
            string(row.AltMotSeqStart),
            string(row.RefMaxScore),
            string(row.RefStart),
            string(row.RefStop),
            row.RefStrand,
            string(row.RefSumScore),
            string(row.RefTotalMotifs),
            string(row.RefTotalMax),
            string(row.AltMaxScore),
            string(row.AltStart),
            string(row.AltStop),
            row.AltStrand,
            string(row.AltSumScore),
            string(row.AltTotalMotifs),
            string(row.AltTotalMax),
            string(row.RefPrMax),
            string(row.AltPrMax),
            string(row.LR_RefAlt),
            string(row.PR_RefAlt),
            string(row.RefSeq),
            string(row.AltSeq),
            row.ID,
        ]
        write(io, join(fields, '\t'))
        write(io, '\n')
    end
end

function main(args)
    if length(args) < 4
        error("Usage: julia --project scripts/score_vcf_jaspar.jl <vcf_or_vcfgz> <reference_fasta> <jaspar_transfac_txt> <output_tsv> [minprmax] [mindeltapr] [chunk_size] [split_multiallelic] [drop_symbolic]")
    end

    vcf_file = args[1]
    fasta_file = args[2]
    motif_file = args[3]
    out_file = args[4]

    minprmax = length(args) >= 5 ? parse(Float64, args[5]) : 0.85
    mindeltapr = length(args) >= 6 ? parse(Float64, args[6]) : 0.10
    chunk_size = length(args) >= 7 ? parse(Int, args[7]) : 5000
    split_multiallelic = length(args) >= 8 ? parse_bool_arg(args[8]) : true
    drop_symbolic = length(args) >= 9 ? parse_bool_arg(args[9]) : true

    println("Loading motifs from TRANSFAC: ", motif_file)
    _, motifs = loadtransfac(motif_file)
    println("Loaded motifs: ", length(motifs))
    println("Scoring with minprmax=", minprmax, " mindeltapr=", mindeltapr, " chunk_size=", chunk_size)

    open(vcf_file, "r") do _
    end

    input_io = endswith(lowercase(vcf_file), ".gz") ? open(`gzip -cd $(vcf_file)`, "r") : open(vcf_file, "r")

    total_variants = 0
    total_output_rows = 0
    chunk_index = 0

    open(out_file, "w") do out
        write(out, "MotifName\tMotifID\tRefMotSeqStart\tAltMotSeqStart\tRefMaxScore\tRefStart\tRefStop\tRefStrand\tRefSumScore\tRefTotalMotifs\tRefTotalMax\tAltMaxScore\tAltStart\tAltStop\tAltStrand\tAltSumScore\tAltTotalMotifs\tAltTotalMax\tRefPrMax\tAltPrMax\tLR_RefAlt\tPR_RefAlt\tRefSeq\tAltSeq\tID\n")

        while !eof(input_io)
            vartable = parse_vcf_chunk!(input_io, chunk_size, split_multiallelic, drop_symbolic)
            isempty(vartable) && continue

            chunk_index += 1
            total_variants += nrow(vartable)

            seqtable = loadrefseqs(vartable, fasta_file)
            seqtable[!, :ID] = vartable.ID

            scores = motifscanall(seqtable, motifs; minprmax=minprmax, mindeltapr=mindeltapr)
            total_output_rows += nrow(scores)

            write_scores!(out, scores)

            println("chunk=", chunk_index, " variants=", nrow(vartable), " cumulative_variants=", total_variants, " chunk_rows=", nrow(scores), " cumulative_rows=", total_output_rows)
        end
    end

    close(input_io)
    println("Done. total_variants=", total_variants, " total_output_rows=", total_output_rows, " output=", out_file)
end

main(ARGS)
