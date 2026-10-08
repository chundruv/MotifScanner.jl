## Calibrated motif scoring of VCF variants (JASPAR or HOCOMOCO TRANSFAC files).
##
## julia --project -t 8 scripts/score_vcf_motifs.jl \
##     --vcf dnvs.vcf.gz --fasta GRCh38.fa --motifs H14CORE_transfac_format.txt.gz \
##     --out scores.tsv [--taxa vertebrates] [--tfs expressed_tfs.txt] \
##     [--pthresh 1e-4] [--maxpval 1e-4] [--mindlog10p 0] [--minprmax -Inf] [--mindeltapr 0] \
##     [--flank 40] [--pseudocount 0.01] [--resolution 0.01] [--chunk 5000] \
##     [--split-multiallelic true] [--drop-symbolic true]
##
## Outputs
##   <out>                one row per variant x motif passing the filters (see ?MotifScanner.scanmots)
##   <out>.motifs.tsv     motif metadata + calibration (threshold score, PrMax equivalent, best k-mer p)
##   <out>.variants.tsv   every parsed variant: RefMatch flag and number of rows written, giving
##                        the denominator for burden tests (variants without rows were scored, no hit)
##
## Defaults: --taxa vertebrates is applied only when the file has tax_group annotations (JASPAR);
## --maxpval defaults to --pthresh, so a row is written only if at least one allele has a site.

using MotifScanner
using DataFrames

const DEFAULTS = Dict(
    "taxa" => "vertebrates", "tfs" => "", "pthresh" => "1e-4", "maxpval" => "",
    "mindlog10p" => "0", "minprmax" => "-Inf", "mindeltapr" => "0", "flank" => "40",
    "pseudocount" => "0.01", "resolution" => "0.01", "chunk" => "5000",
    "split-multiallelic" => "true", "drop-symbolic" => "true")
const REQUIRED = ("vcf", "fasta", "motifs", "out")

parse_bool_arg(x::AbstractString) = lowercase(x) in ("1", "true", "t", "yes", "y")

function parse_args(args)
    opts = copy(DEFAULTS)
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") || error("unexpected argument $a")
        k = a[3:end]
        (haskey(opts, k) || k in REQUIRED) || error("unknown option --$k")
        i + 1 <= length(args) || error("missing value for --$k")
        opts[k] = args[i + 1]
        i += 2
    end
    for k in REQUIRED
        haskey(opts, k) || error("--$k is required")
    end
    isempty(opts["maxpval"]) && (opts["maxpval"] = opts["pthresh"])
    opts
end

function parse_vcf_chunk!(io, chunk_size::Int, split_multiallelic::Bool, drop_symbolic::Bool)
    chroms = String[]; starts = Int[]; stops = Int[]; ids = String[]; refs = String[]; alts = String[]
    while length(starts) < chunk_size && !eof(io)
        line = readline(io)
        (isempty(line) || startswith(line, "#")) && continue
        fields = split(line, '\t')
        length(fields) < 5 && continue
        chrom = fields[1]
        pos = tryparse(Int, fields[2])
        isnothing(pos) && continue
        ref = uppercase(fields[4])
        row_alts = split_multiallelic ? split(fields[5], ',') : [fields[5]]
        for alt_raw in row_alts
            alt = uppercase(strip(alt_raw))
            (alt == "." || isempty(alt)) && continue
            drop_symbolic && (occursin('<', alt) || occursin('>', alt) || alt == "*") && continue
            occursin(r"^[ACGTN]+$", ref) && occursin(r"^[ACGTN]+$", alt) || continue
            # one ID per ALT allele so multiallelic records stay distinguishable
            row_id = fields[3] == "." ? string(chrom, ":", pos, ":", ref, ">", alt) : String(fields[3])
            push!(chroms, chrom); push!(starts, pos); push!(stops, pos + length(ref) - 1)
            push!(ids, row_id); push!(refs, ref); push!(alts, alt)
        end
    end
    DataFrame(chrom=chroms, start=starts, stop=stops, ID=ids, ref=refs, alt=alts)
end

function write_tsv!(io, df::DataFrame; header=false)
    header && (write(io, join(names(df), '\t')); write(io, '\n'))
    for row in eachrow(df)
        write(io, join((ismissing(x) ? "NA" : string(x) for x in row), '\t'))
        write(io, '\n')
    end
end

function main(args)
    o = parse_args(args)
    pthresh = parse(Float64, o["pthresh"])
    background = MotifScanner.background_human_zero_markov()
    taxa = isempty(o["taxa"]) || lowercase(o["taxa"]) == "all" ? nothing : split(o["taxa"], ',')
    tfs = isempty(o["tfs"]) ? nothing : filter(!isempty, strip.(readlines(o["tfs"])))

    println("Loading motifs: ", o["motifs"])
    meta, motifs = loadtransfac(o["motifs"], parse(Float64, o["pseudocount"]), background; taxa=taxa, tfs=tfs)
    isempty(motifs) && error("no motifs left after filtering")
    motifs = calibrate(motifs; background=background, resolution=parse(Float64, o["resolution"]), pthresh=pthresh)
    cal = motifcalibration_table(motifs)
    open(o["out"] * ".motifs.tsv", "w") do io
        write_tsv!(io, hcat(cal, select(meta, Not([:MotifName, :MotifID, :Length])); makeunique=true); header=true)
    end
    nuncallable = sum(cal.BestKmerP .> pthresh)
    nuncallable > 0 && println(stderr, "NOTE: ", nuncallable, " motifs cannot reach p <= ", pthresh,
        " (best k-mer p > pthresh; mostly short motifs). They are still scored (p-values, ProbNorm) but never called.")
    println("Motifs: ", length(motifs), "  pthresh=", pthresh, "  maxpval=", o["maxpval"],
            "  Julia threads: ", Threads.nthreads())

    split_ma = parse_bool_arg(o["split-multiallelic"]); drop_sym = parse_bool_arg(o["drop-symbolic"])
    chunk = parse(Int, o["chunk"]); flank = parse(Int, o["flank"])
    kw = (pthresh=pthresh, maxpval=parse(Float64, o["maxpval"]), mindlog10p=parse(Float64, o["mindlog10p"]),
          minprmax=parse(Float64, o["minprmax"]), mindeltapr=parse(Float64, o["mindeltapr"]))

    input_io = endswith(lowercase(o["vcf"]), ".gz") ? open(`gzip -cd $(o["vcf"])`, "r") : open(o["vcf"], "r")
    nvar = 0; nrows = 0; nmismatch = 0; ichunk = 0
    open(o["out"], "w") do out
        open(o["out"] * ".variants.tsv", "w") do vout
            write(vout, "ID\tchrom\tpos\tref\talt\tRefMatch\tRowsWritten\tCalls_loss\tCalls_gain\n")
            first_chunk = true
            while !eof(input_io)
                vt = parse_vcf_chunk!(input_io, chunk, split_ma, drop_sym)
                isempty(vt) && continue
                ichunk += 1; nvar += nrow(vt)
                st = loadrefseqs(vt, o["fasta"]; mml=flank)
                st[!, :ID] = vt.ID
                nmismatch += sum(.!st.RefMatch)
                sc = motifscanall(st, motifs; kw...)
                write_tsv!(out, sc; header=first_chunk); first_chunk = false
                nrows += nrow(sc)
                cnt = Dict{String,NTuple{3,Int}}()
                for (id, c) in zip(sc.ID, sc.Call)
                    t = get(cnt, id, (0, 0, 0))
                    cnt[id] = (t[1] + 1, t[2] + (c == "loss"), t[3] + (c == "gain"))
                end
                for (r, m) in zip(eachrow(vt), st.RefMatch)
                    t = get(cnt, r.ID, (0, 0, 0))
                    write(vout, join((r.ID, r.chrom, r.start, r.ref, r.alt, m, t...), '\t'), '\n')
                end
                println("chunk=", ichunk, " variants=", nrow(vt), " cumulative=", nvar,
                        " ref_mismatch=", nmismatch, " rows=", nrows)
            end
        end
    end
    close(input_io)
    nmismatch > 0 && println(stderr, "WARNING: ", nmismatch, " variants had REF not matching the FASTA and were not scored")
    println("Done. variants=", nvar, " rows=", nrows, " output=", o["out"])
end

main(ARGS)
