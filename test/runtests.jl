using MotifScanner
using Test
using DataFrames
using BioSequences

@testset "MotifScanner.jl" begin
    @testset "loadvcf" begin
        vcf = joinpath(mktempdir(), "vars.vcf")
        open(vcf, "w") do io
            write(io, "##fileformat=VCFv4.2\n")
            write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            write(io, "chr1\t10\trs1\tA\tG\t.\tPASS\t.\n")
            write(io, "chr2\t20\t.\tAT\tA,AG\t.\tPASS\t.\n")
            write(io, "chr3\t30\trs3\tC\t<DEL>\t.\tPASS\t.\n")
        end

        vars = loadvcf(vcf)
        @test nrow(vars) == 3
        @test names(vars) == ["chrom", "start", "stop", "ID", "ref", "alt"]

        @test vars[1, :chrom] == "chr1"
        @test vars[1, :start] == 10
        @test vars[1, :stop] == 10
        @test vars[1, :ID] == "rs1"
        @test vars[1, :ref] == "A"
        @test vars[1, :alt] == "G"

        @test vars[2, :chrom] == "chr2"
        @test vars[2, :start] == 20
        @test vars[2, :stop] == 21
        @test vars[2, :ID] == "chr2:20:AT>A,AG"
        @test vars[2, :ref] == "AT"
        @test vars[2, :alt] == "A"

        @test vars[3, :chrom] == "chr2"
        @test vars[3, :start] == 20
        @test vars[3, :stop] == 21
        @test vars[3, :ref] == "AT"
        @test vars[3, :alt] == "AG"

        vars_no_split = loadvcf(vcf; split_multiallelic=false)
        @test nrow(vars_no_split) == 2
        @test vars_no_split[2, :alt] == "A,AG"

        vars_keep_symbolic = loadvcf(vcf; drop_symbolic=false)
        @test nrow(vars_keep_symbolic) == 4
        @test vars_keep_symbolic[4, :alt] == "<DEL>"
    end

    @testset "normalized score threshold filtering" begin
        pbg = [3.0 3.0 3.0;
               0.0 0.0 0.0;
               0.0 0.0 0.0;
               0.0 0.0 0.0]
        motif = (name="constant", id="const1", pwm=zeros(4, 3), pbg=pbg, pbg_rc=MotifScanner.rcm(pbg), maxscore=9.0)

        refseq = LongDNA{4}("AAAAAA")
        altseq = LongDNA{4}("CCCCCC")
        refind = 3:3
        altind = 3:3

        df_keep = MotifScanner.scanmots(refseq, altseq, refind, altind, [motif]; minprmax=0.99)
        @test nrow(df_keep) == 1

        df_drop = MotifScanner.scanmots(refseq, altseq, refind, altind, [motif]; minprmax=1.01)
        @test nrow(df_drop) == 0

        df_keep_delta = MotifScanner.scanmots(refseq, altseq, refind, altind, [motif]; minprmax=0.5, mindeltapr=0.5)
        @test nrow(df_keep_delta) == 1

        df_drop_delta = MotifScanner.scanmots(refseq, altseq, refind, altind, [motif]; minprmax=0.5, mindeltapr=1.01)
        @test nrow(df_drop_delta) == 0

        seqtable = DataFrame(ID=["v1"], refseq=[refseq], altseq=[altseq], refind=[refind], altind=[altind])
        all_keep = motifscanall(seqtable, [motif]; minprmax=0.99)
        @test nrow(all_keep) == 1

        all_drop = motifscanall(seqtable, [motif]; minprmax=1.01)
        @test nrow(all_drop) == 0
    end

    @testset "loadtransfac: JASPAR and HOCOMOCO formats, filters" begin
        d = mktempdir()
        jf = joinpath(d, "j.txt")
        open(jf, "w") do io
            write(io, "AC MA0001.1\nXX\nID AAA\nXX\nDE MA0001.1 Foo::Bar ; From JASPAR\nPO\tA\tC\tG\tT\n01\t10\t0\t0\t0\n02\t0\t10\t0\t0\n03\t0\t0\t10\t0\nXX\nCC tax_group:vertebrates\nXX\n//\n")
            write(io, "AC MA0002.1\nXX\nID BBB\nXX\nDE MA0002.1 PlantTF ; From JASPAR\nPO\tA\tC\tG\tT\n01\t0\t0\t0\t10\n02\t0\t0\t10\t0\nXX\nCC tax_group:plants\nXX\n//\n")
        end
        meta, mots = loadtransfac(jf; verbose=false)
        @test length(mots) == 2 && meta.MotifID == ["MA0001.1", "MA0002.1"]
        @test size(mots[1].pwm) == (4, 3) && size(mots[2].pwm) == (4, 2)   # metadata aligned with matrices
        meta, mots = loadtransfac(jf; taxa=["vertebrates"], verbose=false)
        @test meta.MotifName == ["Foo::Bar"] && meta.Genes == ["FOO::BAR"]
        meta, mots = loadtransfac(jf; tfs=["bar"], verbose=false)
        @test meta.MotifID == ["MA0001.1"]

        hf = joinpath(d, "h.txt")
        open(hf, "w") do io
            write(io, "VV  HOCOMOCOv14, OCT 2025\nXX\n//\nAC  M10000\nXX\nID  AHR.H14CORE.0.P.B\nXX\nNA  AHR_HUMAN\nXX\nBF  AHR_HUMAN; Species: human\nXX\nP0      A      C      G      T\n01    176    307    392    128      v\n02      3    986      1     13      C\n03      9      2    986      6      G\nXX\n//\n")
        end
        meta, mots = loadtransfac(hf; taxa=["vertebrates"], verbose=false)   # no tax_group -> kept
        @test length(mots) == 1 && meta.MotifID == ["AHR.H14CORE.0.P.B"] && meta.MotifName == ["AHR"]
        @test meta.Source == ["HOCOMOCO"] && meta.Accession == ["M10000"]
        @test isapprox(sum(mots[1].pwm, dims=1), ones(1, 3))
    end

    @testset "calibrate: DP null equals brute-force enumeration" begin
        bg = MotifScanner.background_human_zero_markov()
        pwm = [0.7 0.1 0.1 0.25; 0.1 0.7 0.1 0.25; 0.1 0.1 0.7 0.25; 0.1 0.1 0.1 0.25]
        pbg = log2.(pwm) .- log2.(bg)
        m = calibrate((name="t", id="t", pwm=pwm, pbg=pbg, pbg_rc=MotifScanner.rcm(pbg),
                       maxscore=sum(maximum(pbg, dims=1))); background=bg, pthresh=0.01)
        scores = Int[]; probs = Float64[]
        for k in Iterators.product(1:4, 1:4, 1:4, 1:4)
            push!(scores, sum(m.ipbg[k[j], j] for j in 1:4)); push!(probs, prod(bg[k[j]] for j in 1:4))
        end
        for s in unique(scores)
            @test isapprox(motifpval(m, s), sum(probs[scores .>= s]); atol=1e-12)
        end
        @test motifpval(m, typemin(Int)) == 1.0
        @test sum(probs[scores .>= m.thrint]) <= 0.01 < sum(probs[scores .>= m.thrint - 1]) + 1e-12
    end

    @testset "scanmots: calibrated calls, N exclusion, signs" begin
        bg = MotifScanner.background_human_zero_markov()
        pwm = [0.97 0.01 0.01 0.01 0.01 0.01; 0.01 0.97 0.01 0.01 0.01 0.01;
               0.01 0.01 0.97 0.97 0.01 0.01; 0.01 0.01 0.01 0.01 0.97 0.97]   # ACGGTT
        pbg = log2.(pwm) .- log2.(bg)
        m = calibrate((name="acggtt", id="m1", pwm=pwm, pbg=pbg, pbg_rc=MotifScanner.rcm(pbg),
                       maxscore=sum(maximum(pbg, dims=1))); pthresh=1e-3)   # a 6-mer cannot reach 1e-4
        ref = LongDNA{4}("CCCCCACGGTTCCCCC"); alt = LongDNA{4}("CCCCCACAGTTCCCCC")
        df = MotifScanner.scanmots(ref, alt, 8:8, 8:8, [m])
        @test nrow(df) == 1 && df.Call[1] == "loss"
        @test df.RefPval[1] <= 1e-3 < df.AltPval[1]
        @test m.sf[end] > 1e-4                         # best k-mer p of a 6-mer
        @test df.Log10P_RefAlt[1] > 0 && df.ProbNorm_AltRef[1] < 0 && df.PR_RefAlt[1] > 0
        @test df.RefTotalMotifs[1] == 1 && df.AltTotalMotifs[1] == 0
        df2 = MotifScanner.scanmots(alt, ref, 8:8, 8:8, [m])
        @test df2.Call[1] == "gain"
        # placement containing N is excluded, not scored as neutral
        nref = LongDNA{4}("CCCCCACGNTTCCCCC")
        r = scanwindow(nref, m, MotifScanner.motifscanseq(nref, 8:8, m))
        @test r.nvalid == 1 && r.start == 1            # only placement 3:8 avoids the N at 9
        allN = LongDNA{4}("CCCCCNNNNNNNNNNCCCCC")
        r2 = scanwindow(allN, m, MotifScanner.motifscanseq(allN, 10:10, m))
        @test r2.nvalid == 0 && r2.intscore == typemin(Int) && motifpval(m, r2.intscore) == 1.0
        # uncalibrated motifs are calibrated on the fly (defaults)
        dfu = MotifScanner.scanmots(ref, alt, 8:8, 8:8,
            [(name="x", id="x", pwm=pwm, pbg=pbg, pbg_rc=MotifScanner.rcm(pbg), maxscore=sum(maximum(pbg, dims=1)))])
        @test nrow(dfu) == 1 && dfu.PThresh[1] == 1e-4
    end

    @testset "loadrefseqs: REF check against FASTA" begin
        d = mktempdir(); fa = joinpath(d, "g.fa")
        seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"   # 60 bp
        open(fa, "w") do io; write(io, ">chrT\n", seq, "\n"); end
        open(fa * ".fai", "w") do io; write(io, "chrT\t60\t6\t60\t61\n"); end
        vt = DataFrame(chrom=["chrT", "chrT", "chrT"], start=[5, 6, 2], stop=[5, 6, 3],
                       ID=["ok", "bad", "near_start"], ref=["A", "G", "CG"], alt=["G", "T", "C"])
        st = loadrefseqs(vt, fa; mml=10)
        @test st.RefMatch == [true, false, true]
        @test st.seqstart[3] == 1                       # window clamped at chromosome start
        @test string(st.refseq[3][st.refind[3]]) == "CG"
    end
end
