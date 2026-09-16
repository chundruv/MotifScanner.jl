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
end
