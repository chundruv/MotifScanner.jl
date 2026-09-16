module MotifScanner

using BioSequences
using StatsBase
using DataFrames

using FASTX
using ProgressMeter

## This is a package for scanning motifs in DNA sequences and plotting sequence motifs.



export loadmeme, loadmemelibrary, loadhomer, scanmotif, scanmotstats, consensus, plotseq, plotletter!, seqlogo!, seqlogo, loadrefseqs, loadvcf, motifscanall, loadtransfac

include("loadmeme.jl")
include("scanning.jl")
# include("plotseqlogo.jl")
# include("seqlogodata.jl")

include("variants.jl")
end
