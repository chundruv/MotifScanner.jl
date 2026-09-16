
## reverse complement
rcm(M) = reverse(reverse(M, dims=1), dims=2)

motifpbg_rc(motif) = hasproperty(motif, :pbg_rc) ? motif.pbg_rc : rcm(motif.pbg)
motifmaxscore(motif) = hasproperty(motif, :maxscore) ? motif.maxscore : sum(maximum(motif.pbg, dims=1))


### loop for forward and reverse motif scanning
## no doubt there are efficiency savings available here
### eg don't allocate rot - just flip indicies!
function scanmotif(seq, mot)
    scanmotif(seq, mot, rcm(mot))
end

function scanmotif(seq, mot, rot)
    n = length(seq)
    m = size(mot, 2)
    fscores = zeros(Float64, n - m + 1)
    rscores = zeros(Float64, n - m + 1)
    for i = 1:(n - m + 1)
        @inbounds @simd for j = 1:m
            if seq[i + j - 1] == DNA_A
                fscores[i] += mot[1, j]
                rscores[i] += rot[1, j]
            elseif seq[i + j - 1] == DNA_C
                fscores[i] += mot[2, j]
                rscores[i] += rot[2, j]
            elseif seq[i + j - 1] == DNA_G
                fscores[i] += mot[3, j]
                rscores[i] += rot[3, j]
            elseif seq[i + j - 1] == DNA_T
                fscores[i] += mot[4, j]
                rscores[i] += rot[4, j]
            end
        end
    end
    fscores, rscores
end

### scan max
function scanmax(seq, motif)
    fs, rs = scanmotif(seq, motif.pbg, motifpbg_rc(motif))
    fm, fi = findmax(fs)
    rm, ri = findmax(rs)
    n = size(motif.pbg, 2)
    if fm > rm
        maxscore = fm
        start = fi
        stop = fi + n - 1
        strand = "+"
    
    else
        maxscore = rm
        start = ri
        stop  = ri + n - 1
        strand = "-"
    end
    maxscore, start, stop, strand
end

### scan sum max
function scansummax(seq, motif, thr=0)
    fs, rs = scanmotif(seq, motif.pbg, motifpbg_rc(motif))
    maxscore = motifmaxscore(motif)
    thrscore = thr * maxscore

    totalmotifs = 0
    totalmax = 0
    sumscore_num = 0.0

    fm = -Inf
    fi = 1
    rm = -Inf
    ri = 1
    @inbounds for i in eachindex(fs, rs)
        f = fs[i]
        r = rs[i]

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
            fi = i
        end
        if r > rm
            rm = r
            ri = i
        end
    end

    sumscore = sumscore_num / maxscore

    n = size(motif.pbg, 2)
    if fm > rm
        maxval = fm
        start = fi
        stop = fi + n - 1
        strand = "+"
    else
        maxval = rm
        start = ri
        stop  = ri + n - 1
        strand = "-"
    end
    maxval, start, stop, strand, sumscore, totalmotifs, totalmax
end



### motif scan stats
function scanmotstats(mot, seq::T, thr=5) where{T}
        
    fs, rs = scanmotif(seq, mot.pbg, motifpbg_rc(mot))
    
    ## neg ecdf
    nec = ecdf([-fs ; -rs])
    
    res = DataFrame(Motif=String[], start=Int[], stop=Int[], score=Float64[], strand=String[], emp_p=Float64[], prmax=Float64[], seq=Vector{T}())
    
    maxscore = motifmaxscore(mot)
    
    n = size(mot.pbg, 2)
     @inbounds for f in eachindex(fs)
         if fs[f] > thr
             push!(res, (mot.name, f, f + n - 1, fs[f], "+", nec(-fs[f]), fs[f]/maxscore, seq[f:f+n-1]))
         end
    end
    
     @inbounds for r in eachindex(rs)
         if rs[r] > thr
             push!(res, (mot.name, r, r + n - 1, rs[r], "-", nec(-rs[r]), rs[r]/maxscore, reverse_complement(seq[r:r+n-1])))
         end
    end
    res
end

function consensus(mot)
    bases = (DNA_A, DNA_C, DNA_G, DNA_T)
    LongDNASeq(mapreduce(c -> bases[argmax(c)], vcat, eachcol(mot.pwm)))
end

function motifind(motname, mots)
    ms = findall(m -> occursin(motname, m.name), mots) 
    if isempty(ms)
       error("$motname not found") 
    elseif length(ms) > 1
        error("Multiple matches for $motname\n$([mots[m].name for m in ms])")
    else
       first(ms)
    end 
    
end

function consensus(motname, mots)
    ms = filter(m -> occursin(motname, m.name), mots) 
    if isempty(ms)
       error("$motname not found") 
    elseif length(ms) > 1
        error("Multiple matches for $motname\n$([m.name for m in ms])")
    else
        consensus(first(ms))
    end
end
