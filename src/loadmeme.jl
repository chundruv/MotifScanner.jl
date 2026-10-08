
background_human_zero_markov()   = [0.2953, 0.2047, 0.2047, 0.2953]
background_xenopus_zero_markov() = [0.2996, 0.2004, 0.2004, 0.2996]

### Load meme into named tuple
function loadmeme(file, w = 1e-2, background=background_human_zero_markov())
   
    motifname = ""
    motifid   = ""
    mode = 1
    data = Vector{Vector{Float64}}()
    for line in eachline(file)
        if mode == 1
            if startswith(line, "MOTIF")
                fields = split(line)
                motifid   = fields[2]
                motifname = fields[3]
            elseif startswith(line, "letter-probability")
                mode = 2
            end
        else
            if startswith(line, "URL")
                mode = 1
                break
            end
            v = parse.(Float64, split(line))
            push!(data, v)
        end
    end
    pwm = reduce(hcat, data)
    pwm .+= w
    pwm ./= sum(pwm, dims=1)
    
    pbg = log2.(pwm) .- log2.(background)
    pbg_rc = rcm(pbg)
    maxscore = sum(maximum(pbg, dims=1))
    (name=motifname, id=motifid, pwm=pwm, pbg=pbg, pbg_rc=pbg_rc, maxscore=maxscore)
end


function loadmemelibrary(file, w = 1e-2, background=background_human_zero_markov())

    motifnames = String[]
    motifids = String[]
    motiflengths = Int[]
    nsites = Int[]
    mode = :notreading
    currentmotif = Vector{Vector{Float64}}()
    pwms = Vector{Matrix{Float64}}()
    
    for line in eachline(file)
        if mode == :notreading
            if startswith(line, "MOTIF")
                
                fields = split(line)
                motifid   = fields[2]
                motifname = fields[3]
                
                push!(motifnames, motifname)
                push!(motifids, motifid)
                
                mode = :reading
            end
        elseif mode == :reading
            
            if startswith(line, "letter-probability")
                currentmotif = Vector{Vector{Float64}}()
                ldata = replace(line, "letter-probability matrix: " => "")
                fields = split.(split(replace(ldata, "= " => "=")), "=")
                
                dd = Dict(k => v for (k, v) in fields)
                push!(motiflengths, parse(Int, dd["w"]))
                push!(nsites, parse(Int, dd["nsites"]))

            elseif startswith(line, "URL")
                mode = :notreading
                pwm = reduce(hcat, currentmotif)
                pwm .+= w
                pwm ./= sum(pwm, dims=1)
                push!(pwms, pwm)
                
            else
                v = parse.(Float64, split(line))
                push!(currentmotif, v)
            end
        end
    end

    motifmeta = DataFrame(MotifName=motifnames, MotifID=motifids, Lengths=motiflengths, NumSites=nsites)
    pbgs = [log2.(pwm) .- log2.(background) for pwm in pwms]
    pbg_rcs = rcm.(pbgs)
    maxscores = [sum(maximum(pbg, dims=1)) for pbg in pbgs]

    motifdata = [(name=motifname, id=id, pwm=pwm, pbg=pbg, pbg_rc=pbg_rc, maxscore=maxscore) for (motifname, id, pwm, pbg, pbg_rc, maxscore) in zip(motifnames, motifids, pwms, pbgs, pbg_rcs, maxscores)]

    motifmeta, motifdata

end

function loadhomer(file, w = 1e-2, background=background_human_zero_markov())

    io = open(file)
    line = readline(io)

    fields = split(line)
    consensus = fields[1][2:end]
    longname = fields[2]
    name = first(split(longname, "/"))
    sct = parse(Float64, fields[3])


    data = Vector{Vector{Float64}}()
    for line in eachline(io)
        v = parse.(Float64, split(line))
        push!(data, v)
    end



    pwm = reduce(hcat, data)
    pwm .+= w
    pwm ./= sum(pwm, dims=1)
    
    pbg = log2.(pwm) .- log2.(background)
    pbg_rc = rcm(pbg)
    maxscore = sum(maximum(pbg, dims=1))
    (name=name, id=longname, pwm=pwm, pbg=pbg, pbg_rc=pbg_rc, maxscore=maxscore, sct=sct, consensus=consensus)
end



### load motifs from a TRANSFAC-format file (JASPAR or HOCOMOCO v14; plain or .gz)
"""
    loadtransfac(file, w=1e-2, background=background_human_zero_markov();
                 taxa=nothing, tfs=nothing, verbose=true)

Read a TRANSFAC matrix file. Handles both JASPAR (`DE <id> <name>`, `PO`, `CC key:value`)
and HOCOMOCO v14 (`ID <motif>`, `NA <TF>_HUMAN`, `P0`, trailing consensus column).
Matrix rows are parsed from columns 2-5 (A C G T); any further column is ignored.

Metadata rows and matrices are kept in file order, so `motifmeta[i, :]` always
describes `motifdata[i]`.

Filters (applied before log-odds are built):
- `taxa`: e.g. `["vertebrates"]`. Applied only if the file carries `tax_group`
  annotations (JASPAR). HOCOMOCO has none; it is human-only, so all motifs are kept.
- `tfs`: collection of TF gene symbols (case-insensitive). A motif is kept if any
  component of its name (dimers split on `::`) is in the set.

Returns `(motifmeta::DataFrame, motifdata::Vector{NamedTuple})` as before.
`motifmeta` gains `Accession`, `Genes` and `Source` columns; CC properties may be
`missing` for motifs that lack them.
"""
function loadtransfac(file, w = 1e-2, background=MotifScanner.background_human_zero_markov();
                      taxa=nothing, tfs=nothing, verbose=true)

    io = endswith(lowercase(file), ".gz") ? open(`gzip -cd $(file)`, "r") : open(file, "r")

    accs = String[]; ids = String[]; names_ = String[]; sources = String[]
    pwms = Vector{Matrix{Float64}}(); nsites = Float64[]
    props = Vector{Dict{String,String}}()

    inrec = false; inmat = false
    ac = ""; idtag = ""; na = ""; deid = ""; dename = ""
    rows = Vector{Vector{Float64}}(); pr = Dict{String,String}()

    function finish!()
        isempty(rows) && return
        if !isempty(deid)                        # JASPAR
            mid, mname, src = deid, dename, "JASPAR"
            if startswith(mname, mid)            # JASPAR 2022+ prefixes the name with the id
                mname = replace(mname, string(mid, ".") => "")
            end
        else                                     # HOCOMOCO (or generic)
            mid = isempty(idtag) ? ac : idtag
            mname = isempty(na) ? first(split(mid, '.')) : replace(na, r"_(HUMAN|MOUSE)$" => "")
            src = occursin("H14", mid) || occursin("H13", mid) || occursin("H12", mid) ? "HOCOMOCO" : "TRANSFAC"
        end
        pwm = reduce(hcat, rows)
        push!(nsites, mean(sum(pwm, dims=1)))
        pwm = pwm ./ sum(pwm, dims=1)
        pwm .+= w
        pwm ./= sum(pwm, dims=1)
        push!(accs, ac); push!(ids, mid); push!(names_, mname); push!(sources, src)
        push!(pwms, pwm); push!(props, pr)
        return
    end

    for line in eachline(io)
        s = strip(line)
        isempty(s) && continue
        tag = length(s) >= 2 ? s[1:2] : s
        if tag == "AC"
            inrec && finish!()
            inrec = true; inmat = false
            ac = strip(s[3:end]); idtag = ""; na = ""; deid = ""; dename = ""
            rows = Vector{Vector{Float64}}(); pr = Dict{String,String}()
        elseif !inrec
            continue
        elseif tag == "P0" || tag == "PO"
            f = split(s)
            @assert f[2:5] == ["A", "C", "G", "T"] "unexpected matrix header: $s"
            inmat = true
        elseif inmat && isdigit(s[1])
            f = split(s)
            push!(rows, parse.(Float64, f[2:5]))
        elseif tag == "XX"
            inmat = false
        elseif tag == "//"
            inmat = false
            finish!(); inrec = false; rows = Vector{Vector{Float64}}()
        elseif tag == "ID"
            idtag = strip(s[3:end])
        elseif tag == "NA"
            na = strip(s[3:end])
        elseif tag == "DE"
            f = split(s)
            deid = f[2]; dename = length(f) >= 3 ? f[3] : f[2]
        elseif tag == "CC"
            ind = findfirst(':', s)
            isnothing(ind) || (pr[strip(s[3:(ind-1)])] = strip(s[(ind+1):end]))
        end
    end
    inrec && finish!()
    close(io)

    genes = [uppercase.(split(replace(n, r"\(var\.\d+\)" => ""), "::")) for n in names_]
    keep = trues(length(pwms))
    if !isnothing(taxa)
        if any(p -> haskey(p, "tax_group"), props)
            tset = Set(lowercase.(taxa))
            keep .&= [lowercase(get(p, "tax_group", "")) in tset for p in props]
        elseif verbose
            println(stderr, "loadtransfac: no tax_group annotation in $(basename(file)); taxa filter not applied")
        end
    end
    if !isnothing(tfs)
        tfset = Set(uppercase.(strip.(collect(tfs))))
        keep .&= [any(g -> g in tfset, gs) for gs in genes]
    end
    verbose && println(stderr, "loadtransfac: kept $(sum(keep)) of $(length(keep)) motifs from $(basename(file))")

    idx = findall(keep)
    motifmeta = DataFrame(MotifName=names_[idx], MotifID=ids[idx], Accession=accs[idx],
                          Genes=[join(g, "::") for g in genes[idx]], Source=sources[idx],
                          NumSites=nsites[idx], Length=size.(pwms[idx], 2))
    for k in sort(unique(reduce(vcat, [collect(keys(p)) for p in props[idx]]; init=String[])))
        motifmeta[!, k] = [get(p, k, missing) for p in props[idx]]
    end
    motifmeta.MotifFam = motiffamily.(motifmeta.MotifName)

    pbgs = [log2.(pwm) .- log2.(background) for pwm in pwms[idx]]
    motifdata = [(name=n, id=i, pwm=pwm, pbg=pbg, pbg_rc=rcm(pbg), maxscore=sum(maximum(pbg, dims=1)))
                 for (n, i, pwm, pbg) in zip(motifmeta.MotifName, motifmeta.MotifID, pwms[idx], pbgs)]
    motifmeta, motifdata
end

function motiffamily(motifname)
    mf = uppercase.(motifname)
    mf = replace(mf, r"\([A-Za-z0-9.]*\)$" => "")
    
    mf = replace(mf, r"[0-9]*$" => "")
    mf = replace(mf, r"-$" => "")
    mf = replace(mf, r"FOX[A-Z]" => "FOX")
end