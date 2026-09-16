
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



### load motifs from JASPAR transfac file
function loadtransfac(file, w = 1e-2, background=MotifScanner.background_human_zero_markov())

    io = open(file) 
    mode = 1
    
    motifname = ""
    motifid = ""
    currentmotif = Vector{Vector{Float64}}()
    pwms = Vector{Matrix{Float64}}()
    nsites = Float64[]
    
    motifdata = DataFrame(MotifName=String[], MotifID=String[], Property=String[], Value=String[])
    
    for line in eachline(io)
        
        fields = split(line)
        
        if mode == 1
            if fields[1] == "AC"
                mode = 2
                currentmotif = Vector{Vector{Float64}}()
            end
        elseif mode == 2
            if fields[1] == "DE"
                motifid = fields[2]
                motifname = fields[3]
                if startswith(motifname, motifid) ## JASPAR 2022 release appears to have altered the format of JASPAR 2020, remove the leading motifid
                    motifname = replace(motifname, string(motifid, ".") => "")
                end
            elseif fields[1] == "PO"
                @assert fields[2] == "A"
                @assert fields[3] == "C"
                @assert fields[4] == "G"
                @assert fields[5] == "T"
                mode = 3
            elseif fields[1] == "CC"
                ind = findfirst(':', line)
                
                cfields = [line[4:(ind-1)], line[(ind+1):end]]
                
                if length(cfields) != 2
                    @show line
                    @show cfields
                    error("")
                end
                push!(motifdata, (motifname, motifid, cfields...))
            elseif fields[1] == "//"
                mode = 1
            end
        elseif mode == 3
            if fields[1] == "XX"
                mode = 2
                pwm = reduce(hcat, currentmotif)
                v = sum(pwm, dims=1)
                push!(nsites, mean(v))
                pwm ./= sum(pwm, dims=1)
                pwm .+= w
                pwm ./= sum(pwm, dims=1)
                push!(pwms, pwm)
            else
                v = parse.(Float64, fields[2:end])
                push!(currentmotif, v)
            end
        end
    end
    close(io)
    motifmeta = unstack(motifdata, :Property, :Value)

    
    motifmeta[!, :NumSites] = nsites
    motifmeta[!, :Length] = size.(pwms, 2)

    @assert iszero(mapreduce(c -> sum(ismissing, c), +, eachcol(motifmeta)))
    dropmissing!(motifmeta)
    
    pbgs = [log2.(pwm) .- log2.(background) for pwm in pwms]
    pbg_rcs = rcm.(pbgs)
    maxscores = [sum(maximum(pbg, dims=1)) for pbg in pbgs]

    motifdata = [(name=motifname, id=id, pwm=pwm, pbg=pbg, pbg_rc=pbg_rc, maxscore=maxscore) for (motifname, id, pwm, pbg, pbg_rc, maxscore) in zip(motifmeta.MotifName, motifmeta.MotifID, pwms, pbgs, pbg_rcs, maxscores)]
    motifmeta.MotifFam = motiffamily.(motifmeta.MotifName);
    motifmeta, motifdata

end

function motiffamily(motifname)
    mf = uppercase.(motifname)
    mf = replace(mf, r"\([A-Za-z0-9.]*\)$" => "")
    
    mf = replace(mf, r"[0-9]*$" => "")
    mf = replace(mf, r"-$" => "")
    mf = replace(mf, r"FOX[A-Z]" => "FOX")
end