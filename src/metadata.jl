export getmetadata, match_ids, getsessions, DATALOCS

DATALOCS = ["/vtad1/data/amorley_merged/",
            "/vtad1/data/rrothaermel_merged/",
            "/vtad1/data/amorley_sf_merged/",
            "/vtad1/data/amorley_mergedLFP/",
            "/vtad1/data/rrothaermel_mergedLFP/",
           "/vtad2/data/ddLab_merged/"]

function finddata(bsnm, datalocs=DATALOCS)
    for fileloc in datalocs
        for fileloc_t in [
				"$fileloc$bsnm.par",
                "$fileloc$bsnm/$bsnm.par",
                "/mnfs/$fileloc$bsnm.par",
                "/mnfs/$fileloc$bsnm/$bsnm.par"]
            isfile(fileloc_t) && return String(split(fileloc_t, '.')[1])
        end
    end
    return ""
end 

"""
function getmetadata(bsnm::String, dataloc::Array{String,1})
    Read in all the metadata for a given recording day. Required files:
        - \$bsnm.par
        - \$bsnm.resofs
        - \$bsnm.desen
        - \$bsnm.desel

        Optional Files:
        - \$bsnm.des
        - \$bsnm.trigchannels
        - \$bsnm.trigchannels2
"""
function getmetadata(bsnm::String, datalocs=DATALOCS; warn_=true, by_tet=true)
    fullbase = finddata(bsnm, datalocs)
 
	fullbase == "" && error("Par file for basename $bsnm not found anywhere!")

    # Read some metadata files
    par = readdlm("$fullbase.par", Any);
    resofs = readdlm("$fullbase.resofs", Int64);
    desen = readlines("$fullbase.desen");
    tetrodes = readlines("$fullbase.desel");

    # check if spikes
    spkflag = isfile("$fullbase.des") & isfile("$fullbase.res")
    if spkflag
        des = readlines("$fullbase.des");
        if by_tet
            des_tet = [readlines("$fullbase.des.$i") for (i,t) in enumerate(tetrodes)]
        end
    end


    # Get some useful stuff from the par file
    numchan = Int64(par[1,1])
    channels = convert(Array{Int64,1},par[4:3+size(tetrodes,1),2])
    fs = 1/(par[2,1]/1000000) #μs -> fs

    # Create structure - is dict the best option here?
    info = Dict(
    "bsnm"         => bsnm,
	"fullbase"     => fullbase,
    "fileloc"      => join(split(fullbase,'/')[1:end-1],'/'),
    "samplerate"   => Int64(par[2,1]),
    "fs"           => fs,
    "trackingfs"   => 20000/512,
    "numses"       => length(desen),
    "numchan"      => numchan,
    "bestchannels" => channels,
    "desen"        => convert(Array{String,1}, desen),
    "tetrodes"     => convert(Array{String,1}, tetrodes),
    "sesnames"     => convert(Array{String,1}, par[5+size(tetrodes,1):end,1]),
    "resofs"       => resofs,
    "seslengths"   => vcat(resofs[1], diff(resofs,dims=1)),
    "whlscale"     => 4.0
    )

    if spkflag
        info["cellIDs"]  = strip.(des)
        if by_tet
            info["cellIDs_by_tet"]  = [strip.(d) for d in des_tet]
            info["tetlist"] = reduce(vcat, [[ti for i in t] for (ti,t) in 
                                    enumerate(info["cellIDs_by_tet"])])
        end
        info["numcells"] = size(des,1)
    end

    if isfile("$fullbase.trigchannels")
        # read in trig channel files
        trigchannel = readdlm("$fullbase.trigchannels")[:,1];
        try
            trigchannel2 = readdlm("$fullbase.trigchannels2")
            trigchannel = cat(trigchannel, trigchannel2, dims=1)
        catch
        end
        # Add custom trigchannel names
        # Read in TrigPulses
        for i in trigchannel
            info[i] = getpulses(info, i)
        end
    else
        warn_ && warn("No .trigchannels file found")
    end

    #return meta
    return info
end

function getsessions(metadata,sess_types)
    desen = metadata["desen"]
    findall(.|([occursin.(s, desen) for s in sess_types]...))
end

function maxtime(metadata::Dict, sessions=:)
    sum([metadata["resofs"][1];diff(metadata["resofs"],dims=1)][sessions])
end

"""
function match_ids(cellIDs::Vector{String},uniqIDs::Vector{String}=regions)
    Get BitArray matching a vector of cell IDs to a vector of regions.
"""
function match_ids(cellIDs,uniqIDs=regions)
    cat([[occursin(id,x) for x in cellIDs] for id in uniqIDs]..., dims=2)
end

"""
function getmetadata(daystoprocess::Array, dataloc::String)
    Deprecated ...
    Create a batch metadata Dict with metatdata from each day in "days2process"
"""
function getmetadata(daystoprocess::Array, dataloc)
    batchmeta = Dict{String,Dict}();
    for (ind,j) in enumerate(daystoprocess)
        day2processfileloc = "$dataloc$j/"
        batchmeta[j] = getmetadata(j, day2processfileloc);
    end
    return batchmeta
end

function filtsortcells(metadata, exclude=["lick"])
    cells2use = !any(metadata["cellIDs"] .== exclude,2)
    cellIDs = metadata["cellIDs"][cells2use]
    return cells2use, cellIDs, sortperm(cellIDs), length(cellIDs)
end

function getpulses(metadata::Dict, pulsename::AbstractString)
    all_pulse=Dict()
    sessions_pulse_files = [joinpath(metadata["fileloc"],"$x.$(pulsename)_pulse") 
    for x in metadata["sesnames"]]
    for (ind,eachses) in enumerate(sessions_pulse_files)
        if stat("$eachses").size == 0; continue ; end
        track_file = readdlm("$eachses",Int64)
        all_pulse[ind] = track_file
    end
    all_pulse
end
