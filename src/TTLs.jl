export loadtrigtimes, getwinaround

function trimpulse(pulse, arr, width)
    dmin,dmax = width, size(arr,1)-width
    pulse[searchsortedfirst(pulse, dmin):(searchsortedlast(pulse, dmax,
                lt=!isless, rev=true))]
end
function trimpulse(pulse, arr, width::Tuple{T,T}) where T
    dmin,dmax = width[1], size(arr,1)-width[2]
    pulse[searchsortedfirst(pulse, dmin):(searchsortedlast(pulse, dmax,
                lt=!isless, rev=true))]
end

function pulseused(trig, arr::AbstractArray, width)
    pulsesused = falses(size(trig))
    dmin,dmax = (width, size(arr,1)-width)
    pulsesused[searchsortedfirst(trig, dmin):(searchsortedlast(trig, dmax,
         lt=!isless, rev=true))] .= true
    return pulsesused
end

function bouts(t::TimeStamps, mintime=1.)
    tdiff = diff(t.timestamps).<(mintime*t.fs)
    onsets = TimeStamps(t.timestamps[[true;.!(tdiff)]], 0, t.maxtime, t.fs)
    offsets = TimeStamps(t.timestamps[[.!(tdiff); true]], 0, t.maxtime, t.fs)
    return onsets, offsets
end

function getwinaround(trig::AbstractArray, arr::AbstractArray, width)
    getwinaround(trig, arr, (width,width))
end

function getwinaround(trig::AbstractArray, arr::AbstractArray,
                      width::Tuple{T,T}) where T
    trig = trimpulse(trig, arr, width)
    aroundtrig = zeros(eltype(arr), sum(width), length(trig), size(arr,2))
    for i in 1:size(arr,2)
        getwinaround!(view(aroundtrig,:,:,i), trig, view(arr,:,i), width)
    end
    aroundtrig
end

function get_range(i, width)
    return (i-width):(i+width-1)
end
function get_range(i, width::Tuple{T,T}) where T
    return (i-width[1]):(i+width[2]-1)
end

function getwinaround!(aroundtrig, trig::Array{Float64,1}, arr, width)
    int_arr = interpolate(arr, BSpline(Linear()))
    for (indi,i) in enumerate(trig)
        for (indj,j) in enumerate(get_range(i, width))
            aroundtrig[indj,indi] = copy(int_arr(j))
        end
    end
end

function getwinaround!(aroundtrig, trig::Array{Int,1}, arr, width)
    #=@inbounds @fastmath =#for (indi,i) in enumerate(trig)
        for (indj,j) in enumerate(get_range(i, width))
            aroundtrig[indj,indi] = copy(arr[j])
        end
    end
end

"""
	readdlmfast(fname, T::Type; kwargs...)
Pre-calculates dims (assuming first line is the same as rest. And sets use_mmap to false)
Can be much faster than readdlm. Only tested for Ints and Floats.
"""
function readdlmfast(fname, T::Type; kwargs...)
    dim1 = countlines(fname)
    (dim1 == 0) && return Array{T,2}(undef,0,2)
    dim2 = open(fname) do f
        size(DelimitedFiles.readdlm_string(readline(f),
			DelimitedFiles.invalid_dlm(Char), T, '\n', true, Dict()),2)
    end
    readdlm(fname, T; dims=(dim1,dim2), use_mmap=true, kwargs...)
end

function loadtrigtimes(trig, metadata, sessions, addoffset=true; dims=nothing, fs=20000,
    header=1)
    bsnm = metadata["bsnm"]
    dataloc = metadata["fileloc"]
    session_lengths = [metadata["resofs"][1]; diff(metadata["resofs"], dims=1)]
    ntrials = zeros(Int,length(sessions))
    use_meta = false #get(metadata,trig,0) != 0
    
    # Read in First Session
    ses = sessions[1]
    trigtimes = use_meta ? copy(metadata[trig][ses]) : readdlmfast("$dataloc/$(bsnm)_$ses.$trig",Int)[header:end,1]
    offset = addoffset ? session_lengths[ses] : 0
    ntrials[1] = size(trigtimes,1) 
  
    # Read in remainder of sessions
    for (ind,ses) in enumerate(sessions[2:end])
        sestrigtime = use_meta ? copy(metadata[trig][ses]) : readdlmfast("$dataloc/$(bsnm)_$ses.$trig",Int)[header:end,1]
        sestrigtime .+= offset
        trigtimes = [trigtimes; sestrigtime]
        offset += addoffset ? session_lengths[ses] : 0
        ntrials[ind+1] = size(sestrigtime,1)
    end
    return TimeStamps(trigtimes[:], 0, sum(session_lengths[sessions]), metadata["fs"])
end

function gettrackingaroundpulses(batchmetadata, basenames, pulsename::AbstractString, sessions_used, pre, post)   
    # for each recording day get all the pulses for the sessions that were used
    trackingaroundpulses = Dict()
    epochedtracking = []
    isfirst = true
    for (j,bsnm) in enumerate(basenames)
        sessions = sessions_used[j,:]
        sessions = sessions[sessions.!=0]
        for i in sessions
            sesname = "$(bsnm)_$i"
            onsets = batchmetadata[bsnm][pulsename][i][:,1]
            tracking = batchmetadata[bsnm]["tracking"][i]

            # for each of these pulses get tracking either side of their onsets
            # +1 second around so that I can cut a second off after upsampling
            temp = get_path(tracking,onsets,pre+1,post+1)
            trackingaroundpulses[sesname] = temp

            # concatenate the pulse times
            if isfirst; epochedtracking = temp; isfirst = false
            else
                epochedtracking = cat(3, epochedtracking, temp)
            end
        end
    end

    epochedtracking, trackingaroundpulses
end

function getmeanpos(tracking,timestamps)
    trk = cat(tracking...)
    licktimestrk = round.(Int,fstimes(timestamps, trk.fs))
    (NaNMath.mean.([trk.xpos[licktimestrk],trk.ypos[licktimestrk]],)...,)::Tuple{Float64,Float64}
end

function getmeanlickpos(metadata,sessions,tracking_ses)
    if isfile("$(metadata["fullbase"]).licking_pulse")
        timestamps = loadtrigtimes("licking_pulse", metadata, sessions)
    elseif isfile("$(metadata["fullbase"]).lickA_pulse")
        timestamps = loadtrigtimes("lickA_pulse", metadata, sessions)
    elseif isfile("$(metadata["fullbase"]).lickB_pulse")
        timestamps = loadtrigtimes("lickB_pulse", metadata, sessions)
    else
        error("Licking Pulse File not found")
    end
    getmeanpos(tracking_ses,timestamps)
end
