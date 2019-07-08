export resample, get_tracking_around

function load(trk::Type{Tracking}, metadata, sessions)
    loadtracking(metadata, sessions)
end

import DSP.resample

function resample2_(trk::Tracking, to_fs)
    Tracking(DSP.resample(trk.xpos, to_fs/trk.fs),
    DSP.resample(trk.ypos, to_fs/trk.fs),
    DSP.resample(trk.speed, to_fs/trk.fs),
    trk.smoothing, trk.pixelpercm, to_fs)
end

## this f**ing resample function doesn't give the right number of samples when
## resampling factor is greater than 2 grrr...
function resample_(trk::Tracking, to_fs)
    exprs = []
    for trk_series in Symbol[:xpos, :ypos, :speed]
        push!(exprs, :(DSP.resample(trk.$trk_series, $to_fs/fs(trk))))
    end
    Tracking(eval.(exprs)..., trk.smoothing, trk.pixelpercm, to_fs)
end

function resample(trk1::Tracking, to_fs)
    while !(trk1.fs ≈ to_fs)
        if to_fs/fs(trk1) > 1.5
            trk1 = resample2_(trk1, 1.5 * fs(trk1))
        else
            trk1 = resample2_(trk1, to_fs)
        end
    end
    trk1
end

"""
function fill_tracking!(track_file, interpolate=true)
    Either interpolate over -1 values (where tracking has been misdetected)

    Or set values to NaN if `interpolate = false`
        """
        function fill_tracking!(track_file, interpolate=true)
            f(x)=(x > 0)
            # set values before tracking started to NaN
            for i in 1:(findfirst(f,track_file[:,1])-1)
                track_file[i,:] = [NaN, NaN]
            end

            # loop over all rows in tracking file
            # from the first non -1 onwards
            for i in findfirst(f,track_file[:,1]):size(track_file,1)
                if track_file[i,1] == -1.0 # if mis/not detected at that point
                    if interpolate
                        # search for next position where detection successful 
                        nextrealind = findnext(f, track_file[:,1], i)
                        if nextrealind==nothing; break; end # no more successful locations

                        interpolate!(view(track_file,:,1), i, nextrealind) # x-coord
                        interpolate!(view(track_file,:,2), i, nextrealind) # y-coord

                    else # if -1s should be set to missing values instead
                    track_file[i,:] = [NaN, NaN]
                end
            end
        end

        # set values after tracking finished to NaN
        for i in size(track_file,1):-1:1
            if track_file[i,1] < 0
                track_file[i,:] = [NaN, NaN]
            else
                break
            end
        end
    end

"""
interpolate!(track_file, ind, nextrealind)
Linearly interpolate column (col) of vector from ind to nextrealind.
"""
function interpolate!(track_file, ind, nextrealind)
    prevreal = track_file[ind-1]
    nextreal = track_file[nextrealind]
    # gradient of line to use
    jump = (nextreal-prevreal)/(nextrealind-ind+1)
    if jump == 0 # otherwise will throw error on range creation
        setindex!(track_file, repeat([prevreal], nextrealind-ind), ind:nextrealind-1)
    else 
        setindex!(track_file, prevreal:jump:nextreal, ind-1:nextrealind)
    end
end

function getspeed(tracking::Array{Float64,2}, pixelpercm, fs)
    speed = zeros(size(tracking,1))
    for ind in 1:size(tracking,1)
        xpos₁,ypos₁ = (tracking[ind,1],tracking[ind,2]) 
        isnan(xpos₁) && (speed[ind]=NaN; continue)
        xpos₂,ypos₂ = (tracking[ind-1,1], tracking[ind-1,2])
        isnan(xpos₂) && (speed[ind]=NaN; continue)
        d = √((xpos₂-xpos₁)^2+(ypos₂-ypos₁)^2) / pixelpercm
        speed[ind] = d * fs
    end
    return speed
end

"""
function loadtracking(metadata::Dict, interpolate=true)
Get the tracking for each session in a recording day.
    """
function loadtracking(metadata::Dict, sessions=:, interpolate=true, smoothing=10)
    sessions_whl = ["$x.whl" for x in metadata["sesnames"][sessions]]
    fileloc = metadata["fileloc"]
    all_tracking = Vector{Tracking}(undef,length(sessions_whl))
    for (ind,eachses) in enumerate(sessions_whl)
        track_file = readdlmfast("$fileloc/$eachses", Float64)
        fill_tracking!(track_file, interpolate) # deals with of -1s
        # exclude nans from filter
        nans = .!(isnan.(track_file[:,1]))
        tofilt = view(track_file, findfirst(nans):findlast(nans), :)
        imfilter!(tofilt, tofilt, KernelFactors.IIRGaussian([smoothing;0]))
        speed = getspeed(track_file, metadata["whlscale"], metadata["trackingfs"])
        # create Tracking object
        all_tracking[ind] = Tracking(track_file[:,1], track_file[:,2],
                                     speed,
                                     smoothing,
                                     metadata["whlscale"],
                                     metadata["trackingfs"])
    end
    all_tracking
end

"""
function gettrackbatch!(batchmetadata)
    Adds another entry to each dict in batch metadata called tracking
"""
function gettrackbatch!(batchmetadata)
    for (k,v) in batchmetadata
        v["tracking"] = gettrack(v)
    end
end

"""
    findinloc{T<:Real}(xy, loc::Array{Tuple{T,T},2})
Return the first index where both x & y positions are within the
limits specified by `loc`.
"""
function findinloc(xy::Array{T,2}, loc::Array{Tuple{T,T},2}) where T<:Real
    findfirst(all(between.(xy, loc),2))
end

function get_tracking_around(trig::TimeStamps,metadata,sessions,tracking,
                             n_samples;
                             fields=[:xpos :ypos])
    trig_fs              = fstimes(trig,tracking.fs);
    tracking_around      = getwinaround(trig_fs,
                                        reduce(hcat,
                                               getfield.([tracking], fields)),
                                        n_samples)
    return tracking_around
end
function get_tracking_around(trig::String,metadata,sessions,tracking,n_samples;
                             fields=[:xpos :ypos])
    trig_timestamps      = loadtrigtimes("$(trig)_pulse",metadata,sessions);
    trig_fs              = fstimes(trig_timestamps,tracking.fs);
    tracking_around      = getwinaround(trig_fs, reduce(hcat, getfield.([tracking], fields)), n_samples)
    return tracking_around
end
