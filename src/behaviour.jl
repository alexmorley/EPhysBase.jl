"""
    function get_path(tracking::Array, onsets::Array, length::Int)
Gets tracking from onsets for length in seconds
Returns an array of paths
"""
function get_path(tracking::Array, onsets::Array, length::Int; verbose=false)
    duration = floor(Int, (20000/512) * length)
    path_after_onset = zeros(duration+1,2,size(onsets,1))
    for (ind, samples) in enumerate(onsets)
        onset = floor(Int,(samples/512))
        try
            path_after_onset[:,:,ind] = tracking[onset:onset+duration,:]
        catch
            verbose && @show ind, onset, size(tracking,1)
        end
        #plot(path_after_onset[:,1,ind], path_after_onset[:,2,ind])
    end
    path_after_onset
end

function get_path(tracking::Array, onsets::Array, pre::Int, post::Int)
    pre = floor(Int, (20000/512) * pre)
    post = floor(Int, (20000/512) * post)
    path_after_onset = zeros(pre+post+1,2,size(onsets,1))
    for (ind, samples) in enumerate(onsets)
        onset = floor(Int,(samples/512))
        path_after_onset[:,:,ind] = tracking[onset-pre:onset+post,:]
        #plot(path_after_onset[:,1,ind], path_after_onset[:,2,ind])
    end
    path_after_onset
end

"""
    function is_tracking(tracking, ls::Int, rs::Int, bm::Int, tp::Int)
true for each value within box defined by ls,rs,bm,tp. false otherwise
"""
function is_tracking(tracking, ls::Int, rs::Int, bm::Int, tp::Int)
    is_loc = falses(size(tracking,1))
    for ind in 1:size(tracking,1)
        x = tracking[ind,1]
        y = tracking[ind,2]
        if x>ls && x<rs && y>bm && y<tp
            is_loc[ind]=true
        end
    end
    is_loc
end 

##helper function for is_tracking to visualise size of box
function plot_box(ls::Real, rs::Real, bm::Real, tp::Real, col="b")

    x1 = ls; y1 = bm
    x2 = rs; y2 = bm
    lx = linspace(x1,x2,100); ly = linspace(y1,y2,100)

    x3 = ls; y3 = tp
    x4 = rs; y4 = tp
    lx3 = linspace(x3,x4,100); ly3 = linspace(y3,y4,100)

    x5 = rs; y5 = bm
    x6 = rs; y6 = tp
    lx5 = linspace(x5,x6,100); ly5 = linspace(y5,y6,100)

    x7 = ls; y7 = bm
    x8 = ls; y8 = tp
    lx7 = linspace(x7,x8,100); ly7 = linspace(y7,y8,100)

    plot(lx, ly, col)
    plot(lx3, ly3, col)
    plot(lx5, ly5, col)
    plot(lx7, ly7, col)
end

function time2loc(tracking::Array, onsets::Array, epochlen::Int, location::Array{Int})
    tracking_epochs = get_path(tracking, onsets, epochlen);
    answer = zeros(size(tracking_epochs,3))
    for ind in 1:size(tracking_epochs,3)
        in_area = is_tracking(tracking_epochs[:,:,ind], location...)
        if sum(in_area) == 0; 
            answer[ind]=NaN
        else
            time = [x[1] for x in collect(enumerate(in_area))[in_area]][1];
            answer[ind]=time
        end
    end
    answer*(512/20000)
end

function choose_session(session_type::AbstractString, metadata)
    choosen_ses = falses(size(metadata["desen"],1))
    for ind in 1:size(metadata["desen"],1)
        choosen_ses[ind] = occursin(session_type, metadata["desen"][ind])
    end
    choosen_ses
end

"""
    function multises_time2loc(options::Dict)
takes dictionary object 
    e.g. options = Dict(
        \"session_type\" => \"sq flick\",
        \"tracking\" => tracking,
        \"onsets\" => toneSucrose,
        \"epochlen\" => 20,
        \"location\" => [100,135,110,160],
        \"metadata\" => metadata
    )
where tracking & toneSucrose are generated using gmeta.getpulses
"""
function multises_time2loc(options::Dict)
    sessions2analyse = choose_session(options["session_type"], options["metadata"])
    answer=Dict()
    for (ind,i) in enumerate(sessions2analyse)
        if i == true
            subopts = (options["tracking"][ind],options["onsets"][ind][:,1], options["epochlen"], options["location"])
            answer[ind] = time2loc(subopts...)
        end
    end
    answer
end


function getspeed(tracking, pixelspercm=4, trackingsamplerate=(512/20000), smooth=3)
    distancebetweetpoints = zeros(size(tracking,1), 1)
    if smooth > 0
        tracking = smooth(tracking, 1, smooth, 2)
    end
    
    for ind in 1:(size(tracking,1)-1)
        distancebetweetpoints[ind+1] = norm(tracking[ind,:]-tracking[ind+1,:])
    end
    
    distancebetweetpoints[1] = distancebetweetpoints[2]
    distanceincm = distancebetweetpoints / pixelspercm
    speed = distanceincm/trackingsamplerate
    return speed
end


function getspeedalltrials(epochedtracking; pixelspercm=(190/50), trackingsamplerate=(512/20000), smooth=11)
    epochedspeed = zeros(size(epochedtracking,1), size(epochedtracking,3))
    for i in 1:size(epochedtracking,3)
        epochedspeed[:,i] = getspeed(epochedtracking[:,:,i], pixelspercm, trackingsamplerate, smooth)
    end
    return epochedspeed
end


function trialPerformance(basenames, sessions_used, batchmetadata;
    location::Array{Int,1}=[100,135,110,160], epochlen::Int=25)
    trialSuc, trialQuin = Array{Float64}(), Array{Float64}()
    for (ind, bsnm) in enumerate(basenames)
        metadata = batchmetadata[bsnm]
        sessions = sessions_used[ind,:]
        for session in sessions
            tracking = metadata["tracking"][session]
            onsetsSuc = metadata["toneSucrose"][session][:,1]
            onsetsQuin = metadata["toneQuinine"][session][:,1]
 
            results = similar(onsetsSuc)
            results = time2loc(tracking, onsetsSuc, epochlen, location)
            trialSuc = cat(1, trialSuc, results)
            results = time2loc(tracking, onsetsQuin, epochlen, location)
            trialQuin =  cat(1, trialQuin, results)
        end
    end
    return(trialSuc[2:end], trialQuin[2:end])
end

function percentInside(trials)
    sum(!isnan(trials))/size(trials,1)
end

nanmean1D(trials) = mean(trials[!isnan(trials)])

function persessionperformance(bsnms, sessions, batchmetadata, epochlen=15)
    per_session_perf = zeros(length(sessions),2)
	for (ind,session) in enumerate(sessions)
        day = ind <= size(sessions,1) ? bsnms[ind] : bsnms[ind-size(sessions,1)]
        tempsuc,tempquin = trialPerformance([day],
            [session], batchmetadata, epochlen=epochlen)
        per_session_perf[ind,1] = percentInside(tempsuc)
        per_session_perf[ind,2] = percentInside(tempquin)
    end
    return per_session_perf
end
