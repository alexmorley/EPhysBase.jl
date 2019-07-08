export load
# extends functions (so far from DSP.jl) to work with EPhysBase types
export getfreqrange, lpfilt, notchfilt, interp_range

function load(T::Type{LFP}, metadata, sessions::Array{Int,1}, fs=1250.)
    lfp = load.([LFP],[metadata],sessions)
    cat(lfp...)
end

function load(T::Type{LFP}, metadata, session::Int, fs=1250.; opts...)
    if fs == 1250.
        ext = "eeg"
        chans = metadata["bestchannels"].+1
        nchans = metadata["numchan"]
    elseif fs == 5000.
        ext = "eegh"
        nchans = length(metadata["desel"])
        chans = 1:nchans
    end

    fileloc = metadata["fileloc"]; bsnm = metadata["bsnm"]
    path = "$fileloc/$(bsnm)_$session.$ext"
    return loadlfp(path, chans, nchans, fs, metadata["tetrodes"])
end

function loadlfp(path::AbstractString, chans, nchans, fs, ids=string.(1:12))
    nsamples = floor(Int,stat(path).size/(nchans*2))
    eeg = Mmap.mmap(path, Array{Int16,2}, (nchans,nsamples), grow=false)
    return LFP([Signal(eeg[chan,:], fs) for chan in chans], ids, fs)
end

import Base.getindex
getindex(lfp::LFP, i::Int) = lfp.lfp[i]
function getindex(lfp::LFP, i::Array{I,1}) where I<:Int
    LFP(getindex.([lfp], i), lfp.ids[i], lfp.fs)
end

import DSP.Filters.resample
function resample(t::Signal, to_fs::Real)
    out = resample(EPhys.signal(t), to_fs/fs(t))
    Signal(out, to_fs)
end

function resample(t::LFP, to_fs::Real)
    LFP(resample.(t.lfp, [to_fs]), t.ids, to_fs)
end

## extend some of the functions in DSP.jl to take Signal and LFP type objects
for pfunc in [:periodogram, :welch_pgram, :mt_pgram, :spectrogram, :stft]
    @eval begin
        import DSP.Periodograms.$pfunc
        function $pfunc(s::Signal; fs=s.fs, kwargs=())
            $pfunc(s.signal; fs=fs, kwargs...)
        end
        function $pfunc(s::LFP; kwargs=())
            $pfunc.(s.lfp; fs=s.fs, kwargs...)
        end
    end    
end

function getfreqrange(p::DSP.Periodograms.Periodogram, freqs::Tuple{<:Real,<:Real})
    inds = freqs[1] .< p.freq .< freqs[2]
    return p.freq[inds], p.power[inds]
end


## add some filter convinience funcs
function lpfilt(lpcutoff, s::Signal)
    responsetype = Lowpass(lpcutoff, fs=s.fs)
    designmethod = Butterworth(4)
    sigfilt = filt(digitalfilter(responsetype, designmethod),s.signal)
    return Signal(sigfilt,s.fs)
end

function notchfilt(lpcutoff, s::Signal)
    notchfilters = SecondOrderSections([iirnotch(x, 0.5,
            fs=s.fs) for x in 50:50],1)
    sigfilt = filt(notchfilters, s.signal)
    return Signal(sigfilt,s.fs)
end

notchfilt(s::Signal) = notchfilt(350, s)
lpfilt(s::Signal)    = lpfilt(350, s)

notchfilt(l::LFP)    = LFP(notchfilt.(350, l.lfp), l.ids, l.fs)
lpfilt(l::LFP)       = LFP(   lpfilt.(350, l.lfp), l.ids, l.fs)

####
function interp_range(y, range::Union{BitArray{1},Array{Bool,1}}, ησ::Real = 3)
    range2len = sum(range)
    range2pre,range2post = falses(range),falses(range)
    range2pre[(findfirst(range)-range2len):findfirst(range)-1] .= true
    range2post[(findlast(range)+1):findlast(range)+range2len] .= true
    range2 = (range2pre, range2post)
    
    y2 = copy(y)
    interp_range!(y2, y, range, range2, ησ)
    return y2
end

function interp_range(y, range, range2, ησ::Real = 3)
    y2 = copy(y)
    interp_range!(y2, y, range, range2, ησ)
    return y2
end

function interp_range!(y2, y, range, range2, ησ = 3)
    N = sum(range)
    
    ## construct a basic model of pre and post using mean, std & unform noise
    rangepre,rangepost = range2
    η1,η2 = std(y[rangepre]), std(y[rangepost])
    η = ησ*(rand(N)-0.5).*linspace(η1, η2, N)
    μ1,μ2 = mean(y[rangepre]), mean(y[rangepost])
    μ = linspace(μ1, μ2, N)
    
    ## fill in the range with values constructed from the model
    rep = μ .+ η;
    y2[range] = rep
    return y2
end
