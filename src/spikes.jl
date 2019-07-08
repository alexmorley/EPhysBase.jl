export SpikeTimes, SpikeMatrix, BinnedSpikes
export loadspikes, rates, filter

export ratefilt
export findspikesaround, inst_rate, inst_sterr, addjitter
export binZ, binZ3 

# Functions for loading in spikes
loadspikes(T::Type, metadata, sessions, args...) = T(loadresclu(metadata,sessions)..., args...)

function loadspikes(Ts::Vector{T}, metadata, sessions, args...) where T<:Type
    spikes = loadresclu(metadata,sessions)
    return [Tc(spikes...) for Tc in Ts]
end

### Spike Times ###
struct SpikeTimes{T,I<:Integer}
    spiketimes::Vector{TimeStamps{T}}
    spikeIDs::Vector{I}
end

load(T::Type{SpikeTimes}, metadata, sessions, args...) = loadspikes(T, metadata, sessions, args...)

mintime(spikes::SpikeTimes) = 0
maxtime(spikes::SpikeTimes) = spikes.spiketimes[1].maxtime
fs(spikes::SpikeTimes) = spikes.spiketimes[1].fs
rates(spikes::SpikeTimes) = (length.(spikes.spiketimes)./maxtime(spikes))*fs(spikes)

import Base.filter
function filter(I::AbstractArray, spiketimes::SpikeTimes)
    inds = findall(I)
    SpikeTimes(spiketimes[:][inds], spiketimes.spikeIDs[inds])
end

function filter(f::Function, spiketimes::SpikeTimes)
    inds = findall(f,spiketimes[:])
    SpikeTimes(spiketimes[:][inds], spiketimes.spikeIDs[inds])
end
rate(ts) = ts.fs*(length(ts.timestamps)/ts.maxtime)
ratefilter(ts,lo,hi) = lo < rate(ts) < hi
badid(id) = startswith(id, "u") | startswith(id, "l")
idfilter(ids) = .!(badid.(ids))

function SpikeTimes(allspiketimes::TimeStamps, allspikeIDs::TimeStamps, nIDs)
    spikeIDs = collect(2:nIDs+1) 
    spiketimes = Vector{TimeStamps{eltype(allspikeIDs.timestamps)}}(undef,nIDs)

    for (i,id) in enumerate(spikeIDs)
        spiketimes[i] = TimeStamps(allspiketimes.timestamps[allspikeIDs.timestamps.==id],
                                   0,
                                   allspiketimes.maxtime,
                                   allspiketimes.fs)
    end 
    return SpikeTimes(spiketimes, spikeIDs)
end


function loadresclu(metadata, sessions, suffix="")
    allspiketimes = loadtrigtimes("res$suffix", metadata, sessions)
    allspikeIDs   = loadtrigtimes("clu$suffix", metadata, sessions, false, header=2)
    nIDs          = maximum(allspikeIDs.timestamps) - 1 #metadata["cellIDs"])
    return allspiketimes, allspikeIDs, nIDs
end

### Spike Matrix ###
struct SpikeMatrix{N} # <: TimeSeries{T,N} 
    spikematrix::BitArray{N}
    fs::Float64
end 

SpikeMatrix(a,b,c) = SpikeMatrix(a,b)

function SpikeMatrix(allspiketimes::TimeStamps{T}, allspikeIDs::TimeStamps{T}, sr::Float64=0.001) where T
    binsize = sr * allspiketimes.fs
    nbins = cld(allspiketimes.maxtime, binsize)
    spkmatrix = BitArray{2}(Int(nbins), maximum(unique(allspikeIDs.timestamps)))
    spikes2mat!(spkmatrix,allspiketimes.timestamps,allspikeIDs.timestamps,binsize)
    return SpikeMatrix(spkmatrix, 1/sr)
end

function spikes2mat!(spkmatrix::BitArray{2}, spiketimes::Array{T,1}, spikeIDs::Array{Int,1}, binsize::Float64) where T<:Real
    for i in 1:length(spiketimes)
        spkmatrix[Int(cld(spiketimes[i],binsize)),spikeIDs[i]] = true
    end
end

function spikes2mat(spikes::SpikeTimes, binsize::Float64, T::Type)
    spkmatrix = zeros(T, SpikeTimes.spiketimes.maxtime, nclu, )
    spikes2mat!(spikes, binsize, spkmatrix)
    spkmatrix
end

### Binned Spike Matrix ###
mutable struct BinnedSpikes
    count::Array{Int,2}
    binsize::Float64
end

BinnedSpikes(res::TimeStamps, clu::TimeStamps, nids::Int, binsize::Int) = begin
    T = maximum(res.timestamps)
    edges = binsize:binsize:(T+binsize)
    out = zeros(Int,length(edges),nids+1)

    edg = 1
    for (ind,i) in enumerate(res.timestamps)
        if i < edges[edg]
            out[edg,clu.timestamps[ind]]+=1
        else
            edg += 1
        end
    end
    BinnedSpikes(out[:,2:end],binsize/res.fs)
end

function Base.getindex(S::SpikeTimes, id::Int64)
    ind = findfirst(S.spikeIDs.==id)
    return ind === nothing ?  error("No SpikeID $id") : S.spiketimes[ind]
end
function Base.getindex(S::SpikeTimes, ids::Array{Int64,1})
    match = [x in S.spikeIDs for x in ids]
    if all(match)
        return [S.spiketimes[id] for id in ids]
    else
        wrongids = ids[.!(match)]
        error("No SpikeID $wrongids")
    end
end

Base.getindex(S::SpikeTimes, ::Colon) = S.spiketimes

"""
	ratefilt(cellIDs::Array{<:String,1}, spiketimes::SpikeTimes, minrate::Float64, maxrate::Float64)
Only get (non-lick) cells that fire between two `minrate` and `maxrate` (in Hz).
"""
function ratefilt(cellIDs::Array{<:AbstractString,1}, spiketimes::SpikeTimes, minrate::Float64, maxrate::Float64)
    out = [false; cellIDs.!="lick"]
    fs_spk = spiketimes.spiketimes[2].fs
    fs_max = spiketimes.spiketimes[2].maxtime
    minspk,maxspk = [minrate, maxrate] .* (fs_max/fs_spk)
    out[2:end] .= out[2:end] .& (minspk .< length.(spiketimes.spiketimes) .< maxspk)
    out::Array{Bool,1}
end

## this can easily be devectorized for better perf
function shuffle(rng, spikes::SpikeTimes)
    alltimes = merge(spikes.spiketimes...)
    lengths = length.(spikes.spiketimes)
    ranges = range.([1;cumsum(lengths)[1:end-1].+1],lengths)
    temptimes = shuffle(rng,alltimes.timestamps)
    spiketimes = TimeStamps.(sort.(getindex.([temptimes], ranges)),
                            [mintime(spikes)],
                            [maxtime(spikes)],
                            [fs(spikes)])
    SpikeTimes(spiketimes, spikes.spikeIDs)
end

shuffle(spikes::SpikeTimes) =  shuffle(GLOBAL_RNG, spikes)

"""
    addjitter(allspkt::SpikeTimes, jitter::AbstractArray)
Shift each spiketime in `allspkt` by a random element from `jitter`.
"""
function addjitter(allspkt::SpikeTimes, jitter::AbstractArray)
    SpikeTimes([st+rand(jitter,length(st)) for st in allspkt.spiketimes], allspkt.spikeIDs)
end

function findspikesaround(spiketimes,trigt,pre::Real,post::Real)
    spikesaround = Array{Float64,1}[]
    f=1::Int
    @inbounds for t in trigt
        spikes = view(spiketimes,f:length(spiketimes))
        firstsp::Int = searchsortedfirst(spikes, t-pre)-1
        lastsp = searchsortedfirst(spikes, t+post)-2
        push!(spikesaround, spiketimes[firstsp+f:lastsp+f].-t)
        f = firstsp
    end
    return spikesaround
end

function findspikesaround(spiketimes::TimeStamps,trigt::TimeStamps,pre::Real,post::Real, fs=1.)
    TimeStamps.(findspikesaround(fstimes(spiketimes,fs), fstimes(trigt,fs), pre, post),
                floor(Int,-pre),
                ceil(Int,post),
                fs)
end

function inst_rate(ts::TimeStamps, fs=1000.::Float64,
        bandwidth_range=(0,2)::Tuple, f=kde)
    t = mintime(ts):(1/fs):(maxtime(ts)/fs(ts))
    times = realtimes(ts)
    d = f(times, t).density
    return d.*(rate(ts)/mean(d))
end

function inst_rate(tss::Array{TimeStamps{T},1}, fs=1000.::Float64) where T
    inst_rate(merge(tss...))./length(tss)
end

function inst_sterr(tss::Array{TimeStamps{T},1}, fs=1000.::Float64) where T
    std(reduce(hcat,inst_rate.(tss, fs)),2)[:]/sqrt(length(tss))
end

id(x,args...) = x
function binZ(spiketimes::SpikeTimes, binsize=0.02, t=zscore)
    binnedspiketimes = fit.([Histogram], times.(spiketimes[:]),
        [0:binsize*fs(spiketimes):maxtime(spiketimes)], closed=:right)
    Z = t(cat(getfield.(binnedspiketimes,[:weights])...,dims=2),1) # faster than reduce(hcat,...) here
    nan2zero!(Z)
    return Z
end

function binZ3(spiketimes,t,l=-1.0:.01:1.)
    z_ = slpxcorr.(spiketimes[:],[t],[l])
    Z = Float64.(permutedims(cat(3,z_...),(3,1,2)))
    Z .= (Z .- mean(Z,(2,3)))./std(Z,(2,3))
    nan2zero!(Z)
    return Z
end

"""
    function pxcorr(x::AbstractVector, y::AbstractVector, bins::AbstractVector)
Binned cross-correlation of event times
x and y are monotonically increasing vectors of event times, bins is
a monotonically increasing vector that specifies time bins for the
cross-correlation function.

function pxcorr(x::TimeStamps, y::TimeStamps, bins::AbstractVector)
"""
function pxcorr(x::TimeStamps, y::TimeStamps, bins::AbstractVector; auto=false)
	pxcorr(realtimes(x), realtimes(y), bins, auto=auto)::Array{Int64,1}
end

function slpxcorr(x::TimeStamps, y::TimeStamps, bins::AbstractVector)
    slpxcorr(realtimes(x), realtimes(y), bins)::Array{Int64,2}
end

function pxcorr!(c::AbstractVector, x::TimeStamps, y::TimeStamps, bins::AbstractVector; auto=false)
    pxcorr!(c, x.timestamps/x.fs, y.timestamps/y.fs, bins, auto)
end

function pxcorr(x::Array{T,1}, y::T, bins::T; auto=false) where T<:AbstractVector
    counts = zeros(Int, length(bins)-1, length(x))
    for (ind, i) in enumerate(x)
        pxcorr!(view(counts,:,ind), i, y, bins, auto)
    end
    return counts
end

function pxcorr(x::T, y::T, bins::T;auto=false) where T<:AbstractVector
    counts = zeros(Int, length(bins)-1)
    pxcorr!(counts, x, y, bins, auto)
    return counts
end

function pxcorr!(counts::X, x::T, y::T, bins::T, auto::Bool) where {T<:AbstractVector,X<:AbstractVector}
    nx = length(x)
    minbin = minimum(bins)
    maxbin = maximum(bins)
    nbins = length(bins)
    
    for yspk in y
        xlo = searchsortedfirst(x, yspk+minbin, Base.Sort.Forward)
        xhi = searchsortedlast(x, yspk+maxbin, xlo, nx, Base.Sort.Forward)
        xlo = xlo < 1 ? 1 : xlo
        xhi = xhi > nx ? nx : xhi
        xlo <= xhi
        for i = xlo:xhi
            yspkdiff = x[i] - yspk
            ((yspkdiff == 0) & auto) && continue
            bin = searchsortedfirst(bins, yspkdiff)
            bin = bin > nbins ? nbins : bin
            (bin == 1) && continue
            counts[bin-1] += 1
        end
    end
end

function slpxcorr!(counts::AbstractArray,
    x::AbstractVector, y::AbstractVector, bins::AbstractVector)
    nx = length(x)
    minbin = minimum(bins)
    maxbin = maximum(bins)
    nbins = length(bins)
    for (yind,yspk) in enumerate(y)
        xlo = searchsortedfirst(x, yspk+minbin, Base.Sort.Forward)
        xhi = searchsortedlast(x, yspk+maxbin, xlo, nx, Base.Sort.Forward)
        xlo = xlo < 1 ? 1 : xlo
        xhi = xhi > nx ? nx : xhi
        for i = xlo:xhi
            yspkdiff = x[i] - yspk
            bin = searchsortedfirst(bins, yspkdiff)
            bin = bin > nbins ? nbins : bin
            bin == 1 && continue
            counts[bin-1,yind] += 1
        end
    end
end
    
function slpxcorr(x::AbstractVector, y::AbstractVector, bins::AbstractVector)
    counts = zeros(Int, length(bins)-1, length(y))
    slpxcorr!(counts, x, y, bins)
    counts
end

function slpxcorr!(counts::Array{eltype(T),3}, x::Array{T,1}, y::T, bins::AbstractVector) where T<:AbstractVector
    for (indx,dx) in enumerate(x)
        slpxcorr!(view(counts,:,:,indx), dx, y, bins)
    end
end

function slpxcorr(x::Array{T,1}, y::T, bins::AbstractVector) where {T<:AbstractVector}
    counts = zeros(Int, length(bins), length(y), length(x))
    slpxcorr!(counts, x, y, bins)
    counts::Array{Int64,3}
end

function slpxcorr(x::Array{T,1}, y::Array{T,1}, bins::AbstractVector) where {T<:AbstractVector}
    counts = zeros(Int, length(bins), length(y[1]), length(x), length(y))
    for (indy,dy) in enumerate(y)
        for (indx,dx) in enumerate(x)
            slpxcorr!(view(counts,:,:,indx,indy), dx, dy, bins)
        end
    end
    counts
end

"""
	getspiketimesarray(spikematrix)
Get vector of spike times (vectors) from a spike matrix (bit/bool array)
"""
getspiketimesarray(spikematrix) = [find(spikematrix[:,x]) for x in 1:size(spikematrix,2)]

function pxcorr!(cout::Array{Int64,3}, st::Array{Array{Int64,1},1}, bins::AbstractVector, auto::Bool)
    for i in 1:length(st); for j in 1:length(st)
            pxcorr!(view(cout,:,i,j), st[i], st[j], bins, auto)
    end; end
end

"""
Calculate cross-correlegrams between all cells pairs
    Returns lags X cells X cells Array
"""
function pxcorr(st::Array{Array{Int64,1},1}, bins::AbstractVector;auto=false)
    out = zeros(Int, length(bins), length(st), length(st))
    pxcorr!(out, st, bins, auto)
    return out
end

function getallcrosscors(spiketimes::Array{Array{Int64,1},1},
    bins::AbstractVector, normalise::Bool)
    bins=[bins[1]-diff(bins)[1]; collect(bins)]
    allcrosscors = pxcorr(spiketimes,bins)
    normalise && (allcrosscors = znorm(allcrosscors,1);)
    allcrosscors[:,:,:]
end

"""
	function getallcrosscors(spiketimes, bins, normalise)
Gets all the crosscorrelations between list of timestamps.
"""
function getallcrosscors(spikes::SpikeTimes{T,I}, bins, normalise, cells2use=:) where {T,I<:Int}
    return getallcrosscors(getfield.(spikes.spiketimes[cells2use], [:timestamps]),
    bins, normalise)
end
        

"""
	function selectcrosscors!(allcrosscors, usediag, just_lower,
    excludecors(modified))
Selects which cross correlations to keep.
`usediag = true` will keep autocorrelations
`justlower = true` will only use unidirectional crosscorrelations
"""
function selectcrosscors!(allcrosscors, usediag, just_lower, 
    excludecors=falses(allcrosscors[1,:]))
    size(allcrosscors,2)
    excludecors[isnan(allcrosscors[1,:])] = true
    if any(isnan(allcrosscors))
        warn("Crosscors $(sum(isnan(allcrosscors[1,:]))) have NaN values,
            adding to exclude list")
            #... setting to 0")
        #allcrosscors[isnan(allcrosscors)] = 0.
    end
    !usediag && (excludecors[diagind(size(allcrosscors,2,3)...)]=true)
	allcrosscors = just_lower ? mapslices(tril, allcrosscors,(2,3)) : allcrosscors
    excludecors[sum(allcrosscors,1).==0.]=true
    !excludecors
end

# Functions that calculate Cell Spiking Properties
### NB all of these functions operate on a time X cell spikematrix
### spikematrix::Array(Int32, 2) for a single rec session/day


function binspikematrix(spikematrix, binsize)
    edges = 0:binsize:size(spikematrix,1)+binsize
    cluIDs = 1:size(spikematrix,2)
    cat(2,map(x->fit(Histogram, find(spikematrix[:,x]), edges).weights, cluIDs)...)
end

function inint(x::Float64, intv::Array{Float64,2})
	@inbounds @simd for i in 1:size(intv,1)
		if x > intv[i,1] && x < intv[i,2]
			return true
		end
	end
	return false
end

function res_filt(res::Array{T,1}, intv::Array{T,2}) where T<:Real
	isinint::Array{Bool,1} = map(x-> inint(x,intv), res)
end

## Basic Functions

function firingrate(spikematrix, binsize=20)
	toHz = 20000/binsize
	rate = (sum(spikematrix,1)/size(spikematrix,1))*toHz
	return(squeeze(rate,1))
end

"""
    Get spiketimes from spike matrix within a window
"""
getspikes(spikematrix,window) = map(x->find(view(spikematrix,window,x)), 1:size(spikematrix,2))

