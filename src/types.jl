export fs,realtimes,fstimes,times,maxtime,mintime,merge,merge_labels
export TimeStamps, Signal, LFP, Tracking, TimeSeries

### Type Layout 

abstract type TimeSeries{T} end

struct Signal{T} <: TimeSeries{T} 
    signal::Vector{T} 
    fs::Float64 
end 

struct LFP{T} <: TimeSeries{T}
    lfp::Vector{Signal{T}}
    ids::Vector{String}
    fs::Float64
end

# this should be a vector of points (or tuples)
# then signal(Tracking) would return that vector
struct Tracking{T} <: TimeSeries{T}
    xpos::Vector{T}
    ypos::Vector{T}
    speed::Vector{T}
    smoothing::Int64
    pixelpercm::Float64
    fs::Float64
end

import Base.length
fs(t::Tracking) = t.fs
length(t::Tracking) = length(t.xpos)

"""
    xy(t::Tracking)
Get x and y locations in a Time by Dim Array
"""
xy(t::Tracking) = cat(getfield.(t, [:xpos; :ypos])...,dims=2)
between(x, lim::Tuple) = (x > lim[1]) & (x < lim[2])

import Base.cat
function cat(s::Vararg{Signal{T}}) where T
    @assert all(fs.(s).==fs(s[1]))
    Signal(cat([x.signal for x in s]...,dims=1),s[1].fs)
end

function cat(s::LFP{T}...) where T
    n = length(s[1].lfp)
    s_out = Vector{Signal{T}}(undef, n)
    for i in 1:n
        s_out[i] = cat([x.lfp[i] for x in s]...)
    end
    LFP(s_out,s[1].ids,s[1].fs)
end

function cat(t::Tracking...)
    constfields = [:smoothing,:pixelpercm, :fs]
    checkfield_eq.([t], constfields)
    constvalues = getfield.([t[1]], constfields)
    catfields   = [:xpos :ypos :speed]
    tmparr      = getfield.(t,catfields)
    vectors     = [cat(tmparr[:,x]...,dims=1) for x in eachindex(catfields)]
    trk_tup     = (vectors..., constvalues...)
    return Tracking(trk_tup...)
end

fs(s::TimeSeries) = s.fs
signal(s::Signal) = s.signal
signal(lfp::LFP) = signal.(lfp.lfp)

abstract type PointProcess{T,N} end

# notes 
# maybe rename to stamp matrix * binnedstamp... 
#   
struct TimeStamps{T}
    timestamps::Array{T,1}
    mintime::Int
    maxtime::Int
    fs::Float64
end

length(ts::TimeStamps) = length(ts.timestamps)

import Base.-, Base.+
function -(ts::TimeStamps,d::T) where T<:Real
    TimeStamps(ts.timestamps .- d*ts.fs, ts.mintime, ts.maxtime, ts.fs)
end
function +(ts::TimeStamps,d::T) where T<:Real
    TimeStamps(ts.timestamps .+ d*ts.fs, ts.mintime, ts.maxtime, ts.fs)
end

function -(ts::TimeStamps,d::Array{T,1}) where T<:Real
    TimeStamps(ts.timestamps .- (d.*ts.fs), ts.mintime, ts.maxtime, ts.fs)
end
function +(ts::TimeStamps,d::Array{T,1}) where T<:Real
    TimeStamps(ts.timestamps .+ (d.*ts.fs), ts.mintime, ts.maxtime, ts.fs)
end


fs(ts::TimeStamps) = ts.fs
realtimes(ts::TimeStamps) = ts.timestamps./ts.fs
fstimes(ts::TimeStamps, fs) = ts.timestamps*(fs/ts.fs)
times(ts::TimeStamps) = ts.timestamps
maxtime(ts::TimeStamps) = ts.maxtime
mintime(ts::TimeStamps) = ts.mintime
rate(ts::TimeStamps) = fs(ts)*(length(ts)/maxtime(ts))

fs(x::Array{T,1}) where T = all(fs.(x).==fs.(x)) ? fs(x[1]) : error()

import Base.merge
function merge(x::TimeStamps...)
    all(fs.(x).==fs.(x)) || error()
    all(maxtime.(x).==maxtime.(x)) || error()
    TimeStamps(sort(cat(times.(x)...,dims=1)), mintime(x[1]), maxtime(x[1]), fs(x[1]))
end

function merge_labels(x::TimeStamps...)
    all(fs.(x).==fs.(x)) || error()
    all(maxtime.(x).==maxtime.(x)) || error()
    labels = [i*ones(Int,j) for (i,j) in enumerate(length.(x))]
    cat(labels..., dims=1)[sortperm(cat(times.(x)...,dims=1))]
end

struct Assembly{P} 
    detectionmethod::Function 
    projectionmatrix::Array{P} 
end 

struct AssemblyActivations{T} <: TimeSeries{T} 
    signal::Signal{T} 
    assembly::Assembly{T} 
end 

abstract type XYGrid{T} end

struct Map{T} <: XYGrid{T} 
    binedges::Array{T,2} 
    binvales::Array{T,2} 
    binsizetocm::Float64 
end

## utils
checkfield_eq(x, fieldname::Symbol) = begin
    @assert all(getfield.(x,fieldname).==getfield(x[1],fieldname))
end
