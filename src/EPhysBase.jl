module EPhysBase

# export common function definitions
export load

using StatsBase,
Images,
utils,
Interpolations,
NaNMath,
DSP,
Mmap,
Distributed,
Random,
DelimitedFiles,
LinearAlgebra

import Random.shuffle
GLOBAL_RNG = Random.GLOBAL_RNG

include("types.jl")

include("metadata.jl")

include("lfp.jl")

include("spikes.jl")

include("tracking.jl")

include("behaviour.jl")

include("TTLs.jl")

end
