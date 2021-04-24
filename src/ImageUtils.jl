module ImageUtils

using ImageMagick
using Reexport
using ColorVectorSpace
using AxisArrays
const axes = Base.axes
using Random
using Statistics
using LinearAlgebra
using Interpolations
using Unitful
using FFTW
import TestImages

@reexport using Colors
@reexport using FixedPointNumbers
@reexport using Images

export shepp_logan

include("Arrays.jl")
include("Coloring.jl")
include("Slicing.jl")
include("Transformation.jl")
include("PerfusionMaps.jl")
include("FileExport.jl")
include("Permutation.jl")
include("AnalyzeNifti.jl")

shepp_logan(args...) = Float64.(TestImages.shepp_logan(args...))

end # module
