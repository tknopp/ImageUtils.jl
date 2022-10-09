module ImageUtils

using AxisArrays
using ColorSchemes
using ColorVectorSpace
using FFTW
using ImageMagick
using Interpolations
using LinearAlgebra
using Random
using Reexport
using Roots
using Statistics
using TestImages
using Unitful
using NIfTI

@reexport using Colors
@reexport using FixedPointNumbers
@reexport using Images

const axes = Base.axes

export shepp_logan

include("Arrays.jl")
include("Coloring.jl")
include("Slicing.jl")
include("Transformation.jl")
include("PerfusionMaps.jl")
include("FileExport.jl")
include("Permutation.jl")
include("AnalyzeNifti.jl")
include("Nifti.jl") # TODO: Proper integration

shepp_logan(args...) = Float64.(TestImages.shepp_logan(args...))

end # module
