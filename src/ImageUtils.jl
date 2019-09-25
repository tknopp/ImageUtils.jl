module ImageUtils

using ImageMagick
using Reexport
#using FileIO
#import FileIO.save
using ColorVectorSpace
using AxisArrays
const axes = Base.axes
using Random
using Statistics
using LinearAlgebra
using Interpolations
using ImageCore
using Unitful
using FFTW

@reexport using Colors
@reexport using FixedPointNumbers
@reexport using Images
@reexport using ImageMetadata
@reexport using ImageAxes

include("Arrays.jl")
include("Coloring.jl")
include("Slicing.jl")
include("Transformation.jl")
include("PerfusionMaps.jl")
include("FileExport.jl")
include("Permutation.jl")
include("AnalyzeNifti.jl")

end # module
