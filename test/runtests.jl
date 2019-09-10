using ImageUtils
using Test
#using FileIO
using Random
using ImageMagick

mkpath("img/")
mkpath("correct/")

macro testImg(filename)
  return :(
    im1 = ImageMagick.load(joinpath("img", $filename));
    im2 = ImageMagick.load(joinpath("correct", $filename));
    @test im1 == im2
    )
end

include("Coloring.jl")
include("Slicing.jl")
include("Transformation.jl")
