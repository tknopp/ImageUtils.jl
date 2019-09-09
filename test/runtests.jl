using ImageUtils
using Test
using FileIO
using Random

mkpath("img/")
mkpath("correct/")

macro testImg(filename)
  return :(
    im1 = load(joinpath("img", $filename));
    im2 = load(joinpath("correct", $filename));
    @test im1 == im2
    )
end

include("Coloring.jl")

