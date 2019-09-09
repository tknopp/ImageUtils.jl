@testset "Testing General submodule" begin

@test ColoringParams(0, 1, 0) == ColoringParams(0, 1, "gray")

Random.seed!(1)

N = 64 
I = Images.shepp_logan(N,N)

IC = colorize(I, ColoringParams(0,1,"gray"), 0.0, 1.0)
@test typeof(IC) == Array{RGBA{Normed{UInt8,8}},2}

exportImage("img/coloring1.png", IC)
@testImg("coloring1.png")

IC = colorize(I, ColoringParams(0,1,"viridis"), 0.0, 1.0)
exportImage("img/coloring2.png", IC)
@testImg("coloring2.png")

exportImage("img/coloring3.png", I, vmin=0.0, vmax=0.5, colormap="blue", normalize=true)
@testImg("coloring3.png")

# Movie

Itemp = randn(30,30,100)
exportMovie("img/movie.gif", Itemp, vmin=0.0, vmax=1.0, colormap="viridis", normalize=true)
@testImg("movie.gif")

# Overlay

O = zeros(N,N)
O[10:30,10:30] .= 1.0

OC = colorize(O, ColoringParams(0,1,"red"), 0.0, 1.0)

over = overlay(IC, OC, cmap("red")[1])
exportImage("img/coloring4.png", over)
@testImg("coloring4.png")


over = blend(IC, OC)
exportImage("img/coloring5.png", over)
@testImg("coloring5.png")


I = reshape(Images.shepp_logan(N,N),1,N,N)
I = cat(I,I,dims=1)





end
