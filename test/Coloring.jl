@testset "Testing General submodule" begin

@test ColoringParams(0, 1, 0) == ColoringParams(0, 1, "gray")


N = 64 
I = Images.shepp_logan(N,N)

IC = colorize(I, ColoringParams(0,1,"gray"), 0.0, 1.0)
@test typeof(IC) == Array{RGBA{Normed{UInt8,8}},2}

exportImage("img/coloring1.png", IC)

IC = colorize(I, ColoringParams(0,1,"viridis"), 0.0, 1.0)
exportImage("img/coloring2.png", IC)

exportImage("img/coloring3.png", I, vmin=0.0, vmax=0.5, colormap="blue", normalize=true)


Itemp = randn(30,30,100)
exportMovie("img/movie.gif", Itemp, vmin=0.0, vmax=1.0, colormap="viridis", normalize=true)

I = reshape(Images.shepp_logan(N,N),1,N,N)
I = cat(I,I,dims=1)





end
