@testset "Testing Coloring Submodule" begin

@test ColoringParams(0, 1, 0) == ColoringParams(0, 1, "gray")

Random.seed!(1)

N = 64
I = shepp_logan(N)

IC = colorize(I, ColoringParams(0,1,"gray"), 0.0, 1.0)
@test typeof(IC) == Array{RGBA{Normed{UInt8,8}},2}

exportImage("img/coloring1.png", IC)
@testImg("coloring1.png")

IC = colorize(I, ColoringParams(0,1,"viridis"), 0.0, 1.0)
exportImage("img/coloring2.png", IC)
@testImg("coloring2.png")

exportImage("img/coloring3.png", I, vmin=0.0, vmax=0.5, colormap="vangogh", normalize=true)
@testImg("coloring3.png")

# Movie

Itemp = Array(reshape(range(0,1,length=90000),30,30,100))
exportMovie("img/movie.gif", Itemp, vmin=0.0, vmax=1.0, colormap="viridis", normalize=true)
@testImg("movie.gif")

exportMovies("img/movie.gif", [Itemp,Itemp,Itemp])
@testImg("movie_xy.gif")
@testImg("movie_xz.gif")
@testImg("movie_yz.gif")


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


I = reshape(shepp_logan(N),1,N,N)
O = reshape(O,1,N,N)
I = cat(I,O,dims=1)

OC = colorize(I, [ColoringParams(0,1,"gray"), ColoringParams(0,1,"red"),],
     [0.0,0.0], [1.0,1.0])

exportImage("img/coloring6.png", OC)
@testImg("coloring6.png")


## Complex coloring

# create 2D colorbar
hues = range(0,2π,length=50)
saturations = 0:0.05:1.0
HS = vec(tuple.(hues, saturations'))
D = kron(saturations,(exp.(im*hues))')
# complex coloring
Y = complexColoring(abs.(D), angle.(D))
exportImage("img/coloring7.png", Y)
@testImg("coloring7.png")

# normalizeGray
colormap = ImageUtils.ColorSchemes.phase
@test ImageUtils.normalizeGray(colormap[1]) == RGB{Float64}(0.8979049716100672,0.6795060499249219,0.28645928277939586)
@test @test_logs (:warn,"Normalization of gray value failed, g ≈ 0.056 is used. Use g >= 0.056 for a consistent normalization.") ImageUtils.normalizeGray(colormap[1],0.05) == RGB{Float64}(0.18396207457891003,0.0,0.0)
@test @test_logs (:warn,"Normalization of gray value failed, g ≈ 0.925 is used. Use g <= 0.925 for a consistent normalization.") ImageUtils.normalizeGray(colormap[1],0.95) == RGB{Float64}(1.0,0.9603124734843839,0.5548057892335283)

end
