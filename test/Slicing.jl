@testset "Testing Slicing submodule" begin

I = ones(1,3,3,3,3)

IA = makeAxisArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256, 1.0 )
IT = TransformedArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256, [0.0,0.0,1.0])

for U in [IA,IT]
  @test sliceTimeDim(U, 3) == ones(1,3,3,3)
  @test sliceTimeDim(U, "MIP") == ones(1,3,3,3)
end

IAL = sliceTimeDim(IA, 3)
ITL = sliceTimeDim(IT, 3)

for U in [IAL,ITL]
  x,y,z = threeSlices(U,"MIP")
  @test x == ones(1,3,3)
  @test y == ones(1,3,3)
  @test z == ones(1,3,3)
  x,y,z = threeSlices(U,[1,1,1])
  @test x == ones(1,3,3)
  @test y == ones(1,3,3)
  @test z == ones(1,3,3)
end


I = ones(3,3,3)

IAL = makeAxisArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256 )
ITL = TransformedArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256, [0.0,0.0,1.0])

for U in [IAL,ITL]
  x,y,z = threeSlices(U,"MIP")
  @test x == ones(3,3)
  @test y == ones(3,3)
  @test z == ones(3,3)
  x,y,z = threeSlices(U,[1,1,1])
  @test x == ones(3,3)
  @test y == ones(3,3)
  @test z == ones(3,3)
end

I = ones(1,3,3,3,3)

IA = makeAxisArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256, 1.0 )
IT = TransformedArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256, [0.0,0.0,1.0])
IAL = sliceTimeDim(IA, 3)
ITL = sliceTimeDim(IT, 3)

params = Dict{Symbol,Any}()
params[:sliceX] = 1
params[:sliceY] = 1
params[:sliceZ] = 1
params[:spatialMIP] = true
params[:blendChannels] = true
getColoredSlices(ITL, nothing, nothing, [ColoringParams(0,1,"gray")], [0.0], [1.0], params)


end
