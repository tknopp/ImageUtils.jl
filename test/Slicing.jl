@testset "Testing Slicing Submodule" begin

I = ones(1,3,3,3,3)

IA = makeAxisArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256, 1.0 )


@test sliceTimeDim(IA, 3) == ones(1,3,3,3)
@test sliceTimeDim(IA, "MIP") == ones(1,3,3,3)

IAL = sliceTimeDim(IA, 3)

x,y,z = threeSlices(IAL,"MIP")
@test x == ones(1,3,3)
@test y == ones(1,3,3)
@test z == ones(1,3,3)
x,y,z = threeSlices(IAL,[1,1,1])
@test x == ones(1,3,3)
@test y == ones(1,3,3)
@test z == ones(1,3,3)


I = ones(3,3,3)

IAL = makeAxisArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256 )

x,y,z = threeSlices(IAL,"MIP")
@test x == ones(3,3)
@test y == ones(3,3)
@test z == ones(3,3)
x,y,z = threeSlices(IAL,[1,1,1])
@test x == ones(3,3)
@test y == ones(3,3)
@test z == ones(3,3)

I = ones(1,3,3,3,3)

IA = makeAxisArray(I, [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256, 1.0 )
IAL = sliceTimeDim(IA, 3)

params = Dict{Symbol,Any}()
params[:sliceX] = 1
params[:sliceY] = 1
params[:sliceZ] = 1
params[:spatialMIP] = true
params[:blendChannels] = true
params[:complexBlending] = false
params[:activeChannel] = 1

xx,yy,zz = getColoredSlices(IAL, nothing, [ColoringParams(0,1,"gray")], [0.0], [1.0], params)
@show typeof(xx)


xx,yy,zz = getColoredSlicesMovie(IA, nothing, [ColoringParams(0,1,"gray")], params)
@show typeof(xx)

end
