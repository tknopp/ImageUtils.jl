@testset "Testing Transformation Submodule" begin

N = 256
I = reshape(Images.shepp_logan(N,N),N,N,1)


IInterp = interpolateToGrid(I, [1.0,1.0,1.0]./256, [0.0,0.0,0.0], [0.0,0.0,1.0],
                            [2.0,2.0,2.0], [0.4,0.4,0.0], [400,400,1])
exportImage("img/trafo1.png", IInterp[:,:,1])
@testImg("trafo1.png")

IA = makeAxisArray(reshape(I,N,N,1), [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256 )
IInterp2 = interpolateToGrid(IA, [0.0,0.0,1.0], [2.0,2.0,2.0], [0.4,0.4,0.0], [400,400,1])
exportImage("img/trafo2.png", IInterp2[:,:,1])
@testImg("trafo2.png")

IM = ImageMeta(IA, Dict{String,Any}("rotation"=>[0.0,0.0,1.0]))
IInterp3 = interpolateToGrid(IM, [2.0,2.0,2.0], [0.4,0.4,0.0], [400,400,1])
exportImage("img/trafo3.png", IInterp3[:,:,1])
@testImg("trafo3.png")

IT = TransformedArray(reshape(I,N,N,1), [1.0,1.0,1.0]./256, -[0.5,0.5,0.0].+[1.0,1.0,0.0]./256, [0.0,0.0,1.0])
IInterp4 = interpolateToGrid(IT, [2.0,2.0,2.0], [0.4,0.4,0.0], [400,400,1])
exportImage("img/trafo4.png", IInterp4[:,:,1]) #TODO
@testImg("trafo4.png")


end
