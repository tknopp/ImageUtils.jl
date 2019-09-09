@testset "Testing General submodule" begin

N = 256
I = reshape(Images.shepp_logan(N,N),N,N,1)

IInterp = interpolateToGrid(I, [1.0,1.0,1.0], [0.0,0.0,0.0], [1.3,0.0,0.0], 
                           [2.0,2.0,2.0], [0.4,0.4,0.4], [400,400,1])
                           
exportImage("img/trafo1.png", IInterp)
#@testImg("trafo1.png")
                           
end
