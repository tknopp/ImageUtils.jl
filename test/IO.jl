using ImageUtils.AxisArrays
using ImageUtils.Unitful
using ImageUtils.NIfTI
@testset "IO tests" begin

  a = zeros(Float32,10,20,30,1,2,3)

  N = size(a)
  pixelspacing_q = (1,1,1).*u"mm"
  offset=(0,0,0).*u"mm"
  a = AxisArray(a, Axis{:color}(1:N[1]),
		 Axis{:x}(range(offset[1],step=pixelspacing_q[1],length=N[2])),
		 Axis{:y}(range(offset[2],step=pixelspacing_q[2],length=N[3])),
		 Axis{:z}(range(offset[3],step=pixelspacing_q[3],length=N[4])))

  filename_nii = "test.nii"
  saveImage(filename_nii, a)
  b = loadImage(filename_nii)
  @test a == b

  # nifti version <= 0.6.0
  ni = niread(filename_nii)
  pixelspacing = voxel_size(ni.header)
  c = makeAxisArray(ni.raw, [pixelspacing...])
  @test a == c
end