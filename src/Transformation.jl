export interpolateToGrid, interpolateToCommonGrid, interpolateToRefImage,
       indexFromBGToFG


function interpolateToGrid(image::TransformedArray{T,3}, fovOut::Vector{Float64},
               offsetOut::Vector{Float64}, gridSizeOut::Vector{Int64}; kargs...) where T

  pixspacing = collect(pixelspacing(image))
  offset = collect(imcenter(image))

  pixelspacingOut = fovOut ./ gridSizeOut

  imOut = interpolateToGrid(image.data, pixspacing, offset, image.rot, fovOut, offsetOut, gridSizeOut)

  offset = (offsetOut .- 0.5.*fovOut .+ 0.5.*pixelspacingOut)
  imOutAxis = TransformedArray(imOut, pixelspacingOut, offset)
  return imOutAxis
end


function interpolateToGrid(image::ImageMeta{T,3}, fovOut::Vector{Float64},
                offsetOut::Vector{Float64}, gridSizeOut::Vector{Int64}; kargs...) where T

  rot = get(image, "rotation", [0.0,0.0,0.0])
  imOutAxis = interpolateToGrid(image.data, rot, fovOut, offsetOut, gridSizeOut)
  return copyproperties(image, imOutAxis)
end

function interpolateToGrid(image::AxisArray{T,3}, rot::Vector{Float64}, fovOut::Vector{Float64},
                offsetOut::Vector{Float64}, gridSizeOut::Vector{Int64}; kargs...) where T

  #pixspacing = collect(converttometer(pixelspacing(image)))
  #offset = collect(converttometer(imcenter(image)))
  pixspacing = collect(pixelspacing(image))
  offset = collect(imcenter(image))

  pixelspacingOut = fovOut ./ gridSizeOut

  imOut = interpolateToGrid(image.data, pixspacing, offset, rot, fovOut, offsetOut, gridSizeOut)

  #offset = (offsetOut .- 0.5.*fovOut .+ 0.5.*pixelspacingOut) .* 1000 .* u"mm"
  #imOutAxis = AxisArray(imOut, (:x,:y,:z),tuple((pixelspacingOut .* 1000 .* u"mm")...),tuple(offset...))
  offset = (offsetOut .- 0.5.*fovOut .+ 0.5.*pixelspacingOut)
  imOutAxis = AxisArray(imOut, (:x,:y,:z),tuple((pixelspacingOut)...),tuple(offset...))


  return imOutAxis
end

function interpolateToGrid(image::Array{T,3}, pixspacing::Vector{Float64}, offset::Vector{Float64},
                           rot::Vector{Float64}, fovOut::Vector{Float64}, offsetOut::Vector{Float64},
                           gridSizeOut::Vector{Int64}; interpDegree=1) where T

  imOut = zeros(eltype(image), gridSizeOut[1], gridSizeOut[2], gridSizeOut[3])

  if interpDegree == 1
    interpType = BSpline(Linear())
  elseif interpDegree == 2
    interpType = BSpline(Quadratic(Reflect()))
  else
    interpType = BSpline(Cubic(Reflect()))
  end

  tmp = ( size(image,1)==1 ? NoInterp() : interpType,
          size(image,2)==1 ? NoInterp() : interpType,
          size(image,3)==1 ? NoInterp() : interpType )

  # Interpolations.jl
  imInterp = extrapolate(interpolate(image, tmp), zero(eltype(image)))


  imSize = [size(image,1),size(image,2),size(image,3)]
  imFov = imSize.*pixspacing


  offsetPixel = (-fovOut./2 + imFov./2 + offsetOut - offset)./ pixspacing
  #diffSpacing = (fovOut ./ gridSizeOut) ./ pixspacing
  diffSpacing = (fovOut ./ pixspacing ) ./ gridSizeOut
  pixelspacingOut = fovOut ./ gridSizeOut

  r = rot
  Rx = [1 0 0; 0 cos(r[1]) -sin(r[1]); 0 sin(r[1]) cos(r[1])]
  Ry = [cos(r[2]) 0 sin(r[2]); 0 1 0; -sin(r[2]) 0 cos(r[2])]
  Rz = [cos(r[3]) -sin(r[3]) 0; sin(r[3]) cos(r[3]) 0; 0 0 1]

  R = Rz*Ry*Rx

  ax = 0.5+gridSizeOut[1]/2
  ay = 0.5+gridSizeOut[2]/2
  az = 0.5+gridSizeOut[3]/2

  bx = R[:,1]*diffSpacing[1]
  by = R[:,2]*diffSpacing[2]
  bz = R[:,3]*diffSpacing[3]

  cx = gridSizeOut[1]/2*diffSpacing[1]+0.5+offsetPixel[1]
  cy = gridSizeOut[2]/2*diffSpacing[2]+0.5+offsetPixel[2]
  cz = gridSizeOut[3]/2*diffSpacing[3]+0.5+offsetPixel[3]

  _innerInterpolation!(imOut, imInterp, imSize, gridSizeOut, ax, ay, az, bx, by, bz, cx, cy, cz)

  return imOut
end

function _innerInterpolation!(imOut, imInterp, imSize, gridSizeOut, ax, ay, az, bx, by, bz, cx, cy, cz)
  for iz=1:gridSizeOut[3]
    for iy=1:gridSizeOut[2]
      v = (iz-az) .* bz + (iy-ay) .* by - ax .* bx
      for ix=1:gridSizeOut[1]
         for d=1:3
           v[d] += bx[d]
         end
         x,y,z = v[1]+cx, v[2]+cy, v[3]+cz
         imOut[ix,iy,iz] = imInterp(x, y, z)
      end
    end
  end
  return
end



### interpolateToRefImage ###


function interpolateToRefImage(refImage::TransformedArray{T,3}, dynImage::TransformedArray{T,3}; offset=nothing,kargs...) where T

  refVoxelSize = collect(pixelspacing(refImage))
  refSize = collect(size(refImage))[1:3]
  refFov = collect(size(refImage))[1:3].*refVoxelSize
  offset==nothing ? refOffset = collect(imcenter(refImage)) : refOffset=offset

  UInterp = interpolateToGrid(dynImage,refFov,refOffset,refSize; kargs...)

  return UInterp
end

function interpolateToRefImage(refImage::TransformedArray{T,3}, dynImage::TransformedArray{T,4}; offset=nothing,kargs...) where T
  #refVoxelSize = collect(converttometer(pixelspacing(refImage)))
  refVoxelSize = collect(pixelspacing(refImage))
  refSize = [size(refImage)...][1:3]
  refFov = [size(refImage)...][1:3].*refVoxelSize
  #offset==nothing ? refOffset = collect(converttometer(imcenter(refImage))) : refOffset=offset
  offset==nothing ? refOffset = collect(imcenter(refImage)) : refOffset=offset

  U = zeros(eltype(dynImage),refSize...,size(dynImage,4))
  for l=1:size(dynImage,4)
    im = getindex(dynImage,:,:,:,l)
    U[:,:,:,l] = interpolateToGrid(im,refFov,refOffset,refSize; kargs...)
  end

  UInterpAxis = AxisArray(U,AxisArrays.axes(refImage)..., AxisArrays.axes(dynImage)[4])
  UInterp = copyproperties(dynImage,UInterpAxis) #This is not entirely correct. dynImage has a different pixelspacing

  return UInterp
end


function interpolateToRefImage(dataBG::TransformedArray{T,3}, data::TransformedArray{T,4}, params::Dict) where T
     
  data.offset = [data.offset[1]+params[:transX]+params[:transBGX],
                    data.offset[1]+params[:transY]+params[:transBGY],
                    data.offset[1]+params[:transZ]+params[:transBGZ]]

  data.rot = [params[:rotX]+params[:rotBGX],
                     params[:rotY]+params[:rotBGY],
                     params[:rotZ]+params[:rotBGZ]]

  out = TransformedArray(zeros(T,size(data,1),size(dataBG,1),
                         size(dataBG,2), size(dataBG,3)))
  for c=1:size(out,1)
    out[c,:,:,:] = interpolateToRefImage(dataBG, data[c,:,:,:], interpDegree=1)
  end

  return out
end

# This version maps the bg on itself
function interpolateToRefImage(dataBG, params::Dict)
  dataBG_ = deepcopy(dataBG)
  #dataBG_.offset = changeCenter(dataBG__, [params[:transBGX], params[:transBGY], params[:transBGZ]])
  dataBG_.offset = [dataBG_.offset[1]+params[:transBGX], 
                    dataBG_.offset[2]+params[:transBGY], 
                    dataBG_.offset[3]+params[:transBGZ]]


  dataBG_.rot = [params[:rotBGX], params[:rotBGY], params[:rotBGZ]]

  dataBGInterp = interpolateToRefImage(dataBG, dataBG_, interpDegree=1)

  return dataBGInterp
end







# difficult to dispatch currently
function interpolateToRefImageAllFr(dataBG, im, params::Dict)

  L = size(im[1],4)

  data_all = zeros(eltype(im[1]),size(dataBG,1),size(dataBG,2),size(dataBG,3),L)

  data_ = [getindex(d,:,:,:,1) for d in im]
  dataInterp = interpolateToRefImage(dataBG, data_, params)
  data_all[:,:,:,1] = dataInterp[1]


  for l=2:L
    data_ = [getindex(d,:,:,:,l) for d in im]
    dataInterp = interpolateToRefImage(dataBG, data_, params)
    data_all[:,:,:,l] = dataInterp[1]
  end

  dataAllAxis = AxisArray(data_all,AxisArrays.axes(dataInterp[1])..., AxisArrays.axes(im[1])[4])
  return copyproperties(im[1], dataAllAxis)
end

















function indexFromBGToFG(refImage, data, params::Dict; offset=nothing)

  image = data[1,:,:,:]

  #image.offset = changeCenter(image_, [params[:transX]+params[:transBGX],
  #                  params[:transY]+params[:transBGY],
  #                  params[:transZ]+params[:transBGZ]])
                    
  image.offset = [image.offset[1]+params[:transX]+params[:transBGX],
                    image.offset[1]+params[:transY]+params[:transBGY],
                    image.offset[1]+params[:transZ]+params[:transBGZ]]

  image.rot = [params[:rotX]+params[:rotBGX],
                     params[:rotY]+params[:rotBGY],
                     params[:rotZ]+params[:rotBGZ]]

  refVoxelSize = collect(pixelspacing(refImage))
  refSize = [size(refImage)...][1:3]
  refFov = [size(refImage)...][1:3].*refVoxelSize
  offset==nothing ? refOffset = collect(imcenter(refImage)) : refOffset=offset

  fovOut = refFov
  offsetOut = refOffset
  gridSizeOut = refSize

  imVoxelSize = collect(pixelspacing(image))
  imSize = [size(image,1),size(image,2),size(image,3)]
  @show imSize size(imSize) imVoxelSize size(imVoxelSize)
  imFov = imSize.*imVoxelSize
  imOffset = collect(imcenter(image))
  rotation = image.rot #get(image, "rotation", [0.0,0.0,0.0])


  offsetPixel = (-fovOut./2 + imFov./2 + offsetOut-imOffset)./ imVoxelSize
  #diffSpacing = (fovOut ./ gridSizeOut) ./ imVoxelSize
  diffSpacing = (fovOut ./ imVoxelSize ) ./ gridSizeOut

  r = rotation
  Rx = [1 0 0; 0 cos(r[1]) -sin(r[1]); 0 sin(r[1]) cos(r[1])]
  Ry = [cos(r[2]) 0 sin(r[2]); 0 1 0; -sin(r[2]) 0 cos(r[2])]
  Rz = [cos(r[3]) -sin(r[3]) 0; sin(r[3]) cos(r[3]) 0; 0 0 1]

  R = Rz*Ry*Rx

  ax = 0.5+gridSizeOut[1]/2
  ay = 0.5+gridSizeOut[2]/2
  az = 0.5+gridSizeOut[3]/2

  bx = R[:,1]*diffSpacing[1]
  by = R[:,2]*diffSpacing[2]
  bz = R[:,3]*diffSpacing[3]

  cx = gridSizeOut[1]/2*diffSpacing[1]+0.5+offsetPixel[1]
  cy = gridSizeOut[2]/2*diffSpacing[2]+0.5+offsetPixel[2]
  cz = gridSizeOut[3]/2*diffSpacing[3]+0.5+offsetPixel[3]

  ix, iy, iz = params[:sliceX],params[:sliceY],params[:sliceZ]


  v = (iz-az) .* bz + (iy-ay) .* by - ax .* bx

  for d=1:3
    v[d] += bx[d].*ix
  end
  x,y,z = v[1]+cx, v[2]+cy, v[3]+cz

  @debug "" x y z

  if x>=1 && y >=1 && z>= 1 && x<= imSize[1] && y<= imSize[2] && z<= imSize[3]
    return (round(Int64,x),round(Int64,y),round(Int64,z))
  else
    return (0,0,0)
  end

end
