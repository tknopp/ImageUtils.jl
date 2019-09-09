using Interpolations

export translation, minVoxelSize, joinedFOV, joinedGridSize,
       interpolateToGrid, interpolateToCommonGrid, interpolateToRefImage,
       indexFromBGToFG

function translation(im0, im1)
  N = size(im0)

  f0 = fft(im0)
  f1 = fft(im1)
  ir = abs(ifft(f0 .* conj(f1) ./ (abs(f0) .* abs(f1))))
  t = [ind2sub(N, indmax(ir))...] .- 1
  for d=1:length(t)
    if t[d] > N[d] / 2
      t[d] -= N[d]
    end
  end
  t
end





function minVoxelSize(images...)
  minVox = Float64[1e10,1e10,1e10]

  for im in images
    currVoxSize = collect(converttometer(pixelspacing(im)))
    for i=1:3
      minVox[i] = min(minVox[i], currVoxSize[i])
    end
  end
  minVox
end

function joinedFOV(images...)

  fovMin = Float64[1e9,1e9,1e9]
  fovMax = Float64[-1e9,-1e9,-1e9]

  for im in images
    currOffset = collect(converttometer(imcenter(im)))
    for i=1:3
      currFov = size(im,i).* converttometer(pixelspacing(im))[i]
      fovMin[i] = min(fovMin[i], currOffset[i] - currFov/2)
      fovMax[i] = max(fovMax[i], currOffset[i] + currFov/2)
    end
  end

  fov = fovMax .- fovMin
  offset = (fovMax .+ fovMin) ./ 2
  return fov, offset
end

function joinedGridSize(images...)
  fov, offset = joinedFOV(images...)
  minVoxSize = minVoxelSize(images...)
  return ceil(Int, fov ./ minVoxSize)
end

function interpolateToGrid(image,fovOut, offsetOut, gridSizeOut, R )
  interpType = BSpline(Linear())
  imOut = zeros(eltype(image), gridSizeOut[1],gridSizeOut[2],gridSizeOut[3])
  imInterp = extrapolate(interpolate(data(image), interpType),zero(eltype(image)))
  imVoxelSize = collect(converttometer(pixelspacing(image)))
  imSize = [size(image,1),size(image,2),size(image,3)]
  imFov = imSize.*imVoxelSize
  imOffset = collect(converttometer(imcenter(image)))
  offsetPixel = (-fovOut./2 + imFov./2 + offsetOut-imOffset)./ imVoxelSize
  diffSpacing = (fovOut ./ imVoxelSize ) ./ gridSizeOut
  pixelspacingOut = fovOut ./ gridSizeOut
  @debug "" imVoxelSize imSize imFov imOffset offsetPixel diffSpacing pixelspacingOut
  display(R)
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

  offset = (offsetOut .- 0.5.*fovOut .+ 0.5.*pixelspacingOut) .* 1000 .* u"mm"
  imOutAxis = AxisArray(imOut, (:x,:y,:z),tuple((pixelspacingOut .* 1000 .* u"mm")...),tuple(offset...))

  imOut = copyproperties(image, imOutAxis) # create an image object
  return imOut

end

function interpolateToGrid(image, fovOut, offsetOut, gridSizeOut; interpDegree=1)
  imOut = zeros(eltype(image), gridSizeOut[1],gridSizeOut[2],gridSizeOut[3])

  if interpDegree == 1
    interpType = BSpline(Linear())
  elseif interpDegree == 2
    interpType = BSpline(Quadratic(Reflect()))
  else
    interpType = BSpline(Cubic(Reflect()))
  end

  # Interpolations.jl
  imInterp = extrapolate(interpolate(data(image), interpType),zero(eltype(image)))

  imVoxelSize = collect(converttometer(pixelspacing(image)))
  imSize = [size(image,1),size(image,2),size(image,3)]
  imFov = imSize.*imVoxelSize
  imOffset = collect(converttometer(imcenter(image)))
  rotation = get(image, "rotation", [0.0,0.0,0.0])


  offsetPixel = (-fovOut./2 + imFov./2 + offsetOut-imOffset)./ imVoxelSize
  #diffSpacing = (fovOut ./ gridSizeOut) ./ imVoxelSize
  diffSpacing = (fovOut ./ imVoxelSize ) ./ gridSizeOut
  pixelspacingOut = fovOut ./ gridSizeOut

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

  _innerInterpolation!(imOut, imInterp, imSize, gridSizeOut, ax, ay, az, bx, by, bz, cx, cy, cz)

  offset = (offsetOut .- 0.5.*fovOut .+ 0.5.*pixelspacingOut) .* 1000 .* u"mm"
  imOutAxis = AxisArray(imOut, (:x,:y,:z),tuple((pixelspacingOut .* 1000 .* u"mm")...),tuple(offset...))

  imOut = copyproperties(image, imOutAxis) # create an image object
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
end

function interpolateToCommonGrid(refImage, dynImage; kargs...)
  error("TODO: Port to new image infrastructure")
  if ndims(dynImage) == 4
    buf = getindex(dynImage,:,:,:,1)
  else
    buf = dynImage
  end
  # getindexim is like [] but the result is an image object, which is essential
  # in order to preserve the metadata of the image

  fovJoint, offsetJoint = joinedFOV( refImage, buf)
  sizeJoint = joinedGridSize( refImage, buf )

  Uref = interpolateToGrid(refImage,fovJoint,offsetJoint,sizeJoint; kargs...)

  if ndims(dynImage) == 3
    UInterp =  interpolateToGrid(dynImage,fovJoint,offsetJoint,sizeJoint; kargs...)
  else

    U = zeros(eltype(dynImage),sizeJoint...,size(dynImage,4))
    for l=1:size(dynImage,4)
      im = getindex(dynImage,:,:,:,l)
      U[:,:,:,l] = interpolateToGrid(im,fovJoint,offsetJoint,sizeJoint; kargs...)
    end
    UInterp = copyproperties(dynImage,U) #This is not entirely correct. dynImage has a different pixelspacing
    UInterp["timedim"] = 4
    UInterp["pixelspacing"] = fovJoint ./ sizeJoint
  end
  return Uref, UInterp
end

function interpolateToRefImage(background, foreground, R, t; offset=nothing)
  foreground = changeCenter(foreground, t)
  refVoxelSize = collect(converttometer(pixelspacing(background)))
  refSize = [size(background)...][1:3]
  refFov = [size(background)...][1:3].*refVoxelSize
  offset==nothing ? refOffset = collect(converttometer(imcenter(background))) : refOffset=offset
  if ndims(foreground) == 3
    @debug "interpTorefImage" refVoxelSize refSize refFov refOffset
    UInterp =  interpolateToGrid(foreground,refFov,refOffset,refSize,R)
    return UInterp
  end
end

function interpolateToRefImage(refImage, dynImage; offset=nothing,kargs...)
  refVoxelSize = collect(converttometer(pixelspacing(refImage)))
  refSize = [size(refImage)...][1:3]
  refFov = [size(refImage)...][1:3].*refVoxelSize
  offset==nothing ? refOffset = collect(converttometer(imcenter(refImage))) : refOffset=offset

  if ndims(dynImage) == 3
    @debug "interpTorefImage" refVoxelSize refSize refFov refOffset
    UInterp =  interpolateToGrid(dynImage,refFov,refOffset,refSize; kargs...)
  else

    U = zeros(eltype(dynImage),refSize...,size(dynImage,4))
    for l=1:size(dynImage,4)
      im = getindex(dynImage,:,:,:,l)
      U[:,:,:,l] = interpolateToGrid(im,refFov,refOffset,refSize; kargs...)
    end

    UInterpAxis = AxisArray(U,AxisArrays.axes(refImage)..., AxisArrays.axes(dynImage)[4])
    UInterp = copyproperties(dynImage,UInterpAxis) #This is not entirely correct. dynImage has a different pixelspacing
  end
  return UInterp
end

function changeCenter(im::ImageMeta, newcenter)
  newRanges_ = collect(axisvalues(im)[1:3])
  newRanges = [newRanges_[d] .+ newcenter[d]*u"m" .- imcenter(im)[d] for d=1:3]

  if timeaxis(data(im)) == nothing
    im_ = AxisArray(im.data, Axis{:x}(newRanges[1]), Axis{:y}(newRanges[2]), Axis{:z}(newRanges[3]))
  else
    im_ = AxisArray(im.data, Axis{:x}(newRanges[1]), Axis{:y}(newRanges[2]), Axis{:z}(newRanges[3]), timeaxis(im.data))
  end
  return copyproperties(im,im_)
end

function interpolateToRefImage(dataBG, data, params::Dict)

  data_ = Any[]
  for d in data
    d_ = changeCenter(d, [params[:transX]+params[:transBGX],
                    params[:transY]+params[:transBGY],
                    params[:transZ]+params[:transBGZ]])

    d_["rotation"] = [params[:rotX]+params[:rotBGX],
                     params[:rotY]+params[:rotBGY],
                     params[:rotZ]+params[:rotBGZ]]
    push!(data_,d_)
  end

  data_ = [interpolateToRefImage(dataBG,d,interpDegree=1) for d in data_]

  return data_
end

# This version maps the bg on itself
function interpolateToRefImage(dataBG, params::Dict)
  dataBG__ = deepcopy(dataBG)
  dataBG_ = changeCenter(dataBG__, [params[:transBGX], params[:transBGY], params[:transBGZ]])

  dataBG_["rotation"] = [params[:rotBGX], params[:rotBGY], params[:rotBGZ]]

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

  image_ = data[1]

  image = changeCenter(image_, [params[:transX]+params[:transBGX],
                    params[:transY]+params[:transBGY],
                    params[:transZ]+params[:transBGZ]])

  image["rotation"] = [params[:rotX]+params[:rotBGX],
                     params[:rotY]+params[:rotBGY],
                     params[:rotZ]+params[:rotBGZ]]

  refVoxelSize = collect(converttometer(pixelspacing(refImage)))
  refSize = [size(refImage)...][1:3]
  refFov = [size(refImage)...][1:3].*refVoxelSize
  offset==nothing ? refOffset = collect(converttometer(imcenter(refImage))) : refOffset=offset

  fovOut = refFov
  offsetOut = refOffset
  gridSizeOut = refSize

  imVoxelSize = collect(converttometer(pixelspacing(image)))
  imSize = [size(image,1),size(image,2),size(image,3)]
  imFov = imSize.*imVoxelSize
  imOffset = collect(converttometer(imcenter(image)))
  rotation = get(image, "rotation", [0.0,0.0,0.0])


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
