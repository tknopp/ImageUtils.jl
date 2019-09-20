export sliceColorDim, sliceTimeDim, sliceSpatialDim, threeSlices,
       getColoredSlices, getColoredSlicesMovie

### sliceColorDim ###

sliceColorDim(A::AbstractArray{T,3}, proj::Integer) where T = A[proj,:,:]
sliceColorDim(A::AbstractArray{T,4}, proj::Integer) where T = A[proj,:,:,:]

### sliceTimeDim ###

function sliceTimeDim(A::AbstractArray{T,D}, proj::String) where {T,D}
  if  proj == "MIP"
    return dropdims(maximum(A, dims=D), dims=D)
  elseif proj == "TTP"
    error("not yet implemented!")
    #data_ = [maximum(d, 4) for d in m.data]
    #data_ = [  mip(d,4) for d in m.data]
  #elseif params[:frameProj] == 2 && ndims(m.data[1]) == 4
  #data_ = [timetopeak(d, alpha=params[:TTPThresh], alpha2=params[:TTPThresh]) for d in m.data]
  end
end

function sliceTimeDim(A::AbstractArray{T,5}, proj::Integer) where T
  return getindex(A,:,:,:,:,proj)
end

function sliceTimeDim(A::AbstractArray{T,4}, proj::Integer) where T
  return getindex(A,:,:,:,proj)
end

### sliceSpatialDim ###

function sliceSpatialDim(A::AbstractArray{T,3}, proj::Union{Tuple,Vector}, dim::Integer) where T
  if dim == 1
    return A[proj[1],:,:]
  elseif dim == 2
    return A[:,proj[2],:]
  else
    return A[:,:,proj[3]]
  end
end

function sliceSpatialDim(A::AbstractArray{T,4}, proj::Union{Tuple,Vector}, dim::Integer) where T
  if dim == 1
    return A[:,proj[1],:,:]
  elseif dim == 2
    return A[:,:,proj[2],:]
  else
    return A[:,:,:,proj[3]]
  end
end

function sliceSpatialDim(A::AbstractArray{T,D}, proj::String, dim::Integer) where {T,D}
  d = (D == 4) ? dim + 1 : dim
  if proj == "MIP"
    return dropdims(maximum(A, dims=d),dims=d)
  elseif proj == "SUM"
    return dropdims(sum(A, d),dims=d)
  else
    error("proj=$proj not available!")
  end
end

### threeSlices ###

function threeSlices(data::AbstractArray{T,3}, proj) where T
  zx = reverse(permutedims(sliceSpatialDim(data, proj, 2),(2,1)),dims=2)
  zy = permutedims(sliceSpatialDim(data, proj, 1), (2,1))
  xy = reverse(sliceSpatialDim(data, proj, 3), dims=1)
  return zx, zy, xy
end

function threeSlices(data::AbstractArray{T,4}, proj) where T
  zx = sliceSpatialDim(data, proj, 2)
  zy = sliceSpatialDim(data, proj, 1)
  xy = sliceSpatialDim(data, proj, 3)
  return _changeDims(zx, zy, xy)
end

function _changeDims(zx, zy, xy)
  zx_ = reverse(permutedims(zx,(1,3,2)),dims=3)
  zy_ = permutedims(zy, (1,3,2))
  xy_ = reverse(xy, dims=2)
  return zx_, zy_, xy_
end

### getColoredSlices ###

function getColoredSlices(data::AbstractArray{T,4}, dataBG, coloring,
                          minval, maxval, params) where {T}
  slices = (params[:sliceX],params[:sliceY],params[:sliceZ])
  proj = params[:spatialMIP] ? "MIP" : slices
  
  return getColoredSlices(data, dataBG, coloring, minval, maxval, proj, params[:blendChannels], params[:activeChannel], params[:coloringBG], get(params,:hideFG, false), get(params,:hideBG, false),
      get(params,:translucentBlending, false), params[:spatialMIPBG] ? "MIP" : slices)                          
end

function getColoredSlices(data::AbstractArray{T,4}, dataBG, coloring,
                          minval, maxval, proj, blendChannels, 
                          activeChannels, coloringBG, hideFG, hideBG, 
                          translucentBlending, projBG) where {T}

  zx, zy, xy = threeSlices(data, proj)

  #zx = data[:,slices[1],:,:]
  #zy = data[:,:,slices[2],:]
  #xy = data[:,:,:,slices[3]]

  cdata_zx = colorize(zx,coloring,minval,maxval,blendChannels,activeChannels)
  cdata_zy = colorize(zy,coloring,minval,maxval,blendChannels,activeChannels)
  cdata_xy = colorize(xy,coloring,minval,maxval,blendChannels,activeChannels)

  if !isempty(dataBG)
    zxBG, zyBG, xyBG = threeSlices(dataBG, projBG)  #no MIP for BG data!!!

    minval,maxval = extrema(dataBG)
    cdataBG_zx = colorize(zxBG,coloringBG,minval,maxval)
    cdataBG_zy = colorize(zyBG,coloringBG,minval,maxval)
    cdataBG_xy = colorize(xyBG,coloringBG,minval,maxval)

    if !hideFG && !hideBG
      if translucentBlending
        cdata_zx = blend(cdataBG_zx, cdata_zx)
        cdata_zy = blend(cdataBG_zy, cdata_zy)
        cdata_xy = blend(cdataBG_xy, cdata_xy)
      else
        #all_maps = existing_cmaps()
        colormap = cmap( coloring[1].cmap)

        cdata_zx = overlay(cdataBG_zx, cdata_zx, first(colormap))
        cdata_zy = overlay(cdataBG_zy, cdata_zy, first(colormap))
        cdata_xy = overlay(cdataBG_xy, cdata_xy, first(colormap))
      end
    elseif hideFG
      cdata_zx = cdataBG_zx
      cdata_zy = cdataBG_zy
      cdata_xy = cdataBG_xy
    end
  end

  return cdata_zx, cdata_zy, cdata_xy
end

export getInterpAndColoredSlices

function getInterpAndColoredSlices(data, dataBG, coloring, minval,maxval, params)
  if dataBG != nothing
    dataFG = interpolateToRefImage(dataBG, data, params)
    dataBG = interpolateToRefImage(dataBG, params)
  else
    dataFG = data
    dataBG = nothing
  end

  cdata_zx, cdata_zy, cdata_xy = getColoredSlices(dataFG, dataBG, coloring, minval, maxval, params)
  return cdata_zx, cdata_zy, cdata_xy
end

function getColoredSlicesMovie(data, dataBG, coloring, params)
    L = size(data[1],4)

    maxval = [maximum(d) for d in data]
    minval = [minimum(d) for d in data]

    data_ = [getindex(d,:,:,:,1) for d in data]
    cdata_zx, cdata_zy, cdata_xy = getInterpAndColoredSlices(data_, dataBG, coloring, minval,maxval, params)

    cdata_xy_all = AxisArray(zeros(eltype(cdata_xy),size(cdata_xy,1),size(cdata_xy,2),L),AxisArrays.axes(cdata_xy)..., AxisArrays.axes(data[1])[4])
    cdata_zx_all = AxisArray(zeros(eltype(cdata_zx),size(cdata_zx,1),size(cdata_zx,2),L),AxisArrays.axes(cdata_zx)..., AxisArrays.axes(data[1])[4])
    cdata_zy_all = AxisArray(zeros(eltype(cdata_zy),size(cdata_zy,1),size(cdata_zy,2),L),AxisArrays.axes(cdata_zy)..., AxisArrays.axes(data[1])[4])

    cdata_xy_all[:,:,1] = cdata_xy
    cdata_zy_all[:,:,1] = cdata_zy
    cdata_zx_all[:,:,1] = cdata_zx

    for l=2:L
      data_ = [getindex(d,:,:,:,l) for d in data]

      cdata_zx, cdata_zy, cdata_xy = getInterpAndColoredSlices(data_, dataBG, coloring, minval, maxval, params)

      cdata_xy_all[:,:,l] = cdata_xy
      cdata_zy_all[:,:,l] = cdata_zy
      cdata_zx_all[:,:,l] = cdata_zx
    end
    return ImageMeta[copyproperties(data[1],cdata_xy_all), copyproperties(data[1],cdata_zx_all), copyproperties(data[1],cdata_zy_all)]

end
