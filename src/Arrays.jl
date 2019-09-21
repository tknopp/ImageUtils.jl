export makeAxisArray, converttometer, imcenter

converttometer(x::NTuple{D,T}) where {D,T<:Quantity} = ustrip.(uconvert.(u"m",x))
converttometer(x::NTuple{D,T}) where {D,T<:AbstractFloat} = x

imcenter(img::AxisArray) = map(x->(0.5*(last(x)+first(x))), 
                               ImageAxes.filter_space_axes(AxisArrays.axes(img), axisvalues(img)))
imcenter(img::ImageMeta) = imcenter(data(img))

function makeAxisArray(array::Array{T,5}, pixelspacing, offset, dt) where T
  N = size(array)
  im = AxisArray(array, Axis{:color}(1:N[1]),
		 Axis{:x}(range(offset[1],step=pixelspacing[1],length=N[2])),
		 Axis{:y}(range(offset[2],step=pixelspacing[2],length=N[3])),
		 Axis{:z}(range(offset[3],step=pixelspacing[3],length=N[4])),
		 Axis{:time}(range(0*unit(dt),step=dt,length=N[5])))
  return im
end

function makeAxisArray(array::Array{T,4}, pixelspacing, offset) where T
  N = size(array)
  im = AxisArray(array, Axis{:color}(1:N[1]),
		 Axis{:x}(range(offset[1],step=pixelspacing[1],length=N[2])),
		 Axis{:y}(range(offset[2],step=pixelspacing[2],length=N[3])),
		 Axis{:z}(range(offset[3],step=pixelspacing[3],length=N[4])))
  return im
end

function makeAxisArray(array::Array{T,3}, pixelspacing, offset) where T
  N = size(array)
  im = AxisArray(array,
		 Axis{:x}(range(offset[1],step=pixelspacing[1],length=N[1])),
		 Axis{:y}(range(offset[2],step=pixelspacing[2],length=N[2])),
		 Axis{:z}(range(offset[3],step=pixelspacing[3],length=N[3])))
  return im
end

