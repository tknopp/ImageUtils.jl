export makeAxisArray, converttometer, imcenter, TransformedArray

converttometer(x) = ustrip.(uconvert.(u"m",x))
imcenter(img::AxisArray) = map(x->(0.5*(last(x)+first(x))), ImageAxes.filter_space_axes(AxisArrays.axes(img), axisvalues(img)))
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



mutable struct TransformedArray{T,N} <: AbstractArray{T,N}
  data::Array{T,N}
  spacing::Vector{Float64}
  offset::Vector{Float64}
  rot::Vector{Float64}


  function TransformedArray(data::AbstractArray{T,N},
   spacing::Vector{Float64}, offset::Vector{Float64}=zeros(N), rot::Vector{Float64}=zeros(N)) where {T,N}
  new{T,N}(data, spacing, offset, rot)
  end

end


Base.size(A::TransformedArray) = size(A.data)
Base.axes(A::TransformedArray) = axes(A.data)
datatype(::Type{TransformedArray{T,N}}) where {T,N} = Array{T,N}
Base.IndexStyle(::Type{M}) where {M<:TransformedArray} = IndexStyle(datatype(M))

# getindex and setindex!

@inline function Base.getindex(img::TransformedArray{T,1}, i::Int) where T
  @boundscheck checkbounds(img.data, i)
  @inbounds ret = img.data[i]
  ret
end
@inline function Base.getindex(img::TransformedArray, i::Int)
  @boundscheck checkbounds(img.data, i)
  @inbounds ret = img.data[i]
  ret
end
@inline function Base.getindex(img::TransformedArray{T,N}, I::Vararg{Int,N}) where {T,N}
  @boundscheck checkbounds(img.data, I...)
  @inbounds ret = img.data[I...]
  ret
end

function Base.getindex(img::TransformedArray{T,N}, c, x::Integer, y::Union{Colon,AbstractRange},
                       z::Union{Colon,AbstractRange}, rest...) where {T,N}
  ret = TransformedArray(img.data[c,x,y,z,rest...],
                         img.spacing[2:3], img.offset[2:3], img.rot) # what to do with rot?
  ret
end

function Base.getindex(img::TransformedArray{T,N}, c, x::Union{Colon,AbstractRange}, y::Integer,
                       z::Union{Colon,AbstractRange}, rest...) where {T,N}
  ret = TransformedArray(img.data[c,x,y,z,rest...],
                         img.spacing[1:2:3], img.offset[1:2:3], img.rot) # what to do with rot?
  ret
end

function Base.getindex(img::TransformedArray{T,N}, c, x::Union{Colon,AbstractRange},
                       y::Union{Colon,AbstractRange},
                       z::Integer, rest...) where {T,N}
  ret = TransformedArray(img.data[c,x,y,z,rest...],
                         img.spacing[1:2], img.offset[1:2], img.rot) # what to do with rot?
  ret
end

function Base.getindex(img::TransformedArray{T,3}, x::Integer, y::Union{Colon,AbstractRange},
                       z::Union{Colon,AbstractRange}, rest...) where {T}
  ret = TransformedArray(img.data[x,y,z,rest...],
                         img.spacing[2:3], img.offset[2:3], img.rot) # what to do with rot?
  ret
end

function Base.getindex(img::TransformedArray{T,3}, x::Union{Colon,AbstractRange}, y::Integer,
                       z::Union{Colon,AbstractRange}, rest...) where {T}
  ret = TransformedArray(img.data[x,y,z,rest...],
                         img.spacing[1:2:3], img.offset[1:2:3], img.rot) # what to do with rot?
  ret
end

function Base.getindex(img::TransformedArray{T,3}, x::Union{Colon,AbstractRange},
                       y::Union{Colon,AbstractRange},
                       z::Integer, rest...) where {T}
  ret = TransformedArray(img.data[x,y,z,rest...],
                         img.spacing[1:2], img.offset[1:2], img.rot) # what to do with rot?
  ret
end


@inline function Base.setindex!(img::TransformedArray{T,1}, val, i::Int) where T
  @boundscheck checkbounds(img.data, i)
  @inbounds img.data[i] = val
  val
end
@inline function Base.setindex!(img::TransformedArray, val, i::Int)
  @boundscheck checkbounds(img.data, i)
  @inbounds img.data[i] = val
  val
end
@inline function Base.setindex!(img::TransformedArray{T,N}, val, I::Vararg{Int,N}) where {T,N}
  @boundscheck checkbounds(img.data, I...)
  @inbounds img.data[I...] = val
  val
end


#Base.copy(img::ImageMeta) = ImageMeta(copy(img.data), deepcopy(img.properties))

#Base.convert(::Type{ImageMeta}, A::ImageMeta) = A
#Base.convert(::Type{ImageMeta}, A::AbstractArray) = ImageMeta(A)
#Base.convert(::Type{ImageMeta{T}}, A::ImageMeta{T}) where {T} = A
#Base.convert(::Type{ImageMeta{T}}, A::ImageMeta) where {T} = shareproperties(A, convert(Array{T}, A.data))
#Base.convert(::Type{ImageMeta{T}}, A::AbstractArray) where {T} = ImageMeta(convert(Array{T}, A))

# similar
Base.similar(img::TransformedArray, ::Type{T}, shape::Dims) where {T} =
  TransformedArray(similar(img.data, T, shape), deepcopy(img.rot), deepcopy(img.spacing), deepcopy(img.offset))

ImageCore.pixelspacing(img::TransformedArray) = img.spacing

imcenter(img::TransformedArray) = img.offset .+ 0.5.*img.spacing.*collect(size(img.data)) .- 0.5.*img.spacing
