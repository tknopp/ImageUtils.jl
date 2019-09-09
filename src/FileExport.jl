export exportImage, exportMovie, exportMovies

exportImage(filename, im::ImageMeta; kargs...) = exportImage(filename, im.data; kargs...)

function exportImage(filename, im::AbstractMatrix{T}; vmin=0.0, vmax=1.0,
                              colormap="gray", normalize=true, kargs...) where {T<:Real}
  imC = colorize( im, vmin, vmax, cmap(colormap), normalize=normalize )
  exportImage(filename, imC; kargs...)
end

function exportImage(filename, im::AbstractMatrix{T}; pixelResizeFactor=1) where {T<:Colorant}
  file, ext = splitext(filename)

  if pixelResizeFactor > 1
    im = repeat(im, inner=[pixelResizeFactor,pixelResizeFactor])
  end

  minPxSpacing = minimum(pixelspacing(im))
  newSize = ceil.(Int64, collect(pixelspacing(im)) / minPxSpacing .* collect(size(im)) )

  dataResized = Images.imresize(im,(newSize[1],newSize[2]))
  rgbdata = convert(Array{RGB},dataResized)
  rgbdata = reverse(rgbdata,dims=1)
  save(filename, rgbdata)
end

function exportImage(filename, data::Vector{T}; kargs...) where {T<:AbstractMatrix}
    file, ext = splitext(filename)
    exportImage(file*"_xy.png", data[1]; kargs...)
    exportImage(file*"_xz.png", data[2]; kargs...)
    exportImage(file*"_yz.png", data[3]; kargs...)
end


### export movies ###

exportMovie(filename, im::ImageMeta; kargs...) = exportMovie(filename, im.data; kargs...)

function exportMovie(filename, im::AbstractArray{T,3}; vmin=0.0, vmax=1.0,
                              colormap="gray", normalize=true, kargs...) where {T<:Real}
  imC = colorize( im, vmin, vmax, cmap(colormap), normalize=normalize )
  exportMovie(filename, imC; kargs...)
end

function exportMovie(filename, im::AbstractArray{T,3}; pixelResizeFactor=1) where {T<:Colorant}
  file, ext = splitext(filename)
  filename_ = file*".gif"

  if pixelResizeFactor > 1
    data = repeat(data,inner=[pixelResizeFactor,pixelResizeFactor,1])
  end

  minPxSpacing = minimum(pixelspacing(im))
  newSize = ceil.(Int64, collect(pixelspacing(im)[1:2]) / minPxSpacing .* collect(size(im)[1:2]) )

  datai = similar(im,(newSize[1],newSize[2],size(im,3)))

  for l=1:size(im,3)
    datai[:,:,l] = Images.imresize(im[:,:,l],(newSize[1],newSize[2]))
  end

  rgbdata = convert(Array{RGB},datai)
  rgbdata = reverse(rgbdata,dims=1)
  @info "saving $filename"
  save(filename, rgbdata)
end


function exportMovies(filename, data::Vector; kargs...)
  file, ext = splitext(filename)
  exportMovie(file*"_xy.gif", data[1]; kargs...)
  exportMovie(file*"_xz.gif", data[2]; kargs...)
  exportMovie(file*"_yz.gif", data[3]; kargs...)
end


