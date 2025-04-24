export ColoringParams, colorize, blend, reorderColor, existing_cmaps, cmap,
       overlay, linearDoge, complexColoring, RGBAGradient, important_cmaps

struct ColoringParams
  cmin::Float64
  cmax::Float64
  cmap::String
end

function ColoringParams(cmin::Real, cmax::Real, i::Int)
  all_maps = existing_cmaps()
  colormap = all_maps[i+1]
  return ColoringParams(cmin, cmax, colormap)
end

#####################
# Weighting function
#####################

"""
    weighted_color_mean(w1, c1, c2)
Returns the color `w1*c1 + (1-w1)*c2` that is the weighted mean of `c1` and
`c2`, where `c1` has a weight 0 ≤ `w1` ≤ 1.
"""
function weighted_color_mean(w1::Real, c1::Colorant, c2::Colorant)
    weight1 = convert(promote_type(eltype(c1), eltype(c2)),w1)
    weight2 = weight1 >= 0 && weight1 <= 1 ? oftype(weight1,1-weight1) : throw(DomainError())
    mapc((x,y)->weight1*x+weight2*y, c1, c2)
end

#####################
# Colorizing an image
#####################

function getBounds( data, coloring, minval, maxval)
  rangeval = maxval - minval
  lowerBound = minval + coloring.cmin * rangeval
  higherBound = minval + coloring.cmax * rangeval
  return lowerBound, higherBound
end

function colorize(data::AbstractArray{T}, coloring::Vector{ColoringParams},
                  minval::Vector, maxval::Vector, blendChannels, activeChannel, complexBlending) where {T<:Real}
  if blendChannels
    if complexBlending
      if size(data,1) != 2
        error("Data needs to be of length 2")
      end
      lowerBound, higherBound = getBounds(data, coloring[1], minval[1], maxval[1]) # scale cmax
      return complexColoring(sliceColorDim(data,1), sliceColorDim(data,2), amax=higherBound)
    else
      return colorize(data, coloring, minval, maxval)
    end
  else
    return colorize(data, coloring, minval, maxval, activeChannel)
  end
end

function colorize(data::AbstractArray{T}, coloring::Vector{ColoringParams},
                 minval::Vector, maxval::Vector) where {T<:Real}
  I = [colorize(sliceColorDim(data,i), coloring[i], minval[i], maxval[i]) for i = 1:size(data,1)]
  return linearDodge(I)
end

function colorize(data::AbstractArray{T}, coloring::Vector{ColoringParams},
                 minval, maxval, chan::Int) where {T<:Real}
  return colorize(sliceColorDim(data,chan), coloring[chan], minval[chan], maxval[chan])
end

function colorize(data::AbstractArray{T}, c::ColoringParams,
                  minval::Real, maxval::Real) where {T<:Real}

  lowerBound, higherBound = getBounds(data, c, minval, maxval)

  cdata = colorize( data, lowerBound, higherBound, cmap(c.cmap),
                    normalize=false)

  return cdata
end


"""
    coloredimage = colorize(inputimage,wmin,wmax,cmap)

This function uses the windowing parameters `w1` and `w2` together with the
color map `cmap` to colorize the `inputimage`. Using the keywordargument
`normalize` the data will be scaled such that all data lies in [0,1].
"""
function colorize(inputimage::AbstractArray, wmin::T, wmax::T, cmap;
                  normalize=true) where {T<:Number}
  if normalize
    image = copy(inputimage)
    # minimum and maximum value of image
    (d_min,d_max)=extrema(image)
    #windowing c -> α ∈ [0,1]		
    if d_max != d_min
      image = (d -> (d-d_min) / (d_max-d_min)).(image)
    elseif d_max != 0
      image = image ./ d_max
    end
  else
    image = inputimage
  end

  return _colorize(image, wmin, wmax, cmap)
end

function _colorize(image::AbstractArray,wmin,wmax,cmap::Vector{T}) where {T<:Colorant}
  C = similar(image, T)
  _colorize!(image,C,wmin,wmax,cmap)
  return C
end

function _colorize!(image,C,wmin,wmax,cmap)
  for l=1:length(image)
    C[l] = _colorize(image[l],wmin,wmax,cmap)
  end
  C
end

function _colorize(x::Number,wmin,wmax,cmap)
  #scaling to discrete color
  L = length(cmap)
  y = (x-wmin) / (wmax-wmin)*(L-1) + 1

  if y <= 1 || wmax == wmin
    c = first(cmap)
  elseif y >= L
    c = last(cmap)
  else
    yidx = floor(Int, y)
    c = weighted_color_mean(y-yidx,cmap[yidx+1],cmap[yidx])
  end
  c
end

##################
#  Linear dodge  #
##################
"""
	  imageC = linearDodge(imageA,imageB)

Blends imageA and imageB using the linear dodge mode, which sums pixelwise with
subsequent clamping.

This function may be used for images of the same size: size(imageA)==size(imageB)
or for a static imageA and a time series of static images imageB.
"""
function linearDodge(images::Vector)
	out = zero(images[1])
	for imag in images
		out .+= imag
	end
	return clamp01!(out)
end

linearDodge(images...) = linearDodge([images...])

function complexColoring(C::Array{T,2}; colormap=ColorSchemes.phase, normalizeG=true, g=0.7, amax=maximum(abs.(C))) where T<:Complex 
  return complexColoring(abs.(C),angle.(C),colormap=colormap, normalizeG=normalizeG, g=g, amax=amax)
end

function complexColoring(amp, phase; colormap=ColorSchemes.phase, normalizeG::Bool=true, g=0.7, amax=maximum(amp))
    if normalizeG
        colormap = ColorScheme(normalizeGray.(colormap,g))
    end
    # normalize amplitude
    I = amp./amax
    clamp01!(I) # using desaturation colormap: make sure to map intensity between 0 and 1
    rawImage = similar(I,RGBA{N0f8})
    for n=1:length(I)
      rawImage[n] = I[n]*get(colormap,phase[n],(-π,Float64(π)))
    end
    return convert.(RGBA{N0f8},rawImage)
end

function overlay(imageBG::AbstractArray{T,D}, imageFG::AbstractArray{U,D},
                  translucentColor) where {T,U,D}
  C = similar(imageFG, RGBA{N0f8})
  overlay!(imageBG, imageFG, C, 0, translucentColor)
  return C
end

# This version operates on a 3D bg and a 4D fg
function overlay(imageBG::AbstractArray{T,3}, imageFG::AbstractArray{U,4},
                  translucentColor) where {T,U}
  C = similar(imageFG, RGBA{N0f8})
  offset = 0
  for l=1:size(C,4)
    overlay!(imageBG, imageFG, C, offset, translucentColor)
    offset += length(imageBG)
  end
  return C
end

function overlay!(imageBG, imageFG, C, offset, translucentColor)
  for l=1:length(imageBG)
    IFG = imageFG[l+offset]
    IBG = imageBG[l]
    C[l+offset] = ( red(IFG) == red(translucentColor) &&
                    green(IFG) == green(translucentColor) &&
		                blue(IFG) == blue(translucentColor)) ? IBG : IFG
  end
  return C
end

##################
# Alpha blending
##################

# This version works on the RGB level and alpha
# is associated to the foreground color.
function blend(bg::C, fg::C, alpha) where {C<:Color}
  #alpha = oftype(eltype(C),alpha)
  output = alpha*fg + (1-alpha)*bg
end

# Type stable?
function blend(bg::C, fg::C) where {C<:TransparentColor}
  α = alpha(bg)
  β  = alpha(fg) * (one(α)-α)
  return mapc((x,y)->α*x + β*y, bg, fg)
end

# this should be replaced by weighted_color_mean(w1, c1, c2) as soon as
# weighted_color_mean supports all colors since that is what it is
function blend(bg::Color, fg::TransparentColor)
  output = blend(bg, color(fg), alpha(fg))
end

# This version operates on two equally dimensioned images
function blend(imageBG::AbstractArray{T,D}, imageFG::AbstractArray{U,D}) where {T,U,D}
  C = similar(imageFG, RGBA{N0f8})
  blend!(imageBG, imageFG, C)
  return C
end

# This version operates on a 3D bg and a 4D fg
function blend(imageBG::AbstractArray{T,3}, imageFG::AbstractArray{U,4}) where {T,U}
  C = similar(imageFG, RGBA{N0f8})
  offset = 0
  for l=1:size(C,4)
    blend!(imageBG, imageFG, C, offset)
    offset += length(imageBG)
  end
  return C
end

function blend!(imageBG, imageFG, C, offset=0)
  for l=1:length(imageBG)
    C[l+offset] = blend(imageBG[l],imageFG[l+offset])
  end
  return C
end

##################
# reorder color
##################

"""
	  reorderColor(image)

separates the color channels of `imgage`.
"""
function reorderColor(image::AbstractArray{T}) where {T<:Colorant}
  imageOut = similar(image, UInt8, size(image)..., 3)
  for l=1:length(image)
    imageOut[l] = reinterpret(UInt8,red(image[l]))
    imageOut[l+length(image)] = reinterpret(UInt8,green(image[l]))
    imageOut[l+2*length(image)] = reinterpret(UInt8,blue(image[l]))
  end
  imageOut
end

function reorderColor(image::AbstractArray{UInt8})
  imageOut = similar(image, RGB{N0f8}, size(image)[1:end-1]...)
  for l=1:length(imageOut)
    imageOut[l] = RGB{N0f8}(reinterpret(N0f8,image[l]),reinterpret(Ufixed8,image[l+length(imageOut)]),reinterpret(N0f8,image[l+2*length(imageOut)]))
  end
  imageOut
end

##################
# Color maps
##################

"""
    `normalizeGray(colormap,g=0.45)`

Transform a color lightness value in Lab space with the aim to output a color
of similar color but with the specified gray value.
"""
function normalizeGray(c::T,g=0.7) where {T<:Colorant}
    cluv = Lab(c)
    # test range
    minG = Gray.(Lab(0,cluv.a,cluv.b)).val
    maxG = Gray.(Lab(100,cluv.a,cluv.b)).val
    if g < minG
      g = minG
      @warn "Normalization of gray value failed, g ≈ $(ceil(minG,digits=3)) is used. Use g >= $(ceil(minG,digits=3)) for a consistent normalization."
    elseif g > maxG
      g = maxG
      @warn "Normalization of gray value failed, g ≈ $(floor(maxG,digits=3)) is used. Use g <= $(floor(maxG,digits=3)) for a consistent normalization."
    end

    l = find_zero(l->Gray.(Lab(l,cluv.a,cluv.b)).val-g, (0, 100))

    return T(Lab(l,cluv.a,cluv.b))
end

"""
    RGBAGradient(colorm,α0,α1)

Maps a linear alpha-value Gradient from α0 to α1 on a colormap.
"""
function RGBAGradient(colorm::Vector{C},α0,α1) where {C<:Colorant}
  n = length(colorm)
  colorm = convert.(RGB{N0f8},colorm)
  αs = range(α0,α1,length=n)
  [RGBA(colorm[i],αs[i]) for i=1:n]
end

"""
    cmap(color)

Creates a colormap from a single `color`, with a color gradient from
transparent black to the opaque `color`.
"""
function cmap(color::T) where {T<:Colorant}
	c1 = RGBA{N0f8}(0,0,0,0)
	c2 = convert(RGBA{N0f8},color)
	return RGBA{N0f8}[c1,c2]
end

"""
    cmap(colors)

Maps a vector of tuples (length 3 or 4) to `RGBA{N0f8}`.
"""
cmap(colors::Vector{T}) where {T<:Tuple} = map(x->RGBA{N0f8}(x...),colors)

"""
    cmap(str::AbstractString)

Chooses colormap determined by `str` from `existing_cmaps()` and returns the colormap with a linear alpha gradient from 0 to 1.
"""
function cmap(str::AbstractString)::Vector{RGBA{Normed{UInt8,8}}}
  if str == "gray"
    return cmap(parse(Colorant,"white"))
  elseif str in ["blue" "green" "red"]
    return cmap(parse(Colorant,str))
  else
    colorm = Symbol(str)
    if colorm in keys(colorschemes)
      return RGBAGradient(colorschemes[colorm].colors,0,1)
    else
      @warn "$colorm is not existing in ColorSchemes.jl. Colorscheme grays is used instead."
      return RGBAGradient(colorschemes[:grays].colors,0,1)
    end
  end
end

function existing_cmaps()
  vcat(["gray", "blue", "green", "red"],String.(collect(keys(ColorSchemes.colorschemes))))
end

function important_cmaps()
  ["gray", "blue", "green", "red", "grays", "viridis", "plasma", "magma", "davos", "lajolla", "delta", "cork", "vik", "phase"]
end


