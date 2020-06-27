export loaddata_analyze, savedata_analyze, read_analyze_data, read_analyze_hdr

#####################
#  Analyze export
#####################

baremodule ElementType
  FLOAT32 = 16
  INT16 = 4
  RGBA = 2304
  RGB = 128
end

function _bitpix(elemType)
  if elemType == ElementType.FLOAT32
    return 32
  elseif elemType == ElementType.INT16
    return 16
  elseif elemType == ElementType.RGBA
    return 32
  elseif elemType == ElementType.RGB
    return 24
  else
    error("elemType ", elemType, " not yet supported")
  end
end

function savedata_analyze(filename::AbstractString, image::ImageMeta{T};
                          kargs...) where {T<:ColorTypes.TransparentColor}
  rgbimage = convert(RGB{N0f8},image)
  savedata_analyze{T<:TransparentColor}(filename, rgbimage; kargs...)
end

function savedata_analyze(filename::AbstractString, im::ImageMeta; verbose = false, saveas16bit=false, perm=nothing, permRGBData=false)
  filenamebase, ext = splitext(filename)

  im[:permRGBData] = permRGBData

  if get(im,"dim",3) < 3
    warn("voxel size will be wrong for 1D/2D data")
  end
  T = eltype(im)
  if T <: AbstractFloat
    c = arraydata(im)
    im[:cmin] = minimum(c)
    im[:cmax] = maximum(c)
    if saveas16bit
        im[:scalingFactor] = 2^14 / maximum(c)
        c *= im[:scalingFactor]
        c = round(Int16,c)
        im[:datatype] = ElementType.INT16
    else
        im[:scalingFactor] = 1
        if eltype(c) != Float32
          #c = map(Float32,c) #this does not work for franzi
          c = AxisArray(map(f, c.data), c.axes...)
        end
        im[:datatype] = ElementType.FLOAT32
    end
    im[:size] = [size(c,i) for i=1:4]
    im[:ndims] = ndims(c)
  else
      c = arraydata(im).data
      im[:cmin] = 0
      im[:cmax] = 1
      im[:scalingFactor] = 1
    if T <: RGB
      # This is what is returned by colorize or blend
      im[:datatype] = ElementType.RGB
      im[:size] = [size(c,i) for i=1:4]
      im[:ndims] = ndims(c)

    elseif T == UInt8
      # This is what is returned by fusedVolumes
      im[:datatype] = ElementType.RGB
      c = reorderColor(c)
      im[:size] = [size(c,i) for i=1:4]
      im[:ndims] = ndims(c)

    else
      error("Colored data has to be passed as ...")
    end
    if permRGBData

      if ndims(c) == 3
        c = reinterpret(UInt8, c, (3, size(c,1), size(c,2), size(c,3)))
        c = permutedims(c,[2,3,1,4])
      else
        c = reinterpret(UInt8, c, (3, size(c,1), size(c,2), size(c,3), size(c,4)))
        c = permutedims(c,[2,3,1,4,5])
      end
    end
  end

  header = properties(im)
  header[:pixelspacing] = collect(ustrip.(pixelspacing(im)))
  header[:offset] = collect(imcenter(im))


  if ext == ".nii"
    im[:nii_file] = true
    im[:vox_offset] = 352

    open(filenamebase*".nii","w") do fd
      write_analyze_hdr(fd, header)
      write(fd, Int32(0))
      write(fd, c)
    end
  else
    im[:nii_file] = false
    im[:vox_offset] = 0
    # write data
    open(filenamebase*".img","w") do fd
      write(fd, c)
    end

    # write Header
    open(filenamebase*".hdr","w") do fd
      write_analyze_hdr(fd, header)
    end
  end

end

function loaddata_analyze(filename::AbstractString)
  filenamebase, ext = splitext(filename)

  header = nothing
  c = nothing

  if ext == ".nii"
    c = open(filenamebase*".nii","r") do fd
      header = read_analyze_hdr(fd)
      read(fd,Int32)
      read_analyze_data(fd,header)
    end
  else
    open(filenamebase*".hdr","r") do fd
      header = read_analyze_hdr(fd)
    end

    c = open(filenamebase*".img","r") do fd
      read_analyze_data(fd,header)
    end
  end

  pixspacing = header[:pixelspacing]*1000u"mm"
  offset = header[:offset]*1000u"mm" .- 0.5.*pixspacing.*header[:size][1:3] .+ 0.5.*pixspacing
  if ndims(c) == 4
    im = AxisArray(c, (:x,:y,:z,:time), tuple(pixspacing...,one(typeof(header[:pixelspacing][1]))*u"s"),
                                   tuple(offset..., zero(typeof(header[:offset][1]))*u"s"))
  else
    im = AxisArray(c, (:x,:y,:z), tuple(pixspacing...), tuple(offset... ))
  end

  delete!(header, "pixelspacing")
  delete!(header, "timedim")

  imMeta = ImageMeta(im, header)

  return imMeta
end

function read_analyze_data(fd,header)
  if header[:datatype] == ElementType.RGB
    # colored data
    s = header[:size]
    if header[:permRGBData] == false
      c = read!(fd, Array{RGB{N0f8}}(undef, header[:size]...))
    else
      if header[:dim] == 3
        @info "3D data"
        c = read!(fd, Array{Int8}(undef, s[1], s[2], 3, s[3]))
	#c = read(fd, UInt8, s[1], s[2], 3, s[3])
        c = permutedims(c,[3,1,2,4])
        c = reinterpret(RGB{N0f8}, c, (size(c,2),size(c,3),size(c,4)))
      else
        @info "4D data"
        c = read!(fd, Array{Int8}(undef, s[1], s[2], 3, s[3], s[4]))
        #c = read(fd, UInt8, s[1], s[2], 3, s[3], s[4])
        c = permutedims(c,[3,1,2,4,5])
        c = reinterpret(RGB{N0f8}, c, (size(c,2),size(c,3),size(c,4),size(c,5)))
      end
    end
  else
    if header[:datatype] == ElementType.FLOAT32
      c = read!(fd, Array{Float32}(undef, header[:size]...))
    else
      c = map(Float32, read!(fd, Array{Int16}(undef, header[:size]...)) )
    end

    if header[:scalingFactor] > 0
      c *= (1/header[:scalingFactor])
    end
  end
  return c
end

function write_analyze_hdr(fd, header)

    s = header[:pixelspacing].*1000

    write(fd, Int32(348)) #sizeof_hdr::Int32

    data_type = "MPI"
    write(fd, data_type)
    write(fd, zeros(Int8,10-length(data_type))) #char data_type[10];

    write(fd, zeros(Int8,18)) #char db_name[18];
    write(fd, Int32(16384)) #  int extents;
    write(fd, Int16(0)) #short int session_error;
    write(fd, 'r') #char regular;
    write(fd, '0') #char hkey_un0;

    write(fd, Int16(header[:ndims])) #short int dim[8];
    write(fd, Int16(header[:size][1]))
    write(fd, Int16(header[:size][2]))
    write(fd, Int16(header[:size][3]))
    write(fd, Int16(header[:size][4]))
    write(fd, ones(Int16,3))


    write(fd, "mm  ") #vox_units	uchar[4]		spatial units of measure for a voxel
    write(fd, "mmol(Fe)") # cal_units	uchar[8]	name of the calibration unit


    write(fd, zeros(Int16,1)) # unused1	Int16


    write(fd, Int16(header[:datatype])) #short int datatype; -> Float32
    write(fd, Int16(_bitpix(header[:datatype]))) #short int bitpix;
    write(fd, Int16(0)) #short int dim_un0;

    write(fd, Float32(0)) # number of spatial dims
    write(fd, Float32(s[1]))
    write(fd, Float32(s[2]))
    write(fd, Float32(s[3]))
    write(fd, zeros(Float32,4)) # float pixdim[8];


    write(fd, Float32(header[:vox_offset])) #float vox_offset;
    #float funused1;
    #float funused2;

    off = ustrip.(uconvert.(Unitful.mm,header[:offset]))
    write(fd, Float32(off[1]))
    write(fd, Float32(off[2]))
    write(fd, Float32(off[3])) #float funused3;
    write(fd, Float32(header[:scalingFactor])) #float cal_max;  HERE WE STORE THE SCALING FACTOR!
    write(fd, Float32(0)) #float cal_min;
    write(fd, Float32(header[:permRGBData])) #float compressed;
    write(fd, Float32(0)) #float verified;
    #write(fd, round(Int32,header[:cmin]))
    #write(fd, round(Int32,header[:cmax]))#int glmax,glmin;
    write(fd, zero(Int32))
    write(fd, Int32(typemax(Int32)))#int glmax,glmin;

    description = Vector{UInt8}(header[:experimentDescription])
    description = length(description) > 9 ? description[1:9] : description
    @debug "" description
    write(fd, description)
    write(fd, zeros(Int8,80-length(description))) #char descrip[80];

    write(fd, zeros(Int8,24)) #char aux_file[24];
    write(fd, Int8(0)) #char orient;
    write(fd, zeros(Int8,10)) #char originator[10];
    write(fd, zeros(Int8,10)) #char generated[10];

    scannum = Vector{UInt8}(string(header[:experimentNumber]))
    scannum = length(scannum) > 9 ? scannum[1:9] : scannum
    @debug "" scannum
    write(fd, scannum)
    write(fd, zeros(Int8,10-length(scannum))) #char scannum[10];


    subjName = Vector{UInt8}(header[:experimentSubject])
    subjName = length(subjName) > 9 ? subjName[1:9] : subjName
    @debug "" subjName
    write(fd, subjName)
    write(fd, zeros(Int8,10-length(subjName))) #char patient_id[10];

    dateStr, timeStr = split("$(header[:acqStartTime])","T")
    dateStr = prod(split(dateStr,"-"))
    timeStr = split(timeStr,".")[1]
    timeStr = prod(split(timeStr,":"))

    date = Vector{UInt8}(dateStr)
    date = length(date) > 9 ? date[1:9] : date
    write(fd, date)
    write(fd, zeros(Int8,10-length(date))) #char exp_date[10];

    time = Vector{UInt8}(timeStr)
    time = length(time) > 9 ? time[1:9] : time
    write(fd, time)
    write(fd, zeros(Int8,10-sizeof(time))) #char exp_time[10];

    write(fd, zeros(Int8,3)) #char hist_un0[3];
    write(fd, Int32(0)) #int views
    write(fd, Int32(0)) #  int vols_added;
    write(fd, Int32(0)) #  int start_field;
    write(fd, Int32(0)) #  int field_skip;
    write(fd, Int32(0))
    write(fd, Int32(0)) #  int omax, omin;
    write(fd, Int32(0)) #  int smax, smin;
    if !header[:nii_file]
      write(fd, UInt32(0x6E693100))
    else
      write(fd, UInt32(0x6E2B3100))
    end
end

function read_analyze_hdr(fd)
    header = Dict{Symbol,Any}()

    read(fd, Int32) #sizeof_hdr::Int32

    data_type = String( read!(fd, Array{UInt8}(undef, 10)) ) #char data_type[10];
    header[:data_type] = data_type

    read!(fd, Array{Int8}(undef, 18)) #char db_name[18];
    read(fd, Int32) #  int extents;
    read(fd, Int16) #short int session_error;
    read(fd, Int8) #char regular;
    read(fd, Int8) #char hkey_un0;

    dim = read(fd, Int16)
    header[:dim] = dim
    sizes = read!(fd, Array{Int16}(undef, 7)) #short int dim[8];
    header[:size] = sizes[1:dim]

    read!(fd, Array{Int8}(undef, 4)) #vox_units	uchar[4]		spatial units of measure for a voxel
    read!(fd, Array{Int8}(undef, 8)) # cal_units	uchar[8]									name of the calibration unit


    read(fd, Int16) # unused1	Int16

    header[:datatype] = read(fd, Int16) #short int datatype; -> Float32
    read(fd, Int16) #short int bitpix;
    read(fd, Int16) #short int dim_un0;

    read(fd, Float32)
    voxelSize = read!(fd, Array{Float32}(undef, 7)) # float pixdim[8];
    header[:pixelspacing] = voxelSize[1:3] ./ 1000 # this is restricted to 3D...

    if header[:dim] == 4
      header[:timedim] = 4
    else
      header[:timedim] = 0 #This is likely not always correct
    end

    read(fd, Float32) #float vox_offset;
    #float funused1;
    #float funused2;
    off = read!(fd, Array{Float32}(undef,3)) #float funused3;
    header[:offset] = off
    if header[:data_type][1:3]!="MPI"
      header[:offset][:] .= 0.0
    end
    header[:scalingFactor] = read(fd, Float32) #float cal_max; #HERE WE READ THE SCALING FACTOR
    read(fd, Float32) #float cal_min;
    header[:permRGBData] = (read(fd, Float32) != 0.0) #float compressed;
    read(fd, Float32) #float verified;
    cmin = read(fd, Int32)
    header[:cmin] = cmin
    cmax = read(fd, Int32)
    header[:cmax] = cmax #int glmax,glmin;

    descr = String( read!(fd, Array{UInt8}(undef, 80)) ) #char descrip[80];
    header[:experimentDescription] = descr

    read!(fd, Array{UInt8}(undef, 24)) #char aux_file[24];
    read(fd, UInt8) #char orient;
    read!(fd, Array{UInt8}(undef, 10)) #char originator[10];
    read!(fd, Array{UInt8}(undef, 10)) #char generated[10];

    scannum =  String( read!(fd, Array{UInt8}(undef, 10)) ) #char scannum[10];
    header[:experimentNumber] = scannum

    patient_id =  String( read!(fd, Array{UInt8}(undef, 10)) ) #char patient_id[10];
    header[:experimentSubject] = patient_id

    exp_date =  String( read!(fd, Array{UInt8}(undef, 10)) ) #char exp_date[10];
    header[:date] = exp_date

    exp_time =  String( read!(fd, Array{UInt8}(undef, 10)) ) #char exp_time[10];
    header[:time] = exp_time


    read!(fd, Array{UInt8}(undef, 3)) #char hist_un0[3];
    read(fd, Int32) #int views
    read(fd, Int32) #  int vols_added;
    read(fd, Int32) #  int start_field;
    read(fd, Int32) #  int field_skip;
    read(fd, Int32)
    read(fd, Int32) #  int omax, omin;
    read(fd, Int32)
    read(fd, Int32) #  int smax, smin;

  return header
end
