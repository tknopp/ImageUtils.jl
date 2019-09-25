export TTP, MTT, CBV, CBF, CBF_maxgrad, MTT_fromCBF, lowpasshann, hann, slidingwindowavg

"""
    hann(L::Int)
Hann window of length L, returns array of length L with entries from [0,1]

"""
function hann(L::Int)
  v = zeros(L);
  for i = 1:L
    v[i] = 0.5*(1-cos(2*pi*(i)/(L-1)))
  end
  return v
end

"""
    lowpasshann(tseries[;windowsize=12])
Applies low pass Hann filter pixelwise on the Fourier transformed temporal data. Returns back tranformed array with data
"""
function lowpasshann(tseries::Array;windowsize=12)

  (px,py,pz,pt) = size(tseries)

  hant=zeros(1,1,1,pt)
  
  hant[1,1,1,round(Int,pt/2)-round(Int,round(Int,windowsize))+1:round(Int,pt/2)]=hann(round(Int,round(Int,windowsize)))
  hant[1,1,1,round(Int,pt/2)+1:round(Int,pt/2)+round(Int,round(Int,windowsize))]=hann(round(Int,round(Int,windowsize)))
  hant[1,1,1,round(Int,pt/2)-round(Int,round(Int,windowsize/2))+1:round(Int,pt/2)+round(Int,round(Int,windowsize/2))]=ones(round(Int,round(Int,windowsize)),1)

  tdata=fft(tseries,4)
  tdata=abs.((ifft(broadcast(*,tdata,ifftshift((hant),4)),4)))

  return tdata

end

function lowpasshann(tseries::ImageMeta;windowsize=12)
  odata = data(tseries)
  return lowpasshann(odata,windowsize=12)
end

"""
    slidingwindowavg(tdata::Array[;filterWidth=20])
Applies a sliding window average on the temporal dimension (fourth dimension)
"""
function slidingwindowavg(tdata::Array;filterWidth=20)
	if (filterWidth > size(tdata,4))
	   filterWidth=floor(Int,size(tdata,4)/2)
           @warn "filterWidth was chosen too large. Reduced to floor(Int,size(tdata,4)/2"
	end
  ldata=tdata;
  for t=1:size(tdata,4)
      if t<size(tdata,4)-filterWidth
          ldata[:,:,:,t]=mean(tdata[:,:,:,t:t+filterWidth],dims=4)
      else
          ldata[:,:,:,t]=mean(tdata[:,:,:,end-filterWidth:end],dims=4)
      end
  end
  return ldata
end
"""
    CBF_maxgrad(tseries::ImageMeta[;positionArtery=[1,1,1], alpha=0.4, alpha2=2,windowsize=6])
Calculates the maximum gradient pixelwise along temporal dimension for pixels inside the mask and divides by maximum concentration at the position of the artery. Alternative: CBF which calculates the cerebral blood flow from the ratio CBV/MTT
"""
function CBF_maxgrad(tseries::ImageMeta;positionArtery=[1,1,1], alpha=0.4, alpha2=2,windowsize=6)
  return CBF_maxgrad(tseries.data;positionArtery=positionArtery, alpha=alpha, alpha2=alpha2,windowsize=windowsize)
end

"""
    CBF_maxgrad(tseries::Array[;positionArtery=[1,1,1], alpha=0.4, alpha2=2,windowsize=6])
Calculates the maximum gradient pixelwise along temporal dimension for pixels inside the mask and divides by maximum concentration at the position of the artery. Alternative: CBF which calculates the cerebral blood flow from the ratio CBV/MTT
"""
function CBF_maxgrad(tseries::Array;positionArtery=[1,1,1], alpha=0.4, alpha2=2,windowsize=6)
  mask = generateMaskFromMIP(tseries[:,:,:,:],alpha, alpha2)
  tmaxgrad, maxgrad = maximumGradient(tseries, mask)
  
  tseries = lowpasshann(tseries,windowsize=windowsize)
  println(positionArtery)
  CBF = zeros(size(maxgrad))
  
  for i in CartesianIndices(maxgrad)
    if (mask[i]>0)
      CBF[i] = maxgrad[i]/maximum(tseries[positionArtery[1],positionArtery[2],positionArtery[3],:])	     
    else
      CBF[i] = 0
    end
  end
  return CBF

end
"""
    CBV(tseries[; positionArtery=[1,1,1], alpha = 0.4, alpha2=0.4,windowsize=12])
Calculates the pixelwise sum along the temporal dimension and normalizes to sum at position of artery
"""
function CBV(tseries; positionArtery=[1,1,1], alpha = 0.4, alpha2=0.4,windowsize=12)
  mask = generateMaskFromMIP(tseries[:,:,:,:],alpha, alpha2)   
  tseries = lowpasshann(tseries,windowsize=windowsize)
  dimension = length(size(tseries))	
  return sum(tseries,dims=dimension)/sum(vec(tseries[positionArtery[1],positionArtery[2],positionArtery[3],:])).*mask
end

"""
    firstMoment(tseries::Array)
Calculates the first moment (sum(t*c(t))) pixelwise along temporal dimension (4)
"""
function firstMoment(tseries::Array)

  (px,py,pz,pt) = size(tseries)
  time = 1:pt
   firstMom = zeros(px,py,pz,1)
   for ib = 1:px
     for jb = 1:py
       for kb = 1:pz
	       firstMom[ib,jb,kb,1]=sum(tseries[ib,jb,kb,:].*time)
       end
     end
   end
  
  return firstMom
end
"""
    generateMaskFromMIP(tseries,alpha,alpha2[;windowsize=12])
Not all pixels in the DF FoV experience signal uptake. Mask covers only regions with signal increase during time
"""
function generateMaskFromMIP(tseries,alpha,alpha2;windowsize=12)
  tseries = lowpasshann(tseries,windowsize=windowsize)	
  dataMIP = maximum(tseries,dims=4)
  mask = zeros(size(dataMIP))
  mask2 = zeros(size(dataMIP))
  mask2[dataMIP .> maximum(dataMIP)*alpha2] .= 1
  mask[dataMIP .> mean(tseries[:,:,:,1:5],dims=4)*alpha ] .= 1
 
  return mask.*mask2
end
"""
    MTT(tseries::ImageMeta[; alpha = 0.4, alpha2 =0.4])
Calculates the mean transit time pixelwise from the first moment MTT =(sum(c(t)*t))/sum(c(t)) along the temporal dimension(4)
"""
function MTT(tseries::ImageMeta; alpha = 0.4, alpha2 =0.4,DFPeriodInS=0.02154)
  return MTT(tseries.data; alpha = alpha, alpha2 =alpha2,DFPeriodInS=0.02154)  
end
"""
    MTT(tseries::Array[; alpha = 0.4, alpha2 =0.4])
Calculates the mean transit time pixelwise from the first moment MTT =(sum(c(t)*t))/sum(c(t)) along the temporal dimension(4)
"""
function MTT(tseries::Array; alpha = 0.4, alpha2 =0.4,DFPeriodInS=0.02154)

  mask = generateMaskFromMIP(tseries,alpha, alpha2)	
  firstMom = firstMoment(tseries)
  IntegratedConcentration = sum(tseries,dims=4)
  mtt = firstMom.*mask./IntegratedConcentration
  for i in CartesianIndices(mtt)
	  if (IntegratedConcentration[i] == 0)
                mtt[i] = 0
          end
  end
  return mtt.*DFPeriodInS
end

"""
    TTP(tseries::ImageMeta[;alpha=0.4,alpha2=0.4,windowsize=12])
Calculates the time in number of repetition times (Multiply with T_R of method to get time in s)
"""
function TTP(tseries::ImageMeta; alpha=0.4,alpha2=0.4,windowsize=12,DFPeriodInS=0.02154)
  dataFG = maximum(tseries, dims=4)
  dataFG.data = TTP(tseries.data,alpha=alpha,alpha2=alpha2,windowsize=windowsize)
  return dataFG
end

function TTP(tseries::Array; alpha=0.4,alpha2=0.4,windowsize=12,DFPeriodInS=0.02154)

  dataFG = maximum(tseries, dims=4)
  dataMIP = deepcopy(dataFG)

  ldata = lowpasshann(tseries,windowsize=windowsize)
  #ldata = slidingwindowavg(ldata)
  (px,py,pz,pt) = size(tseries)

   maxdat = zeros(px,py,pz,1)
   ttp = zeros(px,py,pz,1)
   for ib = 1:px
     for jb = 1:py
       for kb = 1:pz
         buf = vec(ldata[ib,jb,kb,:])
         (maxdat[ib,jb,kb,1],ttp[ib,jb,kb,1])=findmax(buf)
       end
     end
   end
  mask=generateMaskFromMIP(tseries[:,:,:,:],alpha,alpha2) 
  dataFG = float(ttp).*mask
  return dataFG.*DFPeriodInS

end

function findFirstEntrySmallerThreshHigherInt(tseries,thresh,minInt)
  comparison = tseries .< thresh
  comparison[1:minInt] .= false
  positionInArray = findmax(comparison)[2]
  if (findmax(comparison)[1] == 0)
    positionInArray = length(tseries)
  end 
  return positionInArray
end

function maximumGradient(data,mask)
  maxgrad,tmaxgrad=maximumGradient(data)
  return tmaxgrad.*mask,maxgrad.*mask
end

function maximumGradient(data)
  ldata = lowpasshann(data,windowsize=6)
  
  (px,py,pz,pt) = size(data)
  normfactors = mean(ldata[:,:,:,:],dims=4)
  normdata = broadcast(-,ldata,normfactors)
  
  tempdiff = zeros(size(normdata)[1],size(normdata)[2],size(normdata)[3],size(normdata)[4]-1)
  for i = 1:size(tempdiff)[4]
    tempdiff[:,:,:,i] = normdata[:,:,:,i+1]-normdata[:,:,:,i]
  end

  tmaxgrad = zeros(px,py,pz,1)
  maxgrad= zeros(px,py,pz,1)
  for ib = 1:px
    for jb = 1:py
      for kb = 1:pz
        buf = vec(tempdiff[ib,jb,kb,:])
        (maxgrad[ib,jb,kb,1],tmaxgrad[ib,jb,kb,1])=findmax(buf)
      end
    end
  end

  return maxgrad,tmaxgrad
end

"""
    findBeginningEndingBolus(tseries;alpha=0.4, alpha2=0.4)
Determines the beginning of the bolus (bolust1) as the time of maximum derivative and the end (bolust2) as the first time the concentration falls short of concentration(bolust1)
"""
function findBeginningEndingBolus(tseries;alpha=0.4, alpha2=0.4)

  ldata = lowpasshann(tseries,windowsize=6)
  
  (px,py,pz,pt) = size(tseries)
  normfactors = mean(ldata[:,:,:,:],dims=4)
  normdata = broadcast(-,ldata,normfactors)
  mask = generateMaskFromMIP(tseries[:,:,:,:],alpha,alpha2)
  maxgrad,tmaxgrad=maximumGradient(tseries)
  
 
  bolust1=minimum(tmaxgrad.*mask+(1 .-mask).*maximum(tmaxgrad))

  tBolusEnd = zeros(px,py,pz,1)
   
  for ib = 1:px
    for jb = 1:py
      for kb = 1:pz
        concentrationthresh = normdata[ib,jb,kb,Int(tmaxgrad[ib,jb,kb,1])]
        tBolusEnd[ib,jb,kb,1] = findFirstEntrySmallerThreshHigherInt(normdata[ib,jb,kb,:],concentrationthresh,Int(tmaxgrad[ib,jb,kb,1]))        
      end
    end
  end
  bolust2=maximum(tBolusEnd.*mask)
   
  return bolust1,bolust2

end


"""
    MTT_fromCBF(tseries::ImageMeta[;positionArtery=[1,1,1], alpha=0.4, alpha2=2])
Calculates the mean transit time (MTT) pixelwise from the ratio CBV/CBF_maxgrad. The resulting values seem to be more reasonable from the function MTT.
"""
function MTT_fromCBF(tseries::ImageMeta;positionArtery=[1,1,1], alpha=0.4, alpha2=2,DFPeriodInS=0.02154)
  return  MTT_fromCBF(tseries.data;positionArtery=positionArtery, alpha=alpha, alpha2=alpha2,DFPeriodInS=0.02154)
end


"""
    MTT_fromCBF(tseries::Array[;positionArtery=[1,1,1], alpha=0.4, alpha2=2])
Calculates the mean transit time (MTT) pixelwise from the ratio CBV/CBF_maxgrad. The resulting values seem to be more reasonable from the function MTT.
"""
function MTT_fromCBF(tseries::Array;positionArtery=[1,1,1], alpha=0.4, alpha2=2,DFPeriodInS=0.02154)
  cbf = CBF_maxgrad(tseries,positionArtery=positionArtery,alpha=alpha,alpha2=alpha2)
  cbv = CBV(tseries, positionArtery=positionArtery, alpha = alpha, alpha2=alpha2)
  mask = generateMaskFromMIP(tseries[:,:,:,:],alpha,alpha2)
  return cbv./((cbf.*mask)+(1 .-mask)).*mask.*DFPeriodInS
end

"""
    CBF(tseries, alpha=0.4, alpha2 = 0.4)
Calculates the cerebral blood flow from the ratio CBV/MTT. Alternative: CBF_maxgrad as the normalized maximum gradient of the bolus. (CBF_maxgrad seems to give more realistic values)
"""
function CBF(tseries; positionArtery=[1,1,1], alpha = 0.4, alpha2=0.4,windowsize=12)
  cbv = CBV(tseries,positionArtery=positionArtery,alpha=alpha,alpha2=alpha2,windowsize=windowsize)
  mtt = MTT(tseries,alpha=alpha,alpha2=alpha2)
  cbf =zeros(size(cbv))
  mask = generateMaskFromMIP(tseries[:,:,:,:],alpha,alpha2)
  cbf = cbv.*mask./(mtt.*mask+(1 .-mask))
  
  return cbf
end
