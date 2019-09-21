export applyPermutions, applyPermutionsRev
export permuteCombinations, permuteCombinationsName, flippings, flippingsName, convertMode2PermFlip

function convertMode2PermFlip(mode)
  if mode == 0
    perm = [2,1,3]
    flip = [1,2,3]
  elseif mode == 1
    perm = [2,1,3]
    flip = [2]
  elseif mode == 2
    perm = [2,3,1]
    flip = []
  else
    perm = [3,1,2]
    flip = [1,2]
  end
  return perm,flip
end

function applyPermutions(im, mode)
  perm, flip = convertMode2PermFlip(mode)

  im_ = permutedims(im,perm)
  for i in flip
    im_ = copyproperties(im, reverse(data(im_),dims=i))
  end
  return im_
end

function applyPermutionsRev(im, mode)
  perm, flip = convertMode2PermFlip(mode)
  im_ = im
  for i in flip
    im_ = copyproperties(im,reverse(data(im_),dims=i))
  end
  im_ = permutedims(im_,perm)
  return im_
end

function applyPermutions(im, perm, flip)
  im_ = permutedims(im,perm)
  for i in flip
    im_ = copyproperties(im, reverse(data(im_),dims=i))
  end
  return im_
end

function applyPermutionsRev(im, perm, flip)
  im_ = im
  for i in flip
    im_ = copyproperties(im,reverse(data(im_),dims=i))
  end
  im_ = permutedims(im_,perm)
  return im_
end

function permuteCombinations()
  perms = [
  [1,2,3],
  [1,3,2],
  [2,1,3],
  [2,3,1],
  [3,1,2],
  [3,2,1]]
 return perms
end

function permuteCombinationsName()
  permsName = [
  "[1,2,3]",
  "[1,3,2]",
  "[2,1,3] Isotropic 3D/Coronal MRI",
  "[2,3,1] Sagittal MRI",
  "[3,1,2] Transversal MRI",
  "[3,2,1]"]
 return permsName
end

function flippings()
  f = [
  Int64[],
  [1],
  [2],
  [3],
  [1,2],
  [1,3],
  [2,3],
  [1,2,3]]
 return f
end

function flippingsName()
  f =[
  "[] Transversal MRI",
  "[1]",
  "[2] Coronal MRI",
  "[3]",
  "[1,2] Sagittal MRI",
  "[1,3]",
  "[2,3]",
  "[1,2,3] Isotropic 3D"]
 return f
end
