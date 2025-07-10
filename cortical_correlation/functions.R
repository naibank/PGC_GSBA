
#---------------------------------------------------
# Function to rotate brain regions on sphere: Spin Perm: Get Perm-IDs for Select LH ROIs 

rotate.parcellation_LHonly = function(coord.l,nrot=10000) {
  
  # check that coordinate dimensions are correct
  if (dim(coord.l)[1]==3) {
    print('transposing coordinates to be of dimension nROI x 3')
    coord.l = t(coord.l)
  }
  
  
  nroi.l = dim(coord.l)[1]   # n(regions) in the left hemisphere
  nroi = nroi.l
  
  perm.id = array(0,dim=c(nroi,nrot)); # initialise output array
  r = 0; c = 0; # count successful (r) and unsuccessful (c) iterations
  
  # UPDATED 16/10/2019 - set up updated permutation scheme 
  I1 = diag(3); I1[1,1] = -1;
  # main loop -  use of "while" is to ensure any rotation that maps to itself is excluded (this is rare, but can happen)
  while (r < nrot) {
    
    # UPDATED 16/10/2019
    A = matrix(rnorm(9, mean = 0, sd = 1), nrow = 3, ncol = 3)
    qrdec = qr(A)       # QR decomposition
    TL = qr.Q(qrdec)    # Q matrix
    temp = qr.R(qrdec)  # R matrix
    TL = TL%*%diag(sign(diag(temp)))
    if (det(TL)<0) {
      TL[,1] = -TL[,1]
    }
    # reflect across the Y-Z plane for right hemisphere
    TR = I1 %*% TL %*% I1;
    coord.l.rot = coord.l %*% TL; # transformed (rotated) left coordinates
    
    # after rotation, find "best" match between rotated and unrotated coordinates
    # first, calculate distance between initial coordinates and rotated ones
    dist.l = array(0,dim=c(nroi.l,nroi.l));
    
    # UPDATED 5/9/2019 - change of rotated variable name to "coord.l/r.rot" (from coord.l/r.rot.xyz)
    for (i in 1:nroi.l) { # left
      for (j in 1:nroi.l) {
        dist.l[i,j] = sqrt( sum( (coord.l[i,]-coord.l.rot[j,])^2 ) )
      }
    }
    
    # LEFT
    # calculate distances, proceed in order of "most distant minimum"
    # -> for each unrotated region find closest rotated region (minimum), then assign the most distant pair (maximum of the minima), 
    # as this region is the hardest to match and would only become harder as other regions are assigned
    temp.dist.l = dist.l
    rot.l = c(); ref.l = c();
    #tba.r = tba.c = 1:nroi.l # rows and columns that are yet "to be assigned"
    for (i in 1:nroi.l) {
      # max(min) (described above)
      ref.ix = which( rowMins(temp.dist.l,na.rm=T) == max(rowMins(temp.dist.l,na.rm=T),na.rm=T) )   # "furthest" row
      rot.ix = which( temp.dist.l[ref.ix,] == min(temp.dist.l[ref.ix,],na.rm=T) ) # closest region
      
      # # alternative option: mean of row - take the closest match for unrotated region that is on average furthest from rotated regions
      # ref.ix = which(nanmean(temp.dist.l,2)==nanmax(nanmean(temp.dist.l,2)))    # "furthest" row
      # rot.ix = which(temp.dist.l(ref.ix,:)==nanmin(temp.dist.l(ref.ix,:)))      # closest region    
      ref.l = c(ref.l,ref.ix) # store reference and rotated indices
      rot.l = c(rot.l,rot.ix)
      temp.dist.l[,rot.ix] = NA # set temporary column indices to NaN, to be disregarded in next iteration
      temp.dist.l[ref.ix,] = 0 # because in the above form of the code, R doesn't deal well with whole rows and columns of NaN, set row to low value (which won't matter as furthest rows are assigned first)
      #temp.dist.l[,rot.ix] = NA # set temporary indices to NaN, to be disregarded in next iteration
      #temp.dist.l[ref.ix,] = NA
    }
    
    # mapping is x->y
    # collate vectors from both hemispheres + sort mapping according to "reference" vector
    ref.lr = c(ref.l); rot.lr = c(rot.l);
    b = sort(ref.lr,index.return=T); 
    ref.lr.sort = ref.lr[b$ix]; rot.lr.sort = rot.lr[b$ix];
    
    # verify that permutation worked (output should be vector with values 1:nroi = 1:(nroi_l+nroi_r))
    if (!all(sort(rot.lr.sort,decreasing=F)==c(1:nroi))) {
      #save.image('~/Desktop/perm_error.RData')
      browser("permutation error")
    }
    
    # verify that permutation does not map to itself
    if (!all(rot.lr.sort==c(1:nroi))) {
      r = r+1
      perm.id[,r] = rot.lr.sort # if it doesn't, store it
    } else {
      c = c+1
      print(paste('map to itself n. ',toString(c),sep=''))
    }
    
    # track progress
    #if (r%%10==0) print(paste('permutation ',toString(r),' of ',toString(nrot),sep=''))
    
  }
  
  return(perm.id)
  
}

#---------------------------------------------------
# Get Perm-IDs for a given set of ROI-indices for Glasser
fGet_perm_ids_input_ROI_indices_GlasserLH = function(input_LH_ROI_indices,nIterNull = 1000){
  
  # input: input_LH_ROI_indices -> indices of select ROIs for LH-Glasser (range within 1 to 180)
  # Read Brain region sphere coordinates for Glasser HCP brain parcellation
  sphere_HCP <- read.table("sphere_HCP.txt", quote="\"", comment.char="")
  
  #coord.l = as.matrix(sphere_HCP[c(1:180),])
  coord.l = as.matrix(sphere_HCP[input_LH_ROI_indices,])
  
  # call function to get permuted IDs
  perm.id.LH.select =rotate.parcellation_LHonly(coord.l,nrot=nIterNull) 
  
  return(perm.id.LH.select)
}


#------------------------------------------------------------------
# Function for running spin Perm test 
# REF: # Frantisek Vasa, fv247@cam.ac.uk, June 2017 - July 2018
#perm.sphere.p = function(x,y,perm.id,corr.type='spearman') {

perm.sphere.p_return_corr_array = function(x,y,perm.id,corr.type='pearson'){
  
  nroi = dim(perm.id)[1]  # number of regions
  nperm = dim(perm.id)[2] # number of permutations
  
  rho.emp = cor(x,y,method=corr.type, use="complete.obs")  # empirical correlation
  
  # permutation of measures
  x.perm = y.perm = array(NA,dim=c(nroi,nperm))
  for (r in 1:nperm) {
    for (i in 1:nroi) {
      x.perm[i,r] = x[perm.id[i,r]]
      y.perm[i,r] = y[perm.id[i,r]]
    }
  }
  
  # correlation to unpermuted measures
  rho.null.xy = rho.null.yx = vector(length=nperm)
  for (r in 1:nperm) {
    rho.null.xy[r] = cor(x.perm[,r],y,method=corr.type, use="complete.obs")
    rho.null.yx[r] = cor(y.perm[,r],x,method=corr.type, use="complete.obs")
  }
  
  # p-value definition depends on the sign of the empirical correlation
  if (rho.emp>0) {
    p.perm.xy = sum(rho.null.xy>rho.emp)/nperm
    p.perm.yx = sum(rho.null.yx>rho.emp)/nperm
  } else { 
    p.perm.xy = sum(rho.null.xy<rho.emp)/nperm
    p.perm.yx = sum(rho.null.yx<rho.emp)/nperm
  } 
  
  # return average p-value
  pval_avg <- (p.perm.xy+p.perm.yx)/2
  
  # check if p-value is 0
  if(pval_avg == 0){
    pval_avg = 1/nperm
  }
  
  res_list <- list(rho.emp,pval_avg,rho.null.xy,rho.null.yx)
  
  return(res_list)
  
}

#------------------------------------------------------------------
# Function to Compute correlation with two given maps: subset of ROIs
fGet_Corr_pval_SpinTest_two_maps_selectROIs_CorrArray <- function(array_es_1,array_es_2,select_LH_ROI_indices,nIterNull,corr.type='pearson'){
  
  # Get permuted regional IDs based on spherical rotations of Glasser brain atlas
  #nIterNull = 1000  # 10000
  df_perm_ids_LH = fGet_perm_ids_input_ROI_indices_GlasserLH(select_LH_ROI_indices,nIterNull)
  
  temp_cor = cor(array_es_1,array_es_2,method=corr.type, use="complete.obs")  # empirical correlation
  
  # spin test
  res_list <- perm.sphere.p_return_corr_array(array_es_1,array_es_2,df_perm_ids_LH,corr.type=corr.type) 
  
  temp_pvalue = res_list[[2]]
  rho.null.xy = res_list[[3]]
  rho.null.yx = res_list[[4]]
  
  res_cor_pval <- list(temp_cor,temp_pvalue,rho.null.xy,rho.null.yx)
  
  return(res_cor_pval)
}
