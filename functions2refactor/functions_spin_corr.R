library(corrplot)
library(Hmisc)

get_68perm.id <- function(nrot=10000,seed=666){
  lh_centroid_data=read.csv('F:/Google Drive/post-doc/p_factor/summary_stat/allenbrain/spin_test/lh_centroid.csv')
  rh_centroid_data=read.csv('F:/Google Drive/post-doc/p_factor/summary_stat/allenbrain/spin_test/rh_centroid.csv')
  lr_centroid=rbind(lh_centroid_data,rh_centroid_data)
  colnames(lr_centroid)=c('label','centroid1','centroid2','centroid3')
  centroid.l=lh_centroid_data %>%  select(starts_with('centroid'))
  centroid.r=rh_centroid_data %>%  select(starts_with('centroid'))  
  perm.id=rotate.parcellation(as.matrix(centroid.l),as.matrix(centroid.r),nrot=nrot)
  
  return(perm.id)
}


caculate_spin_p <- function(df2cor,perm.id){
  p.spin <-matrix(nrow = ncol(df2cor), ncol = ncol(df2cor))
  for (i in (1:ncol(df2cor))) {
    for (j in (i:ncol(df2cor))){
      if (i==j){p.spin[i,j]=NA}
      else{
        p.spin[i,j] <- perm.sphere.p(df2cor[,i],df2cor[,j],perm.id,corr.type='pearson')}
    }
  }
  #       [,1] [,2] [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
  # [1,]    0    0    0 0.005 0.005 0.120 0.055 0.000 0.010 0.000
  # [2,]   NA    0    0 0.000 0.005 0.485 0.000 0.000 0.010 0.000
  # [3,]   NA   NA    0 0.005 0.000 0.180 0.010 0.000 0.000 0.000
  # [4,]   NA   NA   NA 0.000 0.005 0.000 0.555 0.150 0.455 0.255
  # [5,]   NA   NA   NA    NA 0.000 0.320 0.025 0.000 0.000 0.005
  # [6,]   NA   NA   NA    NA    NA 0.000 0.090 0.015 0.000 0.000
  # [7,]   NA   NA   NA    NA    NA    NA 0.000 0.190 0.005 0.030
  # [8,]   NA   NA   NA    NA    NA    NA    NA 0.000 0.000 0.000
  # [9,]   NA   NA   NA    NA    NA    NA    NA    NA 0.000 0.000
  # [10,]   NA   NA   NA    NA    NA    NA    NA    NA    NA 0.000
  
  
  # copy upper to lower
  p.spin[lower.tri(p.spin)]=t(p.spin)[lower.tri(p.spin)]
  #       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
  # [1,] 0.000 0.000 0.000 0.005 0.005 0.120 0.055 0.000 0.010 0.000
  # [2,] 0.000 0.000 0.000 0.000 0.005 0.485 0.000 0.000 0.010 0.000
  # [3,] 0.000 0.000 0.000 0.005 0.000 0.180 0.010 0.000 0.000 0.000
  # [4,] 0.005 0.000 0.005 0.000 0.005 0.000 0.555 0.150 0.455 0.255
  # [5,] 0.005 0.005 0.000 0.005 0.000 0.320 0.025 0.000 0.000 0.005
  # [6,] 0.120 0.485 0.180 0.000 0.320 0.000 0.090 0.015 0.000 0.000
  # [7,] 0.055 0.000 0.010 0.555 0.025 0.090 0.000 0.190 0.005 0.030
  # [8,] 0.000 0.000 0.000 0.150 0.000 0.015 0.190 0.000 0.000 0.000
  # [9,] 0.010 0.010 0.000 0.455 0.000 0.000 0.005 0.000 0.000 0.000
  # [10,] 0.000 0.000 0.000 0.255 0.005 0.000 0.030 0.000 0.000 0.000
  
  nam=dimnames(df2cor)[[2]]
  dimnames(p.spin) =list(nam,nam)
  return(p.spin)
}


 # plot correlation matrix and show signficant correlation with spin_p<0.05
cor_plot <- function(cor.list,...){
  # type: upper or lower 
  col2use=colorRampPalette(c("steelblue1","white","firebrick1"))
  corrplot(cor.list$r, method="color", col=col2use(200),  
           addCoef.col = "black", 
           tl.col="black", tl.srt=45, 
           p.mat = cor.list$p.spin, sig.level = 0.05, insig = "blank", 
           diag=FALSE ,
           number.cex=0.8,cl.ratio = 0.3, addgrid.col='grey',...
  )
}


# correlation results to save
cor_t2save <- function(cor.list){
  r2save=cor.list$r
  p2save=cor.list$P
  p.spin2save=cor.list$p.spin
  
  r2save[lower.tri(r2save)] <- NA
  p2save[lower.tri(p2save)] <- NA
  p.spin2save[lower.tri(p.spin2save)] <- NA
  t2save=rbind(r2save,p2save,p.spin2save)
  return(t2save)
}


# shuffle loading according to perm.id
build_null_brain_df <- function(perm.id,brain_df){
  brain_data=as.matrix(brain_df)
  perm.n=dim(perm.id)[2]
  all_null_brain=list()
  for (idx in c(1:perm.n)){
    tmp_brain_data=brain_data
    idx_name=sprintf('null_brain_%d',idx)
    perm.idx=perm.id[,idx]
    null_brain=tmp_brain_data[perm.idx]
    all_null_brain[[idx_name]]=null_brain
  }
  null_brain_df=as.data.frame(all_null_brain)
  return(null_brain_df)
}


# p is calculated as the quantile of the true r in the null distribution
# note: a better way is to +1 so that p will not be 0.

caculate_p_from_null <- function(pos_null_dist,neg_null_dist, r_val){
  if (r_val>0){
    p_val=sum(pos_null_dist>r_val)/length(pos_null_dist)
  }else{
    p_val=sum(neg_null_dist<r_val)/length(neg_null_dist)
  }
  return(p_val)
}

# Function to generate a permutation map from a set of cortical regions of interest to itself, 
# while (approximately) preserving contiguity and hemispheric symmetry.
# The function is based on a rotation of the FreeSurfer projection of coordinates
# of a set of regions of interest on the sphere.
#
# Inputs:
# coord.l       coordinates of left hemisphere regions on the sphere        array of size [n(LH regions) x 3]
# coord.r       coordinates of right hemisphere regions on the sphere       array of size [n(RH regions) x 3]
# nrot          number of rotations (default = 10000)                       scalar
#
# Output:
# perm.id      array of permutations, from set of regions to itself        array of size [n(total regions) x nrot]
#
# required library: matrixStats (for rowMins function)
#
# Frantisek Vasa, fv247@cam.ac.uk, June 2017 - July 2018
#       Updated on 16/10/2019 with permutation scheme that uniformly samples the space of permutations on the sphere
#       See github repo (@frantisekvasa) and references within for details

rotate.parcellation = function(coord.l,coord.r,nrot=10000) {
  
  # check that coordinate dimensions are correct
  if (!all(dim(coord.l)[2]==3,dim(coord.r)[2]==3)) {
    if (all(dim(coord.l)[1]==3,dim(coord.r)[1]==3)) {
      print('transposing coordinates to be of dimension nROI x 3')
      coord.l = t(coord.l)
      coord.r = t(coord.r)
    }
  }
  
  nroi.l = dim(coord.l)[1]   # n(regions) in the left hemisphere
  nroi.r = dim(coord.r)[1]   # n(regions) in the right hemisphere
  nroi = nroi.l+nroi.r       # total n(regions)
  
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
    coord.r.rot = coord.r %*% TR; # transformed (rotated) right coordinates
    
    # OLD PERMUTATION SCHEME - COMMENTED OUT 16/10/2019
    # # choose three angles at random - x, y, z
    # # random angles
    # rm(.Random.seed, envir=globalenv()) # reset random seed
    # ax = 2*pi*runif(1)
    # ay = 2*pi*runif(1)
    # az = 2*pi*runif(1)
    # 
    # ### rotation matrices
    # # left hemisphere
    # rx.l = cbind(c(1,0,0),c(0,cos(ax),sin(ax)),c(0,-sin(ax),cos(ax)))
    # ry.l = cbind(c(cos(ay),0,-sin(ay)),c(0,1,0),c(sin(ay),0,cos(ay)))
    # rz.l = cbind(c(cos(az),sin(az),0),c(-sin(az),cos(az),0),c(0,0,1))
    # # right hemisphere - same magnitude of rotation, but signs for y and z axes are flipped to retain symmetry
    # rx.r = cbind(c(1,0,0),c(0,cos(ax),sin(ax)),c(0,-sin(ax),cos(ax)))
    # ry.r = cbind(c(cos(-ay),0,-sin(-ay)),c(0,1,0),c(sin(-ay),0,cos(-ay)))
    # rz.r = cbind(c(cos(-az),sin(-az),0),c(-sin(-az),cos(-az),0),c(0,0,1))
    # 
    # # perform rotation (mutiply coordinates by rotation matrices, for n-by-3 matrix)
    # # left hemisphere
    # coord.l.rot.x = coord.l %*% rx.l
    # coord.l.rot.xy = coord.l.rot.x %*% ry.l
    # coord.l.rot.xyz = coord.l.rot.xy %*% rz.l
    # # right hemisphere
    # coord.r.rot.x = coord.r %*% rx.r
    # coord.r.rot.xy = coord.r.rot.x %*% ry.r
    # coord.r.rot.xyz = coord.r.rot.xy %*% rz.r
    
    # after rotation, find "best" match between rotated and unrotated coordinates
    # first, calculate distance between initial coordinates and rotated ones
    dist.l = array(0,dim=c(nroi.l,nroi.l));
    dist.r = array(0,dim=c(nroi.r,nroi.r));
    # OLD PERMUTATION SCHEME - COMMENTED OUT 16/10/2019
    # for (i in 1:nroi.l) { # left
    #   for (j in 1:nroi.l) {
    #     dist.l[i,j] = sqrt( sum( (coord.l[i,]-coord.l.rot.xyz[j,])^2 ) )
    #   }
    # }
    # for (i in 1:nroi.r) { # right
    #   for (j in 1:nroi.r) {
    #     dist.r[i,j] = sqrt( sum( (coord.r[i,]-coord.r.rot.xyz[j,])^2 ) )
    #   }
    # }
    # UPDATED 5/9/2019 - change of rotated variable name to "coord.l/r.rot" (from coord.l/r.rot.xyz)
    for (i in 1:nroi.l) { # left
      for (j in 1:nroi.l) {
        dist.l[i,j] = sqrt( sum( (coord.l[i,]-coord.l.rot[j,])^2 ) )
      }
    }
    for (i in 1:nroi.r) { # right
      for (j in 1:nroi.r) {
        dist.r[i,j] = sqrt( sum( (coord.r[i,]-coord.r.rot[j,])^2 ) )
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
    
    # RIGHT
    # calculate distances, proceed in order of "most distant minimum"
    # -> for each unrotated region find closest rotated region (minimum), then assign the most distant pair (maximum of the minima), 
    # as this region is the hardest to match and would only become harder as other regions are assigned
    temp.dist.r = dist.r;
    rot.r = c(); ref.r = c();
    for (i in 1:nroi.r) {
      # max(min) (described above)
      ref.ix = which( rowMins(temp.dist.r,na.rm=T) == max(rowMins(temp.dist.r,na.rm=T),na.rm=T) )   # "furthest" row
      rot.ix = which( temp.dist.r[ref.ix,] == min(temp.dist.r[ref.ix,],na.rm=T) )             # closest region
      
      # # alternative option: mean of row - take the closest match for unrotated region that is on average furthest from rotated regions
      # ref.ix = which(nanmean(temp.dist.r,2)==nanmax(nanmean(temp.dist.r,2)))    # "furthest" row
      # rot.ix = which(temp.dist.r(ref.ix,:)==nanmin(temp.dist.r(ref.ix,:)))      # closest region
      ref.r = c(ref.r,ref.ix) # store reference and rotated indices
      rot.r = c(rot.r,rot.ix)
      temp.dist.r[,rot.ix] = NA # set temporary column indices to NaN, to be disregarded in next iteration
      temp.dist.r[ref.ix,] = 0 # because in the above form of the code, R doesn't deal well with whole rows and columns of NaN, set row to low value (which won't matter as furthest rows are assigned first)
      #temp.dist.l[,rot.ix] = NA # set temporary indices to NaN, to be disregarded in next iteration
      #temp.dist.l[ref.ix,] = NA
    }
    
    # mapping is x->y
    # collate vectors from both hemispheres + sort mapping according to "reference" vector
    ref.lr = c(ref.l,nroi.l+ref.r); rot.lr = c(rot.l,nroi.l+rot.r);
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
    if (r%%10==0) print(paste('permutation ',toString(r),' of ',toString(nrot),sep=''))
    
  }
  
  return(perm.id)
  
}




# Function to generate a p-value for the spatial correlation between two parcellated cortical surface maps, 
# using a set of spherical permutations of regions of interest (which can be generated using the function "rotate_parcellation").
# The function performs the permutation in both directions; i.e.: by permute both measures, 
# before correlating each permuted measure to the unpermuted version of the other measure
#
# Inputs:
# x                 one of two maps to be correlated                                                                    vector
# y                 second of two maps to be correlated                                                                 vector
# perm_id           array of permutations, from set of regions to itself (as generated by "rotate_parcellation")        array of size [n(total regions) x nrot]
# corr_type         type of correlation                                                                                 "spearman" (default) or "pearson"
#
# Output:
# p_perm            permutation p-value
#
# Frantisek Vasa, fv247@cam.ac.uk, June 2017 - July 2018

perm.sphere.p = function(x,y,perm.id,corr.type='spearman') {
  
  nroi = dim(perm.id)[1]  # number of regions
  nperm = dim(perm.id)[2] # number of permutations
  
  rho.emp = cor(x,y,method=corr.type)  # empirical correlation
  
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
    rho.null.xy[r] = cor(x.perm[,r],y,method=corr.type)
    rho.null.yx[r] = cor(y.perm[,r],x,method=corr.type)
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
  return((p.perm.xy+p.perm.yx)/2)
  
}