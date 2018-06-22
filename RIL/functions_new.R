check_onemap = function(file, nround){
  # Reading .raw file
  data = read.table(file, skip = 3) # Only considering Onemap files
  data = as.matrix(data)
  nmarkers = nrow(data)
  
  check_twopts = rf_estimate(data, nmarkers, nround)
  
  return(check_twopts)
  
}

cross_AB = function(data){
  data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
  n = ncol(data)
  t = table(as.data.frame(t(data)))
  if (length(t) == 4){
  n1 = t[1,1]
  n2 = t[1,2]
  n3 = t[2,1]
  n4 = t[2,2]
  r = (n2+n3)/(2*(n1+n4))
  } else if (length(t) == 2){
    n1 = t[1,1]
    n3 = t[2,1]
    r = n3/(2*n1)
  } else {r=0}
  return(r)
}

# Estimating recombination frequency between all markers
rf_estimate = function(data, nmarkers, nround = 6){
  data_p = as.matrix(data)
  
  # Creating list for results
  CC <- matrix(NA,nmarkers,nmarkers)
  rownames(CC) <- colnames(CC) <- substr(as.vector(data_p[,1]),2,4)
  analysis.check <- list()
  analysis.check$'CC' <- CC
  analysis.check$'CR' <- CC
  analysis.check$'RC' <- CC
  analysis.check$'RR' <- CC
  
  for(i in 1:(nrow(data_p)-1)){
    for (j in (i+1):nrow(data_p)){
#      if (data_p[i,2] == "A" & data_p[j,2] == "A"){
        outp = cross_AB(data_p[c(i,j),-c(1,2)])
        analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp#[1,1]
#        analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
#        analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
#        analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]    
#      }
   }
  }
  analysis.check[] <- lapply(analysis.check,round,nround)
  return(analysis.check)
}