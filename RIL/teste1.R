# Reading .raw file
data_p = read.table("onemap_example_riself.raw", skip = 3)
data = as.matrix(data_p)
nmarkers = nrow(data)

cross_1 = function(data){
  data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
  n = ncol(data)
  t = table(as.data.frame(t(data)))
  n1 = t[1,1]
  n2 = t[1,2]
  n3 = t[2,1]
  n4 = t[2,2]
  r1 = (n2 + n3)/n
  r2 = (n1 + n4)/n
  r = matrix(c(r1,r2,r1,r2),2,2)
  return(r)
}

# Estimating recombination frequency between all markers
rf_estimate = function(data_p, nmarkers, nround = 6){
  source("cross_types.R")
  data_p = as.matrix(data_p)
  
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
      if (data_p[i,2] == "A" & data_p[j,2] == "A"){
        outp = cross_1(data_p[c(i,j),-c(1,2)])
        analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
        analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
        analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
        analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]    
      }
      else if (data_p[i,2] == "A" & data_p[j,2] == "B1"){
        outp = cross_13a(data_p[c(i,j),-c(1,2)])
        analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
        analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
        analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
        analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]    
      }
    }
  }
  analysis.check[] <- lapply(analysis.check,round,nround)
  return(analysis.check)
}