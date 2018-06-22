library(onemap)
check_onemap = function(threshold = 0.015, nround = 2){
# Load data
load("check_onemap.RData")

# Recombination frequency
twopts = rf_2pts(example_out, LOD = suggest_lod(example_out), max.rf = 0.5)

# Round values
twopts$analysis[] <- lapply(twopts$analysis, round, nround)
results[] <- lapply(results, round, nround)

# Removing first NA value (A1xA2) -> Origin: Onemap function error
results[[1]][2,1] <- 0.25
results[[2]][2,1] <- 0.25
results[[3]][2,1] <- 0.25
results[[4]][2,1] <- 0.25

# Transform NA values to zero
results[[1]][which(is.na(results[[1]]))] <- 0
results[[2]][which(is.na(results[[2]]))] <- 0
results[[3]][which(is.na(results[[3]]))] <- 0
results[[4]][which(is.na(results[[4]]))] <- 0

# Calculating difference between estimated values
CC11 = results[[1]][] - twopts$analysis$CC
CC22 = results[[2]][] - twopts$analysis$CR
CC33 = results[[3]][] - twopts$analysis$RC
CC44 = results[[4]][] - twopts$analysis$RR

for(i in 1:81){
  for(j in (i+1):82){
    CC11[i,j] = CC11[j,i]
    if (abs(CC11[j,i]) > threshold) {warning('Different RF values detected.')}
    indice1 = which((abs(CC11[]) > threshold), arr.ind = T)
      }
}

for(i in 1:81){
  for(j in (i+1):82){
    CC22[i,j] = CC22[j,i]
    if (abs(CC22[j,i]) > threshold) {warning('Different RF values detected.')}
    indice2 = which((abs(CC22[]) > threshold), arr.ind = T)
  }
}

for(i in 1:81){
  for(j in (i+1):82){
    CC33[i,j] = CC33[j,i]
    if (abs(CC33[j,i]) > threshold) {warning('Different RF values detected.')}
    indice3 = which((abs(CC33[]) > threshold), arr.ind = T)
  }
}

for(i in 1:81){
  for(j in (i+1):82){
    CC44[i,j] = CC44[j,i]
    if (abs(CC44[j,i]) > threshold) {warning('Different RF values detected.')}
    indice4 = which((abs(CC44[]) > threshold), arr.ind = T)
  }
}
index = list(indice1, indice2, indice3, indice4)
index[[2]]
return(warnings(), index)
}

results[[4]][2,6]
twopts$analysis$RR[6,2]

a1 = c(results[[1]][1,6],results[[2]][1,6],results[[3]][1,6],results[[4]][1,6])
a2 = c(twopts$analysis$CC[6,1],twopts$analysis$CR[6,1],twopts$analysis$RC[6,1],twopts$analysis$RR[6,1])
match(a1,a2)

b1 = c(results[[1]][2,6],results[[2]][2,6],results[[3]][2,6],results[[4]][2,6])
b2 = c(twopts$analysis$CC[6,2],twopts$analysis$CR[6,2],twopts$analysis$RC[6,2],twopts$analysis$RR[6,2])
match(b1,b2)

c1 = c(results[[1]][5,6],results[[2]][5,6],results[[3]][5,6],results[[4]][5,6])
c2 = c(twopts$analysis$CC[6,5],twopts$analysis$CR[6,5],twopts$analysis$RC[6,5],twopts$analysis$RR[6,5])
match(b1,b2)

CC44[22,7]

results[[1]][4,77]
twopts$analysis$CC[77,4]
