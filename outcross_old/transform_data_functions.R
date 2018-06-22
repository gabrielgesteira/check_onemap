# Randomly creates (%) of missing data
random_datalost = function(data, missing = 0.1){
data_t = data[,-c(1,2)]
data_t[sample(1:length(data_t), round(missing*length(data_t)), replace = FALSE)] <- "-"
data = cbind(data[,c(1,2)],data_t)
return(data)
}

# Test for presence of all marker types (IMPORTANT FOR OTHER FUNCTIONS)
test_mtypes = function(data) {
data_p = as.data.frame(data)
lvls = matrix(levels(data_p[,2]))
labels_test = matrix(c("A.1", "A.2", "A.3", "A.4", "B1.5", "B2.6", "B3.7", "C.8", "D1.10", "D1.11", "D1.12", "D1.13", "D1.9", "D2.14", "D2.15", "D2.16", "D2.17", "D2.18"),18,1)
res = matrix(NA, 18, 1)
for(i in 1:nrow(lvls)){
  res[i] = isTRUE(lvls[i] == labels_test[i])
}
return(res)
}

# Transform data to Onemap/MapMaker pattern
transform_data = function(data, nmarkers, individuals) {
dat = data[,-c(1:5)]
test = matrix(NA, nmarkers, individuals)
k=1
for (i in 1:ncol(test)){
  test[,i] = paste0(dat[,k],dat[,k+1])
  k = k+2
}
test = tolower(test)
data = as.matrix(data)
test = cbind(data[,1], test)
return(test)
}

# Randomize marker types
randomize_marker_type = function(data, nmarkers, individuals, guide) {
  markers = data[,1]
  types = matrix(NA, nmarkers, 1)
  dat = as.data.frame(t(data[,-1]))
  #guide = sample(1:18, nmarkers, replace = T)
  for(i in 1:ncol(dat)){
  if  (guide[i] == 1) {
  levels(dat[,i]) = c("ac", "ad", "bc", "bd") #A.1
  types[i,] = paste0("A.1")}
    else if (guide[i] == 2) {
  levels(dat[,i]) = c("a", "ac", "ba", "bc") #A.2
  types[i,] = paste0("A.2")}
    else if (guide[i] == 3) {
  levels(dat[,i]) = c("ac", "a", "bc", "b") #A.3
  types[i,] = paste0("A.3")}
    else if (guide[i] == 4) {
  levels(dat[,i]) = c("ab", "a", "b", "o") #A.4
  types[i,] = paste0("A.4")}
    else if (guide[i] == 5) {
  levels(dat[,i]) = c("ab", "a", "a", "b") #B1.5
  types[i,] = paste0("B1.5")}
    else if (guide[i] == 6) {
  levels(dat[,i]) = c("ab", "a", "a", "b") #B2.6
  types[i,] = paste0("B2.6")}
    else if (guide[i] == 7) {
  levels(dat[,i]) = c("a", "ab", "ab", "b") #B3.7
  types[i,] = paste0("B3.7")}
    else if (guide[i] == 8) {
  levels(dat[,i]) = c("a", "a", "a", "o") #C.8
  types[i,] = paste0("C.8")}
    else if (guide[i] == 9) {
  levels(dat[,i]) = c("ac", "ac", "bc", "bc") #D1.9
  types[i,] = paste0("D1.9")}
    else if (guide[i] == 10) {
  levels(dat[,i]) = c("a", "a", "ab", "ab") #D1.10
  types[i,] = paste0("D1.10")}
    else if (guide[i] == 11) {
  levels(dat[,i]) = c("a", "a", "b", "b") #D1.11
  types[i,] = paste0("D1.11")}
    else if (guide[i] == 12) {
  levels(dat[,i]) = c("ab", "ab", "a", "a") #D1.12
  types[i,] = paste0("D1.12")}
    else if (guide[i] == 13) {
  levels(dat[,i]) = c("a", "a", "o", "o") #D1.13
  types[i,] = paste0("D1.13")}
    else if (guide[i] == 14) {
  levels(dat[,i]) = c("ac", "ac", "bc", "bc") #D2.14
  types[i,] = paste0("D2.14")}
    else if (guide[i] == 15) {
  levels(dat[,i]) = c("a", "a", "ab", "ab") #D2.15
  types[i,] = paste0("D2.15")}
    else if (guide[i] == 16) {
  levels(dat[,i]) = c("a", "a", "b", "b") #D2.16
  types[i,] = paste0("D2.16")}
    else if (guide[i] == 17) {
  levels(dat[,i]) = c("ab", "ab", "a", "a") #D2.17
  types[i,] = paste0("D2.17")}
    else if (guide[i] == 18) {
  levels(dat[,i]) = c("a", "a", "o", "o") #D2.18
  types[i,] = paste0("D2.18")}
  }
markers = paste0("*", markers)
data = cbind(markers,types, t(dat))
#colnames(data) = seq(1, (ncol(data)), 1)
fr1 = matrix(c("data", "type", "outcross", rep(NA, (ncol(data)-3))),1, ncol(data))
fr2 = matrix(c(individuals, nmarkers, "0", "0", "0", rep(NA, (ncol(data)-5))), 1, ncol(data))
ind = matrix(NA, 1, ncol(data))
ind[,1:(ncol(ind)-2)] = seq(1, individuals, 1)
data = rbind(fr1, fr2, ind, data)
write.table(data, file = "data_output.raw", sep = " ", quote = F, row.names = F, col.names = F, na = "")
data = data[-c(1,2,3),]
return(data)
}

## Transform to general marker types (A.1, B1.5, B2.6, B3.7, C.8, D1.10, D2.15)
generalize_mtype = function(data_p){
  data_t = as.data.frame(t(data_p[,-c(1,2)]))
  
  for(i in 1:nrow(data_p)){
    if (data_p[i,2] == "A.2"){levels(data_t[,i]) = c("ac", "ad", "bc", "bd") #A.1 labels
    }
    else if (data_p[i,2] == "A.3"){levels(data_t[,i]) = c("ad", "ac", "bd", "bc") #A.1 labels
    }
    else if (data_p[i,2] == "A.4"){levels(data_t[,i]) = c("ad", "ac", "bc", "bd") #A.1 labels
    }
    else if (data_p[i,2] == "D1.9"){levels(data_t[,i]) = c("aa", "ab") #D1.10 labels
    }
    else if (data_p[i,2] == "D1.10"){levels(data_t[,i]) = c("aa", "ab") #D1.10 labels
    }
    else if (data_p[i,2] == "D1.11"){levels(data_t[,i]) = c("aa", "ab") #D1.10 labels
    }
    else if (data_p[i,2] == "D1.12"){levels(data_t[,i]) = c("ab", "aa") #D1.10 labels
    }
    else if (data_p[i,2] == "D1.13"){levels(data_t[,i]) = c("aa", "ab") #D1.10 labels
    }
    else if (data_p[i,2] == "D2.14"){levels(data_t[,i]) = c("aa", "ab") #D2.15 labels
    }
    else if (data_p[i,2] == "D2.15"){levels(data_t[,i]) = c("aa", "ab") #D2.15 labels
    }
    else if (data_p[i,2] == "D2.16"){levels(data_t[,i]) = c("aa", "ab") #D2.15 labels
    }
    else if (data_p[i,2] == "D2.17"){levels(data_t[,i]) = c("ab", "aa") #D2.15 labels
    }
    else if (data_p[i,2] == "D2.18"){levels(data_t[,i]) = c("aa", "ab") #D2.15 labels
    }
  }

# Changing marker labels  
  for(i in 1:nrow(data_p)){
    if (data_p[i,2] == "A.1"){data_p[i,2] = "A"}
    else if (data_p[i,2] == "A.2"){data_p[i,2] = "A"}
    else if (data_p[i,2] == "A.3"){data_p[i,2] = "A"}
    else if (data_p[i,2] == "A.4"){data_p[i,2] = "A"}
    else if (data_p[i,2] == "B1.5"){data_p[i,2] = "B1"}
    else if (data_p[i,2] == "B2.6"){data_p[i,2] = "B2"}
    else if (data_p[i,2] == "B3.7"){data_p[i,2] = "B3"}
    else if (data_p[i,2] == "C.8"){data_p[i,2] = "C8"}
    else if (data_p[i,2] == "D1.9"){data_p[i,2] = "D1"}
    else if (data_p[i,2] == "D1.10"){data_p[i,2] = "D1"}
    else if (data_p[i,2] == "D1.11"){data_p[i,2] = "D1"}
    else if (data_p[i,2] == "D1.12"){data_p[i,2] = "D1"}
    else if (data_p[i,2] == "D1.13"){data_p[i,2] = "D1"}
    else if (data_p[i,2] == "D2.14"){data_p[i,2] = "D2"}
    else if (data_p[i,2] == "D2.15"){data_p[i,2] = "D2"}
    else if (data_p[i,2] == "D2.16"){data_p[i,2] = "D2"}
    else if (data_p[i,2] == "D2.17"){data_p[i,2] = "D2"}
    else if (data_p[i,2] == "D2.18"){data_p[i,2] = "D2"}
  }

  data_p = cbind(data_p[,c(1,2)], t(data_t))
  #data_p = as.data.frame(data_p)
  #levels(data_p[,2]) = c("A", "A", "A", "A", "B1", "B2", "B3", "C8", "D1", "D1", "D1", "D1", "D1", "D2", "D2", "D2", "D2", "D2")
  return(data_p)
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
    outp = cross_11(data_p[c(i,j),-c(1,2)])
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
    else if (data_p[i,2] == "B1" & data_p[j,2] == "A"){
    outp = cross_13a(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "A" & data_p[j,2] == "B2"){
    outp = cross_13b(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B2" & data_p[j,2] == "A"){
    outp = cross_13b(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "A" & data_p[j,2] == "B3"){
    outp = cross_8(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B3" & data_p[j,2] == "A"){
    outp = cross_8(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "A" & data_p[j,2] == "C8"){
    outp = cross_12(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "C8" & data_p[j,2] == "A"){
    outp = cross_12(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "A" & data_p[j,2] == "D1"){
    outp = cross_3a(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D1" & data_p[j,2] == "A"){
    outp = cross_3a(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "A" & data_p[j,2] == "D2"){
    outp = cross_3b(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D2" & data_p[j,2] == "A"){
    outp = cross_3b(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B1" & data_p[j,2] == "B1"){
    outp = cross_16a(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B1" & data_p[j,2] == "B2"){
    outp = cross_17(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B2" & data_p[j,2] == "B1"){
    outp = cross_17(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B1" & data_p[j,2] == "B3"){
    outp = cross_10a(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B3" & data_p[j,2] == "B1"){
    outp = cross_10a(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B1" & data_p[j,2] == "C8"){
    outp = cross_15a(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "C8" & data_p[j,2] == "B1"){
    outp = cross_15a(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B1" & data_p[j,2] == "D1"){
    outp = cross_5(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D1" & data_p[j,2] == "B1"){
    outp = cross_5(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B1" & data_p[j,2] == "D2"){
    outp = cross_6(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D2" & data_p[j,2] == "B1"){
    outp = cross_6(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B2" & data_p[j,2] == "B2"){
    outp = cross_16b(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B2" & data_p[j,2] == "B3"){
    outp = cross_10b(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B3" & data_p[j,2] == "B2"){
    outp = cross_10b(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B2" & data_p[j,2] == "C8"){
    outp = cross_15b(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "C8" & data_p[j,2] == "B2"){
    outp = cross_15b(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B2" & data_p[j,2] == "D1"){
    outp = cross_6(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D1" & data_p[j,2] == "B2"){
    outp = cross_6(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B2" & data_p[j,2] == "D2"){
    outp = cross_5(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D2" & data_p[j,2] == "B2"){
    outp = cross_5(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B3" & data_p[j,2] == "B3"){
    outp = cross_7(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B3" & data_p[j,2] == "C8"){
    outp = cross_9(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "C8" & data_p[j,2] == "B3"){
    outp = cross_9(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B3" & data_p[j,2] == "D1"){
    outp = cross_2(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D1" & data_p[j,2] == "B3"){
    outp = cross_2(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "B3" & data_p[j,2] == "D2"){
    outp = cross_2(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D2" & data_p[j,2] == "B3"){
    outp = cross_2(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "C8" & data_p[j,2] == "C8"){
    outp = cross_14(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "C8" & data_p[j,2] == "D1"){
    outp = cross_4(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D1" & data_p[j,2] == "C8"){
    outp = cross_4(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "C8" & data_p[j,2] == "D2"){
    outp = cross_4(data_p[c(j,i),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D2" & data_p[j,2] == "C8"){
    outp = cross_4(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D1" & data_p[j,2] == "D1"){
    outp = cross_1(data_p[c(i,j),-c(1,2)])
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    }
    else if (data_p[i,2] == "D1" & data_p[j,2] == "D2"){
    outp = matrix(c(0,0,0,0),2,2)
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    warning("Missing cross information", call. = F)
    }
    else if (data_p[i,2] == "D2" & data_p[j,2] == "D1"){
    outp = matrix(c(0,0,0,0),2,2)
    analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
    analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
    analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
    analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
    warning("Missing cross information", call. = F)
    }
    else if (data_p[i,2] == "D2" & data_p[j,2] == "D2"){
    outp = cross_1(data_p[c(i,j),-c(1,2)])
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
