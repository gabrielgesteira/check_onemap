library(onemap)
check_onemap = function(){
load("check_onemap.RData")
twopts = rf_2pts(example_out, LOD = suggest_lod(example_out), max.rf = 0.5)
twopts$analysis[] <- lapply(twopts$analysis,round,2)
results[[1]][which(is.na(results[[1]]))] <- 0
results[] <- lapply(results,round,2)
CC11 = results[[1]][] - twopts$analysis$CC
summary(CC11)
as = matrix(NA, 82, 82)
for(i in 1:81){
  for(j in (i+1):82){
    CC11[i,j] = CC11[j,i]
    if (abs(CC11[j,i]) > 0.02) {print('Empty')}
    else print('Not Empty')
    indicex = which((abs(CC11[]) > 0.02), arr.ind = T)
    }
}


# Problem: C8 X C8

apply(CC11, 2, function(c)sum(c!=0))

results$CC[2,1] = 2
twopts$analysis$CC[42,15]
results[[1]][42,15]
results[[1]][which(is.na(results[[1]]))] <- 0
lapply(results[1],function(x) ifelse(x == NA,0,x), how = "replace")



# TUTORIAL
example_out = read_onemap(inputfile = "data_out.raw")
example_out
plot(example_out, all = F)
test_segregation_of_a_marker(example_out, 1)
segreg_test = test_segregation(example_out)
print(segreg_test)
plot(segreg_test)
select_segreg(segreg_test, distorted = F)
select_segreg(segreg_test, distorted = T)
LOD_sug = suggest_lod(example_out)
twopts = rf_2pts(example_out, LOD = LOD_sug, max.rf = 0.5)
print(twopts, c("A1", "A4"))
sequence = make_seq(twopts, "all")
gr = group(sequence)
g1 = make_seq(gr, 1)
g2 = make_seq(gr, 2)
g1o = order_seq(g1)
g2o = order_seq(g2)
X11()
rf_graph_table(make_seq(g1o),inter=FALSE)
