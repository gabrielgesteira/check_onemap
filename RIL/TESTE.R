#### CHECK ONEMAP ####
source("functions_new.R")

file = "onemap_example_riself.raw"
nround=6

teste = check_onemap(file, nround = 6)
teste$CC


library(onemap)
example_out = read_onemap(inputfile = "data_out.raw")

# Saving RData
save(example_out, results, file="check_onemap.RData")
