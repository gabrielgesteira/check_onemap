#### CHECK ONEMAP ####
data = read.table("example_out_genotypes.dat", h=T)
source("transform_data_functions.R")
set.seed(12345)

# Transform data
data = transform_data(data, 82, 1000)

# Randomize marker types (Onemap/MapMaker pattern)
marker_type_seq = sample(c(1:18), 82, replace = T)
data_p = randomize_marker_type(data, 82, 1000, marker_type_seq)

# Test for presence of all marker types
test_mtypes(data_p)

# Reading generated .raw file
data_p = read.table("data_out.raw", skip = 3)
data_p = as.matrix(data_p)

## Transform to general marker types (A.1, B1.5, B2.6, B3.7, C.8, D1.10, D2.15)
data_p = generalize_mtype(data_p)

# Randomly creates % lost data
data_p = random_datalost(data_p)

# Counting recombinants
results = rf_estimate(data_p, 82, 4)
results[[1]][2,3]

# Writing in Onemap format
fr1 = matrix(c("data", "type", "outcross", rep(NA, (ncol(data_p)-3))),1, ncol(data_p))
fr2 = matrix(c(1000, 82, "0", "0", "0", rep(NA, (ncol(data_p)-5))), 1, ncol(data_p))
ind = matrix(NA, 1, ncol(data_p))
ind[,1:(ncol(ind)-2)] = seq(1, 1000, 1)
data = rbind(fr1, fr2, ind, data_p)
write.table(data, file = "data_out.raw", sep = " ", quote = F, row.names = F, col.names = F, na = "")

library(onemap)
example_out = read_onemap(inputfile = "data_out.raw")

# Saving RData
save(example_out, results, file="check_onemap.RData")
