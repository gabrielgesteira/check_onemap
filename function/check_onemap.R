################################
###       CHECK ONEMAP        ##
################################
### Required package: onemap ###
################################

# Calling function
source("functions.R")

# Declaring file (containing only markers and in Onemap format), cross type and round number
ril = check_onemap(file = "onemap_example_riself.raw", cross = "ril", nround = 4)
View(ril$DIFF)

bc = check_onemap(file = "onemap_example_bc.raw", cross = "bc", nround = 4)
View(bc$DIFF)

outc = check_onemap(file = "data_out.raw", cross = "outcross", nround = 4)
View(outc$DIFF)
View(outc$CC)
