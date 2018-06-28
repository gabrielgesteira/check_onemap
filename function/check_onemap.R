################################
###       CHECK ONEMAP       ###
################################
### Required package: onemap ###
################################

# Calling function
source("functions.R")

## RIL example
# Declaring file (containing only markers and in Onemap format), cross type and round number
ril = check_onemap(file = "onemap_example_riself.raw", cross = "ril", nround = 4)

# Viewing results:
View(ril$CC)
View(ril$DIFF)

## BC example
# Declaring file (containing only markers and in Onemap format), cross type and round number
bc = check_onemap(file = "onemap_example_bc.raw", cross = "bc", nround = 4)

# Viewing results:
View(bc$CC)
View(bc$DIFF)

# OUTCROSS example
# Declaring file (containing only markers and in Onemap format), cross type and round number
outc = check_onemap(file = "data_out.raw", cross = "outcross", nround = 4)

# Viewing results:
View(outc$CC)
View(outc$CR)
View(outc$RC)
View(outc$RR)
View(outc$DIFFCC)
View(outc$DIFFCR)
View(outc$DIFFRC)
View(outc$DIFFRR)
