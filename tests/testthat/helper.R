#####
# Prepare data sets for test
##

## eco data
data("eco")

## rats data
rats <- read.table("rats.raw", header = TRUE)
rats <- rats[!is.na(rats$response), ]
rats$subject <- factor(rats$subject)
rats$group <- factor(rats$group)
rats$treatment <- factor(ifelse(rats$control, 1, ifelse(rats$low, 2, 3)),
                         labels = c("control", "low", "high"))
set.seed(42)
rats$noise1 <- rnorm(nrow(rats))
rats$noise2 <- sample(c("grey", "white"), size = nrow(rats), replace = TRUE)

