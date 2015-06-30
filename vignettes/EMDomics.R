## ------------------------------------------------------------------------
exp_data <- rnorm(100)
names(exp_data) <- paste("sample", 1:100)

groupA.labels <- rep("A",50)
groupB.labels <- rep("B",20)
groupC.labels <- rep("C",30)

labels <- c(groupA.labels, groupB.labels, groupC.labels)
names(labels) <- names(exp_data)

## ------------------------------------------------------------------------
library(ggplot2)

df <- as.data.frame(exp_data)
df$group[1:50] <- "A"
df$group[51:70] <- "B"
df$group[71:100] <- "C"

ggplot(df, aes(exp_data, fill=group)) + geom_density(alpha=0.5)

## ------------------------------------------------------------------------
library(EMDomics)
calculate_emd_gene(exp_data, labels, names(exp_data))

## ------------------------------------------------------------------------
exp_data2 <- exp_data
mod_vec <- sample(c(2,-2), 50, replace=TRUE)
exp_data2[1:50] <- exp_data2[1:50] + mod_vec

