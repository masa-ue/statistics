###rm(list=ls())
library(gdata)
library(np)

aaa <- read.table("data.txt")
zzz <- as.numeric(aaa[,5])-2
print(mean(zzz[zzz>-0.5]))
print(sqrt(var(zzz[zzz>-0.5]))*qnorm(0.95)/sqrt(784))