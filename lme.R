library(nlme)

lm_alt <- lme(value ~ 1 + WT, random = ~ 1 | samples, data = data, method="ML") #test
lm_null <- lme(value ~ 1, random = ~ 1 | samples, data = data, method="ML") #control
lr <- 2*(as.numeric(logLik(lm_alt)) - as.numeric(logLik(lm_null))) #lr
pchisq(lr, df=1, lower.tail=F) #chisq
