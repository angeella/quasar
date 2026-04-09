
##############################Real data analysis#####################

#Load packages
library(MASS)
library(quasar)

#Load db
data("birthwt")

#Specify quantiles levels
tau <- c(1:9)/10

#Some basic preprocessing
birthwt$race <- as.factor(birthwt$race)

#fit the model
mod <- rq(bwt ~ age + lwt + race +smoke + ht + ui + ptl + ftv,
          data=birthwt, tau=tau)

#estimate rank score tests
res <- closedTesting(mod, X = "smoke")

#Adjusted pvalues (bonferroni)
p.adjust(res$p.value, "bonferroni")

#Adjusted pvalues (holm-bonferroni)
p.adjust(res$p.value, "holm")

#Adjusted pvalues (closed-testing)
res$p.value.adjusted
