##########################################################################
##########################################################################

library(nortest)
#library(lmtest)

x = read.table("~/Dropbox/Tese/Resultados/série3/cor/cegs_1_2.txt")
Y = read.table("~/Dropbox/Tese/Resultados/série3/cor/contiguity.txt")
y = log10(Y)
#y = Y

#normality test (Lilliefors (Kolmogorov-Smirnov)):
normX = lillie.test(x$V1)
normY = lillie.test(y$V1)
#X:
normX$p.value
#Y:
normY$p.value

Cor = cor(x$V1, y$V1)
linRegr = lm(y$V1~x$V1)
sum = summary(linRegr)

#Pearson correlation coefficient:
Cor

#Coefficient of determination:
sum$r.squared

#Intercept:
int = sum$coefficients[1,"Estimate"]
int

#Slope:
slo = sum$coefficients[2,"Estimate"]
slo

#P-value of T-test for slope = 0:
sum$coefficients[2,"Pr(>|t|)"]

#scatterplot and regression line
jpeg('plot1.jpg', width = 800, height = 480, quality = 500)

plot(y$V1 ~ x$V1, pch = 16, xlab = "Number of CEGs (group A)", ylab = "Log of contiguity", main = "")
title(main = list("Contiguity vs CEGs (group A)", cex = 2, col = "black", font = 6))
abline(int, slo, lty = 2, lwd = 2, col = "black")
#mtext(bquote(y == .(int) + .(slo)*x), adj=1, padj=0)

dev.off()

#plot residuals
jpeg('plot2.jpg', width = 800, height = 480, quality = 500)

plot(fitted(linRegr), residuals(linRegr), pch = 16, xlab = "Predicted contiguity (log)", ylab = "Residuals", main = "")
abline(h = 0, lty = 2, lwd = 2, col = "black")
title(main = list("Residuals vs predicted contiguity", cex = 2, col = "black", font = 6))

dev.off()

#plot histogram
jpeg('plot3.jpg', width = 800, height = 480, quality = 500)

hist(residuals(linRegr), col="darkgray", main = "", xlab =  "Residuals", ylab = "Frequency")
title(main = list("Residuals histogram", cex = 2, col = "black", font = 6))

dev.off()

#normality test (Lilliefors (Kolmogorov-Smirnov)):
normRes = lillie.test(residuals(linRegr))
normRes$p.value
#####################################################
Cor
sum$r.squared
sum$coefficients[2,"Pr(>|t|)"]
residuals = 10^fitted(linRegr) - 10^y$V1
RMSE = sqrt(mean((residuals)^2))
#RMSE = sqrt(mean((residuals(linRegr))^2))
RMSE