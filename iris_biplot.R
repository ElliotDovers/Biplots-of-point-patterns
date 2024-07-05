# example biplot from iris dataset
library(corrplot)
library(plotrix)

plot.res <- 500

pca <- prcomp(iris[ , -5], scale. = T) # prcomp auto centers/scales the data
Sigma <- cor(scale(iris[,-5]))

# could equivalently use SVD
tmp <- svd(scale(iris[ , -5]))
G <- tmp$u[,1:2]
H <- tmp$v[,1:2] %*% t(diag(tmp$d[1:2]))

G <- pca$x[,1:2]
H <- pca$rotation[,1:2] %*% diag(pca$sdev[1:2])

# get the estimated correlation
Sigma_hat <- H %*% t(H)
# set up nicer names
nice.names <- c("Sepal L.", "Sepal W.", "Petal L.", "Petal W.")
dimnames(Sigma) <- list(nice.names, nice.names)
dimnames(Sigma_hat) <- list(nice.names, nice.names)

# rotate the data?
G <- pca$x[,1:2]
H <- pca$rotation[,1:2] %*% diag(pca$sdev[1:2])
theta <- -pi/6 + pi
R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
G <- G %*% R
H <- H %*% R

png(filename = paste0(getwd(), "/eg_biplot.png"), width = 6.2 * plot.res, height = 5 * plot.res, res = plot.res)
# layout(matrix(c(1,1,1,2,3,4), nrow = 2, byrow = T), widths = rep(1/3, 3), heights = c(0.5,0.5))
par(mfrow = c(2,2))
par(mar = c(2.1, 2.1, 3.1, 1))
alpha <- max(apply(G, 1, norm, type = "2")) / max(apply(H, 1, norm, type = "2"))
loads <- H * alpha
plot(G, col = iris$Species, xlim = range(c(G[,1], loads[,1])), asp = 1, pch = as.numeric(iris$Species), cex = 0.75,
     ylim = range(c(G[,2], loads[,2])), xlab = "Principle Component 1", ylab = "Principle Component 2")
mtext(expression(bold("A: Iris Dataset Biplot")), side = 3, line = 1.4, adj = 0, at = -6.5)
mtext(expression(bold(paste(X%~~%U, "\u039b"^T, " (points U, lines \u039b)"))), side = 3, line = 0.15, adj = 0, cex = 0.8, at = -6.5)
arrows(x0 = rep(0, nrow(loads)), y0 = rep(0, nrow(loads)), x1 = loads[,1], y1 = loads[,2], length = 0.0, angle = 30, lwd = 0.5)
text(x = loads[,1], y = loads[,2] + c(0.15,0.15,-0.15,0.15), labels = nice.names, cex = 0.75)
draw.arc(x = 0, y = 0, radius = 2.5, angle1 = sin(loads[2, 2] / loads[2, 1]), angle2 = (pi/2) - sin(loads[1, 2] / loads[1, 1]) + (sin(loads[1, 2] / loads[1, 1])/3.2), col = "royalblue", lwd = 1.5)
text(x = loads[1,1] + 4.5, y = loads[1,2] + 0.2, labels = expression("Cov("~X[1]~","~X[2]~") "%~~%-0.09), cex = 1, col = "royalblue")

par(mar = c(2.1, 2.1, 2.1, 1))
corrplot(Sigma_hat, type = "lower", method = "square", tl.srt = 45, diag = F, tl.col = "black", tl.cex = 0.75, addCoef.col = c("royalblue", rep("black", 5)), cl.pos = "n", col = "white", mar = c(0, 0, 2.2, 0))
mtext(expression(bold(paste("B: Estimated Covariance, \u039b\u039b"^T))), side = 3, line = 0.4, adj = 0)

par(mar = c(3.5, 3.5, 2.5, 0.5))
with(iris, plot(Sepal.Width, Sepal.Length, col = Species, cex = 0.65, pch = as.numeric(Species), xlab = "", ylab = ""))
mtext(expression(bold(paste("C: Raw Iris Data (", X[1], " and ", X[2], ")"))), side = 3, line = 0.65, adj = 0, at = 1.325)
mtext(expression(paste(X[2], ": Sepal W.")), side = 1, line = 2.5, cex = 0.6)
mtext(expression(paste(X[1], ": Sepal L.")), side = 2, line = 2.5, cex = 0.6)
with(iris, plot(Petal.Width, Petal.Length, col = Species, cex = 0.65, pch = as.numeric(Species), xlab = "", ylab = ""))
mtext(expression(bold(paste("D: Raw Iris Data (", X[3], " and ", X[4], ")"))), side = 3, line = 0.65, adj = 0, at = -0.625)
legend("bottomright", legend = c("", "", levels(iris$Species)), bty = "n", pch = c(NA, NA, 1:3), col = c(NA, NA, 1:3))
mtext(expression(paste(X[4], ": Petal W.")), side = 1, line = 2.5, cex = 0.6)
mtext(expression(paste(X[3], ": Petal L.")), side = 2, line = 2.5, cex = 0.6)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# for talk:

png(filename = paste0(getwd(), "/talk_plot1A.png"), width = 6.2 * plot.res, height = 5 * plot.res, res = plot.res)
par(mar = c(3, 3, 0, 0))
alpha <- max(apply(G, 1, norm, type = "2")) / max(apply(H, 1, norm, type = "2"))
loads <- H * alpha
plot(G, col = iris$Species, xlim = range(c(G[,1], loads[,1])), asp = 1, pch = as.numeric(iris$Species), cex = 0.75,
     ylim = range(c(G[,2], loads[,2])), xlab = "", ylab = "")
# mtext(expression(bold("A: Iris Dataset Biplot")), side = 3, line = 1.4, adj = 0, at = -6.5)
# mtext(expression(bold(paste("Biplot of ", X%~~%"(", g[1], ",...,", g[n], ")(", h[1], ",...,", h[m], ")"^"T", " where we plot points ", g[i], " and lines ", h[j]))), side = 3, line = 0.15, adj = 0, cex = 0.8, at = -4.5)
mtext("Principal Component 1", side = 1, line = 2, cex = 0.8)
mtext("Principal Component 2", side = 2, line = 2.2, cex = 0.8)
arrows(x0 = rep(0, nrow(loads)), y0 = rep(0, nrow(loads)), x1 = loads[,1], y1 = loads[,2], length = 0.0, angle = 30, lwd = 0.5)
text(x = loads[,1], y = loads[,2] + c(0.15,0.15,-0.15,0.15), labels = nice.names, cex = 0.75)
# draw.arc(x = 0, y = 0, radius = 2.5, angle1 = sin(loads[2, 2] / loads[2, 1]), angle2 = (pi/2) - sin(loads[1, 2] / loads[1, 1]) + (sin(loads[1, 2] / loads[1, 1])/3.2), col = "royalblue", lwd = 1.5)
# text(x = loads[1,1] + 4.5, y = loads[1,2] + 0.2, labels = expression("Cov("~X[1]~","~X[2]~") "%~~%-0.09), cex = 1, col = "royalblue")
legend("bottomleft", legend = c("", "", levels(iris$Species)), bty = "n", pch = c(NA, NA, 1:3), col = c(NA, NA, 1:3))
dev.off()

X <- iris[,-5]
dimnames(X) <- list(paste("obs.", 1:nrow(X)), c("Sepal L.", "Sepal W.", "Petal L.", "Petal W."))
U <- as.data.frame(G)
dimnames(U) <- list(paste("obs.", 1:nrow(U)), c("PC1", "PC2"))
Lambda <- as.data.frame(H)
dimnames(Lambda) <- list(c("Sepal L.", "Sepal W.", "Petal L.", "Petal W."), c("PC1", "PC2"))

png(filename = paste0(getwd(), "/talk_plot1B.png"), width = 4 * plot.res, height = 4 * plot.res, res = plot.res)
corrplot(Sigma_hat, type = "lower", method = "square", tl.srt = 45, diag = F, tl.col = "black", tl.cex = 0.75, addCoef.col = rep("black", 6), cl.pos = "n", number.cex = 2.5, col = "white", mar = c(0, 0, 0, 0))
mtext(expression(paste(Sigma%~~%Lambda, Lambda^T)), side = 3, line = 0, adj = 0, at = 2, cex = 3)
dev.off()

png(filename = paste0(getwd(), "/talk_plot1C.png"), width = 4 * plot.res, height = 3 * plot.res, res = plot.res)
par(mar = c(3.2, 3.2, 2, 0))
with(iris, plot(Sepal.Width, Sepal.Length, col = Species, cex = 0.65, pch = as.numeric(Species), xlab = "", ylab = ""))
mtext(expression(bold(paste("Raw Iris Data (", X[1], " and ", X[2], ")"))), side = 3, line = 0.65, adj = 0, at = 1.5)
mtext(expression(paste(X[2], ": Sepal W.")), side = 1, line = 2, cex = 0.8)
mtext(expression(paste(X[1], ": Sepal L.")), side = 2, line = 2, cex = 0.8)
dev.off()

png(filename = paste0(getwd(), "/talk_plot1D.png"), width = 4 * plot.res, height = 3 * plot.res, res = plot.res)
par(mar = c(3.2, 3.2, 2, 0))
with(iris, plot(Petal.Width, Petal.Length, col = Species, cex = 0.65, pch = as.numeric(Species), xlab = "", ylab = ""))
mtext(expression(bold(paste("Raw Iris Data (", X[3], " and ", X[4], ")"))), side = 3, line = 0.65, adj = 0, at = 0)
# legend("bottomright", legend = c("", "", levels(iris$Species)), bty = "n", pch = c(NA, NA, 1:3), col = c(NA, NA, 1:3))
mtext(expression(paste(X[4], ": Petal W.")), side = 1, line = 2, cex = 0.8)
mtext(expression(paste(X[3], ": Petal L.")), side = 2, line = 2, cex = 0.8)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()