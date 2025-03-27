

# Required packages
install.packages(c("MASS", "geometry", "depth"))
library(MASS)
library(geometry)
library(depth)

# --------------------------------------
# Utility Functions
# --------------------------------------

#' Trimming function based on uniform distribution
#' @param x Numeric vector
#' @return Trimmed values between 0.05 and 0.95
trimf <- function(x) {
  punif(x, min = 0.05, max = 0.95)
}

#' Identity function
#' @param x Numeric vector
#' @return Unchanged input
identity <- function(x) {
  x
}

#' Square function
#' @param x Numeric vector
#' @return Squared values
square <- function(x) {
  x^2
}

#' Square root function
#' @param x Numeric vector
#' @return Square root of values
sroot <- function(x) {
  sqrt(x)
}

#' Sigmoid function
#' @param x Numeric vector
#' @return Sigmoid-transformed values
sigmoid <- function(x) {
  1 / (1 + (x / (1 - x))^(-2))
}

# --------------------------------------
# Statistical Functions
# --------------------------------------

#' Distorted M-Quantile Function
#' @param x Numeric vector
#' @param alfa Quantile level (default: 0.5)
#' @param r Power parameter (default: 2)
#' @param g Distortion function (default: identity)
#' @param tol Convergence tolerance (default: exp(-20))
#' @return Computed M-quantile
dist_mquantile <- function(x, alfa = 0.5, r = 2, g = identity, tol = exp(-20)) {
  y <- sort(x)
  n <- length(y)
  w <- numeric(n)
  gtilde <- function(u) 1 - g(1 - u)
  
  for (i in seq_len(n)) {
    w[i] <- gtilde(i / n) - gtilde((i - 1) / n)
  }
  
  a <- min(x)
  b <- max(x)
  diff <- 1
  
  while (abs(diff) > tol) {
    m <- (a + b) / 2
    y1 <- y[y < m]
    y2 <- y[y > m]
    w1 <- w[y < m]
    w2 <- w[y > m]
    
    z <- (m - y1)^(r - 1)
    h <- (y2 - m)^(r - 1)
    K2 <- (1 - alfa) * sum(z * w1)
    K1 <- alfa * sum(h * w2)
    diff <- K1 - K2
    
    if (diff < 0) {
      b <- m
    } else {
      a <- m
    }
  }
  
  m
}

#' Robust Expectile Function
#' @param x Numeric vector
#' @param probs Probability levels (default: seq(0, 1, 0.25))
#' @param lo Lower trimming function (default: trimf)
#' @param up Upper trimming function (default: trimf)
#' @param dec Decimal places for rounding (default: 16)
#' @return Vector of expectiles
rexpectile <- function(x, probs = seq(0, 1, 0.25), lo = trimf, up = trimf, dec = 16) {
  if (!is.vector(x)) stop("Observations must be in vector form.")
  if (any(probs < 0 | probs > 1)) stop("Probs must be between 0 and 1.")
  
  x <- sort(x)
  n <- length(x)
  steps <- seq(1, 0, length.out = n + 1)
  lw <- lo(steps[-(n + 1)]) - lo(steps[-1])
  uw <- up(steps[-(n + 1)]) - up(steps[-1])
  
  e <- mean(x)
  ee <- numeric(length(probs))
  g <- max(abs(x)) * 1e-6
  
  for (k in seq_along(probs)) {
    p <- probs[k]
    if (p == 0) {
      ee[k] <- min(x, na.rm = TRUE)
    } else if (p == 1) {
      ee[k] <- max(x, na.rm = TRUE)
    } else {
      for (it in 1:20) {
        w <- ifelse(x < e, (1 - p) * lw, p * uw)
        enew <- sum(w * x) / sum(w)
        de <- max(abs(enew - e))
        e <- enew
        if (de < g) break
      }
      ee[k] <- e
    }
  }
  
  names(ee) <- probs
  round(ee, dec)
}

#' Robust Quantile Function
#' @param x Numeric vector
#' @param probs Probability levels (default: seq(0, 1, 0.25))
#' @param trim Trimming proportion (default: 0.1)
#' @return Vector of quantiles
rquant <- function(x, probs = seq(0, 1, 0.25), trim = 0.1) {
  x <- sort(x)
  n <- length(x)
  x <- x[floor(n * trim / 2 + 1):(n - floor(n * trim / 2))]
  n <- length(x)
  x[floor(n * probs) + 1]
}

#' Least Squares M-Quantile Regression
#' @param x Predictor vector
#' @param y Response vector
#' @param r Power parameter (default: 2)
#' @param halfspace Use halfspace method (default: FALSE)
#' @param probs Quantile level (default: NULL)
#' @param gtilde Distortion function (default: identity)
#' @param trim Trimming proportion (default: 0)
#' @return Coefficients (intercept, slope)
lsreg <- function(x, y, r = 2, halfspace = FALSE, probs, gtilde = identity, trim = 0) {
  if (halfspace) {
    rquant_obj <- function(a, x, y, probs, trim) {
      res <- y - a * x - rquant(y - a * x, probs = probs, trim = trim)
      q1 <- quantile(res, trim / 2, na.rm = TRUE)
      q2 <- quantile(res, 1 - trim / 2, na.rm = TRUE)
      res <- res[res >= q1 & res <= q2]
      sum(abs((1 - probs) * pmax(-res, 0) + probs * pmax(res, 0)))
    }
    a1 <- optimize(rquant_obj, x = x, y = y, probs = probs, trim = trim, 
                   lower = -1000, upper = 1000, maximum = FALSE)$minimum
    a0 <- rquant(y - a1 * x, probs = probs, trim = trim)
  } else {
    mq_obj <- function(a, x, y, r, probs, gtilde, trim) {
      res <- y - a * x - dist_mquantile(y - a * x, alfa = probs, g = gtilde)
      q1 <- quantile(res, trim / 2)
      q2 <- quantile(res, 1 - trim / 2)
      res <- res[res >= q1 & res <= q2]
      sum((1 - probs) * pmax(-res, 0)^r + probs * pmax(res, 0)^r)
    }
    a1 <- optimize(mq_obj, x = x, y = y, r = r, probs = probs, gtilde = gtilde, 
                   trim = trim, lower = -1000, upper = 1000, maximum = FALSE)$minimum
    a0 <- dist_mquantile(y - a1 * x, alfa = probs, g = gtilde)
  }
  c(a0, a1)
}

# --------------------------------------
# Extreme Points and Regions
# --------------------------------------

#' Robust Expectile Extreme Points
#' @param data Numeric matrix (2 columns)
#' @param nlines Number of lines (default: 1000)
#' @param alpha Quantile level (default: 0.1)
#' @param y Reference point
#' @return Matrix of extreme points
rextreme.points <- function(data, nlines = 1000, alpha = 0.1, y) {
  p <- matrix(nrow = nlines, ncol = 3)
  angle <- seq(0, 2 * pi, length.out = nlines)
  for (i in seq_len(nlines)) {
    p[i, ] <- c(cos(angle[i]), sin(angle[i]), 
                -rexpectile(data[, 1] * cos(angle[i]) + data[, 2] * sin(angle[i]), probs = 1 - alpha))
  }
  pts <- halfspacen(p, y)
  pts <- pts[chull(pts), ]
  rbind(pts, pts[1, ])
}

#' Expectile Extreme Points (No Trimming)
#' @param data Numeric matrix (2 columns)
#' @param nlines Number of lines (default: 1000)
#' @param alpha Quantile level (default: 0.1)
#' @param y Reference point
#' @return Matrix of extreme points
extreme.points <- function(data, nlines = 1000, alpha = 0.1, y) {
  p <- matrix(nrow = nlines, ncol = 3)
  angle <- seq(0, 2 * pi, length.out = nlines)
  for (i in seq_len(nlines)) {
    p[i, ] <- c(cos(angle[i]), sin(angle[i]), 
                -rexpectile(data[, 1] * cos(angle[i]) + data[, 2] * sin(angle[i]), 
                            probs = 1 - alpha, up = identity, lo = identity))
  }
  pts <- halfspacen(p, y)
  pts <- pts[chull(pts), ]
  rbind(pts, pts[1, ])
}

#' Quantile Extreme Points (No Trimming)
#' @param data Numeric matrix (2 columns)
#' @param nlines Number of lines (default: 1000)
#' @param alpha Quantile level (default: 0.1)
#' @param y Reference point
#' @return Matrix of extreme points
qextreme.points <- function(data, nlines = 1000, alpha = 0.1, y) {
  p <- matrix(nrow = nlines, ncol = 3)
  angle <- seq(0, 2 * pi, length.out = nlines)
  for (i in seq_len(nlines)) {
    p[i, ] <- c(cos(angle[i]), sin(angle[i]), 
                -rquant(data[, 1] * cos(angle[i]) + data[, 2] * sin(angle[i]), 
                        probs = 1 - alpha, trim = 0))
  }
  pts <- halfspacen(p, y)
  pts <- pts[chull(pts), ]
  rbind(pts, pts[1, ])
}

#' Robust Quantile Extreme Points
#' @param data Numeric matrix (2 columns)
#' @param nlines Number of lines (default: 1000)
#' @param alpha Quantile level (default: 0.1)
#' @param y Reference point
#' @return Matrix of extreme points
rqextreme.points <- function(data, nlines = 1000, alpha = 0.1, y) {
  p <- matrix(nrow = nlines, ncol = 3)
  angle <- seq(0, 2 * pi, length.out = nlines)
  for (i in seq_len(nlines)) {
    p[i, ] <- c(cos(angle[i]), sin(angle[i]), 
                -rquant(data[, 1] * cos(angle[i]) + data[, 2] * sin(angle[i]), 
                        probs = 1 - alpha, trim = 0.1))
  }
  pts <- halfspacen(p, y)
  pts <- pts[chull(pts), ]
  rbind(pts, pts[1, ])
}

#' M-Quantile Extreme Points
#' @param data Numeric matrix (2 columns)
#' @param nlines Number of lines (default: 1000)
#' @param r Power parameter (default: 1.5)
#' @param g Distortion function (default: identity)
#' @param alpha Quantile level (default: 0.1)
#' @param y Reference point
#' @return Matrix of extreme points
mqextreme.points <- function(data, nlines = 1000, r = 1.5, g = identity, alpha = 0.1, y) {
  p <- matrix(nrow = nlines, ncol = 3)
  angle <- seq(0, 2 * pi, length.out = nlines)
  for (i in seq_len(nlines)) {
    p[i, ] <- c(cos(angle[i]), sin(angle[i]), 
                -dist_mquantile(data[, 1] * cos(angle[i]) + data[, 2] * sin(angle[i]), 
                                r = r, g = g, alfa = 1 - alpha))
  }
  pts <- halfspacen(p, y)
  pts <- pts[chull(pts), ]
  rbind(pts, pts[1, ])
}

#' Robust M-Quantile Extreme Points
#' @param data Numeric matrix (2 columns)
#' @param nlines Number of lines (default: 1000)
#' @param r Power parameter (default: 1.5)
#' @param alpha Quantile level (default: 0.1)
#' @param y Reference point
#' @return Matrix of extreme points
rmqextreme.points <- function(data, nlines = 1000, r = 1.5, alpha = 0.1, y) {
  p <- matrix(nrow = nlines, ncol = 3)
  angle <- seq(0, 2 * pi, length.out = nlines)
  for (i in seq_len(nlines)) {
    p[i, ] <- c(cos(angle[i]), sin(angle[i]), 
                -dist_mquantile(data[, 1] * cos(angle[i]) + data[, 2] * sin(angle[i]), 
                                r = r, g = trimf, alfa = 1 - alpha))
  }
  pts <- halfspacen(p, y)
  pts <- pts[chull(pts), ]
  rbind(pts, pts[1, ])
}

#' Conditional M-Quantile Regression Regions
#' @param x Predictor vector
#' @param Y Response matrix (2 columns)
#' @param x0 Evaluation point
#' @param r Power parameter (default: 2)
#' @param nsupp Number of support lines (default: 100)
#' @param probs Quantile level (default: 0.05)
#' @param trim Trimming proportion (default: 0)
#' @param show.lines Show support lines (default: FALSE)
conditional.regions <- function(x, Y, x0, r = 2, nsupp = 100, probs = 0.05, trim = 0, show.lines = FALSE) {
  x0 <- rbind(1, x0)
  G.inf <- matrix(ncol = 3, nrow = nsupp)
  
  for (t in seq_len(nsupp)) {
    u <- c(cos(t * pi / nsupp), sin(t * pi / nsupp))
    Yu <- Y %*% u
    a <- lsreg(x = x, y = Yu, r = r, probs = probs, trim = trim)
    G.inf[t, ] <- c(u, as.numeric(a %*% x0))
  }
  
  intercept.inf <- G.inf[, 3] / G.inf[, 2]
  slope.inf <- -G.inf[, 1] / G.inf[, 2]
  
  G.sup <- matrix(ncol = 3, nrow = nsupp)
  for (t in seq_len(nsupp)) {
    u <- c(cos(t * pi / nsupp), sin(t * pi / nsupp))
    Yu <- Y %*% u
    a <- lsreg(x = x, y = Yu, r = r, probs = 1 - probs, trim = trim)
    G.sup[t, ] <- c(u, as.numeric(a %*% x0))
  }
  
  intercept.sup <- G.sup[, 3] / G.sup[, 2]
  slope.sup <- -G.sup[, 1] / G.sup[, 2]
  
  G <- rbind(G.inf, G.sup)
  m <- choose(nsupp - 1, 2)
  E.inf <- matrix(0, m, 4)
  index <- 1
  
  for (i in 1:(nsupp - 2)) {
    for (j in (i + 1):(nsupp - 1)) {
      A <- matrix(c(slope.inf[i], slope.inf[j], -1, -1), ncol = 2)
      B <- matrix(c(intercept.inf[i], intercept.inf[j]), ncol = 1)
      E.inf[index, ] <- c(t(solve(A, -B)), i, j)
      index <- index + 1
    }
  }
  
  E.sup <- matrix(0, m, 5)
  index <- 1
  for (i in 1:(nsupp - 2)) {
    for (j in (i + 1):(nsupp - 1)) {
      A <- matrix(c(slope.sup[i], slope.sup[j], -1, -1), ncol = 2)
      B <- matrix(c(intercept.sup[i], intercept.sup[j]), ncol = 1)
      E.sup[index, ] <- c(t(solve(A, -B)), i, j, 1)
      index <- index + 1
    }
  }
  
  E.inf <- cbind(E.inf, rep(1, nrow(E.inf)))
  tol <- 1e-10
  
  for (i in 1:(nsupp - 2)) {
    for (j in 1:nrow(E.inf)) {
      if (E.inf[j, 2] - slope.inf[i] * E.inf[j, 1] - intercept.inf[i] < -tol) {
        E.inf[j, 5] <- -1
      }
    }
  }
  
  E.sup <- rbind(E.inf[E.inf[, 5] == 1, ], E.sup)
  for (i in 1:(nsupp - 1)) {
    for (j in 1:nrow(E.sup)) {
      if (E.sup[j, 2] - slope.sup[i] * E.sup[j, 1] - intercept.sup[i] > tol) {
        E.sup[j, 5] <- -1
      }
    }
  }
  
  F1 <- E.sup[E.sup[, 5] == 1, ]
  for (i in 1:(nsupp - 1)) {
    for (j in 1:nrow(F1)) {
      if (F1[j, 2] - slope.inf[i] * F1[j, 1] - intercept.inf[i] < -tol) {
        F1[j, 5] <- -1
      }
    }
  }
  
  E <- F1[F1[, 5] == 1, 1:2]
  polyx <- E[, 1]
  polyy <- E[, 2]
  polygon(polyx, polyy, col = NA, border = TRUE)
  hpts <- chull(E)
  hpts <- c(hpts, hpts[1])
  lines(E[hpts, ], lwd = 1)
  
  if (show.lines) {
    for (i in seq_len(nrow(G))) {
      abline(a = G[i, 3] / G[i, 2], b = -G[i, 1] / G[i, 2])
    }
  }
}

# --------------------------------------
# Example Usage and Visualization
# --------------------------------------

# Generate sample data
set.seed(1)
data.ok <- mvrnorm(190, mu = c(10, 20), Sigma = matrix(c(5, 4, 4, 4), ncol = 2, byrow = TRUE))
data.ko <- mvrnorm(10, mu = c(16, 16), Sigma = matrix(c(1, 0, 0, 1), ncol = 2, byrow = TRUE))
data <- rbind(data.ok, data.ko)

n1 <- 100
n2 <- 25
n3 <- 75
y1 <- mvrnorm(n1, mu = c(-1, 1.5), Sigma = matrix(c(1, 0.85, 0.85, 1), ncol = 2, byrow = TRUE))
y2 <- mvrnorm(n2, mu = c(0, 0), Sigma = matrix(c(1, 0.85, 0.85, 1), ncol = 2, byrow = TRUE))
y3 <- mvrnorm(n3, mu = c(1, -1.5), Sigma = matrix(c(1, 0.85, 0.85, 1), ncol = 2, byrow = TRUE))
Y <- rbind(y1, y2, y3)
x <- c(rep(-1, n1), rep(0, n2), rep(1, n3))

# Visualization 1: M-Regions
pdf("mregions.pdf", height = 10, width = 6)
par(mfrow = c(2, 2))

plot(data, pch = 3, cex = 0.3, col = "grey", xlab = "", ylab = "", 
     main = "Half-space regions 10% trimming")
for (alpha in c(0.005, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4)) {
  layer <- rqextreme.points(data, nlines = 1000, alpha = alpha, y = c(10, 20))
  lines(layer, lty = 1)
}

plot(data, pch = 3, cex = 0.3, col = "grey", xlab = "", ylab = "", 
     main = "Half-space regions")
for (alpha in c(0.005, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4)) {
  layer <- qextreme.points(data, nlines = 1000, alpha = alpha, y = c(10, 20))
  lines(layer, lty = 1)
}

plot(data, pch = 3, cex = 0.3, col = "grey", xlab = "", ylab = "", 
     main = "M-quantile regions, r=1.5, 10% trimming")
for (alpha in c(0.005, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4)) {
  layer <- rmqextreme.points(data, nlines = 1000, alpha = alpha, y = c(10, 20))
  lines(layer, lty = 1)
}

plot(data, pch = 3, cex = 0.3, col = "grey", xlab = "", ylab = "", 
     main = "M-quantile regions, r=1.5")
for (alpha in c(0.005, 0.02, 0.05, 0.1, 0.2, 0.3)) {
  layer <- mqextreme.points(data, nlines = 1000, alpha = alpha, y = c(10.3, 19.8))
  lines(layer, lty = 1)
}
layer <- mqextreme.points(data, nlines = 1000, alpha = 0.4, y = c(10.1, 19.8))
lines(layer, lty = 1)

dev.off()

# Visualization 2: Conditional Regions 
pdf("conditional_regions.pdf")  
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(y2, pch = 2, cex = 0.5, col = "darkgrey", xaxt = "n", yaxt = "n", bty = "n", 
     xlab = "", ylab = "", xlim = range(Y[, 1]), ylim = range(Y[, 2]), main = " ")
points(y1, pch = 1, cex = 0.5, col = "darkgrey")
points(y3, pch = 3, cex = 0.5, col = "darkgrey")

r0 <- 2
trim <- 0
conditional.regions(x = x, Y = Y, x0 = -1, r = r0, nsupp = 50, probs = 0.05, 
                    trim = trim, show.lines = F)

for (i in c(0.05, 0.15, 0.35, 0.45)) {
  conditional.regions(x = x, Y = Y, x0 = -1, r = r0, nsupp = 100, probs = i, 
                      trim = trim, show.lines = F)
  conditional.regions(x = x, Y = Y, x0 = 0, r = r0, nsupp = 100, probs = i, 
                      trim = trim, show.lines = F)
  conditional.regions(x = x, Y = Y, x0 = 1, r = r0, nsupp = 100, probs = i, 
                      trim = trim, show.lines = FALSE)
}
dev.off()  