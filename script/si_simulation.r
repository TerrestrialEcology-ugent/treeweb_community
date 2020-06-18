#############
# R script to #
# run simulation #
# of community tightness #
# an community distance metric #
##################################

# load libraries
library(ade4)
library(ggplot2)
library(gridExtra)
library(Matrix)

# now expand with positive or negative correlation
n_obs <- 53
n_col1 <- 10
n_col2 <- 10
n_col3 <- 10
n_col4 <- 10

# very tight communities
M <- matrix(rbeta(n_col1 * n_col2 * n_col3 * n_col4, 5, 1), 
            ncol = n_col1 + n_col2 + n_col3 + n_col4,
            nrow = n_col1 + n_col2 + n_col3 + n_col4)
M <- M * ifelse(runif(length(M)) < 0.5, -1, 1) # add randomly negative corr
gdata::lowerTriangle(M) <- gdata::upperTriangle(M, byrow=TRUE)
diag(M) <- 1
# turn into PD
M_pd <- nearPD(M, corr = TRUE)$mat
M_pm <- melt(as.matrix(M_pd))
M_pm$type <- "tight"
M_pm$value[M_pm$value == 1] <- NA

# Cholesky decomposition
L = chol(M_pd)
nvars = dim(L)[1]
r = t(L) %*% matrix(rnorm((n_col1 + n_col2 + n_col3 + n_col4)*n_obs), 
                    nrow=n_col1+n_col2+n_col3+n_col4,
                    ncol=n_obs)
r = t(r)

comm1 <- as.matrix(r)[,1:n_col1]
comm2 <- as.matrix(r)[,(n_col1 + 1):(n_col1 + n_col2)]
comm3 <- as.matrix(r)[,(n_col1 + n_col2 + 1):(n_col1 + n_col2 + n_col3)]
comm4 <- as.matrix(r)[,(n_col1 + n_col2 + n_col3 + 1):(n_col1 + n_col2 + n_col3 + n_col4)]



# run mcoa
d1 <- dudi.pca(comm1, scannf = FALSE)
d2 <- dudi.pca(comm2, scannf = FALSE)
d3 <- dudi.pca(comm3, scannf = FALSE)
d4 <- dudi.pca(comm4, scannf = FALSE)

ktt <- ktab.list.dudi(list(d1, d2, d3, d4))

mm <- mcoa(ktt, scannf = FALSE, nf = 2)

# get the dist for each plots
vhb_xy <- mm$Tl1
vhb_xy$id_plot <- rep(1:53, 4)
vhb_xy$ref_x <- mm$SynVar[rep(1:53, 4),1]
vhb_xy$ref_y <- mm$SynVar[rep(1:53, 4),2]
vhb_xy$d <- with(vhb_xy, sqrt((Axis1-ref_x) ** 2 + (Axis2-ref_y) ** 2))
vhb_xy$type <- "tight - w. negative"
vhb_xy %>%
  group_by(type, id_plot) %>%
  summarise(d = sum(d)) -> vhb_dd


# now with loose communities
M <- matrix(rbeta(n_col1 * n_col2 * n_col3 * n_col4, 1, 6), 
            ncol = n_col1 + n_col2 + n_col3 + n_col4,
            nrow = n_col1 + n_col2 + n_col3 + n_col4)
gdata::lowerTriangle(M) <- gdata::upperTriangle(M, byrow=TRUE)
M <- M * ifelse(runif(length(M)) < 0.5, -1, 1) # add randomly negative corr
diag(M) <- 1
# turn into PD
M_pd <- nearPD(M, corr = TRUE)$mat
# to plot the correlation matrix
M_pm2 <- melt(as.matrix(M_pd))
M_pm2$type <- "loose"
M_pm2$value[M_pm2$value == 1] <- NA

# Cholesky decomposition
L = chol(M_pd)
nvars = dim(L)[1]
r = t(L) %*% matrix(rnorm((n_col1 + n_col2 + n_col3 + n_col4)*n_obs), 
                    nrow=n_col1+n_col2+n_col3+n_col4,
                    ncol=n_obs)
r = t(r)

comm1 <- as.matrix(r)[,1:n_col1]
comm2 <- as.matrix(r)[,(n_col1 + 1):(n_col1 + n_col2)]
comm3 <- as.matrix(r)[,(n_col1 + n_col2 + 1):(n_col1 + n_col2 + n_col3)]
comm4 <- as.matrix(r)[,(n_col1 + n_col2 + n_col3 + 1):(n_col1 + n_col2 + n_col3 + n_col4)]



# run mcoa
d1 <- dudi.pca(comm1, scannf = FALSE)
d2 <- dudi.pca(comm2, scannf = FALSE)
d3 <- dudi.pca(comm3, scannf = FALSE)
d4 <- dudi.pca(comm4, scannf = FALSE)

ktt <- ktab.list.dudi(list(d1, d2, d3, d4))

mm <- mcoa(ktt, scannf = FALSE, nf = 2)

# get the dist for each plots
vhb_xy <- mm$Tl1
vhb_xy$id_plot <- rep(1:53, 4)
vhb_xy$ref_x <- mm$SynVar[rep(1:53, 4),1]
vhb_xy$ref_y <- mm$SynVar[rep(1:53, 4),2]
vhb_xy$d <- with(vhb_xy, sqrt((Axis1-ref_x) ** 2 + (Axis2-ref_y) ** 2))
vhb_xy$type <- "loose - w. negative"
vhb_xy %>%
  group_by(type, id_plot) %>%
  summarise(d = sum(d)) -> vhb_dd2

# now tight within but loose across
M <- matrix(rbeta(n_col1 * n_col2 * n_col3 * n_col4, 1, 6), 
            ncol = n_col1 + n_col2 + n_col3 + n_col4,
            nrow = n_col1 + n_col2 + n_col3 + n_col4)
gdata::lowerTriangle(M) <- gdata::upperTriangle(M, byrow=TRUE)
M <- M * ifelse(runif(length(M)) < 0.5, -1, 1) # add randomly negative corr
# go through the different within community corr
gdata::lowerTriangle(M[1:10, 1:10]) <- rbeta(45, 5, 1)
M[1:10, 1:10] <- M[1:10, 1:10] * ifelse(runif(100) < 0.5, -1, 1)
gdata::upperTriangle(M[1:10, 1:10]) <- gdata::lowerTriangle(M[1:10, 1:10], byrow = TRUE)
gdata::lowerTriangle(M[11:20, 11:20]) <- rbeta(45, 5, 1)
M[11:20, 11:20] <- M[11:20, 11:20] * ifelse(runif(100) < 0.5, -1, 1)
gdata::upperTriangle(M[11:20, 11:20]) <- gdata::lowerTriangle(M[11:20, 11:20], byrow = TRUE)
gdata::lowerTriangle(M[21:30, 21:30]) <- rbeta(45, 5, 1)
M[21:30, 21:30] <- M[21:30, 21:30] * ifelse(runif(100) < 0.5, -1, 1)
gdata::upperTriangle(M[21:30, 21:30]) <- gdata::lowerTriangle(M[21:30, 21:30], byrow = TRUE)
gdata::lowerTriangle(M[31:40, 31:40]) <- rbeta(45, 5, 1)
M[31:40, 31:40] <- M[31:40, 31:40] * ifelse(runif(100) < 0.5, -1, 1)
gdata::upperTriangle(M[31:40, 31:40]) <- gdata::lowerTriangle(M[31:40, 31:40], byrow = TRUE)
diag(M) <- 1

M_pd <- nearPD(M, corr = TRUE)$mat
# to plot the correlation matrix
M_pm3 <- melt(as.matrix(M_pd))
M_pm3$type <- "tight within\nloose between"
M_pm3$value[M_pm3$value == 1] <- NA

# Cholesky decomposition
L = chol(M_pd)
nvars = dim(L)[1]
r = t(L) %*% matrix(rnorm((n_col1 + n_col2 + n_col3 + n_col4)*n_obs), 
                    nrow=n_col1+n_col2+n_col3+n_col4,
                    ncol=n_obs)
r = t(r)

comm1 <- as.matrix(r)[,1:n_col1]
comm2 <- as.matrix(r)[,(n_col1 + 1):(n_col1 + n_col2)]
comm3 <- as.matrix(r)[,(n_col1 + n_col2 + 1):(n_col1 + n_col2 + n_col3)]
comm4 <- as.matrix(r)[,(n_col1 + n_col2 + n_col3 + 1):(n_col1 + n_col2 + n_col3 + n_col4)]



# run mcoa
d1 <- dudi.pca(comm1, scannf = FALSE)
d2 <- dudi.pca(comm2, scannf = FALSE)
d3 <- dudi.pca(comm3, scannf = FALSE)
d4 <- dudi.pca(comm4, scannf = FALSE)

ktt <- ktab.list.dudi(list(d1, d2, d3, d4))

mm <- mcoa(ktt, scannf = FALSE, nf = 2)

# get the dist for each plots
vhb_xy <- mm$Tl1
vhb_xy$id_plot <- rep(1:53, 4)
vhb_xy$ref_x <- mm$SynVar[rep(1:53, 4),1]
vhb_xy$ref_y <- mm$SynVar[rep(1:53, 4),2]
vhb_xy$d <- with(vhb_xy, sqrt((Axis1-ref_x) ** 2 + (Axis2-ref_y) ** 2))
vhb_xy$type <- "tight within\nloose between - w. negative"
vhb_xy %>%
  group_by(type, id_plot) %>%
  summarise(d = sum(d)) -> vhb_dd3

# a plot comparing the correlation
M_a <- rbind(M_pm, M_pm2, M_pm3)
gg_1 <- ggplot(M_a, aes(x=Var1, y = Var2, fill = value)) +
  geom_raster() +
  facet_wrap(~ type) +
  scale_fill_gradient2(low = "red", high = "blue") +
  labs(x = "", y = "",
       title = "Simulated correlation structure between 4 communities") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())

# a plot comparing the d-metric
vhb_a <- rbind(vhb_dd, vhb_dd2, vhb_dd3)
gg_2 <- ggplot(vhb_a, aes(x=type, y=d))+
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.5) +
  stat_summary(fun.data = "mean_cl_boot", color = "red") +
  labs(x = "Type of correlation structure",
       y = "Community distance index",
       title = "Effect of correlation structure between communities\non community distance index")
# nice!

gg_a <- grid.arrange(gg_1, gg_2, nrow = 2)
ggsave("~/PostDoc_Ghent/Synthesis_3/figures/simulation_res_wneg.png", gg_a,
       width = 8, height = 9, units = "in")