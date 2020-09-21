#######################
# R script for the Treeweb       #
# cross-community analysis       #
# analysis reported in main text #
#######################################

# load libraries
library(spdep)
library(ade4)
library(adespatial)
library(plyr)
library(vegan)
library(dplyr)
library(tidyr)
library(DHARMa)
library(ggplot2)
library(gridExtra)

# set working directory
setwd("~/PostDoc_Ghent/Synthesis_3/treeweb_community/")

# load plot data
plot_info <- read.csv("data/synthesis_expldata_raw.csv", stringsAsFactors = FALSE)
plot_xy <- read.csv("data/plot_xy.csv")

# load community data
comms <- c("bird", "bat", "spider", "opiliones", "carabid", "isopods", "diplopod", "herb", "vegetation")
comm_l <- lapply(comms, function(x) read.csv(paste0("data/treeweb_", x, "_format.csv")))
names(comm_l) <- comms

# remove the species with no record in the bird table
comm_l$bird <- comm_l$bird[,-which(colSums(comm_l$bird) == 0)]
# divide herbivore by sampling effort
comm_l$herb <- cbind(comm_l$herb[,1], comm_l$herb[,-1] / plot_info$specrich)
# remove genus level id in isopods
comm_l$isopods <- comm_l$isopods[,-ncol(comm_l$isopods)]

# define spatial neighbors based on distance and fragment ID
xx <- dnearneigh(as.matrix(plot_xy[,c("X", "Y")]), d1 = 0, d2 = 450)
xx[[45]] <- as.integer(c(36, 37)) # plot 45 is in the same fragment than plots 36 and 37

# distance between neighbors
dist_xx <- nbdists(xx, as.matrix(plot_xy[,c("X", "Y")]))
# weight the neighbors by their distance relative to the maximum distance
fdist <- lapply(dist_xx, function(x) 1 - x/max(dist(as.matrix(plot_xy[,c("X", "Y")]))))

# generate listw
listw_xx <- nb2listw(xx, glist = fdist, zero.policy = TRUE)

######### Part 1: co-inertia analysis #######

############## define helper functions
### this first function compute the percent
### of variation explained by the first two
### co-inertia axis, arguments are:
# @pair: a character string of the pair of taxa, like: "bird_bat"
# @nf_1, nf_2, nf_coia : an integer, the number of axis to keep in the PCAs and for the coinertia
# @std : a logical, whether to hellinger-standardize the community data
get_percent_explained <- function(pair, nf_1 = 2, nf_2 = 2, nf_coia = 2, std = TRUE){
  # get names of the two taxa
  taxa_1 <- strsplit(pair, "_")[[1]][1]
  taxa_2 <- strsplit(pair, "_")[[1]][2]
  # get community matrices from the list of communities
  comm1 <- comm_l[[taxa_1]][,-1]
  comm2 <- comm_l[[taxa_2]][,-1]
  # if needed standardize the matrices
  if(std){
    comm1 <- decostand(comm1, method = "hellinger")
    comm2 <- decostand(comm2, method = "hellinger")
  }
  # run the PCA
  dudi_1 <- dudi.pca(comm1, scannf = FALSE, nf = nf_1)
  dudi_2 <- dudi.pca(comm2, scannf = FALSE, nf = nf_2)
  # run the co-inertia
  coia <- coinertia(dudi_1, dudi_2, scannf = FALSE, nf = nf_coia)
  # % of explained dev in the first two axis
  eigs <- coia$eig
  percent <- sum(eigs[1:2]) / sum(eigs)
  
  return(percent)
  
}

### a function to get the RV metric together with p-value generated via spatial null model
### arguments to the function are:
# @pair: a character, the pair of taxa like: "bird_bat"
# @nf_1, nf_2, nf_coia: an integer, the number of axis to keep for the PCAs and the co-inertia
# @std: a logical, whether to hellinger-standardize the community
# @repet: the number of spatially null communities to generate
get_coia_rv <- function(pair, nf_1 = 2, nf_2 = 2, nf_coia = 2, std = TRUE, nrepet = 99){
  
  taxa_1 <- strsplit(pair, "_")[[1]][1]
  taxa_2 <- strsplit(pair, "_")[[1]][2]
  
  comm1 <- comm_l[[taxa_1]][,-1]
  comm2 <- comm_l[[taxa_2]][,-1]
  
  if(std){
    comm1 <- decostand(comm1, method = "hellinger")
    comm2 <- decostand(comm2, method = "hellinger")
  }
  
  dudi_1 <- dudi.pca(comm1, scannf = FALSE, nf = nf_1)
  dudi_2 <- dudi.pca(comm2, scannf = FALSE, nf = nf_2)
  
  coia <- coinertia(dudi_1, dudi_2, scannf = FALSE, nf = nf_coia)
  # get a vector of RV values from null communities
  rv_null <- get_null_coia_rv(pair, nf_1, nf_2, nf_coia, nrepet) 
  # compute the p value, which is the proportion of simulation with RV values larger than the observed one
  pval_null <- 1 - (sum(coia$RV >= rv_null) / nrepet)
  
  out <- data.frame(rv = coia$RV, pval = pval_null, col = ifelse(pval_null < 0.05, "red", "black"))
  
  return(out)
}

# function to get nrepet null RV values based on moran spectral randomization
get_null_coia_rv <- function(pair, nf_1 = 2, nf_2 = 2, nf_coia = 2, nrepet = 99){
  
  taxa_1 <- strsplit(pair, "_")[[1]][1]
  taxa_2 <- strsplit(pair, "_")[[1]][2]
  
  comm1 <- comm_l[[taxa_1]][,-1]
  comm2 <- comm_l[[taxa_2]][,-1]
  # compute moran spectral randomization on the community based on listw object
  comm_n1 <- msr(comm1, listw_xx, method = "singleton", nrepet = nrepet)
  comm_n2 <- msr(comm2, listw_xx, method = "singleton", nrepet = nrepet)
  # run the PCA and co-inertia analysis on these MSR generated communities
  dudi_n1 <- lapply(comm_n1, function(x) dudi.pca(x, scannf = FALSE, nf = 2))
  dudi_n2 <- lapply(comm_n2, function(x) dudi.pca(x, scannf = FALSE, nf = 2))
  
  coia_n <- mapply(function(x, y) coinertia(x, y, scannf = FALSE, nf = 2), x = dudi_n1, y = dudi_n2, SIMPLIFY = FALSE)
  # output the RV value
  rv_null <- ldply(coia_n, function(x) x$RV)
  
  return(rv_null$V1)
}

########## end of helper function definition

# all potential pairs of communities
pairs <- outer(comms, comms, paste, sep = "_")
pairs_u <- pairs[upper.tri(pairs)]

# define the trophic levels
tl <- data.frame(taxa = comms,
                 tl = c(4, 4, 3, 3, 3, 2, 2, 2, 1))

# put all pairs in a data frame together with their trophic distance
pairs_d <- data.frame(pairs = pairs_u, stringsAsFactors = FALSE)
pairs_d <- separate(pairs_d, "pairs", c("from", "to"), sep = "_", remove = FALSE)

pairs_d <- merge(pairs_d, tl, by.x = "from", by.y = "taxa", sort = FALSE)
pairs_d <- merge(pairs_d, tl, by.x = "to", by.y = "taxa", sort=FALSE)
pairs_d$dist <- with(pairs_d, sqrt((tl.x - tl.y)**2))


## get percent explained variation of the first two axis in the coinertia
tt <- sapply(pairs_u, function(x) get_percent_explained(x))

# RV value between each comms
rv_all <- sapply(pairs_u, function(x) get_coia_rv(x, nrepet = 9999), simplify = FALSE) # this can take time to run depending on the number of repetitions
# some data wraggling to put te results in a plotable format
rv_all <- rbind.fill(rv_all)
rv_all$pair <- pairs_u
rv_all <- merge(rv_all, pairs_d, by.x = "pair", by.y = "pairs")
rv_all <- arrange(rv_all, dist, pair)
# some re-coding for nicer plots
rv_all$to <- recode(rv_all$to, bat = "bats", herb = "herbivores", diplopod = "millipedes",
                    carabid = "carabids", opiliones = "harvestmen", spider = "spiders",
                    isopods = "woodlice", vegetation = "vegetation")

rv_all$from <- recode(rv_all$from, bat = "bats", herb = "herbivores", diplopod = "millipedes",
                      carabid = "carabids", opiliones = "harvestmen", spider = "spiders",
                      isopods = "woodlice", vegetation = "vegetation")

rv_all <- unite(rv_all, "pair", c(to, from), sep =":")
rv_all$pair <- factor(rv_all$pair)
rv_all$x <- 1:36

# the plot (figure X in the manuscript)
gg_rv <- ggplot(rv_all, aes(x=x, y=rv, fill=col)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_brewer(type = "qual",palette = "Set1",
                    labels = c("no sig. relation", "sig. relation"),
                    name = "RV permutation\ntest:") +
  geom_vline(xintercept = c(7.5, 25.5, 34.5), linetype = "dashed") +
  annotate("text", 2.75, y = 0.6, label = "Trophic\ndistance: 0") +
  annotate("text", 18.75, y = 0.6, label = "Trophic\ndistance: 1") +
  annotate("text", 27.75, y = 0.6, label = "Trophic\ndistance: 2") +
  annotate("text", 35.5, 0.6, label = "3") +
  scale_x_continuous(breaks = 1:36, labels = rv_all$pair, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.67), expand = c(0, 0)) +
  coord_flip() +
  theme(axis.text.x = element_text(hjust = 0),
        axis.text.y = element_text(vjust = 0.25)) +
  labs(x = "Pairwise community comparison",
       y = "RV value") +
  theme_bw()

ggsave("figures/coin_rv.png", gg_rv)

################# Part 2: multiple co-inertia analysis ########

# an helper function to run the mcoa
# the arguments are:
# @web: a character, the web of interaction to quantify like "bird_bat_spider"
# @nf: an integer, the number of axis to keep
# @option: a character, the weighing used, see ?mcoa for the different options
# @std: a logical, whether to hellinger-standardize the data
mcoa_run <- function(web, nf = 2, option = "inertia", std = TRUE){
  comm_name <- strsplit(web, "_")[[1]]
  
  comms <- lapply(comm_name, function(x) comm_l[[x]])
  
  if(std){
    comms <- lapply(comms, function(comm) decostand(comm[,-1], "hellinger"))
  }
  
  # grab number of columns
  nb_col <- sapply(comms, function(x) ncol(x))
  
  # put all together
  comms_all <- do.call(cbind, comms)
  # create ktab object
  kt <- ktab.data.frame(comms_all, nb_col,
                        tabnames = comm_name)
  
  # run mcoa
  m <- mcoa(kt, option = option, nf = nf, scannf = FALSE)
  
  return(m)
}

## first web of interaction: all groups except for the bats:
web_1 <- "bird_spider_opiliones_carabid_isopods_diplopod_herb_vegetation"
m_1 <- mcoa_run(web_1)

## compute the distance of each group from the plot synthetic scores
plot_info$x <- plot_xy$X
plot_info$y <- plot_xy$Y
vhb_xy <- m_1$Tl1
vhb_xy$id_plot <- rep(1:53, length(web_1))
vhb_xy$ref_x <- m_1$SynVar[rep(1:53, length(web_1)),1]
vhb_xy$ref_y <- m_1$SynVar[rep(1:53, length(web_1)),2]
vhb_xy$d <- with(vhb_xy, sqrt((Axis1-ref_x) ** 2 + (Axis2-ref_y) ** 2))
vhb_xy %>%
  group_by(id_plot) %>%
  summarise(d = sum(d)) %>%
  left_join(plot_info[,c(1,13,14,15, 16, 17)], by = "id_plot") -> vhb_dd

# fit a glm to this
m_vhb <- glm(d ~ fragm_std * speccomb, vhb_dd, family = Gamma(link="log"), contrasts = list(speccomb = "contr.sum"))

## look at residuals
ss <- simulateResiduals(m_vhb)
plot(ss) # good
testSpatialAutocorrelation(ss, x = plot_info$x, y = plot_info$y) # no spatial autocorrelation

## anova
anova(m_vhb, test = "Chisq")

## plot the composition effects
### a new data frame to derive the model predictions
newdat <- expand.grid(speccomb = unique(plot_info$speccomb), fragm_std = 0)

# augment per diversification path
newdat2 <- newdat[c(2, 5, 3, 4, 7, 5, 1, 4, 6, 3, 1, 4),]
newdat2$tree <- rep(c("fsyl", "qrob", "qrub"), each = 4)
# compute the model predictions
newdat2$pred <- predict(m_vhb, newdata = newdat2, type = "response")
newdat2$se <- predict(m_vhb, newdata = newdat2, type = "response", se.fit = TRUE)$se.fit

# some refinments for ordering the species somposition as we want them
newdat2$speccomb <- factor(newdat2$speccomb, levels = c("fsyl","qrob","qrub","fsyl_qrob","fsyl_qrub","qrob_qrub","all"))
vhb_dd$speccomb <- factor(vhb_dd$speccomb, levels = c("fsyl","qrob","qrub","fsyl_qrob","fsyl_qrub","qrob_qrub","all"))

# the plots
gg_syl <- ggplot(subset(vhb_dd, speccomb %in% c("fsyl", "fsyl_qrob", "fsyl_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "fsyl"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "fsyl"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("beech", "beech +\nped. oak", "beech +\nred oak", "all")) +
  labs(x = "", y = "", title = "(a)") +
  ylim(c(2.25, 8.55))

gg_rob <- ggplot(subset(vhb_dd, speccomb %in% c("qrob", "fsyl_qrob", "qrob_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "qrob"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "qrob"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("ped. oak", "beech +\nped. oak", "ped. oak +\nred oak", "all")) +
  labs(x = "", y = "", title = "(b)") +
  ylim(c(2.25, 8.55))

gg_rub <- ggplot(subset(vhb_dd, speccomb %in% c("qrub", "fsyl_qrub", "qrob_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "qrub"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "qrub"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "(c)") +
  scale_x_discrete(labels = c("red oak", "beech +\nred oak", "ped. oak +\nred oak", "all")) +
  ylim(c(2.25, 8.55))


gg_all <- grid.arrange(gg_syl, gg_rob, gg_rub, nrow = 3, bottom = "Tree species composition", left = "Community distance")

ggsave("figures/mcoa_composition_2.png", gg_all, width = 4, height = 8)

# supporting information graph displaying
# the distance of all taxa from the plot centroid scores
# this will be done separately for all diversification pathways

## for fsyl paths
vhb_syl <- subset(vhb_xy, id_plot %in% id_fsyl)
vhb_syl <- arrange(vhb_syl, id_plot)

vhb_syl <- merge(vhb_syl, plot_info, by = "id_plot")
vhb_syl$speccomb <- factor(vhb_syl$speccomb, levels = c("fsyl", "fsyl_qrob", "fsyl_qrub", "all"))


gg_syl <- ggplot(vhb_syl, aes(x=ref_x, y=ref_y, group=id_plot, color = fragm_std)) +
  geom_point(aes(x = Axis1, y = Axis2), size = 0.5) +
  geom_segment(aes(xend=Axis1, yend=Axis2)) +
  geom_point() +
  facet_wrap(~ speccomb, ncol = 4) +
  scale_color_viridis(name = "Fragm.\nintensity", option = "C", 
                      limits = c(-1.7, 2.45)) +
  labs(x="", y = "", title = "(a)") +
  xlim(c(-1.4, 0.5)) +
  ylim(c(-2.3,2.5))

## now for qrob
vhb_rob <- subset(vhb_xy, id_plot %in% id_qrob)
vhb_rob <- arrange(vhb_rob, id_plot)

vhb_rob <- merge(vhb_rob, plot_info, by = "id_plot")
vhb_rob$speccomb <- factor(vhb_rob$speccomb, levels = c("qrob", "fsyl_qrob", "qrob_qrub", "all"))


gg_rob <- ggplot(vhb_rob, aes(x=ref_x, y=ref_y, group=id_plot, color = fragm_std)) +
  geom_point(aes(x = Axis1, y = Axis2), size = 0.5) +
  geom_segment(aes(xend=Axis1, yend=Axis2)) +
  geom_point() +
  facet_wrap(~ speccomb, ncol = 4) +
  scale_color_viridis(name = "Fragm.\nintensity", option = "C",
                      limits = c(-1.7, 2.45)) +
  labs(x="", y = "", title = "(b)")  +
  xlim(c(-1.4, 0.5)) +
  ylim(c(-2.3,2.5))

## now for qrub
vhb_rub <- subset(vhb_xy, id_plot %in% id_qrub)
vhb_rub <- arrange(vhb_rub, id_plot)

vhb_rub <- merge(vhb_rub, plot_info, by = "id_plot")
vhb_rub$speccomb <- factor(vhb_rub$speccomb, levels = c("qrub", "fsyl_qrub", "qrob_qrub", "all"))


gg_rub <- ggplot(vhb_rub, aes(x=ref_x, y=ref_y, group=id_plot, color = fragm_std)) +
  geom_point(aes(x = Axis1, y = Axis2), size = 0.5) +
  geom_segment(aes(xend=Axis1, yend=Axis2)) +
  geom_point() +
  facet_wrap(~ speccomb, ncol = 4) +
  scale_color_viridis(name = "Fragm.\nintensity", option = "C",
                      limits = c(-1.7, 2.45)) +
  labs(x="", y = "", title = "(c)")  +
  xlim(c(-1.4, 0.5)) +
  ylim(c(-2.3,2.5))

# tweaking to get common legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(gg_rob)

# the plot
gg_all <- grid.arrange(arrangeGrob(gg_syl + theme(legend.position = "none"),
                                   gg_rob + theme(legend.position = "none"),
                                   gg_rub + theme(legend.position = "none"), nrow = 3),
                       mylegend, ncol = 2, widths = c(9, 2),
                       bottom = "MCOA axis 1", left = "MCOA axis 2")

ggsave("figures/mcoa_plot.png", gg_all, width = 8, height = 8)



########## end of the main script ############