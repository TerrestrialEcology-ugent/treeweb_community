#######################
# R script for the Treeweb       #
# cross-community analysis       #
# SI analysis                    #
#######################################

# load libraries
library(vegan)
library(plyr)
library(reshape2)
library(boot)
library(dplyr)
library(UpSetR)
library(DHARMa)
library(ade4)
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

########### Part 1: RDA per trophic group ###########

# an helper function to run the RDA
# arguments are:
# @comm: a character, the name of the taxa
# @interaction: a logical, whether the include an interaction term between species composition and fragmentation
# @decostand: a logical, whether to hellinger-standardize the community data
# @species: a character, the tree species in focus
# @...: further argument to pass to the RDA, see ?rda

fit_rda <- function(comm, interaction = TRUE, decostand = TRUE, species = "qrob",...){
  
  print(species)
  # grab the plot index for the focal tree species
  ids <- eval(as.name(paste0("id_",species)))
  # if needed standardize
  if(decostand){
    comm <- decostand(comm, method = "hellinger", MARGIN = 1)
  }
  
  # compute PCNM
  rs <- rowSums(comm) / sum(comm) # the site weights
  pcnmw <- pcnm(dist(plot_xy[,c("X", "Y")]), w = rs) # PCNM based on plot coordinates
  # add the PCNM scores to the dataset
  dat <- cbind(plot_info[,c("speccomb", "fragm_std")], scores(pcnmw))
  # NA for opiliones, fix by putting PCNM to 0
  dat[is.na(dat)] <- 0
  # fit the RDA with or without interaction
  if(interaction){
    rr <- rda(comm[ids,] ~ . + speccomb : fragm_std, dat[ids,], ...)
  } else {
    rr <- rda(comm[ids,] ~ ., dat[ids,],...)
  }
  
  return(rr)
}

### fit the RDA
## plot id for the specific tree composition changes
id_qrob <- which(plot_info$speccomb %in% c("qrob", "fsyl_qrob","qrob_qrub", "all"))
id_qrub <- which(plot_info$speccomb %in% c("qrub", "fsyl_qrub","qrob_qrub", "all"))
id_fsyl <- which(plot_info$speccomb %in% c("fsyl", "fsyl_qrob","fsyl_qrub", "all"))
## the different tree species code
spp <- c("qrob", "qrub", "fsyl")

# go through the different communities
# with interaction
rda_l <- lapply(comm_l, function(comm) sapply(spp, function(sp) fit_rda(comm[,-1], interaction = TRUE, decostand = TRUE, species = sp), 
                                              simplify = FALSE))
### analyze the results
## first table of effects
# grab the RDA results for the Table
## an helper function that grab p-values of the effects and R-square per taxa
eff_rda <- function(taxa = "bird"){
  rr <- rda_l[[taxa]]
  tmp <- ldply(rr, function(x) as.data.frame(anova(x, by = "terms", permutations = how(nperm = 9999))))
  tmp$terms <- rep(c("composition", "fragmentation", paste0("P", 1:14), "interaction", "residual"), 3)
  tmp <- filter(tmp, terms %in% c("composition", "fragmentation", "interaction"))
  tmp$taxa <- taxa
  
  # add R2
  tmp2 <- ldply(rr, function(x) as.data.frame(RsquareAdj(x)$adj.r.squared))
  tmp2$terms <- "R2"
  names(tmp2)[2] <- names(tmp)[5]
  tmp2$taxa <- taxa
  
  tmp <- rbind(tmp2, tmp[,c(1,5,6,7)])
  
  return(tmp)  
}

# apply to all taxa
tabb <- rbind.fill(sapply(comms, function(sp) eff_rda(sp), simplify = FALSE))
# some data wraggling
tabb_w <- dcast(tabb, taxa ~ .id + terms, value.var = "Pr(>F)")
xtable(tabb_w[c(1:3,8, 7,6,4:5,9),c(1,5,2,3,4,9,6,7,8, 13, 10:12)], digits = 3)

## second, effect size of composition / fragmentation effects
# per tree diversification pathways develop indices of tree composition changes
change_qrob <- c("qrob -> qrob_fsyl", "qrob -> qrob_qrub",
                 "qrob_fsyl -> all", "qrob_qrub -> all", "qrob -> all")

change_qrub <- c("qrub -> fsyl_qrub", "qrub -> qrob_qrub",
                 "fsyl_qrub -> all", "qrob_qrub -> all", "qrub -> all")

change_fsyl <- c("fsyl -> qrob_fsyl", "fsyl -> fsyl_qrub",
                 "qrob_fsyl -> all", "fsyl_qrub -> all", "fsyl -> all")

# an helper function to get the distance between centroids of different tree composition
# and also to get the distance of changes due to the fragmentation effect
# the arguments are:
# @rr: an object of class rda, the fitted RDA from vegan
# @tree: a character, the focal tree species, like "qrob"
get_length_rda <- function(rr, tree){
  # define the composition changes
  change_comp <- eval(as.name(paste0("change_", tree)))
  # define the different composition
  change_fragm <- c(grep(tree, unique(plot_info$speccomb), value = TRUE), "all")
  # a new data frame to compute predicted shifts along the gradient of fragmentation
  newdat <- expand.grid(speccomb = change_fragm, fragm_std = c(-1.69, 2.42),
                        PCNM1 = 0, PCNM2 = 0, PCNM3 = 0, PCNM4 = 0,
                        PCNM5 = 0, PCNM6 = 0, PCNM7 = 0, PCNM8 = 0,
                        PCNM9 = 0, PCNM10 = 0, PCNM11 = 0, PCNM12 = 0,
                        PCNM13 = 0, PCNM14 = 0)
  # get the RDA centroid scores
  sc_cn <- scores(rr, choices = c(1,2), display = "cn", scaling = 0)
  # depeding on the focal tree compute the shifts in centroids
  if(tree == "qrob"){
    c1 <- dist(rbind(sc_cn[3,], sc_cn[2,]))
    c2 <- dist(rbind(sc_cn[3,], sc_cn[4,]))
    c3 <- dist(rbind(sc_cn[2,], sc_cn[1,]))
    c4 <- dist(rbind(sc_cn[4,], sc_cn[1,]))
    c5 <- dist(rbind(sc_cn[3,], sc_cn[1,]))
  }
  if(tree == "qrub"){
    c1 <- dist(rbind(sc_cn[4,], sc_cn[2,]))
    c2 <- dist(rbind(sc_cn[4,], sc_cn[3,]))
    c3 <- dist(rbind(sc_cn[2,], sc_cn[1,]))
    c4 <- dist(rbind(sc_cn[3,], sc_cn[1,]))
    c5 <- dist(rbind(sc_cn[4,], sc_cn[1,]))
  }
  if(tree == "fsyl"){
    c1 <- dist(rbind(sc_cn[2,], sc_cn[3,]))
    c2 <- dist(rbind(sc_cn[2,], sc_cn[4,]))
    c3 <- dist(rbind(sc_cn[3,], sc_cn[1,]))
    c4 <- dist(rbind(sc_cn[4,], sc_cn[1,]))
    c5 <- dist(rbind(sc_cn[2,], sc_cn[1,]))
  }
  # compute the shift along the fragmentation gradient
  sc_fragm <- predict(rr, newdata = newdat, type = "lc")
  fragm <- c(dist(sc_fragm[c(1, 5),1:2]), 
             dist(sc_fragm[c(2, 6),1:2]),
             dist(sc_fragm[c(3, 7),1:2]),
             dist(sc_fragm[c(4, 8),1:2])) # if we do not restrict the dist to the first two axis we always get the same dist: 0.76
  
  data.frame(effect = rep(c("composition", "fragmentation"), times = c(5, 4)), change = c(change_comp, change_fragm), length = c(c1, c2, c3, c4, c5, fragm))
  
}

# an helper function to compute the bootstrapped shifts
# arguments are:
# @dat: a data.frame or matrix, the community data
# @id: internal argument for boot
# @N: an integer, the number of species
# @tree: a character, the focal tree species like "qrob"
boot_rda <- function(dat, id, N, tree){
  # the community data
  comm <- dat[,1:N]
  # the predicotrs
  plot_info <- dat[,(N+1):ncol(dat)]
  # standardize the community
  comm <- decostand(comm[,-1], "hellinger", MARGIN = 1)
  # run the RDA with interaction
  rr <- rda(comm[id,] ~ speccomb * fragm_std, plot_info[id,])
  # grab the effect sizes
  out <- get_length_rda(rr, tree)
  return(out$length)
}

# an helper function combining the effect computation with the bootstrapped computation
# the arguments are:
# @taxa: a character, the name of the focal taxa like "bird"
# @tree: a character, the name of the focal tree species like "qrob"
# @R: an integer, the number of bootstrap samples
get_length_df <- function(taxa = "bird", tree = "qrob", R = 10){
  print(taxa) # for debug purposes
  # grab the relevant RDA
  rr <- rda_l[[taxa]]
  # grab the relevant community matrix
  comm <- comm_l[[taxa]]
  # grab the plot index for the tree diversification pathways
  ids <- eval(as.name(paste0("id_", tree)))
  # the species numer
  n_sp <- ncol(comm)
  # put together the community and the predictor data
  boot_dat <- cbind(comm[ids,], plot_info[ids,])
  # compute the effect size
  tab <- get_length_rda(rr[[tree]], tree)
  # compute the bootstrapped effect size
  bb <- boot(data = boot_dat, statistic = boot_rda, R = R, N = n_sp, tree = tree, strata = factor(boot_dat$speccomb))
  # the bootstrapped standard deviation
  tab$std_err <- apply(bb$t, 2, sd)
  
  return(tab)
}

# run this for all taxa and trees
c_rob <- rbind.fill(sapply(comms, function(taxa) get_length_df(taxa = taxa, tree = "qrob", R = 10000), simplify = FALSE))
c_rub <- rbind.fill(sapply(comms, function(taxa) get_length_df(taxa = taxa, tree = "qrub", R = 10000), simplify = FALSE))
c_syl <- rbind.fill(sapply(comms, function(taxa) get_length_df(taxa = taxa, tree = "fsyl", R = 10000), simplify = FALSE))
# put everything together
c_all <- rbind(c_rob, c_rub, c_syl)
# some re-naming for the plots
c_all$tree <- rep(c("qrob", "qrub", "fsyl"), each = 81)
c_all$taxa <- rep(comms, each = 9)
c_all$taxa[c_all$taxa=="herb"] <- "herbivore"
c_all$taxa <- factor(c_all$taxa, 
                     levels = c(comms[1:7], "herbivore", "vegetation"))
c_all$change <- factor(c_all$change, levels = c(change_qrob[1:2], change_qrub[1:2], change_fsyl[1:2], change_qrob[3:4], change_qrub[3], change_qrob[5], change_qrub[5], change_fsyl[5], "qrob", "fsyl", "qrub", "fsyl_qrob", "fsyl_qrub", "qrob_qrub", "all"))

# first plot on the fragmentation effect
c_fragm <- subset(c_all, effect == "fragmentation")
# the position of the bars in the graph
c_fragm$X <- c(rep(c(3, 2, 1, 4), 9), rep(c(8, 7, 6, 9), 9), rep(c(11, 12, 13, 14), 9))

# compute the average per taxa and div paths
c_fragm %>%
  group_by(taxa, tree) %>%
  summarise(avg = mean(length)) -> avg_eff
# end of the dotted average line on the graphs
avg_eff$x_start <- rep(c(10.5, 0.5, 5.5), times = 9)
avg_eff$x_end <- rep(c(14.5, 4.5, 9.5), times = 9)

# nicer name for the facets
facet_n <- c("bird" = "birds", "bat" = "bats", "spider" = "spiders",
                "opiliones" = "harvestmen", "carabid" = "carabids",
                "isopods" = "woodlice", "diplopod" = "millipedes",
                "herbivore" = "herbivores", "vegetation" = "vegetation")

# the plot
gg_f <- ggplot(c_fragm,aes(x=X, y=length, color=tree))+
  geom_bar(stat = "identity", fill = "white")+
  geom_linerange(aes(ymin=length, ymax = length + 2 * std_err))+
  scale_x_continuous(breaks = c(1:4, 6:9, 11:14), labels = c("qrob", "fsyl_qrob", "qrob_qrub", "all", 
                                                             "qrub", "fsyl_qrub", "qrob_qrub", "all",
                                                             "fsyl", "fsyl_qrub", "fsyl_qrob", "all")) +
  facet_wrap(~taxa, labeller = as_labeller(facet_n)) +
  geom_segment(data = avg_eff, aes(x = x_start, xend = x_end, y = avg, yend = avg),
               size = 1, linetype = "dashed") +
  labs(x = "Diversification pathway",
       y = "Strength of fragmentation effect (95% bootstraped CI)"
       ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90))

ggsave("figures/fragmentation_rda.png", gg_f, width = 22.77,
       height = 14.71, units = "cm")

# now the compositional change in a similar graph
c_comp <- subset(c_all, effect == "composition")
c_comp$X <- c(rep(1:5, 9), rep(7:11, 9), rep(13:17, 9))

gg_c <- ggplot(subset(c_comp, change != "fragmentation"), 
               aes(x=X, y=length, group = paste0(change, tree), color = tree)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), fill = "white")+
  geom_linerange(aes(ymin=length, ymax = length + 2* std_err), position = position_dodge(width = 0.9))+
  facet_wrap(~taxa, labeller = as_labeller(facet_n)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9), 
        legend.position = "none") +
  scale_x_continuous(breaks= c(1:5, 7:11, 13:17), labels = subset(c_comp, taxa == "bird")$change) +
  labs(x = "Diversification pathway",
       y = "Strength of diversification effect (95% bootstraped CI)")

ggsave("figures/composition_rda.png", gg_c, width = 22.77,
       height = 14.71, units = "cm")

########### Part 2: compositional turnover across groups, using upset graphs #####

# an helper function formatting the comminty data and outputting the desired plots
upset_comm <- function(comm){
  
  if(names(comm)[1] != "id_plot"){
    names(comm)[1]<- "id_plot"
  }
  
  #format the community
  comm %>%
    pivot_longer(cols = 2:ncol(comm), names_to = "species", values_to = "abundance") %>%
    #rename(id_plot = plot) %>%
    left_join(plot_info, by = "id_plot") %>%
    group_by(speccomb, species) %>%
    summarise(Present = sum(abundance, na.rm = TRUE)) %>%
    mutate(Present = ifelse(Present > 0, 1, 0)) %>%
    pivot_wider(names_from = speccomb, values_from = Present) -> comm_dd
  
  comm_dd <- as.data.frame(comm_dd)
  # one plot per diversification pathways
  u_rob <- upset(comm_dd, sets = c("qrob", "fsyl_qrob", "qrob_qrub", "all"), order.by = "freq",
                 mainbar.y.label = "Number of (shared) species", sets.x.label = "Total species\nnumber",
                 main.bar.color="#1b9e77", matrix.color="#1b9e77",
                 sets.bar.color="#d95f02")
  
  
  u_robc <- cowplot::plot_grid(NULL, u_rob$Main_bar, u_rob$Sizes, u_rob$Matrix,
                               nrow=2, align='hv', rel_heights = c(3,1),
                               rel_widths = c(2,3))
  
  u_rub <- upset(comm_dd, sets = c("qrub", "fsyl_qrub", "qrob_qrub", "all"), order.by = "freq",
                 mainbar.y.label = "Number of (shared) species", sets.x.label = "Total species\nnumber",
                 main.bar.color="#1b9e77", matrix.color="#1b9e77",
                 sets.bar.color="#d95f02")
  
  
  u_rubc <- cowplot::plot_grid(NULL, u_rub$Main_bar, u_rub$Sizes, u_rub$Matrix,
                               nrow=2, align='hv', rel_heights = c(3,1),
                               rel_widths = c(2,3))
  
  u_syl <- upset(comm_dd, sets = c("fsyl", "fsyl_qrub", "fsyl_qrob", "all"), order.by = "freq",
                 mainbar.y.label = "Number of (shared) species", sets.x.label = "Total species\nnumber",
                 main.bar.color="#1b9e77", matrix.color="#1b9e77",
                 sets.bar.color="#d95f02")
  
  
  u_sylc <- cowplot::plot_grid(NULL, u_syl$Main_bar, u_syl$Sizes, u_syl$Matrix,
                               nrow=2, align='hv', rel_heights = c(3,1),
                               rel_widths = c(2,3))
  
  
  gga <- grid.arrange(u_robc, u_rubc, u_sylc, ncol = 3)
  
  
  return(gga)
}

# run this for all taxa
ggsave("figures/upset_bird.png", upset_comm(comm_l$bird), width = 11, height = 5, units = "in")
ggsave("figures/upset_bat.png", upset_comm(comm_l$bat), width = 11, height = 5, units = "in")
ggsave("figures/upset_carabid.png", upset_comm(comm_l$carabid), width = 11, height = 5, units = "in")
ggsave("figures/upset_spider.png", upset_comm(comm_l$spider), width = 11, height = 5, units = "in")
ggsave("figures/upset_isopod.png", upset_comm(comm_l$isopods), width = 11, height = 5, units = "in")
ggsave("figures/upset_diplopod.png", upset_comm(comm_l$diplopod), width = 11, height = 5, units = "in")
ggsave("figures/upset_herbivore.png", upset_comm(comm_l$herb), width = 11, height = 5, units = "in")
ggsave("figures/upset_vegetation.png", upset_comm(comm_l$vegetation), width = 11, height = 5, units = "in")
ggsave("figures/upset_opiliones.png", upset_comm(comm_l$opiliones), width = 11, height = 5, units = "in")


########### Part 3: testing additional community webs ########

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

##full web of interaction: all groups:
web_2 <- "bird_bat_spider_opiliones_carabid_isopods_diplopod_herb_vegetation"
n_taxa <- 9
m_2 <- mcoa_run(web_2)

## compute the distance of each group from the plot synthetic scores
plot_info$x <- plot_xy$X
plot_info$y <- plot_xy$Y
vhb_xy <- m_2$Tl1
vhb_xy$id_plot <- rep(1:53, n_taxa)
vhb_xy$ref_x <- m_2$SynVar[rep(1:53, n_taxa),1]
vhb_xy$ref_y <- m_2$SynVar[rep(1:53, n_taxa),2]
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
anova(m_vhb, test = "Chisq") # species composition effect

## make plot of predicted distance for the different tree composition
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
  ylim(c(3.8, 10.1))

gg_rob <- ggplot(subset(vhb_dd, speccomb %in% c("qrob", "fsyl_qrob", "qrob_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "qrob"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "qrob"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("ped. oak", "beech +\nped. oak", "ped. oak +\nred oak", "all")) +
  labs(x = "", y = "", title = "(b)") +
  ylim(c(3.8, 10.1))

gg_rub <- ggplot(subset(vhb_dd, speccomb %in% c("qrub", "fsyl_qrub", "qrob_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "qrub"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "qrub"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "(c)") +
  scale_x_discrete(labels = c("red oak", "beech +\nred oak", "ped. oak +\nred oak", "all")) +
  ylim(c(3.8, 10.1))


gg_all <- grid.arrange(gg_syl, gg_rob, gg_rub, nrow = 3, bottom = "Tree species composition", left = "Community distance")

ggsave("figures/mcoa_composition_web2.png", gg_all, width = 4, height = 8)


## now only vegetation - bird - carabid - spider
web_3 <- "bird_carabid_spider_vegetation"
n_taxa <- 4
m_3 <- mcoa_run(web_3)

## compute the distance of each group from the plot synthetic scores
plot_info$x <- plot_xy$X
plot_info$y <- plot_xy$Y
vhb_xy <- m_3$Tl1
vhb_xy$id_plot <- rep(1:53, n_taxa)
vhb_xy$ref_x <- m_3$SynVar[rep(1:53, n_taxa),1]
vhb_xy$ref_y <- m_3$SynVar[rep(1:53, n_taxa),2]
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
anova(m_vhb, test = "Chisq") # species composition effect

## make plot of predicted distance for the different tree composition
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
  ylim(c(0.5, 5.2))

gg_rob <- ggplot(subset(vhb_dd, speccomb %in% c("qrob", "fsyl_qrob", "qrob_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "qrob"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "qrob"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("ped. oak", "beech +\nped. oak", "ped. oak +\nred oak", "all")) +
  labs(x = "", y = "", title = "(b)") +
  ylim(c(0.5, 5.2))

gg_rub <- ggplot(subset(vhb_dd, speccomb %in% c("qrub", "fsyl_qrub", "qrob_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "qrub"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "qrub"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "(c)") +
  scale_x_discrete(labels = c("red oak", "beech +\nred oak", "ped. oak +\nred oak", "all")) +
  ylim(c(0.5, 5.2))


gg_all <- grid.arrange(gg_syl, gg_rob, gg_rub, nrow = 3, bottom = "Tree species composition", left = "Community distance")

ggsave("figures/mcoa_composition_web3.png", gg_all, width = 4, height = 8)

##### Part 4: Co-inertia analysis with new groupings ######

# define 5 groups, herbivores, primary producers, predators, detritivores and omnivores
comm_l$primary <- comm_l$vegetation
comm_l$herbivore <- comm_l$herb
comm_l$detritivore <- cbind(comm_l$isopods, comm_l$diplopod[,-1])
comm_l$omnivore <- cbind(comm_l$carabid, comm_l$opiliones[,-1], comm_l$bird[,-1])
comm_l$predator <- cbind(comm_l$bat, comm_l$spider[,-1])

# all potential pairs of communities
comms2 <- c("primary", "herbivore", "detritivore", "omnivore", "predator")
pairs2 <- outer(comms2, comms2, paste, sep = "_")
pairs_u2 <- pairs2[upper.tri(pairs2)]

# define the trophic levels
tl2 <- data.frame(taxa = comms2,
                 tl = c(1, 2, 2, 3, 4))

# put all pairs in a data frame together with their trophic distance
pairs_d2 <- data.frame(pairs = pairs_u2, stringsAsFactors = FALSE)
pairs_d2 <- separate(pairs_d2, "pairs", c("from", "to"), sep = "_", remove = FALSE)

pairs_d2 <- merge(pairs_d2, tl2, by.x = "from", by.y = "taxa", sort = FALSE)
pairs_d2 <- merge(pairs_d2, tl2, by.x = "to", by.y = "taxa", sort=FALSE)
pairs_d2$dist <- with(pairs_d2, sqrt((tl.x - tl.y)**2))


## get percent explained variation of the first two axis in the coinertia
tt2 <- sapply(pairs_u2, function(x) get_percent_explained(x))

# RV value between each comms
rv_all2 <- sapply(pairs_u2, function(x) get_coia_rv(x, nrepet = 9999), simplify = FALSE) # this can take time to run depending on the number of repetitions
# some data wraggling to put te results in a plotable format
rv_all2 <- rbind.fill(rv_all2)
rv_all2$pair <- pairs_u2
rv_all2 <- merge(rv_all2, pairs_d2, by.x = "pair", by.y = "pairs")
rv_all2 <- arrange(rv_all2, dist, pair)

# plot
rv_all2 <- unite(rv_all2, "pair", c(to, from), sep =":")
rv_all2$pair <- factor(rv_all2$pair)
rv_all2$x <- 1:10

# the plot (figure X in the manuscript)
gg_rv2 <- ggplot(rv_all2, aes(x=x, y=rv, fill=col)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_brewer(type = "qual",palette = "Set1",
                    labels = c("no sig. relation", "sig. relation"),
                    name = "RV permutation\ntest:") +
  scale_x_continuous(breaks = 1:10, labels = rv_all2$pair, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.6), expand = c(0, 0)) +
  coord_flip() +
  theme(axis.text.x = element_text(hjust = 0),
        axis.text.y = element_text(vjust = 0.25)) +
  labs(x = "Pairwise community comparison",
       y = "RV value") +
  theme_bw()

ggsave("figures/coin_rv_feeding.png", gg_rv2)


# multiple COIA

## first web of interaction: all groups except for the bats:
web_1 <- "primary_herbivore_detritivore_omnivore_predator"
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
## the plot of composition effects
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
  ylim(c(1.6, 6.2))

gg_rob <- ggplot(subset(vhb_dd, speccomb %in% c("qrob", "fsyl_qrob", "qrob_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "qrob"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "qrob"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("ped. oak", "beech +\nped. oak", "ped. oak +\nred oak", "all")) +
  labs(x = "", y = "", title = "(b)") +
  ylim(c(1.6, 6.2))

gg_rub <- ggplot(subset(vhb_dd, speccomb %in% c("qrub", "fsyl_qrub", "qrob_qrub", "all"))) +
  geom_jitter(aes(x=speccomb, y = d), width = 0.1) +
  geom_point(data=subset(newdat2, tree == "qrub"), aes(x=speccomb, y=pred), color="red", size = 2.5) +
  geom_linerange(data = subset(newdat2, tree == "qrub"), aes(x=speccomb, ymin = pred-2*se, ymax=pred+2*se), color="red") +
  geom_hline(yintercept = exp(coef(m_vhb)[1]), linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "(c)") +
  scale_x_discrete(labels = c("red oak", "beech +\nred oak", "ped. oak +\nred oak", "all")) +
  ylim(c(1.6, 6.2))


gg_all <- grid.arrange(gg_syl, gg_rob, gg_rub, nrow = 3, bottom = "Tree species composition", left = "Community distance")

ggsave("figures/mcoa_composition_feeding.png", gg_all, width = 4, height = 8)

##### Part 5: Fourth corner analysis #######

# explore trait-environment relation via a fourth corner analysis
# load traits data
ff <- list.files("data/", pattern = "_trait")
traits_l <- lapply(ff, function(x) read.csv(paste0("data/", x)))
names(traits_l) <- gsub("_traits.csv", "", ff)

# some refinments, per taxa especiallymaking sure that species order is the
# same between the trait table and the community matrix
traits_l$bat$body_size <- with(traits_l$bat, (body_size_min + body_size_max) / 2)
traits_l$bat$forest_affinity <- factor(traits_l$bat$forest_affinity)
traits_l$bat <- traits_l$bat[,c(1, 5, 4)]
traits_l$bird$forest_affinity <- factor(traits_l$bird$forest_affinity)
traits_l$bird <- traits_l$bird[traits_l$bird$species %in%colnames(comm_l$bird),]
traits_l$carabid <- traits_l$carabid[order(traits_l$carabid$species),]
traits_l$carabid$dispersal <- factor(traits_l$carabid$dispersal)
comm_l$carabid <- comm_l$carabid[,c(1, order(colnames(comm_l$carabid[,-1])) + 1)]
comm_l$spider <- comm_l$spider[,c(1, order(colnames(comm_l$spider[,-1])) + 1)]
traits_l$diplopod <- traits_l$diplopod[order(traits_l$diplopod$species),]
traits_l$diplopod$forest_affinity <- factor(traits_l$diplopod$forest_affinity)
comm_l$diplopod <- comm_l$diplopod[,c(1, order(colnames(comm_l$diplopod[,-1])) + 1)]
traits_l$vegetation <- traits_l$vegetation[which(!is.na(traits_l$vegetation$height) & traits_l$vegetation$forest_affinity != ""),]
for(i in 2:ncol(comm_l$vegetation)){
  colnames(comm_l$vegetation)[i] <- as.character(tt_d$species_scientific[colnames(comm_l$vegetation)[i] == tt_d$n])
}
comm_l$vegetation <- comm_l$vegetation[,c(1, order(colnames(comm_l$vegetation[,-1])) + 1)]
comm_l$vegetation <- comm_l$vegetation[,c(1, which(colnames(comm_l$vegetation) %in% traits_l$vegetation$species))]
names(comm_l)[6] <- "isopod"

## an helper function
do_fourthcorner <- function(comm, ...){
  spe <- comm_l[[comm]][,-1]
  
  trait <- traits_l[[comm]][,2:3]
  
  envs <- plot_info[,c("speccomb", "fragm_std")]
  envs$speccomb <- factor(envs$speccomb)
  
  four <- fourthcorner(envs, spe, trait, ...)
  
  cat(paste0(comm, " is done!\n"))
  
  return(four)
}


fourth_l <- lapply(names(traits_l), function(comm) do_fourthcorner(comm,
                                                                   modeltype = 6,
                                                                   nrepet = 9999))
png("figures/fourth_corner.png", width = 800, height = 800)
par(mfrow = c(3, 3))
for(i in 1:8){
  plot(fourth_l[[i]])
  mtext(names(traits_l)[i], line = 2, adj = 0)
}
dev.off()

############ Part 6: compute sample coverage and diversity metrics

# compute sample coverage
library(iNEXT)

get_coverage <- function(comm){
  
  cc <- comm_l[[comm]]
  # remove plot id
  cc <- cc[,-1]
  # remove plots with no species
  cc <- cc[rowSums(cc) > 0,]
  # remove plots with only one species
  cc <- cc[apply(cc, 1, function(x) sum(x > 0)) != 1,]
  
  # get iNext
  ii <- iNEXT(t(cc), q = 0, datatype = "abundance", knots = 2)
  
  return(mean(ii$DataInfo$SC))
  
}

lapply(comms[-c(2, 9)], function(x) get_coverage(x))


# get plot of tightness vs diversity
library(tidyr)
get_div <- function(x){
  
  rich <- sum(x > 0, na.rm = TRUE)
  p <- x / sum(x, na.rm = TRUE)
  sha <- exp(- sum(p * log(p), na.rm = TRUE))
  simp <- 1 / sum(p ** 2)
  
  out <- data.frame(rich = rich, sha = sha, simp = simp)
  
  return(out)
}


dd <- ldply(comm_l, function(x) {
  out <- adply(x, 1, get_div, .expand = FALSE)
})

dd %>%
  pivot_longer(3:5, names_to = "index") %>%
  rename(id_plot = X1) %>%
  left_join(vhb_dd[,1:2]) -> dd_d

ggplot(dd_d, aes(x=d, y=value)) +
  stat_smooth(method = "lm", se = FALSE) +
  geom_point() +
  facet_grid(.id ~ index) 

dd_d %>%
  group_by(.id, index) %>%
  summarise(R = cor(value, d)) %>%
  arrange(desc(R)) -> cc


