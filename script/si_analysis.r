#######################
# R script for the Treeweb       #
# cross-community analysis       #
# SI analysis                    #
#######################################

# load libraries
library(vegan)


# set working directory
setwd("~/PostDoc_Ghent/Synthesis_3/treeweb_community/")

# load plot data
plot_info <- read.csv("data/synthesis_expldata_raw.csv", sep = " ", stringsAsFactors = FALSE)
plot_xy <- read.csv("data/plot_xy.csv")

# load community data
comms <- c("bird", "bat", "spider", "opiliones", "carabid", "isopods", "diplopod", "herb", "vegetation")
comm_l <- lapply(comms, function(x) read.csv(paste0("data/treeweb_", x, "_format.csv")))
names(comm_l) <- comms

# remove the species with no record in the bird table
comm_l$bird <- comm_l$bird[,-which(colSums(comm_l$bird) == 0)]
# divide herbivore by sampling effort
comm_l$herb <- cbind(comm_l$herb[,1], comm_l$herb[,-1] / plot_info$specrich)

########### Part 1: RDA per trophic group

# an helper function to run the RDA
# arguments are:
# @comm: a character, the name of the taxa
# @interaction: a logical, whether the include an interaction term between species composition and fragmentation
# @std: a logical, whether to hellinger-standardize the community data
# @species: a character, the tree species in focus
# @...: further argument to pass to the RDA, see ?rda

fit_rda <- function(comm, interaction = TRUE, std = TRUE, species = "qrob",...){
  
  print(species)
  ids <- eval(as.name(paste0("id_",species)))
  
  if(std){
    comm <- decostand(comm, method = "hellinger", MARGIN = 1)
  }
  
  # compute PCNM
  rs <- rowSums(comm) / sum(comm)
  pcnmw <- pcnm(dist(plot_xy[,c("X", "Y")]), w = rs)
  
  dat <- cbind(plot_info[,c("speccomb", "fragm_std")], scores(pcnmw))
  # NA for opiliones, fix by putting PCNM to 0
  dat[is.na(dat)] <- 0
  
  if(interaction){
    rr <- rda(comm[ids,] ~ . + speccomb : fragm_std, dat[ids,], ...)
  } else {
    rr <- rda(comm[ids,] ~ ., dat[ids,],...)
  }
  
  return(rr)
}