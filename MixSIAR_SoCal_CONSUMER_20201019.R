### MixSIAR Model for southern California stable isotope data
### Andrew D. Somerville
### Using r package, MixSIar 
### Stock, B. C., Jackson, A. L., Ward, E. J., Parnell, A. C., Phillips, D. L., & Semmens, B. X. (2018). Analyzing mixing systems using a new generation of Bayesian tracer mixing models. PeerJ, 6, e5096.
### Stock, B. C. and B. X. Semmens (2016). MixSIAR GUI User Manual.Version 3.1. https://github.com/brianstock/MixSIAR/. doi:10.5281/zenodo.47719.




library(MixSIAR)




# Load mix data

mix <- load_mix_data(filename="CONSUMER_noSNI.csv",
                     iso_names=c("d13C","d15N"),
                     factors="Prov",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)

# Load source data

source <- load_source_data(filename="SOURCE_20201012.csv",
                           source_factors=NULL,
                           conc_dep=TRUE,
                           data_type="means",
                           mix)

# Load discrimination/TDF data

discr <- load_discr_data(filename="TEF.csv", mix)

# Isospace plot

plot_data(filename="isospace_plot",
          plot_save_pdf=TRUE,
          plot_save_png=FALSE,
          mix,source,discr)

# Calculate standardized convex hull area
#if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr)

################################################################################
# # PRIORS (construct alpha from geographic assumptions [island, coastal, inland])
################################################################################
# Resource importance scored as high = 3, med = 2, low = 1

# Uninformed Prior
Cali.alpha <- c(1,1,1,1)

# Generate alpha hyperparameters scaling sum(alpha)=n.sources
Cali.alpha <- Cali.alpha*length(Cali.alpha)/sum(Cali.alpha)

# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)
Cali.alpha[which(Cali.alpha==0)] <- 0.01

# Plot informative prior
plot_prior(alpha.prior=Cali.alpha,
           source=source,
           plot_save_pdf=FALSE,
           plot_save_png=FALSE,
           filename="prior_plot_Cali_inf")


# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- FALSE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("short" first, then "short")
jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior=1)
#jags.1 <- run_model(run="short", mix, source, discr, model_filename, alpha.prior=1)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.1, mix, source)




###################################################

## OUTPUT PLOTTING ##



library(MASS)
library(R2jags)
library(tidyverse)
library(RColorBrewer)
library(lattice)
library(dplyr)
library(grid)

attach.jags(jags.1)
post.Coastal <- data.frame(Prov = "2 Coastal", Marine.High = p.fac1[,2,1], Marine.Low = p.fac1[,2,2], Plants = p.fac1[,2,3], T.Mammal = p.fac1[,2,4])
post.Inland <- data.frame(Prov = "1 Inland",Marine.High = p.fac1[,1,1], Marine.Low = p.fac1[,1,2], Plants = p.fac1[,1,3], T.Mammal = p.fac1[,1,4])
post.San.Clemente <- data.frame(Prov = "5 San Clemente", Marine.High = p.fac1[,5,1], Marine.Low = p.fac1[,5,2], Plants = p.fac1[,5,3], T.Mammal = p.fac1[,5,4])
post.Santa.Cruz <- data.frame(Prov = "4 Santa Cruz", Marine.High = p.fac1[,4,1], Marine.Low = p.fac1[,4,2], Plants = p.fac1[,4,3], T.Mammal = p.fac1[,4,4])
post.Santa.Rosa <- data.frame(Prov = "3 Santa Rosa", Marine.High = p.fac1[,3,1], Marine.Low = p.fac1[,3,2], Plants = p.fac1[,3,3], T.Mammal = p.fac1[,3,4])


Coastal <- post.Coastal %>% gather(source,value, 2:5)
Inland <-post.Inland %>% gather(source,value, 2:5)
San_Clemente <- post.San.Clemente %>% gather(source,value, 2:5)
Santa_Cruz <- post.Santa.Cruz %>% gather(source,value, 2:5)
Santa_Rosa <- post.Santa.Rosa %>% gather(source,value, 2:5)

all <- rbind(Coastal, Inland, San_Clemente, Santa_Cruz, Santa_Rosa)



########### BOX PLOTS ###########

ggplot(aes(y = value, x = source, fill=source), data = all) + 
  geom_boxplot(outlier.colour = NA) +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_fill_manual(values=c("#225ea8", "#41b6c4", "#7fcdbb", "#ffffd9"))+ 
  xlab("Source") + 
  ylab("Diet proportion") +
  facet_wrap(~Prov, ncol=2) +
  theme(legend.position=c(.80, .15))+
  ggsave(filename = "Cali_BoxPlots_CONSUMER_noSNI_facet.tiff", width=8, height=11)


########## OVERALL VIOLIN PLOTS ###########

ggplot(aes(x = value, y=source, fill=source),data=all)+ 
  geom_violin(trim = FALSE,
        alpha = 0.8)+
  geom_boxplot(alpha=0.3, color="black", width=.1)+
    theme_bw() +
  theme(legend.position="none", )+
  scale_fill_manual(values=c("#225ea8", "#41b6c4", "#7fcdbb", "#ffffd9"))+ 
  facet_wrap(~Prov, ncol=2, scales="free") +
  coord_flip()+
  xlab("Diet Proportions") + 
  ylab("Source") +
  ggsave(filename = "Cali_BoxPlots_CONSUMER_noSNI_violin.tiff", width=8, height=11)


