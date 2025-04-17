library(did)
library(did2s)
library(tidyverse)
library(ggpubr)

# this reads in a suite of functions this script uses predominately to ease the
# the repetive nature of repeating tasks and model fits to different buffer sizes.
source("functions.R")

## script aims -----------------------------------------------------------------

# This script fits two heterogeniety robust DiD models specifically the the
# Callaway and Sant'Anna heterogenirty robust model for dynamic treatment periods
# taken from from Callaway, Brantly and Pedro H.C. Sant'Anna."Difference-in-Differences 
# with Multiple Time Periods." Journal of Econometrics, Vol. 225, No. 2, pp. 200-230, 2021. 
# <https://doi.org/10.1016/j.jeconom.2020.12.001> and a second model from Gardner
#(2022) https://arxiv.org/pdf/2207.05943 that fits a 2-stage DiD model to isolate 
# dynamic treatment effects.

# These models are fit to all "forest" mine clusters classed as having at least
# a third of the pixels in a XX buffer as closed canopy moist tropical forest
# as per the TMF data. Additionally this main analysis assumes you need a cluster
# of at least 5 pixels to be classed as a mine. This is partly to reduce the risk
# of commission errors where 1 or 2 pixels might be erroneously classed as mines.

## Forest loss: Fit dynamic DiD for sequential buffers--------------------------

# read in the 5-pixel mining clusters  30% TMF
load("Data/raw/forest.loss.df/masolele.mines.5pixels.third.forest.1km.buffer.forest.loss.df.RData")
km1 <- mine.buffer.forest.loss.df
load("Data/raw/forest.loss.df/masolele.mines.5pixels.third.forest.2.5km.buffer.forest.loss.df.RData")
km25 <- mine.buffer.forest.loss.df
load("Data/raw/forest.loss.df/masolele.mines.5pixels.third.forest.5km.buffer.forest.loss.df.RData")
km5 <- mine.buffer.forest.loss.df
load("Data/raw/forest.loss.df/masolele.mines.5pixels.third.forest.10km.buffer.forest.loss.df.RData")
km10 <- mine.buffer.forest.loss.df
load("Data/raw/covariates/masolele.mines.final.5pixels.third.forest.final.COVARIATES.RData")
covariates.5pix <- mine.cluster.points.df

mine.forest.loss.5pix.combi <- list("km1" = km1, "km2" = km25, "km5" = km5, "km10" = km10)
length(unique(mine.forest.loss.5pix.combi[[1]]$CLUSTER)) ## 446

# create relative time to mine covariates and add the covariate data.
mine.forest.loss.5pix.combi.did <- lapply(mine.forest.loss.5pix.combi,
                                          did.prep, lead.time = -20, post.time = 20,
                                          type = "loss", covariates = covariates.5pix)

# test.1km <- mine.forest.loss.dist.combi.did[[1]]
# ## visualisation
# ggplot(mine.forest.loss.5pix.combi.did[[1]], aes(rel.year.first, cumulative.forest.loss.prop)) +
#   geom_point() +
#   geom_smooth()
# 
# ggplot(mine.forest.loss.combi.did[[1]],
#        aes(first.loss.year)) +
#   geom_histogram()

# fit Gardner and Callaway and Sant'Anna DiD model
forest.loss.5pix.did <- lapply(mine.forest.loss.5pix.combi.did, fit.dynamic.DiD,
                               yname = "cumulative.forest.loss.perc",
                               xformula = NULL, method = c("csa", "gardner"))

forest.loss.5pix.plts.f <- lapply(forest.loss.5pix.did, plot.DiD2S.f, 
                                  pre.yrs = -5, post.yrs = 10, legend = "none",
                                  method = c("csa", "gardner"),
                                  pre.mine.col = c("#313695", "#abd9e9"),
                                  post.mine.col = c("#fdae61", "#d73027"),
                                  y.label = "Forest loss (% points)")

forest.loss.5pix.plts.s <- lapply(forest.loss.5pix.did, plot.DiD2S.s,
                                  pre.yrs = -5, post.yrs = 10,
                                  method = c("csa", "gardner"),
                                  #pre.col = "#1b7837", post.col = "#762a83",
                                  pre.col = "#4393c3", post.col = "#8c510a",
                                  y.label = "Forest loss (% points)")


bind_rows(forest.loss.5pix.did) %>% 
  write.csv("Outputs/Tables/forest.loss.thirdforest.5pix.mod.results.csv")

## Forest loss: Fit dynamic DiD - for overall 0-5km buffer ---------------------

load("Data/raw/forest.loss.df/masolele.mines.5pixels.third.forest.5km.master.buffer.forest.loss.df.RData")
master.forest.loss.5pix.5km <- mine.buffer.forest.loss.df
mean(master.forest.loss.5pix.5km$cumulative.forest.loss.prop)

# create relative time to mine covariates and add the covariate data.
mine.forest.loss.5pix.5km.master.did <- did.prep(master.forest.loss.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "loss", covariates = covariates.5pix)

ggplot(mine.forest.loss.5pix.5km.master.did, aes(rel.year.first, cumulative.forest.loss.prop,
                                                 colour = treatment)) +
  geom_point() +
  geom_smooth()

# fit Gardner and Callaway and Sant'Anna DiD model
forest.loss.master.5pix.5km.did <- fit.dynamic.DiD(mine.forest.loss.5pix.5km.master.did,
                                                   yname = "cumulative.forest.loss.perc",
                                                   xformula = NULL, method = c("csa", "gardner"))

forest.loss.master.5pix.plts.f <- plot.DiD2S.f(forest.loss.master.5pix.5km.did, 
                                               pre.yrs = -5, post.yrs = 10, legend = "none",
                                               method = c("csa", "gardner"),
                                               pre.mine.col = c("#313695", "#abd9e9"),
                                               post.mine.col = c("#fdae61", "#d73027"),
                                               y.label = "Forest loss (% points)")

forest.loss.master.5pix.plts.s <- plot.DiD2S.s(forest.loss.master.5pix.5km.did,
                                               pre.yrs = -5, post.yrs = 10,
                                               method = c("csa", "gardner"),
                                               #pre.col = "#1b7837", post.col = "#762a83",
                                               pre.col = "#4393c3", post.col = "#8c510a",
                                               y.label = "Forest loss (% points)")

bind_rows(forest.loss.master.5pix.5km.did) %>% 
  write.csv("Outputs/Tables/forest.loss.master.thirdforest.5pix.mod.results.csv")

## Arrange forest loss figure ----------------------------------------------

forest.loss.DiD.5pix.main <- 
  ggarrange(forest.loss.master.5pix.plts.s$`Gardner 2022`,
            ggarrange(forest.loss.5pix.plts.s$km1$`Gardner 2022`,
                      forest.loss.5pix.plts.s$km2$`Gardner 2022` +
                        ylab(""),
                      forest.loss.5pix.plts.s$km5$`Gardner 2022` +
                        ylab(""),
                      forest.loss.5pix.plts.s$km10$`Gardner 2022` +
                        ylab(""),
                      ncol = 4, nrow = 1, labels = c("b", "c", "d", "e"),
                      hjust = 0, font.label = list(size = 9)),
            nrow = 2, heights = c(1, .5),
            labels = c("a", ""), hjust = 0, font.label = list(size = 9))

ggsave(path = "Outputs/Figures/Final/Main", filename = "forest.loss.DiD.5pix.main.png",
       forest.loss.DiD.5pix.main, bg = "white",
       device = "png", width = 17, height = 13, units = "cm")          

## Degradation: Fit dynamic DiD for sequential buffers -------------------------
load("Data/raw/degradation.df/masolele.mines.5pixels.third.forest.1km.buffer.degradation.df.RData")
km1 <- mine.buffer.forest.loss.df
load("Data/raw/degradation.df/masolele.mines.5pixels.third.forest.2.5km.buffer.degradation.df.RData")
km25 <- mine.buffer.forest.loss.df
load("Data/raw/degradation.df/masolele.mines.5pixels.third.forest.5km.buffer.degradation.df.RData")
km5 <- mine.buffer.forest.loss.df
load("Data/raw/degradation.df/masolele.mines.5pixels.third.forest.10km.buffer.degradation.df.RData")
km10 <- mine.buffer.forest.loss.df
load("Data/raw/covariates/masolele.mines.final.5pixels.third.forest.final.COVARIATES.RData")
covariates.5pix <- mine.cluster.points.df

mine.degradation.5pix.combi <- list("km1" = km1, "km2" = km25, "km5" = km5, "km10" = km10)
length(unique(mine.degradation.5pix.combi[[1]]$CLUSTER)) ## 446 clusters

# create relative time to mine covariates and add the covariate data.
mine.degradation.5pix.combi.did <- lapply(mine.degradation.5pix.combi,
                                          did.prep, lead.time = -20, post.time = 20,
                                          type = "degradation", covariates = covariates.5pix)

# fit Gardner and Callaway and Sant'Anna DiD model
degradation.5pix.did <- lapply(mine.degradation.5pix.combi.did, fit.dynamic.DiD,
                               yname = "cumulative.degradation.perc",
                               xformula = NULL, method = c("csa", "gardner"))

degradation.5pix.plts.s <- lapply(degradation.5pix.did, plot.DiD2S.s,
                                  pre.yrs = -5, post.yrs = 10,
                                  method = c("csa", "gardner"),
                                  #pre.col = "#1b7837", post.col = "#762a83",
                                  pre.col = "#4393c3", post.col = "#e08214",
                                  y.label = "Degradation (% points)")

degradation.5pix.plts.f <- lapply(degradation.5pix.did, plot.DiD2S.f, 
                                  pre.yrs = -5, post.yrs = 10, legend = "none",
                                  method = c("csa", "gardner"),
                                  pre.mine.col = c("#313695", "#abd9e9"),
                                  post.mine.col = c("#fdae61", "#d73027"),
                                  y.label = "Degradation (% points)")

bind_rows(degradation.5pix.did) %>% 
  write.csv("Outputs/Tables/degradation.thirdforest.5pix.mod.results.csv")

## Degradation: Fit dynamic DiD - for overall 0-5km buffer----------------------
load("Data/raw/degradation.df/masolele.mines.5pixels.third.forest.5km.master.buffer.degradation.df.RData")
master.degradation.loss.5pix.5km <- mine.buffer.forest.loss.df
mean(master.degradation.loss.5pix.5km$cumulative.degradation.prop)

# create relative time to mine covariates and add the covariate data.
mine.degradation.5pix.5km.master.did <- did.prep(master.degradation.loss.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "degradation", covariates = covariates.5pix)

ggplot(mine.degradation.5pix.5km.master.did, aes(rel.year.first, cumulative.degradation.prop,
                                                 colour = treatment)) +
  geom_point() +
  geom_smooth()

# fit Gardner and Callaway and Sant'Anna DiD model
degradation.master.5pix.5km.did <- fit.dynamic.DiD(mine.degradation.5pix.5km.master.did,
                                                   yname = "cumulative.degradation.perc",
                                                   xformula = NULL, method = c("csa", "gardner"))

degradation.master.5pix.plts.f <- plot.DiD2S.f(degradation.master.5pix.5km.did, 
                                               pre.yrs = -5, post.yrs = 10, legend = "none",
                                               method = c("csa", "gardner"),
                                               pre.mine.col = c("#313695", "#abd9e9"),
                                               post.mine.col = c("#fdae61", "#d73027"),
                                               y.label = "Degradation (% points)")

degradation.master.5pix.plts.s <- plot.DiD2S.s(degradation.master.5pix.5km.did,
                                               pre.yrs = -5, post.yrs = 10,
                                               method = c("csa", "gardner"),
                                               #pre.col = "#1b7837", post.col = "#762a83",
                                               pre.col = "#4393c3", post.col = "#e08214",
                                               y.label = "Degradation (% points)")

bind_rows(degradation.master.5pix.5km.did) %>% 
  write.csv("Outputs/Tables/degradation.master.thirdforest.5pix.mod.results.csv")

## Arrange degradation figure ----------------------------------------------

degradation.DiD.5pix.main <- 
  ggarrange(degradation.master.5pix.plts.s$`Gardner 2022`,
            ggarrange(degradation.5pix.plts.s$km1$`Gardner 2022`,
                      degradation.5pix.plts.s$km2$`Gardner 2022` +
                        ylab(""),
                      degradation.5pix.plts.s$km5$`Gardner 2022` +
                        ylab(""),
                      degradation.5pix.plts.s$km10$`Gardner 2022` +
                        ylab(""),
                      ncol = 4, nrow = 1, labels = c("b", "c", "d", "e"),
                      hjust = 0, font.label = list(size = 9)),
            nrow = 2, heights = c(1, .5),
            labels = c("a", ""), hjust = 0, font.label = list(size = 9))

ggsave(path = "Outputs/Figures/Final/Main", filename = "degradation.DiD.5pix.main.png",
       degradation.DiD.5pix.main, bg = "white",
       device = "png", width = 17, height = 13, units = "cm")          

