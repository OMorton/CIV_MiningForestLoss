library(did)
library(did2s)
library(tidyverse)
library(ggpubr)

# this reads in a suite of functions this script uses predominately to ease the
# the repetive nature of repeating tasks and model fits to different buffer sizes.
source("functions.R")

## SM analysis 1 - refit models with Callaway and Sant'Anna estimator ----------
# on the third forest data at 5 pixel level clusters

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


# fit Callaway and Sant'Anna DiD model
forest.loss.5pix.did <- lapply(mine.forest.loss.5pix.combi.did, fit.dynamic.DiD,
                               yname = "cumulative.forest.loss.perc",
                               xformula = NULL, method = c("csa"))

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

## master 0-5 km buffer
load("Data/raw/forest.loss.df/masolele.mines.5pixels.third.forest.5km.master.buffer.forest.loss.df.RData")
master.forest.loss.5pix.5km <- mine.buffer.forest.loss.df
mean(master.forest.loss.5pix.5km$cumulative.forest.loss.prop)

# create relative time to mine covariates and add the covariate data.
mine.forest.loss.5pix.5km.master.did <- did.prep(master.forest.loss.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "loss", covariates = covariates.5pix)


# fit Callaway and Sant'Anna DiD model
forest.loss.master.5pix.5km.did <- fit.dynamic.DiD(mine.forest.loss.5pix.5km.master.did,
                                                   yname = "cumulative.forest.loss.perc",
                                                   xformula = NULL, method = c("csa"))

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

forest.loss.DiD.5pix.csa <- 
  ggarrange(forest.loss.master.5pix.plts.s$`Callaway & Sant'Anna 2021`,
            ggarrange(forest.loss.5pix.plts.s$km1$`Callaway & Sant'Anna 2021`,
                      forest.loss.5pix.plts.s$km2$`Callaway & Sant'Anna 2021` +
                        ylab(""),
                      forest.loss.5pix.plts.s$km5$`Callaway & Sant'Anna 2021` +
                        ylab(""),
                      forest.loss.5pix.plts.s$km10$`Callaway & Sant'Anna 2021` +
                        ylab(""),
                      ncol = 4, nrow = 1, labels = c("b", "c", "d", "e"),
                      hjust = 0, font.label = list(size = 9)),
            nrow = 2, heights = c(1, .5),
            labels = c("a", ""), hjust = 0, font.label = list(size = 9))

ggsave(path = "Outputs/Figures/Final/SM", filename = "forest.loss.DiD.5pix.csa.png",
       forest.loss.DiD.5pix.csa, bg = "white",
       device = "png", width = 17, height = 13, units = "cm")   

## SM analysis 2 - refit degr models with Callaway and Sant'Anna estimator -----
# on the third forest data at 5 pixel level clusters
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

# fit Callaway and Sant'Anna DiD model
degradation.5pix.did <- lapply(mine.degradation.5pix.combi.did, fit.dynamic.DiD,
                               yname = "cumulative.degradation.perc",
                               xformula = NULL, method = c("csa"))

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

## master 0-5 km
load("Data/raw/degradation.df/masolele.mines.5pixels.third.forest.5km.master.buffer.degradation.df.RData")
master.degradation.loss.5pix.5km <- mine.buffer.forest.loss.df

# create relative time to mine covariates and add the covariate data.
mine.degradation.5pix.5km.master.did <- did.prep(master.degradation.loss.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "degradation", covariates = covariates.5pix)


# fit Callaway and Sant'Anna DiD model
degradation.master.5pix.5km.did <- fit.dynamic.DiD(mine.degradation.5pix.5km.master.did,
                                                   yname = "cumulative.degradation.perc",
                                                   xformula = NULL, method = c("csa"))

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


degradation.DiD.5pix.csa <- 
  ggarrange(degradation.master.5pix.plts.s$`Callaway & Sant'Anna 2021`,
            ggarrange(degradation.5pix.plts.s$km1$`Callaway & Sant'Anna 2021`,
                      degradation.5pix.plts.s$km2$`Callaway & Sant'Anna 2021` +
                        ylab(""),
                      degradation.5pix.plts.s$km5$`Callaway & Sant'Anna 2021` +
                        ylab(""),
                      degradation.5pix.plts.s$km10$`Callaway & Sant'Anna 2021` +
                        ylab(""),
                      ncol = 4, nrow = 1, labels = c("b", "c", "d", "e"),
                      hjust = 0, font.label = list(size = 9)),
            nrow = 2, heights = c(1, .5),
            labels = c("a", ""), hjust = 0, font.label = list(size = 9))

ggsave(path = "Outputs/Figures/Final/SM", filename = "degradation.DiD.5pix.cas.png",
       degradation.DiD.5pix.csa, bg = "white",
       device = "png", width = 17, height = 13, units = "cm")          


## SM analysis 3 - refit forest loss models including covariates ----------------
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

forest.loss.5pix.did.COVARS <- lapply(mine.forest.loss.5pix.combi.did, fit.dynamic.DiD,
                                      yname = "cumulative.forest.loss.perc",
                                      xformula = ~travel.time + pop.density + road.distance +
                                        river.distance + elevation + slope, 
                                      method = c("gardner"))

forest.loss.5pix.covars.plts.f <- lapply(forest.loss.5pix.did.COVARS, plot.DiD2S.f, 
                                         pre.yrs = -5, post.yrs = 8, legend = "none",
                                         method = c("gardner"),
                                         pre.mine.col = c("#313695", "#abd9e9"),
                                         post.mine.col = c("#fdae61", "#d73027"),
                                         y.label = "Forest loss (% points)")

forest.loss.5pix.covars.plts.s <- lapply(forest.loss.5pix.did.COVARS, plot.DiD2S.s,
                                         pre.yrs = -5, post.yrs = 8,
                                         method = c("gardner"),
                                         #pre.col = "#1b7837", post.col = "#762a83",
                                         pre.col = "#4393c3", post.col = "#8c510a",
                                         y.label = "Forest loss (% points)")

# master 0 - 5km
load("Data/raw/forest.loss.df/masolele.mines.5pixels.half.forest.5km.master.buffer.forest.loss.df.RData")
master.forest.loss.5pix.5km <- mine.buffer.forest.loss.df

mine.forest.loss.5pix.5km.master.did <- did.prep(master.forest.loss.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "loss", covariates = covariates.5pix)


forest.loss.master.5pix.did.COVARS <- fit.dynamic.DiD(mine.forest.loss.5pix.5km.master.did,
                                      yname = "cumulative.forest.loss.perc",
                                      xformula = ~travel.time + pop.density + road.distance +
                                        river.distance + elevation + slope, 
                                      method = c("gardner"))

forest.loss.master.5pix.covars.plts.f <- plot.DiD2S.f(forest.loss.master.5pix.did.COVARS, 
                                         pre.yrs = -5, post.yrs = 8, legend = "none",
                                         method = c("gardner"),
                                         pre.mine.col = c("#313695", "#abd9e9"),
                                         post.mine.col = c("#fdae61", "#d73027"),
                                         y.label = "Forest loss (% points)")

forest.loss.master.5pix.covars.plts.s <- plot.DiD2S.s(forest.loss.master.5pix.did.COVARS,
                                         pre.yrs = -5, post.yrs = 8,
                                         method = c("gardner"),
                                        # pre.col = "#1b7837", post.col = "#762a83",
                                        pre.col = "#4393c3", post.col = "#8c510a",
                                         y.label = "Forest loss (% points)")

bind_rows(forest.loss.5pix.did.COVARS) %>% 
  write.csv("Outputs/Tables/forest.loss.thirdforest.5pix.COVARS.mod.results.csv")
bind_rows(forest.loss.master.5pix.did.COVARS) %>% 
  write.csv("Outputs/Tables/forest.loss.master.thirdforest.5pix.COVARS.mod.results.csv")

forest.loss.DiD.5pix.covars.SM <- 
  ggarrange(forest.loss.master.5pix.covars.plts.s$`Gardner 2022`,
            ggarrange(forest.loss.5pix.covars.plts.s$km1$`Gardner 2022`,
                      forest.loss.5pix.covars.plts.s$km2$`Gardner 2022` +
                        ylab(""),
                      forest.loss.5pix.covars.plts.s$km5$`Gardner 2022` +
                        ylab(""),
                      forest.loss.5pix.covars.plts.s$km10$`Gardner 2022` +
                        ylab(""),
                      ncol = 4, nrow = 1, labels = c("b", "c", "d", "e"),
                      hjust = 0, font.label = list(size = 9)),
            nrow = 2, heights = c(1, .5),
            labels = c("a", ""), hjust = 0, font.label = list(size = 9))

ggsave(path = "Outputs/Figures/Final/SM", filename = "forest.loss.DiD.5pix.covars.SM.png",
       forest.loss.DiD.5pix.covars.SM, bg = "white",
       device = "png", width = 17, height = 13, units = "cm") 

## SM analysis 4 - refit degr models including covariates ----------------
# read in the 5-pixel mining clusters  30% TMF
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
length(unique(mine.degradation.5pix.combi[[1]]$CLUSTER)) ## 446

# create relative time to mine covariates and add the covariate data.
mine.degradation.5pix.combi.did <- lapply(mine.degradation.5pix.combi,
                                          did.prep, lead.time = -20, post.time = 20,
                                          type = "degradation", covariates = covariates.5pix)

degradation.5pix.did.COVARS <- lapply(mine.degradation.5pix.combi.did, fit.dynamic.DiD,
                                      yname = "cumulative.degradation.perc",
                                      xformula = ~travel.time + pop.density + road.distance +
                                        river.distance + elevation + slope, 
                                      method = c("gardner"))

degradation.5pix.covars.plts.f <- lapply(degradation.5pix.did.COVARS, plot.DiD2S.f, 
                                         pre.yrs = -5, post.yrs = 8, legend = "none",
                                         method = c("gardner"),
                                         pre.mine.col = c("#313695", "#abd9e9"),
                                         post.mine.col = c("#fdae61", "#d73027"),
                                         y.label = "Degradation (% points)")

degradation.5pix.covars.plts.s <- lapply(degradation.5pix.did.COVARS, plot.DiD2S.s,
                                         pre.yrs = -5, post.yrs = 8,
                                         method = c("gardner"),
                                         #pre.col = "#1b7837", post.col = "#762a83",
                                         pre.col = "#4393c3", post.col = "#e08214",
                                         y.label = "Degradation (% points)")

# master 0 - 5km
load("Data/raw/degradation.df/masolele.mines.5pixels.half.forest.5km.master.buffer.degradation.df.RData")
master.degradation.5pix.5km <- mine.buffer.forest.loss.df

mine.degradation.5pix.5km.master.did <- did.prep(master.degradation.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "degradation", covariates = covariates.5pix)


degradation.master.5pix.did.COVARS <- fit.dynamic.DiD(mine.degradation.5pix.5km.master.did,
                                                      yname = "cumulative.degradation.perc",
                                                      xformula = ~travel.time + pop.density + road.distance +
                                                        river.distance + elevation + slope, 
                                                      method = c("gardner"))

degradation.master.5pix.covars.plts.f <- plot.DiD2S.f(degradation.master.5pix.did.COVARS, 
                                                      pre.yrs = -5, post.yrs = 8, legend = "none",
                                                      method = c("gardner"),
                                                      pre.mine.col = c("#313695", "#abd9e9"),
                                                      post.mine.col = c("#fdae61", "#d73027"),
                                                      y.label = "Degradation (% points)")

degradation.master.5pix.covars.plts.s <- plot.DiD2S.s(degradation.master.5pix.did.COVARS,
                                                      pre.yrs = -5, post.yrs = 8,
                                                      method = c("gardner"),
                                                      # pre.col = "#1b7837", post.col = "#762a83",
                                                      pre.col = "#4393c3", post.col = "#e08214",
                                                      y.label = "Degradation (% points)")

bind_rows(degradation.5pix.did.COVARS) %>% 
  write.csv("Outputs/Tables/degradation.thirdforest.5pix.COVARS.mod.results.csv")
bind_rows(degradation.master.5pix.did.COVARS) %>% 
  write.csv("Outputs/Tables/degradation.master.thirdforest.5pix.COVARS.mod.results.csv")

degradation.DiD.5pix.covars.SM <- 
  ggarrange(degradation.master.5pix.covars.plts.s$`Gardner 2022`,
            ggarrange(degradation.5pix.covars.plts.s$km1$`Gardner 2022`,
                      degradation.5pix.covars.plts.s$km2$`Gardner 2022` +
                        ylab(""),
                      degradation.5pix.covars.plts.s$km5$`Gardner 2022` +
                        ylab(""),
                      degradation.5pix.covars.plts.s$km10$`Gardner 2022` +
                        ylab(""),
                      ncol = 4, nrow = 1, labels = c("b", "c", "d", "e"),
                      hjust = 0, font.label = list(size = 9)),
            nrow = 2, heights = c(1, .5),
            labels = c("a", ""), hjust = 0, font.label = list(size = 9))

ggsave(path = "Outputs/Figures/Final/SM", filename = "degradation.DiD.5pix.covars.SM.png",
       degradation.DiD.5pix.covars.SM, bg = "white",
       device = "png", width = 17, height = 13, units = "cm") 

## SM analysis 5 - refit forest loss with 50% forest cover definition ----------

load("Data/raw/forest.loss.df/masolele.mines.5pixels.half.forest.1km.buffer.forest.loss.df.RData")
km1 <- mine.buffer.forest.loss.df
load("Data/raw/forest.loss.df/masolele.mines.5pixels.half.forest.2.5km.buffer.forest.loss.df.RData")
km25 <- mine.buffer.forest.loss.df
load("Data/raw/forest.loss.df/masolele.mines.5pixels.half.forest.5km.buffer.forest.loss.df.RData")
km5 <- mine.buffer.forest.loss.df
load("Data/raw/forest.loss.df/masolele.mines.5pixels.half.forest.10km.buffer.forest.loss.df.RData")
km10 <- mine.buffer.forest.loss.df
load("Data/raw/covariates/masolele.mines.final.5pixels.half.forest.final.COVARIATES.RData")
covariates.5pix <- mine.cluster.points.df

mine.forest.loss.5pix.combi <- list("km1" = km1, "km2" = km25, "km5" = km5, "km10" = km10)
length(unique(mine.forest.loss.5pix.combi[[1]]$CLUSTER)) ## 446


mine.forest.loss.5pix.combi.did <- lapply(mine.forest.loss.5pix.combi,
                                          did.prep, lead.time = -20, post.time = 20,
                                          type = "loss", covariates = covariates.5pix)

# Fitting 2-stage DiD 
forest.loss.5pix.half.forest.did <- lapply(mine.forest.loss.5pix.combi.did, fit.dynamic.DiD,
                                           yname = "cumulative.forest.loss.perc",
                                           xformula = NULL, method = c("csa", "gardner"))

forest.loss.5pix.plts.f <- lapply(forest.loss.5pix.half.forest.did, plot.DiD2S.f, 
                                  pre.yrs = -5, post.yrs = 10, legend = "none",
                                  method = c("csa", "gardner"),
                                  pre.mine.col = c("#313695", "#abd9e9"),
                                  post.mine.col = c("#fdae61", "#d73027"),
                                  y.label = "Forest loss (% points)")

forest.loss.5pix.plts.s <- lapply(forest.loss.5pix.half.forest.did, plot.DiD2S.s,
                                  pre.yrs = -5, post.yrs = 10,
                                  method = c("csa", "gardner"),
                                  #pre.col = "#1b7837", post.col = "#762a83",
                                  pre.col = "#4393c3", post.col = "#8c510a",
                                  y.label = "Forest loss (% points)")

# write out forest loss results 
bind_rows(forest.loss.5pix.half.forest.did) %>% 
  write.csv("Outputs/Tables/forest.loss.halfforest.5pix.mod.results.csv")

# master 0 - 5km
load("Data/raw/forest.loss.df/masolele.mines.5pixels.half.forest.5km.master.buffer.forest.loss.df.RData")
master.forest.loss.5pix.5km <- mine.buffer.forest.loss.df

mine.forest.loss.5pix.5km.master.did <- did.prep(master.forest.loss.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "loss", covariates = covariates.5pix)

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

forest.loss.DiD.5pix.half.forest.SM <- 
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

ggsave(path = "Outputs/Figures/Final/SM", filename = "forest.loss.DiD.5pix.half.forest.main.png",
       forest.loss.DiD.5pix.half.forest.SM, bg = "white",
       device = "png", width = 17, height = 13, units = "cm") 

## SM analysis 6 - refit forest degr with 50% forest cover definition ----------
load("Data/raw/degradation.df/masolele.mines.5pixels.half.forest.1km.buffer.degradation.df.RData")
km1 <- mine.buffer.forest.loss.df
load("Data/raw/degradation.df/masolele.mines.5pixels.half.forest.2.5km.buffer.degradation.df.RData")
km25 <- mine.buffer.forest.loss.df
load("Data/raw/degradation.df/masolele.mines.5pixels.half.forest.5km.buffer.degradation.df.RData")
km5 <- mine.buffer.forest.loss.df
load("Data/raw/degradation.df/masolele.mines.5pixels.half.forest.10km.buffer.degradation.df.RData")
km10 <- mine.buffer.forest.loss.df
load("Data/raw/covariates/masolele.mines.final.5pixels.half.forest.final.COVARIATES.RData")
covariates.5pix <- mine.cluster.points.df

mine.degradation.5pix.combi <- list("km1" = km1, "km2" = km25, "km5" = km5, "km10" = km10)
length(unique(mine.degradation.5pix.combi[[1]]$CLUSTER)) ## 289


mine.degradation.5pix.combi.did <- lapply(mine.degradation.5pix.combi,
                                          did.prep, lead.time = -20, post.time = 20,
                                          type = "degradation", covariates = covariates.5pix)

# Fitting 2-stage DiD 
degradation.5pix.half.forest.did <- lapply(mine.degradation.5pix.combi.did, fit.dynamic.DiD,
                                           yname = "cumulative.degradation.perc",
                                           xformula = NULL, method = c("csa", "gardner"))

degradation.5pix.plts.f <- lapply(degradation.5pix.half.forest.did, plot.DiD2S.f, 
                                  pre.yrs = -5, post.yrs = 10, legend = "none",
                                  method = c("csa", "gardner"),
                                  pre.mine.col = c("#313695", "#abd9e9"),
                                  post.mine.col = c("#fdae61", "#d73027"),
                                  y.label = "Degradation (% points)")

degradation.5pix.plts.s <- lapply(degradation.5pix.half.forest.did, plot.DiD2S.s,
                                  pre.yrs = -5, post.yrs = 10,
                                  method = c("csa", "gardner"),
                                  #pre.col = "#1b7837", post.col = "#762a83",
                                  pre.col = "#4393c3", post.col = "#e08214",
                                  y.label = "Degradation (% points)")

# write out forest loss results 
bind_rows(degradation.5pix.half.forest.did) %>% 
  write.csv("Outputs/Tables/degradation.halfforest.5pix.mod.results.csv")

# master 0 - 5km
load("Data/raw/degradation.df/masolele.mines.5pixels.half.forest.5km.master.buffer.degradation.df.RData")
master.degradation.5pix.5km <- mine.buffer.forest.loss.df

mine.degradation.5pix.5km.master.did <- did.prep(master.degradation.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "degradation", covariates = covariates.5pix)

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

degradation.DiD.5pix.half.forest.SM <- 
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

ggsave(path = "Outputs/Figures/Final/SM", filename = "degradation.DiD.5pix.half.forest.main.png",
       degradation.DiD.5pix.half.forest.SM, bg = "white",
       device = "png", width = 17, height = 13, units = "cm") 

## SM analysis 7 - Cutting early years -----------------------------------------

## first year of mining data may not necessarily reflect mines actually initiated
## that year, it is simply the first year of data collection. Thus we rerun all 
## models excluding 2001.

load("Data/raw/forest.loss.df/masolele.mines.5pixels.third.forest.5km.master.buffer.forest.loss.df.RData")
master.forest.loss.5pix.5km <- mine.buffer.forest.loss.df
load("Data/raw/covariates/masolele.mines.final.5pixels.third.forest.final.COVARIATES.RData")
covariates.5pix <- mine.cluster.points.df

# create relative time to mine covariates and add the covariate data.
mine.forest.loss.5pix.5km.master.did <- did.prep(master.forest.loss.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "loss", covariates = covariates.5pix)

mine.forest.loss.5pix.5km.master.did %>% group_by(first.loss.year) %>% 
  summarise(length(unique(CLUSTER)))

master.5km.minus2001 <- mine.forest.loss.5pix.5km.master.did %>% filter(first.loss.year != 1)
master.5km.minus2002 <- mine.forest.loss.5pix.5km.master.did %>% filter(first.loss.year >= 3)

# fit Gardner and Callaway and Sant'Anna DiD model
master.5km.minus2001.DiD <- fit.dynamic.DiD(master.5km.minus2001,
                                                   yname = "cumulative.forest.loss.perc",
                                                   xformula = NULL, method = c("csa", "gardner"))


forest.loss.master.5pix2001.plts.f <- plot.DiD2S.f(master.5km.minus2001.DiD, 
                                               pre.yrs = -5, post.yrs = 10, legend = "none",
                                               method = c("csa", "gardner"),
                                               pre.mine.col = c("#4393c3", "#4393c3"),
                                               post.mine.col = c("#8c510a", "#8c510a"),
                                               y.label = "Forest loss (% points)")


# degradation version
load("Data/raw/degradation.df/masolele.mines.5pixels.third.forest.5km.master.buffer.degradation.df.RData")
master.degradation.5pix.5km <- mine.buffer.forest.loss.df

master.degradation.5pix.5km %>% group_by(first.loss.year) %>% 
  summarise(length(unique(CLUSTER)))

# create relative time to mine covariates and add the covariate data.
mine.degradation.5pix.5km.master.did <- did.prep(master.degradation.5pix.5km,
                                                 lead.time = -20, post.time = 20,
                                                 type = "degradation", covariates = covariates.5pix)

master.degradation.5km.minus2001 <- mine.degradation.5pix.5km.master.did %>% filter(first.loss.year != 1)
master.degradation.5km.minus2002 <- mine.degradation.5pix.5km.master.did %>% filter(first.loss.year >= 3)

master.degradation.5km.minus2001.DiD <- fit.dynamic.DiD(master.degradation.5km.minus2001,
                                            yname = "cumulative.degradation.perc",
                                            xformula = NULL, method = c("csa", "gardner"))


degradation.master.5pix2001.plts.f <- plot.DiD2S.f(master.degradation.5km.minus2001.DiD, 
                                                   pre.yrs = -5, post.yrs = 10, legend = "none",
                                                   method = c("csa", "gardner"),
                                                   pre.mine.col = c("#4393c3", "#4393c3"),
                                                   post.mine.col = c("#e08214", "#e08214"),
                                                   y.label = "Degradation (% points)")

 DiD.plt.m2001 <-  ggarrange(forest.loss.master.5pix2001.plts.f$`Gardner 2022`,
                             forest.loss.master.5pix2001.plts.f$`Callaway & Sant'Anna 2021`,
                             degradation.master.5pix2001.plts.f$`Gardner 2022`,
                             degradation.master.5pix2001.plts.f$`Callaway & Sant'Anna 2021`,
                             ncol = 2, nrow = 2, 
                             labels = c("a", "b", "c", "d"))

 ggsave(path = "Outputs/Figures/Final/SM", 
        filename = "DiD.plt.m2001.png",
        DiD.plt.m2001, bg = "white",
        device = "png", width = 20, height = 13, units = "cm") 
 
## SM analysis 8 - incorporate potential spillover effects ---------------------
 ## 5 pix third forest
 load("Data/raw/forest.loss.df/masolele.mines.5pixels.third.forest.5km.master.buffer.forest.loss.df.RData")
 master.forest.loss.5pix.5km <- mine.buffer.forest.loss.df
 
 mine.forest.loss.5pix.5km.master.did <- did.prep(master.forest.loss.5pix.5km,
                                                  lead.time = -20, post.time = 20,
                                                  type = "loss", covariates = covariates.5pix)
 mine.forest.loss.5pix.5km.master.did.spill <- spill.tidy(mine.forest.loss.5pix.5km.master.did)
 
 forest.loss.spill.5pix.5km <- spillover.dynamic.DiD(mine.forest.loss.5pix.5km.master.did.spill, 
                                                     yname = "cumulative.forest.loss.perc")
 
 forest.loss.spillover.thirdforest.5pix.5km.master <- ggplot(filter(forest.loss.spill.5pix.5km, year.since > -7 & year.since < 11), 
                                                             aes(year.since, estimate, colour = effect,
                                                                 shape = effect)) +
   geom_point(position = position_dodge(width = .5)) +
   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, position = position_dodge(width = .5)) +
   geom_vline(xintercept = -1) +
   geom_hline(yintercept = 0, linetype = "dashed") +
   scale_colour_manual(values = c("black", "grey", "#8c510a")) +
   scale_shape_manual(values = c(1, 1, 16)) +
   xlab("Years since mining") +
   ylab("Forest loss (% points)") +
   theme_minimal() +
   theme(legend.title = element_blank(), legend.position = "bottom")
 
 # degradation
 load("Data/raw/degradation.df/masolele.mines.5pixels.third.forest.5km.master.buffer.degradation.df.RData")
 master.degradation.5pix.5km <- mine.buffer.forest.loss.df
 
 mine.degradation.5pix.5km.master.did <- did.prep(master.degradation.5pix.5km,
                                                  lead.time = -20, post.time = 20,
                                                  type = "degradation", covariates = covariates.5pix)
 
 mine.degradation.5pix.5km.master.did.spill <- spill.tidy(mine.degradation.5pix.5km.master.did)
 
 degradation.spill.5pix.5km <- spillover.dynamic.DiD(mine.degradation.5pix.5km.master.did.spill, 
                                                     yname = "cumulative.degradation.perc")
 
 degradation.spillover.thirdforest.5pix.5km.master <- ggplot(filter(degradation.spill.5pix.5km, year.since > -7 & year.since < 11), 
                                                             aes(year.since, estimate, colour = effect,
                                                                 shape = effect)) +
   geom_point(position = position_dodge(width = .5)) +
   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, position = position_dodge(width = .5)) +
   geom_vline(xintercept = -1) +
   geom_hline(yintercept = 0, linetype = "dashed") +
   scale_colour_manual(values = c("black", "grey", "#e08214")) +
   scale_shape_manual(values = c(1, 1, 16)) +
   xlab("Years since mining") +
   ylab("Degradation (% points)") +
   theme_minimal() +
   theme(legend.title = element_blank(), legend.position = "bottom")
 
 
 spillover.arr <- ggarrange(forest.loss.spillover.thirdforest.5pix.5km.master,
                            degradation.spillover.thirdforest.5pix.5km.master,
                            ncol = 1, labels = c("a", "b"), hjust = 0,
                            font.label = list(size = 9))
 
 
 ggsave(path = "Outputs/Figures/Final/SM", filename = "spillover.arr.png",
        spillover.arr, bg = "white",
        device = "png", width = 18, height = 18, units = "cm")
 
 