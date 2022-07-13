library(ggmap)
library(rethinking)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggrepel)
setwd("E:/Bamberg Uni/Typology course/brms_models")



wowa <- readRDS("wowa_recoded.rds")


final_model <- readRDS("upd_prior_fitted/final_training_model_4chains.rds")
stan_data <- readRDS("final_testing_set.rds")
geodistance <- stan_data$geo_dist
post <- extract.samples(final_model)


K <-matrix(0,nrow=24,ncol=24)
for (i in 1:24){
  for (j in 1:24){
    K[i, j] <- median(post$etasq_geo)*
      exp(-median(post$rhosq_geo)*geodistance[i, j]^2)
  }
}
diag(K) <- median(post$etasq_geo) + median(post$delta_geo)
Rho <- round(cov2cor(K), 2)

wowa24 <- wowa %>% 
  group_by(doculect) %>% 
  filter(row_number(doculect) == 1) %>% 
  arrange(doculect) %>% 
  mutate(
    family = case_when(
      affiliation1 == "turk" ~ "Turkic",
      affiliation1 == "kart" ~ "Kartvelian",
      affiliation1 == "iran" ~ "Iranian", 
      affiliation1 == "semi" ~ "Semitic", 
      affiliation1 == "hell" ~ "Hellenic"
    )
    )


wowa24 <- data.table(wowa24)
wowa24[doculect == "northern+ankara"]$longitude <- wowa24[doculect == "northern+ankara"]$longitude + 0.5
wowa24[doculect == "northern+ankara"]$latitude <- wowa24[doculect == "northern+ankara"]$latitude + 0.5


map_wowa <- get_stamenmap(bbox = c(left = 27, bottom = 24, 
                                   right = 64, top = 48), 
                          maptype = "terrain-background", 
                          zoom = 7)

map_points <-
  ggmap(map_wowa)  

for (i in 1:24){
  for (j in 1:24){
    if (i < j){ 
    map_points <- map_points + geom_segment(x = wowa24$longitude[i], 
                                            y = wowa24$latitude[i], 
                                            xend = wowa24$longitude[j], 
                                            yend = wowa24$latitude[j], 
                                            alpha = Rho[i, j]^1.5, 
                                            size = 1)
    }
  }
}

map_points <- map_points +
  geom_point(data = wowa24, aes(x = longitude, y = latitude, color = family), 
             size = 3) + 
  geom_label_repel(data = wowa24, aes(x = longitude, y = latitude, 
                               label = doculect), 
             position = position_jitter(height = ifelse(wowa24$doculect == "cewlig", 1, 0), 
                                        width = ifelse(wowa24$doculect == "cewlig", -1, 0))) +
  scale_color_brewer(palette = "Set1") +
  xlab("") + 
  ylab("") +
  labs(color = "Family") +
  theme(legend.position = "bottom")


ggsave("map_corr.png", width = 12, height = 9.5)

names <- unique(wowa24$doculect)
dist <- c()
corr <- c()


for (i in 1:24){
  for (j in 1:24){
    if (i < j){
      dist <- c(dist, geodistance[i, j])
      corr <- c(corr, Rho[i, j])
    }
  }
}
dist_corr <- data.table(dist, corr)
dist_corr


ggplot(dist_corr) + 
  geom_point(aes(x = dist, y = corr)) +
  xlab("Geographic distance (thousands of km)") +
  ylab("Correlation") +
  theme_classic()

ggsave("dist_corr.png", units = "px", 
       width = 2500, height = 2500)

Rho_copy <- Rho
diag(Rho_copy) <- 0
range(Rho_copy)
##############################################################################

phylodist <- stan_data$phylo
K <- matrix(0, nrow = 24, ncol = 24)
for (i in 1:24){
  for (j in 1:24){
    K[i, j] <- median(post$etasq_phy)*
      exp(-median(post$rhosq_phy)*phylodist[i, j]^2)
    }
  }



diag(K) <- median(post$etasq_phy) + median(post$delta_phy)
Rho_phy <- round(cov2cor(K), 2)

Rho_phy

Rho_phy_copy <- Rho_phy
diag(Rho_phy_copy) <- 0
range(Rho_phy_copy)


dist_phy <- c()
corr_phy <- c()


for (i in 1:24){
  for (j in 1:24){
    if (i < j){
      dist_phy <- c(dist_phy, phylodist[i, j])
      corr_phy <- c(corr_phy, Rho_phy[i, j])
    }
  }
}
dist_corr_phy <- data.table(dist_phy, corr_phy)
dist_corr_phy

ggplot(dist_corr_phy) + 
  geom_point(aes(x = dist_phy, y = corr_phy)) +
  xlab("Phylogenetic distance") +
  ylab("Correlation") +
  ylim(0, 1) +
  theme_minimal()


range(Rho_phy_zero)

