library(readr)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(ggmap)
library(maps)
library(marmap)



gom_bathy <- getNOAA.bathy(
  lon1 = -92,
  lon2 = -82,
  lat1 = 26,
  lat2 = 31,
  resolution = 4
)

bathy_df <- fortify(gom_bathy)

#PART ONE
citation("readr")
citation("tidyverse")
citation("ggplot2")
citation("Polychrome")
citation("ggmap")
citation("maps")
citation("marmap")

DSB_unclean <- read.csv("https://raw.githubusercontent.com/nespar/Biometry_NS/refs/heads/main/DSB_allcruises_moll_datasheet_COPY_csv.csv")
DSB_species <- DSB_unclean |>
  rename(Species_name = Species.Original.ID, 
         N_number = N..number.,
         Date1 = Date,
         SampleID = Station)

DSB_trawl_unclean <-read.csv("https://raw.githubusercontent.com/nespar/Biometry_NS/refs/heads/main/DSB1_DSB2_Trawl_Data.csv")
DSB <- merge(DSB_species, DSB_trawl_unclean, by = "SampleID", all.x = TRUE)


DSB$N_number <- as.numeric(DSB$N_number)
DSB$ML.TL <- as.numeric(DSB$ML.TL) 
DSB$Date1 <- as.character(DSB$Date1)
DSB$Ceph.Hetero.ptero <- tolower(DSB$Ceph.Hetero.ptero)
str(DSB)
sum(is.na(DSB))

#Figure 1. Species Richness
#https://stackoverflow.com/questions/27214282/add-bathymetry-lines-to-ggplot-using-marmap-package-and-getnoaa-bathy 
#https://cran.r-universe.dev/maps/doc/manual.html 

DSB_srich_date <- DSB |>
  group_by(Date1) |>
  mutate(species_rich = n_distinct(Species_name)) |>
  mutate(date_depth = max(MaxSampleDepth)) |>
  arrange(Date1)

US_map <- map_data("usa")

gom_bathy <- getNOAA.bathy(lon1 = -92, lon2 = -82, lat1 = 26, lat2 = 31, resolution = 4)

bathy_df <- fortify(gom_bathy)
Figure_1 <- ggplot() +
  geom_polygon(data = US_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  geom_contour(data = bathy_df, aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -2000, -3000), color = "gray60") +
  geom_point(data = DSB_srich_date, aes(StartLon, StartLat, col = species_rich)) +
  labs(title = "Species Richness Across Sampling Sites in the Gulf of Mexico", 
       caption = "Figure 1. shows the species richness amongst the collected deep-sea pelagic Mollusks across the northern Gulf of Mexico at each sampling site.",
       x = "Longitude", y = "Latitude") +
  coord_cartesian(xlim = c(-92, -82), ylim = c(26, 31)) +
  theme(plot.caption = element_text(hjust = 0.5)) +
  scale_color_viridis_c() +
  theme_minimal()

Figure_1
#Figure 2.
DSB$depth_cat <- cut(DSB$MaxSampleDepth,
  breaks = seq(0, max(DSB$MaxSampleDepth, na.rm = TRUE), by = 100),
  include.lowest = TRUE)

DSB_srich_depth <- DSB |>
  group_by(depth_cat, Ceph.Hetero.ptero) |>
  summarise(species_rich = n_distinct(Species_name)) |>
  filter(!is.na(depth_cat)) |>
  arrange(depth_cat)

Figure_2 <- DSB_srich_depth |>
  ggplot(aes(x = depth_cat, y = species_rich, fill = Ceph.Hetero.ptero)) +
  geom_col(width = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5), plot.caption = element_text(hjust = 0.5)) +
  labs(title = "Species Richness Across Various Depths", x = "Species Richness", y = "Depth Range (m)",
       caption = "Figure 2. shows a bargraph of the total and separate group species richness across depth increments of 100m.")
  
  
Figure_2

#Figure 3
DSB_topspp <- DSB |>
  group_by(Ceph.Hetero.ptero, Species_name) |>
  summarise(species_count = sum(N_number, na.rm = TRUE), .groups = "drop") |>
  group_by(Ceph.Hetero.ptero) |>
  slice_max(species_count, n = 5, with_ties = FALSE)
#determine the top most common species across the 3 groups
#Cephs -- Pterygioteuthis sp, Abralia redfieldi, Cranchia scabra, Japetella diaphana, Abralia sp
#Hetero -- Pterotrachea coronata, Pterotrachae scutata, Pterotrachae sp, Cardiapoda placenta, Carinaria lamarcki
#Ptero -- Clio pyramidata, Clio sp, Diacria trispinosa, Diacria sp, Cavolinia sp

DSB_counts <- DSB |>
  group_by(Species_name, Ceph.Hetero.ptero, MaxSampleDepth) |>
  summarise(species_count = sum(N_number, na.rm = TRUE)) |>
  arrange(desc(species_count))

species_order <- c(
  # Ceph
  "Pterygioteuthis sp", "Abralia redfieldi", "Cranchia scabra", "Japetlla diaphana", "Abralia sp",
  # Hetero 
  "Pterotrachea coronata", "Pterotrachea scutata", "Pterotrachea sp", "Cardiapoda placenta", "Carinaria lamarcki",
  # Ptero
  "Clio sp", "Clio pyramidata", "Diacria trispinosa", "Diacria sp",  "Cavolinia sp")

DSB_top5 <- DSB_counts |>
  filter(Species_name %in% c("Pterygioteuthis sp", "Abralia redfieldi", "Cranchia scabra", "Japetlla diaphana", "Abralia sp",
                             "Pterotrachea coronata", "Pterotrachea scutata", "Pterotrachea sp", "Cardiapoda placenta", "Carinaria lamarcki",
                             "Clio sp", "Clio pyramidata", "Diacria trispinosa", "Diacria sp", "Cavolinia sp")) |>
  mutate(Species_name = factor(Species_name, levels = species_order))

species <- unique(DSB_top5$Species_name)

cols <- setNames(palette36.colors(length(species)), species)


Figure_3 <- DSB_top5 |>
  ggplot(aes(x = MaxSampleDepth, y = species_count, col = Species_name)) +
  geom_point() +
  geom_line() +
  labs(x = "Sample Depth(m)", y = "Species Counts", title = "Species Abundance Across Various Depths",
       caption = "Figure 3. shows the species counts across various depths for the top 5 most abundant species in each group.") +
  scale_color_manual(values = cols) +
  theme(plot.caption = element_text(hjust = 0.5)) +
  facet_wrap(~Ceph.Hetero.ptero, scales = "free")  

Figure_3

#FIX THE COLOR SCHEME OF THIS ONEEEE and title and caption


#Figure 4.
Cephs <- DSB[which(DSB$Ceph.Hetero.ptero == "ceph"), ] 
Heteros <- DSB[which(DSB$Ceph.Hetero.ptero == "hetero"), ] 

Cephs_filter <- Cephs |>
  group_by(Species_name) |>
  filter(n_distinct(MaxSampleDepth) > 2) |>
  ungroup()
#To make the graph less convoluted with tons of species and since we are visualizing how mean length of the species changes across depths only visualize the species that occur at more than 1 depth
#COME BACK TO THIS TO DECIDE HOW YOU WANT TO FILTER IT

ML_depth <- Cephs_filter |>
  group_by(Species_name, MaxSampleDepth) |>
  summarise(mean_ML = mean(ML.TL, na.rm = TRUE))
#to get the means of the mantle length for each species at each depth 

Hetero <- DSB[which(DSB$Ceph.Hetero.ptero == "hetero"), ]
Hetero_filter <- Hetero |>
  group_by(Species_name) |>
  filter(n_distinct(MaxSampleDepth) > 1) |>
  ungroup()

TL_depth <- Hetero_filter |>
  group_by(Species_name, MaxSampleDepth) |>
  summarise(mean_TL = mean(ML.TL, na.rm = TRUE))

species2 <- unique(ML_depth$Species_name)
cols2 <- setNames(palette36.colors(length(species2)), species2)

Figure_4a <- ML_depth |>
  ggplot(aes(x = MaxSampleDepth, y = mean_ML, col = Species_name)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = cols2) +
  theme(plot.caption = element_text(hjust = 0.5)) +
  labs(y = "Mantle Length", x = "Sample Depth(m)", title = "Cephalopod Mantle Length Across Various Depths",
       caption = "Figure 4a. shows how mantle length changes across various depths for cephalopod species. Only species that were collected at more than 2 depths were used.")

Figure_4a

species3 <- unique(TL_depth$Species_name)
cols3 <- setNames(palette36.colors(length(species3)), species3)

Figure_4b <- TL_depth |>
  ggplot(aes(x = MaxSampleDepth, y = mean_TL, col = Species_name)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = cols3) +
  theme(plot.caption = element_text(hjust = 0.5)) +
  labs(y = "Total Length", x = "Sample Depth(m)", title = "Heteropod Length Across Various Depths",
       caption = "Figure 4b. shows how total length changes across various depths for heteropod species. Only species that were collected at more than 1 depth were used.")
  
  
Figure_4b


#PART TWO
twoa <- read_rds("https://github.com/bmaitner/biometry_course/raw/refs/heads/main/midterm/b_2a.RDS")
twob <- read_rds("https://github.com/bmaitner/biometry_course/raw/refs/heads/main/midterm/b_2b.RDS")
twoc <- read_rds("https://github.com/bmaitner/biometry_course/raw/refs/heads/main/midterm/b_2c.RDS")

#2a
str(twoa)
power <- lm(log(Egg.Mass.g) ~ log(Adult.Body.Mass.g), data = twoa)

twoa |>
  ggplot(aes(x = log(Adult.Body.Mass.g), y = log(Egg.Mass.g))) +
  geom_point() +
  geom_smooth(method = "lm", color = "red")
#The log of both axes visualizes the data more clearly, as it spreads out the points a bit more.
twoa |>
  ggplot(aes(x = log10(Egg.Mass.g))) +
  geom_histogram()
mean(twoa$Egg.Mass.g)
median(twoa$Egg.Mass.g)
max(twoa$Egg.Mass.g)

#Based on the scatterplot of the data, a power law would be a good fit. The data seems to increase almost linearly and because the log of the axes visualizes the data the best, the log of a power law ends up being a linear regression.
#Egg Mass (g) increases exponentially with Adult Body Mass (g), then gradually plateaus with an outlier from a very large bird.   
#The original histogram was very skewed right, so the log10 allows for a more normal distribution
#Egg Mass (g) is a continuous variable, with no ceiling or negatives. Additionally, the Egg Mass(g) data is skewed right with a mean of 22.66g and a median of 3.92g.
#Therefore, a lognormal probability distribution is likely a good fit.


#2b
str(twob)
twob|>
  ggplot(aes(x = Total.Fish.Abundance, y = Species.Richness)) +
  geom_point()
  
#The log of the x-axis makes the plot easier to interpret.
mu = mean(twob$Species.Richness, na.rm = TRUE)
vari = var(twob$Species.Richness, na.rm = TRUE)
size <- mu^2 / (vari - mu)

twob|>
  ggplot(aes(x = Species.Richness)) +
  geom_histogram()

#mean = 21.784 and variance = 112.556, so not Poisson
sd(twob$Species.Richness)

#The shape of the plot using the log of the x-axis demonstrates allometric growth that increases with fish abundance. This seems to match a logistic shape. 
#Species Richness is a discrete variable, as it is a count of the species. 
#The histogram of Species Richness does not have any skew and would have a ceiling for the possible species that could exist in that environment.
#Therefore the Negative Binomial distribution would be a good fit for this data. 


#2c
str(twoc)
twoc|>
  ggplot(aes(x = Rugosity, y = Hard.Bottom.Cover.Pct)) +
  geom_point() 


#The shape of this plot fits well with Michaelis-menten as it is saturating and due to the y-axis being a percentage, there is a given asymptote at 100. 
#Hard.Bottom.Cover.Pct is a continuous variable. 
#It is quantified as a percentage of cover, which fits the Beta distribution well.


#PART THREE
sample_size <- 5
a <- 2
b <- 1

nsim <- 400
sdvec <- seq(1, 20, by = 1)
power.sd <- numeric(length(sdvec))
pval <- numeric(nsim)


for(j in 1:length(sdvec)){
  for(i in 1:nsim){
    
    x <- sample(x = 1:20,
                size = sample_size,
                replace = TRUE)
    
    sd <- sdvec[j]
    y_det <- a + b*x
    y <- rnorm(n = length(y_det),
               mean = y_det,
               sd = sd)
    
    m <- lm(y ~ x)
    
    #get p-value
    
    pval[i] <- coef(summary(m))["x","Pr(>|t|)"]
    
  }#end i loop
  power.sd[j] <- sum(pval < 0.05)/nsim
  
}#end j loop

par(mfrow = c(2,2))
plot(power.sd ~ sdvec, main = sample_size, xlab = "Standard Deviation", ylab = "Statistical Power (p-value > 0.05)", ylim = c(0, 1))

sample_size2 <- 20
a <- 2
b <- 1

nsim <- 400
sdvec2 <- seq(1, 20, by = 1)
power.sd2 <- numeric(length(sdvec2))
pval <- numeric(nsim)


for(j in 1:length(sdvec2)){
  for(i in 1:nsim){
    
    x <- sample(x = 1:20,
                size = sample_size2,
                replace = TRUE)
    
    sd <- sdvec2[j]
    y_det <- a + b*x
    y <- rnorm(n = length(y_det),
               mean = y_det,
               sd = sd)
    
    m <- lm(y ~ x)
    
    #get p-value
    
    pval[i] <- coef(summary(m))["x","Pr(>|t|)"]
    
  }#end i loop
  power.sd2[j] <- sum(pval < 0.05)/nsim
  
}#end j loop

plot(power.sd2 ~ sdvec2, main = sample_size2, xlab = "Standard Deviation", ylab = "Statistical Power (p-value > 0.05)", ylim = c(0, 1))

sample_size3 <- 50
a <- 2
b <- 1

nsim <- 400
sdvec3 <- seq(1, 20, by = 1)
power.sd3 <- numeric(length(sdvec3))
pval <- numeric(nsim)


for(j in 1:length(sdvec3)){
  for(i in 1:nsim){
    
    x <- sample(x = 1:20,
                size = sample_size3,
                replace = TRUE)
    
    sd <- sdvec3[j]
    y_det <- a + b*x
    y <- rnorm(n = length(y_det),
               mean = y_det,
               sd = sd)
    
    m <- lm(y ~ x)
    
    #get p-value
    
    pval[i] <- coef(summary(m))["x","Pr(>|t|)"]
    
  }#end i loop
  power.sd3[j] <- sum(pval < 0.05)/nsim
  
}#end j loop

plot(power.sd3 ~ sdvec3, main = sample_size3, xlab = "Standard Deviation", ylab = "Statistical Power (p-value > 0.05)", ylim = c(0, 1))

sample_size4 <- 100
a <- 2
b <- 1

nsim <- 400
sdvec4 <- seq(1, 20, by = 1)
power.sd4 <- numeric(length(sdvec4))
pval <- numeric(nsim)


for(j in 1:length(sdvec4)){
  for(i in 1:nsim){
    
    x <- sample(x = 1:20,
                size = sample_size4,
                replace = TRUE)
    
    sd <- sdvec4[j]
    y_det <- a + b*x
    y <- rnorm(n = length(y_det),
               mean = y_det,
               sd = sd)
    
    m <- lm(y ~ x)
    
    #get p-value
    
    pval[i] <- coef(summary(m))["x","Pr(>|t|)"]
    
  }#end i loop
  power.sd4[j] <- sum(pval < 0.05)/nsim
  
}#end j loop
plot(power.sd4 ~ sdvec4, main = sample_size4, xlab = "Standard Deviation", ylab = "Statistical Power (p-value > 0.05)", ylim = c(0, 1))

#The figures of Statistical Power as a function of Standard Deviations across 4 different sample sizes (5, 20, 50, and 100) demonstrate the relationship between "noise" and sample sizes. 
#Statistical Power increases as sample size increases and as standard deviation decreases.
#These figures indicate that increasing your sample size will produce more reliable results in terms of statistical power, even across varying standard deviations. 





