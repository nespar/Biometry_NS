lint("Assignment_4.R") #OOPS
library(lme4)
library(bbmle)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(marginaleffects)
#This data set is provided by DEEPEND from the Deep Sea Benefits research cruises.
DSB_unclean <- read.csv("https://raw.githubusercontent.com/nespar/Biometry_NS/refs/heads/main/DSB_allcruises_moll_datasheet_COPY_csv.csv")
DSB_species <- DSB_unclean |>
  rename(Species_name = Species.Original.ID, 
         N_number = N..number.,
         Date1 = Date,
         SampleID = Station)

DSB_trawl_unclean <-read.csv("https://raw.githubusercontent.com/nespar/Biometry_NS/refs/heads/main/DSB1_DSB2_Trawl_Data.csv")
DSB <- merge(DSB_species, DSB_trawl_unclean, by = "SampleID", all.x = TRUE) 

DSB <- DSB|>
  rename(Time = TowStartTime_CDT)

DSB$N_number <- as.numeric(DSB$N_number)
DSB$ML.TL <- as.numeric(DSB$ML.TL) 
DSB$Date1 <- as.Date(DSB$Date1, format = "%Y%m%d")
DSB$Ceph.Hetero.ptero <- tolower(DSB$Ceph.Hetero.ptero)
DSB$Ceph.Hetero.ptero <- as.factor(DSB$Ceph.Hetero.ptero)
str(DSB)
sum(is.na(DSB))


DSB_srich <- DSB|>
  group_by(MaxSampleDepth, Ceph.Hetero.ptero) |>
  mutate(species_rich = n_distinct(Species_name)) |>
  filter(!is.na(MaxSampleDepth)) |>
  mutate(sample_year = lubridate::year(Date1)) |>
  ungroup()

mean(DSB_srich$species_rich)
var(DSB_srich$species_rich)
#the mean < var, with a discrete response variable and a clear right skew, negative binomial would be a good choice.

DSB_srich |> ggplot(aes(x = MaxSampleDepth, y = species_rich, color = Ceph.Hetero.ptero)) +
  geom_point() +
  labs(x = "Depth", y = "Species Richness", title = "Species Richness Over Depths Across Cephalopods, Heteropods, and Pteropods")
#The purpose of the Deep Sea Benefits project is to investigate how Mollusk species distribution changes in the Northern Gulf of Mexico across the continental slope.
#The main question that I am looking into is how species richness of Cephalopods, Pteropods, and Heteropods changes with depth in the Northern Gulf of Mexico.

DSB_fit <- MASS::glm.nb(species_rich ~ MaxSampleDepth, data = DSB_srich)

AIC(DSB_fit)
exp(coef(DSB_fit))
exp(confint(DSB_fit))

plot_predictions(DSB_fit, condition = list("MaxSampleDepth")) +
  labs(x = "Sample Depth", y = "Species Richness")

#The most simple model, using Depth as a predictor for species_richness.
#This model would allow us to see that species richness does change with depth across the continental slope.
#The plot shows the correct general shape of the data, but biologically it would not decrease as gradually because of the physical differences of the deep sea zonations. 

DSB_q <- MASS::glm.nb(species_rich ~ (MaxSampleDepth + I(MaxSampleDepth^2)), data = DSB_srich)

AIC(DSB_q)
exp(coef(DSB_q))
exp(confint(DSB_q))

plot_predictions(DSB_q, condition = list("MaxSampleDepth")) +
  labs(x = "Sample Depth", y = "Species Richness")

#To expand this model, I added a quadratic component to the depth predictor because it allows the relationship
#to be curved, which is more representative of species richness in the pelagic ocean.
#This model allows us to have a more biologically accurate representation of the relationship between depth and diversity.
#The plot with the added quadratic component matches the expected trend much more accurately.

DSB_group <- MASS::glm.nb(species_rich ~ Ceph.Hetero.ptero * (MaxSampleDepth + I(MaxSampleDepth^2)), data = DSB_srich)

AIC(DSB_group)
exp(coef(DSB_group))
exp(confint(DSB_group))

plot_predictions(DSB_group, condition = list("MaxSampleDepth", "Ceph.Hetero.ptero")) +
  labs(x = "Sample Depth", y = "Species Richness") 

#In this model, I added a parameter for the Ceph.Hetero.Ptero groups. These groups are entirely different in morphology and ecology, such as trophic levels and abundance, 
#so they can not be analyzed as the same for species richness.
#The coefficients and intervals are much easier to interpret with this model because it allows us to interpret the groups separately
#Additionally, now the model plot looks very much like the sample data, as the difference in trends among the 3 groups is important.

DSB_WC <- MASS::glm.nb(species_rich ~ Ceph.Hetero.ptero * (MaxSampleDepth + I(MaxSampleDepth^2) + End_WaterColumn_Depth), data = DSB_srich)

AIC(DSB_WC)
exp(coef(DSB_WC))
exp(confint(DSB_WC))

plot_predictions(DSB_WC, condition = list("MaxSampleDepth", "Ceph.Hetero.ptero")) +
  labs(x = "Sample Depth", y = "Species Richness") 

#Beyond just the Sample Depth, which is the depth at which samples were collections, the water column depth gives incite to the environment being sampled.
#The actual sample depth may not tell the whole story, so the water column depth could supplement that in terms of the continental slope.\
#The fit of the graph for this model even further improved to resemble the sample data trends.


DSB_year <- MASS::glm.nb(species_rich ~ Ceph.Hetero.ptero * (MaxSampleDepth + I(MaxSampleDepth^2)) + sample_year, data = DSB_srich)

AIC(DSB_year)
exp(coef(DSB_year))
exp(confint(DSB_year))

plot_predictions(DSB_year, condition = list("MaxSampleDepth", "Ceph.Hetero.ptero", "sample_year")) +
  labs(x = "Sample Depth", y = "Species Richness")

#Due to the fact that the research cruises took place across 3 years, I was curious if that had any impact of the species_richness measured.
#With how deep sea sampling is, the specimens that are collected are highly dependent on what you simply come across. 
#This isn't a set difference between years, but rather a result of the relative lack of species abundance and encounters in the deep sea.
#Therefore, the year could have a random effect on the samples collected and impact the model, but there isn't enough variables to add it as a mixed model random effect, so it will be treated as a fixed effect.
#The plots for year don't seem to have a distinct visual difference between them, however the AIC indicates that the model was improved.
#Perhaps with more years, this difference would be better illustrated.

DSB_all <- MASS::glm.nb(species_rich ~ Ceph.Hetero.ptero * (MaxSampleDepth + I(MaxSampleDepth^2)) + End_WaterColumn_Depth + sample_year, data = DSB_srich)
AIC(DSB_all)
exp(coef(DSB_all))
exp(confint(DSB_all))

plot_predictions(DSB_all, condition = list("MaxSampleDepth", "Ceph.Hetero.ptero", "sample_year")) +
  labs(x = "Sample Depth", y = "Species Richness")

#Finally, I was curious if all the parameters combined would give the best representation of the data.
#The sample depth, water column depth, year, and group are all important to the species that are collected.
#I want to see how the effective the model is compared to the others, as this one has so many parameters. 
#With all the variables, the graphs resemble those from the DSB_WC model, with the 2 separate graphs for year.
#However, once again, that year graph does not seem to have any visual difference, but the AIC value is significantly improved.

AICtab(DSB_fit, DSB_q, DSB_group, DSB_year, DSB_WC, DSB_all)
#Based on the AICtab, the DSB_all model has the best fit for the data.
#This makes since because the species richness is known to be strongly correlated between Depth in the pelagic ocean.
#Additionally, it is known that there is distinct differences between diversity among the 3 Mollusk groups samples: Cephalopod, Heteropod, and Pteropod.
#Another component that was added was the Water Column depth because that is more reflective of the environment and can has relevance beyond what depths were sampled. 
#Finally, the year is an important parameter due to the difference among sampling across each sampling effort due to the relative lack of organism abundance of the deep sea.


