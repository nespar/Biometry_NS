# Assignment 3, Part 1: model modification and selection
#
# Three biological questions, each with a starting model already written and
# fitted. Your job, for each one:
#
#   1. Fit the starting model as given and look at its estimates.
#   2. Build alternative models. Each question below lists ideas, but you are
#      free to try others, including adding parameters.
#   3. Compare all of your models for that question with AICtab().
#   4. Say which model (or models) you would keep, and why, in terms of the
#      biology, not only the AIC.
#
# Lecture 18 covers everything you need: the formula version of mle2(), the
# parameters = argument for letting a parameter differ between groups, and the
# rules for a fair comparison.

library(bbmle)

# The data are in the course data package. If you have not installed it yet:
#   install.packages("remotes")
#   remotes::install_github("bmaitner/biometry_course", subdir = "biometrydata")
# If you cannot install it, each section gives a read.csv() or readRDS() line
# that loads the same data from the course repository.

library(biometrydata)


# 1a: Competition between bean beetles ----------------------------------------

  # Seed beetles lay their eggs on beans, and the larvae develop inside them,
  # so beetles sharing a dish of beans compete for food. Each row is one dish:
  # a known number of females of two species placed on 2 g of beans, and the
  # number of offspring of each species that emerged.
  #
  # Here we model the offspring of Callosobruchus maculatus (species 1, "Cmac")
  # when it shares beans with Callosobruchus subinnotatus (species 2, "Csub").
  #
  #   start_females_1   C. maculatus females placed in the dish
  #   start_females_2   C. subinnotatus females placed in the dish
  #   progeny_1         C. maculatus offspring that emerged
  #
  # ?beetle_competition_interspecific describes the experiment.

  beetles <- beetle_competition_interspecific
  # or: beetles <- read.csv("https://github.com/bmaitner/biometry_course/raw/refs/heads/main/data/Beetle_competition/competition_interspecific.csv")

  beetles <- beetles[beetles$pair == "Cmac_Csub" & !is.na(beetles$progeny_1), ]

  # plot if you like: offspring per female against the number of competitors

  plot(beetles$start_females_1 + beetles$start_females_2,
       beetles$progeny_1 / beetles$start_females_1,
       xlab = "females in the dish, both species",
       ylab = "C. maculatus offspring per female")

  # The starting model is a competition model in the Beverton-Holt family:
  #
  #   expected offspring = R * N1 / (1 + (a11 * N1 + a12 * N2)^d)
  #
  #   R     offspring per female with no competition at all
  #   a11   how much each C. maculatus female reduces that (within species)
  #   a12   how much each C. subinnotatus female reduces it (between species)
  #   d     how sharply competition bites as the dish fills up
  #   k     negative binomial size: how variable dishes are
  #
  # Bounds keep every parameter positive, so the optimizer cannot try
  # impossible values. If you remove a parameter, remove it from start, lower
  # and upper too.

  beetle.fit <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 /
                                           (1 + (a11 * start_females_1 +
                                                 a12 * start_females_2)^d),
                                         size = k),
                     start = list(R = 20, a11 = 0.1, a12 = 0.1, d = 1, k = 2),
                     lower = c(R = 0.01, a11 = 0, a12 = 0, d = 0.1, k = 0.01),
                     upper = c(R = 1000, a11 = 10, a12 = 10, d = 10, k = 100),
                     method = "L-BFGS-B",
                     data = beetles)

  coef(beetle.fit)
  AIC(beetle.fit)

  # Questions to build models around:
  #   - Is d needed? With d = 1 (leave it out of the formula) this is the
  #     classic Beverton-Holt model.
  #   - Does C. subinnotatus affect C. maculatus at all? (What value of a12
  #     would mean "no"?)
  #   - Is a C. subinnotatus female as harmful as a C. maculatus female?
  #     (Hint: use one parameter for both.)
  #   - A different shape: the Ricker model,
  #     R * N1 * exp(-a11 * N1 - a12 * N2). Is it nested in the starting model?

  
  beetle.fit2 <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 /
                                           (1 + (a11 * start_females_1 +
                                                   a12 * start_females_2)),
                                         size = k),
                     start = list(R = 20, a11 = 0.1, a12 = 0.1, k = 2),
                     lower = c(R = 0.01, a11 = 0, a12 = 0, k = 0.01),
                     upper = c(R = 1000, a11 = 10, a12 = 10, k = 100),
                     method = "L-BFGS-B",
                     data = beetles)
  
  coef(beetle.fit2)
  AIC(beetle.fit2) 
  #beetle.fit2 without d in the equation is a slightly better fit than the starting model.
  
  beetle.fit3 <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 /
                                           (1 + (a11 * start_females_1 +
                                                   start_females_2)^d),
                                         size = k),
                     start = list(R = 20, a11 = 0.1, d = 1, k = 2),
                     lower = c(R = 0.01, a11 = 0, d = 0.1, k = 0.01),
                     upper = c(R = 1000, a11 = 10, d = 10, k = 100),
                     method = "L-BFGS-B",
                     data = beetles)
  
  coef(beetle.fit3)
  AIC(beetle.fit3)
  #If C. subinnotatus had no effect on C. maculatus, then a12 would not be a parameter in this function because
  #the a12 parameter quantifies the effect of competition from C. subinnotatus. To build a model to test this,
  #I removed the a12 parameter from the model, and it output a higher AIC value, so it's a worse fit. This indicates 
  #that there is an effect from competition by C. subinnotatus.
  

  beetle.fit4 <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 /
                                           (1 + (a11 * start_females_1 +
                                                   a11 * start_females_2)^d),
                                         size = k),
                     start = list(R = 20, a11 = 0.1, d = 1, k = 2),
                     lower = c(R = 0.01, a11 = 0, d = 0.1, k = 0.01),
                     upper = c(R = 1000, a11 = 10, d = 10, k = 100),
                     method = "L-BFGS-B",
                     data = beetles)
  
  coef(beetle.fit4)
  AIC(beetle.fit4)
  
  #I removed the a12 parameter and replaced the reduction of C. sibinnotatus with the parameter for the reduction
  #of C. maculatus, giving them the same parameter. This led to convergence failure, and a worse fit model than the original.
  #Therefore, C. subinnotatus female  is not as as harmful as a C. maculatus female.

  beetle.fit5 <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 * exp(-a11 * start_females_1 - a12 * start_females_2),
                                         size = k),
                     start = list(R = 20, a11 = 0.1, a12 = 0.1, k = 2),
                     lower = c(R = 0.01, a11 = 0, a12 = 0, k = 0.01),
                     upper = c(R = 1000, a11 = 10, a12 = 10, k = 100),
                     method = "L-BFGS-B",
                     data = beetles)
  
  coef(beetle.fit5) 
  AIC(beetle.fit5)
  #The Ricker model has a better fit than the Beverton-Holt model, as the AIC value is lower.
  
  AICtab(beetle.fit, beetle.fit2, beetle.fit3, beetle.fit4, beetle.fit5, base = TRUE)
  #The Ricker model is the best fit model from the alternatives ran. It has the lowest AIC value.
  #This makes sense because the Ricker model accounts for the competition has the parents increase in number,
  #without using the d value seen in the starting model.
  
# 1b: Wing length in hawks and eagles -----------------------------------------

  # Long-distance migrants should benefit from longer wings, which make
  # sustained flight cheaper. Do migratory hawks and eagles (family
  # Accipitridae) have longer wings than residents of the same body mass?
  #
  #   Wing.Length   mm
  #   Mass          g
  #   migration     resident, partial (some populations migrate), or migrant
  #
  # Each row is a species, so 249 species in all. One caution worth a
  # sentence in your answer: related species share traits, so species are not
  # fully independent data points. Accounting for that is beyond this course.

  hawks <- avonet
  # or: hawks <- readRDS(url("https://github.com/bmaitner/biometry_course/raw/refs/heads/main/data/Avonet/AVONET1_BirdLife.rds"))

  hawks <- hawks[hawks$Family1 == "Accipitridae", c("Species1", "Wing.Length", "Mass", "Migration")]
  hawks <- na.omit(hawks)
  hawks$migration <- factor(hawks$Migration, levels = 1:3,
                            labels = c("resident", "partial", "migrant"))

  # plot if you like

  plot(log(hawks$Mass), log(hawks$Wing.Length), col = hawks$migration,
       xlab = "log body mass", ylab = "log wing length")
  legend("topleft", levels(hawks$migration), col = 1:3, pch = 1)

  # The starting model: wing length is lognormal, and its log rises linearly
  # with log mass. That is a power law, wing = e^int * mass^b, the usual way
  # one body measurement scales with another.
  #
  #   int     wing length at a given mass (on the log scale)
  #   b       how fast wing length rises with mass
  #   log_s   the log of the scatter around the line. Fitting the log keeps
  #           the standard deviation positive, as in Lecture 16
  #
  # parameters = lets every one of the three differ among the migration
  # classes, which is deliberately more than the data are likely to support.

  hawk.fit <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                        sdlog = exp(log_s)),
                   start = list(int = 3, b = 0.3, log_s = -2),
                   parameters = list(int ~ migration,
                                     b ~ migration,
                                     log_s ~ migration),
                   data = hawks)

  coef(hawk.fit)
  AIC(hawk.fit)

  # Questions to build models around:
  #   - Do migrants have longer wings for their mass? (Which parameter is that?)
  #   - Does wing length scale with mass differently in migrants?
  #   - Are some migration classes more variable than others?
  #   - Is the distinction between partial migrants and full migrants needed,
  #     or would migrant versus resident do? (You would need a new column.)


  hawk.fit2 <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                        sdlog = exp(log_s)),
                   start = list(int = 3, b = 0.3, log_s = -2),
                   parameters = list(int ~ migration),
                   data = hawks)
  
  coef(hawk.fit2)
  AIC(hawk.fit2)
  #Migrants would have longer wings for their mass if the int parameter is larger in migrants.
  #By only having the int parameter differ with migration status, it gives us a more clear answer
  #and a better fit model.
  
  hawk.fit3 <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                        sdlog = exp(log_s)),
                   start = list(int = 3, b = 0.3, log_s = -2),
                   parameters = list(b ~ migration),
                   data = hawks)
  
  coef(hawk.fit3)
  AIC(hawk.fit3)
  #This question is similar to the hawk one, but using the b parameter only instead to indicate if this gives us a better fit.
  #The fit is better than the starting model, indicating that wing length does scale differently with migrants.
  #It is faster, as the b value is higher for the migrants. 
  
  hawk.fit4 <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                        sdlog = exp(log_s)),
                   start = list(int = 3, b = 0.3, log_s = -2),
                   parameters = list(log_s ~ migration),
                   data = hawks)
  
  coef(hawk.fit4)
  AIC(hawk.fit4)
  #Using the same idea as the other 2 questions but with the log_s parameter, the fit of the model is worse
  #than the original model. 
  
  hawks$migration2 <- factor(
    ifelse(hawks$Migration == 1, "resident", "migrant"),
    levels = c("resident", "migrant")
  )
  
  hawk.fit5 <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                        sdlog = exp(log_s)),
                   start = list(int = 3, b = 0.3, log_s = -2),
                   parameters = list(int ~ migration2),
                   data = hawks)
  
  coef(hawk.fit5)
  AIC(hawk.fit5)
  #With only the distinction being between migrants and residents, this is a better fit model.
  
  AICtab(hawk.fit, hawk.fit2, hawk.fit3, hawk.fit4, hawk.fit5, base = TRUE)
  #The best fit model is hawk.fit5 which only distinguishes migration type as resident vs migrant and
  #only in the int parameter. The partial category likely has lots of variation, that could be insignificant
  #and the other parameters don't need to be broken up for the purpose of this model. 
  
  
# 1c: Reef fish species richness ----------------------------------------------

  # How many fish species does a diver record on a survey of a Florida coral
  # reef, and what explains the differences?
  #
  #   species_richness   number of species recorded on the survey
  #   depth              tens of metres (so 1.5 is 15 m)
  #   visibility         underwater visibility, tens of metres
  #   year               years since 2014
  #   region             Tortugas, Keys, or North (Miami to Martin County)
  #   protected          "yes" inside a marine protected area, otherwise "no"
  #
  # Depth and visibility are in tens of metres because depth is squared in
  # the model, and squaring 30 m gives numbers large enough to stall the
  # optimizer. Rescaling a covariate changes nothing about the biology.

  reef <- reef_fish_surveys
  # or: reef <- read.csv("https://github.com/bmaitner/biometry_course/raw/refs/heads/main/data/Reef_fish/NCRMP_reef_fish_surveys.csv")

  reef$depth <- reef$DEPTH / 10
  reef$visibility <- reef$UNDERWATER_VISIBILITY / 10
  reef$year <- reef$YEAR - 2014
  reef$region <- factor(reef$REGION_DESCRIPTION,
                        levels = c("Dry Tortugas",
                                   "Florida Keys: from Key West to Miami",
                                   "Miami to North Martin County"),
                        labels = c("Tortugas", "Keys", "North"))
  reef$protected <- ifelse(reef$MPA_NAME == "Unprotected", "no", "yes")

  # plot if you like

  plot(reef$DEPTH, reef$species_richness, col = rgb(0, 0, 0, 0.1),
       xlab = "depth (m)", ylab = "species recorded")

  # The starting model: richness is Poisson, its log rises and then falls
  # with depth (the d2 * depth^2 term makes a hump), and it also depends on
  # visibility, year, region and protection.

  reef.fit <- mle2(species_richness ~ dpois(lambda = exp(int +
                                                         d1 * depth +
                                                         d2 * depth^2 +
                                                         v * visibility +
                                                         t * year)),
                   start = list(int = 3, d1 = 0, d2 = 0, v = 0, t = 0),
                   parameters = list(int ~ region + protected),
                   data = reef)

  coef(reef.fit)
  AIC(reef.fit)

  # Questions to build models around:
  #   - Is the Poisson the right distribution? Check before anything else:
  #     what did Lecture 15 use var() / mean() for?
  #   - With 8,560 surveys, AIC may keep nearly every term. If so, which
  #     effects are large enough to matter? Back-transform them: exp() of a
  #     coefficient is the multiplier on expected richness.
  #   - Visibility is a property of the water on the day, not of the reef.
  #     What might a visibility effect be measuring?
  #   - Does protection matter once region is accounted for?
  
  mean(reef$species_richness)
  
  var(reef$species_richness)
  #mean < var, so cant be poissons --> try neg binomial 
  
  reef.fit2 <- mle2(species_richness ~ dnbinom(mu = exp(int +
                                                           d1 * depth +
                                                           d2 * depth^2 +
                                                           v * visibility +
                                                           t * year), size = k ),
                   start = list(int = 3, d1 = 0, d2 = 0, v = 0, t = 0, k = 2),
                   parameters = list(int ~ region + protected),
                   data = reef)
  
  coef(reef.fit2)
  AIC(reef.fit2)
  #negative binomial produces a much better model with a significant decrease in AIC
  #had NaNs produced.

  reef.fit3 <- mle2(species_richness ~ dnbinom(mu = (int +
                                                          d1 * depth +
                                                          d2 * depth^2 +
                                                          v * visibility +
                                                          t * year), size = k ),
                    start = list(int = 3, d1 = 0, d2 = 0, v = 0, t = 0, k = 2),
                    parameters = list(int ~ region + protected),
                    data = reef)
  
  coef(reef.fit3)
  AIC(reef.fit3)
  #Removed the exp to back-transform it, and it had a slightly better fit.
  
  reef.fit4 <- mle2(species_richness ~ dnbinom(mu = (int +
                                                          d1 * depth +
                                                          d2 * depth^2 +
                                                          t * year), size = k ),
                    start = list(int = 3, d1 = 0, d2 = 0, t = 0, k = 2),
                    parameters = list(int ~ region + protected),
                    data = reef)
  
  coef(reef.fit4)
  AIC(reef.fit4)
  #Worse fit without visibility as a parameter. Despite not being a direct measurement of the reef,
  #visibility can also be an important measurement of turbidity and light available to the reef which is key to the nutrients available.
  
  reef.fit5 <- mle2(species_richness ~ dnbinom(mu = (int +
                                                       d1 * depth +
                                                       d2 * depth^2 +
                                                       v * visibility +
                                                       t * year), size = k ),
                    start = list(int = 3, d1 = 0, d2 = 0, v = 0, t = 0, k = 2),
                    parameters = list(int ~ region),
                    data = reef)
  
  coef(reef.fit5)
  AIC(reef.fit5)
  
  #After removing protection as a parameter grouping, the fit was worse. Protection does have to be accounted for
  #because even within the regional sections, there may be sites within protected areas and sites that aren't. 
  
  AICtab(reef.fit, reef.fit2, reef.fit3, reef.fit4, reef.fit5, base = TRUE)
  #Reef.fit3 was the best fit model for these data. This makes sense because all of the parameters originally included
  #play a major role, with visibility impacting light and therefore nutrient availability and with protected regions impacting
  #anthropogenic stressors on the reef. 
  
  
  
#PART 2
library(tidyverse)
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
DSB$Date1 <- as.character(DSB$Date1)
DSB$Ceph.Hetero.ptero <- tolower(DSB$Ceph.Hetero.ptero)
DSB$TowStartTime_CDT <
str(DSB)
sum(is.na(DSB))

DSB_srich <- DSB |>
  group_by(MaxSampleDepth, Ceph.Hetero.ptero, SolarCycle) |>
  summarise(species_rich = n_distinct(Species_name)) |>
  filter(!is.na(MaxSampleDepth)) |>
  arrange(species_rich)

mean(DSB_srich$species_rich)
var(DSB_srich$species_rich)
#the mean < var, with a discrete response variable and a clear right skew, negative binomial would be a good choice.

plot(DSB_srich$species_rich, DSB_srich$MaxSampleDepth, col = as.factor(DSB_srich$Ceph.Hetero.ptero), xlab = "Depth", ylab = "Species Richness")

DSB.fit <- mle2(species_rich ~ dnbinom(mu = exp(int + d1 * MaxSampleDepth + t * ), size = k),
                start = list(int = 1, d1 = 0, k = 2),
                parameters = list(int ~ Ceph.Hetero.ptero),
                data = DSB_srich)
  
coef(DSB.fit)
AIC(DSB.fit)

#try without separating by group
DSB.fit2 <- mle2(species_rich ~ dnbinom(mu = exp(int + d1 * MaxSampleDepth), size = k),
                start = list(int = 1, d1 = 0, k = 2),
                data = DSB_srich)

coef(DSB.fit2)
AIC(DSB.fit2)

#slightly worse than with the int parameter separated by Ceph.Hetero.ptero grouping

#now trying by separating by group and based on Day vs Night because many of these species are Diel Migrating organisms.
DSB.fit3 <- mle2(species_rich ~ dnbinom(mu = exp(int + d1 * MaxSampleDepth), size = k),
                 start = list(int = 1, d1 = 0, k = 2),
                 parameters = list(int ~ Ceph.Hetero.ptero + SolarCycle),
                 data = DSB_srich)

coef(DSB.fit3)
AIC(DSB.fit3)

#The fit is much worse, Solar Cycle doesn't seem to impact the species.richness 

#MAYBE TRY DOING THE BACK-TRANSFORMING THING BUT IDK THAT SOOOOOOO

AICtab(DSB.fit, DSB.fit2, DSB.fit3, base = TRUE)

#Based on the AIC values from the separate models, the original model is the best fit for the data.
#This makes sense because the distribution of species is very different depending on the group.
#However, the solar cycle doesn't seem to have an impact, likely because light is not an important factor in the Deep Sea.