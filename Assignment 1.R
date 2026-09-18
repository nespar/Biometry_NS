#1
#This is data from an ongoing Deepend research project called Deep Sea Benefits. 
#The data has deep sea cephalopod, heteropod, and pteropod data over 3 research cruises over 3 years, looking at distribution changes over the continental slope vs open ocean. 
#Data was obtained directly through Dr. Heather Judkins at USFSP/Deepend, and waiting for the rest of the data, so date will be used to substitute the location data in the mean time. 

#2
DSB_unclean <- read.csv("DSB_allcruises_moll_datasheet_COPY_csv.csv") #entering data into

#3
class(DSB_unclean)
#class is a data.frame
library(tidyverse)

#4
#there are 7 variables which include 6 character variables; Station, Ceph.Hetero.Ptero, Species.Original.ID,
#N..number, ML.TL, Jar.Label; and 1 integer variable; Date.

#5
str(DSB_unclean) #looking to make sure it looks correct with columns and rows
sum(is.na(DSB_unclean)) #checking to see if expected na values are as expected
DSB <- DSB_unclean |>
        rename(Species_name = Species.Original.ID, 
         N_number = N..number.,
         Date1 = Date) #renamed some more complex column names to make simpler
DSB$N_number <- as.numeric(DSB$N_number)
DSB$ML.TL <- as.numeric(DSB$ML.TL) #needed to change these 3 because the str showed it as chr
DSB$Date1 <- as.character(DSB$Date1) #needed to change to the 
DSB$Ceph.Hetero.ptero <- tolower(DSB$Ceph.Hetero.ptero)
str(DSB)
sum(is.na(DSB)) #checking that it worked 

#6 --> did some summary statistics of the ML/TL data further in the code, but had limited numeric variables to work with

Cephs <- DSB[which(DSB$Ceph.Hetero.ptero == "ceph"), ] 
Heteros <- DSB[which(DSB$Ceph.Hetero.ptero == "hetero"), ]
Pteros <- DSB[which(DSB$Ceph.Hetero.ptero == "ptero"), ] #had to subset the 3 general groups, so graphs could be separate

Ceph_counts <- Cephs |>
  group_by(Species_name) |>
  summarise(species_count = sum(N_number, na.rm = TRUE)) |>
  arrange(desc(species_count))
#getting the species_counts, but data has multiple observations for the same species overlapping, so this condenses that 
top_cephs <- Ceph_counts[1:15,] #chose only top 15 because graphs would be too convoluted, not readable 
top_cephs
Hetero_counts <- Heteros |>
  group_by(Species_name) |>
  summarise(species_count = sum(N_number, na.rm = TRUE)) |>
  arrange(desc(species_count))
top_heteros <- Hetero_counts[1:15,]
Ptero_counts <- Pteros |>
  group_by(Species_name) |>
  summarise(species_count = sum(N_number, na.rm = TRUE)) |>
  arrange(desc(species_count))
top_pteros <- Ptero_counts[1:15,]
#repeated the same process amongst all groups

#7.1
par(mfrow = c(1,3)) #showing the 3 graphs for each group next to each other, as 1 figure
barplot(top_cephs$species_count, names.arg = top_cephs$Species_name,
        las = 2,
        xlab = "",
        ylab = "Number of Individuals",
        main = "Cephalopod Species Count",
        cex.names = 1)
#used a barplot to visualize species counts relative to the others within each groups
#plotted top 15 species with the highest counts in a barplot and lowered font of the labels
barplot(top_heteros$species_count, names.arg = top_heteros$Species_name,
        las = 2,
        xlab = "",
        ylab = "Number of Individuals",
        main = "Heteropod Species Count",
        cex.names = 1)
barplot(top_pteros$species_count, names.arg = top_pteros$Species_name,
        las = 2,
        xlab = "",
        ylab = "Number of Individuals",
        main = "Pteropod Species Count",
        cex.names = 1)
par(mfrow = c(1,1))
#repeated the barplot for all 3 groups, then set the plot settings so its back to 1
#Figure 1. shows the 15 most common species in each of the 3 groups across all the DeepSea Benefits research cruises.

#identifying the top 3 for the next figure
head(top_cephs, 3) #Pterygioteuthis sp, Abralia redfieldi, Abralia sp
head(top_heteros, 3) #Pterotrachea coronata, Pterotrachea scutata, Pterotrachea sp
head(top_pteros, 3) #Clio sp, Clio pyramidata, Diacria trispinosa


top3_gen <- DSB |>
  filter(Species_name %in% c("Pterygioteuthis sp", "Abralia redfieldi", "Cranchia scabra",
         "Pterotrachea coronata", "Pterotrachea scutata", "Pterotrachea sp", 
         "Clio sp", "Clio pyramidata", "Diacria trispinosa")) #using only the top 3 species from each group to look at common trends across the research cruises 
Dates_Density <- top3_gen |>
  group_by(Date1, Species_name, Ceph.Hetero.ptero) |>
  summarise(Sp_count_day = sum(N_number, na.rm = TRUE)) |>
  ungroup(Species_name) |>
  mutate(Species_density = ((Sp_count_day/(sum(Sp_count_day, na.rm = TRUE)))) *100) |>
  #REGROUP
  arrange(Date1)
  


library(ggplot2)

#7,2
Dates_Density |>
  ggplot(mapping = aes(x = Date1, y = Species_density, color = Species_name, shape = Ceph.Hetero.ptero)) +
  geom_point() +
  scale_x_discrete(breaks = Dates_Density$Date1) + #dates are not evenly spaced as each cruise was a year apart, so plotted this way because will swap dates for location when that data is recieved. 
  labs(x = "Date", y = "Species Percentage per Day (%)", 
      title = "Species Percentage Through Research Cruise Dates") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) #turns the labels on an angle to make more readable
#Figure 2. shows the Species Percentage per sampling day across the 3 Major Gastropod Groups.

Hetero_species <- Heteros |>
  group_by(Species_name) #groups by species to combine overlapping species from different net samples and research cruise dates

#6
max(Hetero_species$ML.TL)
min(Hetero_species$ML.TL)
mean(Hetero_species$ML.TL)

#7.3
Hetero_species |>
  ggplot(mapping = aes(x = Species_name, y = ML.TL)) +
  geom_boxplot() +
  labs(x = "Species", y = "Mantle Length/Total Length", title = "Mantle Length/Total Length Across Deep-sea Heteropod Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#Figure 3. shows Mantle Length/Total Length of different Deep Sea Heteropod species within the Gulf of America


            



