# Clean the environment
rm(list = ls())

# data: yield data of wheat genotypes sown in 14 field trials. 
# The experimental design was resolvable row-column. 
# Most of the genotypes were replicated twice. 
# A few genotypes were replicated more than twice.

# Load required libraries
# install.packages("pacman") install in case you dont have it
library(pacman)
pacman::p_load(tidyverse, statgenSTA, lme4, lmerTest, emmeans)

# Load the dataset
dat <- read_csv("./data/wheat_Australia_2008_2009_7loc.csv")
dim(dat)  # Display dimensions of the dataset
head(dat)  # Display first few rows of the dataset
glimpse(dat)  # Get a glimpse of the 'rep' variable (2 reps)

# Visualize yield data using a scatter plot
plot(dat$yield)

# Create a Trial Design (TD) object with the phenotypic data
td <- createTD(dat, genotype = "geno", trial = "Trial", loc = "Location",
               year = "Year", repId = "rep", 
               rowCoord = "row", colCoord = "col", 
               trDesign = "res.rowcol",
               trLat = "Lat", trLong = "Lon")

# Visualize the Trial Design object (list)
View(td$WH09)

# Count unique trials
length(unique(dat$Trial))

# Learn more about the createTD functions
?createTD()
unique(dat$Location)

# Plot a map of the trial design using statgenSTA
plot(td, plotType = "map", minLatRange = 50, minLongRange = 50)

# Plot map using ggplot
locations <- dat %>% select(Location, Lat, Lon) %>% distinct()

# Get world map data for reference
australia_map <- map_data("world", region = "Australia")

library(ggrepel)

# Create a ggplot object for the map
static_map <- ggplot() +
  geom_polygon(data = australia_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", col = "black", linewidth = 0.5) +
  geom_point(data = locations, aes(x = Lon, y = Lat), size = 2, color = "red") +
  geom_text_repel(data = locations, aes(x = Lon, y = Lat, label = Location), 
                  size = 2) +
  geom_text(aes(x = 130, y = -20, label = "Author: Luis Fdo. Delgado"), size = 1.5) +
  labs(x =  "Latitude", y = "Longitude") +
  labs(title = "Locations in Austria") +
  theme(plot.title = element_text(face = "bold.italic"),
        plot.subtitle = element_text(face = "italic")) +
  coord_fixed(1.3) # This sets the aspect ratio

# Display the map
print(static_map)

ggsave(paste("images\\map", ".png", sep = "_"),
       plot = static_map, units = "in", dpi = 300, width = 6, height = 5
)

# Plot two specific trial layouts
plot(td, trials = c("KA09", "WH09", "RS09"), showGeno = TRUE, 
     highlight = c("g001", "g005"))

# Boxplot using ggplot
ggplot(dat, aes(x = Trial, y = yield)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) 
   

# Boxplot of yield using statgenSTA
plot(td, plotType = "box", traits = "yield")

# Boxplot using statgenSTA (color trials by location)
plot(td, plotType = "box", colorTrialBy = "loc", traits = "yield")

# Boxplot with trials sorted by ascending yield
plot(td, plotType = "box", colorTrialBy = "loc", orderBy = "ascending", traits = "yield")

# modeling
# Convert factors
dat$rep <- as.factor(dat$rep)
dat$row <- as.factor(dat$row)
dat$col <- as.factor(dat$col)

# Filter data for trial 'WH09' + "RS09"
WH09 <- subset(dat, Trial == "WH09")
glimpse(WH09)

# Fit a mixed model using lme4
m0r <- lmer(yield ~ rep + (1|geno), data = WH09)
summary(m0r)
anova(m0r)

# Print variance components
vc0 <- VarCorr(m0r)
print(vc0, comp = "Variance", digits = 3)

# hheritability
H2 <- 0.0584 / (0.0584 + (0.3005 / 2)) 
H2

# Fit a linear model
m0 <- lm(yield ~ rep + geno, data = WH09)
anova(m0)

# Print genotype means
mean_value = lsmeans(m0, "geno")
mean_value = as.data.frame(mean_value)

# Use the means squared error (MSE) to calculate the standard error of the 
# difference (sed = sqrt(2xMSE/n))
sed_value = sqrt(2 * 0.4137 / 2)

# Fit another mixed model
m1 <- lmer(yield ~ rep + (1|row) + (1|col) + geno, data = WH09)
vc1 <- VarCorr(m1)
print(vc1, comp = "Variance", digits = 3)
anova(m1)

# Fit a mixed model with genotype as a random effect
m2 <- lmer(yield ~ rep + (1|row) + (1|col) + (1|geno), data = WH09)
vc2 <- VarCorr(m2)
print(vc2, comp = "Variance", digits = 3)
anova(m2)
table(WH09$rep)
H2_m2 <- 0.08543 / (0.08543 + (0.18105 / 2))
H2_m2


## Fit series of mixed models with SpATS
m.SpATS <- fitTD(TD = td, traits = "yield", 
                 design = "res.rowcol", 
                 what = c("fixed", "random"), 
                 spatial = TRUE, 
                 engine = "SpATS")

# Extract heritability
extractSTA(m.SpATS, what = "heritability")

# Plot spatial trends for specific trials
plot(m.SpATS, plotType = "spatial", traits = "yield", trials = c("KA09", "WH09"))
plot(m.SpATS, plotType = "spatial", spaTrend = "percentage", 
     traits = "yield", trials = c("KA09", "WH09"))

# Extract BLUEs and BLUPs
B.SpATS <- extractSTA(m.SpATS, what = c("BLUEs", "BLUPs"), keep = "trial")

# Convert results to a data frame and rename columns
BLUEs <- B.SpATS$WH09$BLUEs
colnames(BLUEs)[3] <- "BLUEs_yield"
BLUPs <- B.SpATS$WH09$BLUPs
colnames(BLUPs)[3] <- "BLUPs_yield"
WH09BB <- left_join(BLUPs, BLUEs)
head(WH09BB)

# Plot BLUES vs BLUPs
ggplot(data = WH09BB, mapping = aes(x = BLUEs_yield, y = BLUPs_yield)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color = "red")

# Calculate the range for BLUES and BLUPS
range(WH09BB$BLUEs_yield)
range(WH09BB$BLUPs_yield)

# Extract the BLUEs to a TD object
BLUEs <- STAtoTD(m.SpATS, what = c("BLUEs", "seBLUEs"), keep = "year")
head(BLUEs[[1]])

# Make boxplot and correlation plot for the BLUEs
plot(BLUEs, plotType = "box", traits = "BLUEs_yield", colorTrialBy = "year")
plot(BLUEs, plotType = "cor", traits = "BLUEs_yield")
